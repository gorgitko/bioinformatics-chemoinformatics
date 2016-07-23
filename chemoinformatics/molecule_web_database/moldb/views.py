# TODO: IUPAC naming from Chemdoodle web or some Python module

from django.shortcuts import render, render_to_response
from django.template.loader import render_to_string
from moldb.models import *
from django.template import RequestContext
from .forms import AddMoleculeForm
from django.http import JsonResponse, HttpResponse
from rdkit import Chem
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
import json
from django.utils.datastructures import MultiValueDictKeyError
from django.core.files.temp import NamedTemporaryFile
from wsgiref.util import FileWrapper
from rdkit.Chem import AllChem
import os
import random
from subprocess import Popen
import time
import logging
from django.db.models import Q
from django.db.utils import DataError
from psycopg2 import DataError as DataError_psycopg2
import operator
from functools import reduce
from django_rdkit.models import QMOL, Value
from .settings import *
from .search_paginator import SearchPaginator
from random import randint
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
from django.http import HttpResponseRedirect

logger = logging.getLogger(__name__)

def home(request):
    return render(request,
        'home.html',
        {"molecules_count": Molecule.objects.count()})

@login_required()
def add_molecules(request):
    #request.session["uploaded_mols"] = []
    #request.session["last_added_mols_index"] = 0
    #request.session["upload_finished"] = False

    return render(request,
        'add_molecules.html',
        {"form": AddMoleculeForm()})

def list_molecules(request):
    paginator = SearchPaginator(Molecule,
                                Molecule.objects.all().values_list("id", flat=True),
                                MOLECULES_PER_LIST_PAGE)

    try:
        if "page" in request.GET.keys():
            page = paginator.get_page(request.GET["page"])
        else:
            page = paginator.get_page(1)
    except SearchPaginator.PageDoesNotExist:
        page = None

    return render(request,
        'list_molecules.html',
        {"mols": page, "all_mols_count": paginator.objects_count, "paginator": page, "type": "list"})

def search_molecules(request):
    return render(request,
                  'search_molecules.html',
                  {})

@login_required()
def purchases(request):
    try:
        purchase_progress = Purchase.objects.get(person__id=request.user.id, status__status="in progress")
        purchase_progress.purchased_compounds = []

        for purchase_compound in purchase_progress.purchasecompound_set.all():
            compound = purchase_compound.compound
            compound.purchased_amount = purchase_compound.amount
            purchase_progress.purchased_compounds.append(compound)
    except ObjectDoesNotExist:
        purchase_progress = None

    purchases_finished = Purchase.objects.filter(person__id=request.user.id, status__status="finished")
    purchases_finished_modified = []

    for purchase_finished in purchases_finished:
        purchase_finished_modified = purchase_finished
        purchase_finished_modified.purchased_compounds = []
        for purchase_compound in purchase_finished_modified.purchasecompound_set.all():
            compound = purchase_compound.compound
            compound.purchased_amount = purchase_compound.amount
            purchase_finished_modified.purchased_compounds.append(compound)
        purchases_finished_modified.append(purchase_finished_modified)

    return render(request,
                  'purchases.html',
                  {"purchase_progress": purchase_progress, "purchases_finished": purchases_finished_modified})

@login_required()
def orders(request):
    try:
        order_progress = Order.objects.get(person__id=request.user.id, status__status="in progress")
        order_progress.ordered_compounds = []

        for order_compound in order_progress.ordercompound_set.all():
            compound = order_compound.compound
            compound.ordered_amount = order_compound.amount
            order_progress.ordered_compounds.append(compound)
    except ObjectDoesNotExist:
        order_progress = None

    orders_finished = Order.objects.filter(person__id=request.user.id, status__status="finished")
    orders_finished_modified = []

    for order_finished in orders_finished:
        order_finished_modified = order_finished
        order_finished_modified.ordered_compounds = []
        for order_compound in order_finished_modified.ordercompound_set.all():
            compound = order_compound.compound
            compound.ordered_amount = order_compound.amount
            order_finished_modified.ordered_compounds.append(compound)
        orders_finished_modified.append(order_finished_modified)

    return render(request,
                  'orders.html',
                  {"order_progress": order_progress, "orders_finished": orders_finished_modified})

@login_required()
def logout_user(request):
    logout(request)
    return HttpResponseRedirect('/')

def login_user(request):
    logout(request)

    if request.POST:
        username = request.POST['username']
        password = request.POST['password']

        user = authenticate(username=username, password=password)
        if user is not None and user.is_active:
            login(request, user)
            return HttpResponseRedirect('/')

    return render(request, "login.html")

def test(request):
    return render(request,
                  'temp_bootstrap.html',
                  {"mols": Molecule.objects.all()})

# API

def api_molConverter(request):
    if request.POST:
        format_from = request.POST["format_from"]
        format_to = request.POST["format_to"]
        data = request.POST["data"]

        if format_from == "molfile":
            mol = Chem.MolFromMolBlock(data)
        elif format_from == "smiles":
            mol = Chem.MolFromSmiles(data)
        else:
            return JsonResponse(addMoleculeDictSerialize(False, error="Input format unknown."))

        if mol:
            if format_to == "molfile":
                output = Chem.MolToMolBlock(mol)
            elif format_to == "smiles":
                output = Chem.MolToSmiles(mol)
            else:
                return JsonResponse(addMoleculeDictSerialize(False, error="Output format unknown."))

            if output:
                return JsonResponse({"success": True, "data": output})
        else:
            return JsonResponse(
                addMoleculeDictSerialize(False,
                                         error="Cannot convert from '{}' to '{}': probably invalid data supplied."
                                         .format(format_from, format_to)))

@login_required()
def api_addMolecule(request):
    if request.POST and request.POST["molfile"]:
        return JsonResponseStatus("success", data=[saveMol(molfile=request.POST["molfile"], name=request.POST["name"])])

@login_required()
def api_uploadMolecules(request):
    #request.session["upload_finished"] = False
    #request.session["uploaded_mols"] = []
    #request.session["last_added_mols_index"] = 0

    try:
        filetype = request.POST["filetype"]
    except MultiValueDictKeyError:
        return JsonResponseStatus("error", message="No 'filetype' in POST parameters.")

    try:
        file = request.FILES["file"]
    except MultiValueDictKeyError:
        return JsonResponseStatus("error", message="No file supplied (missing 'file' in FILE parameters).")

    data = []

    try:
        if filetype == "smiles":
            for chunk in file.chunks():
                for line in str(chunk, encoding="utf-8").split("\n"):
                    data.append(saveMol(smiles=line))
                    #request.session["uploaded_mols"].append(saveMol(smiles=line))
        elif filetype == "sdf":
            for mol in AllChem.SDMolSupplier(file.temporary_file_path()):
                data.append(saveMol(rdmol=mol))
                #request.session["uploaded_mols"].append(saveMol(rdmol=mol))
        else:
            return JsonResponseStatus("error", message="Invalid filetype.")
    except Exception as e:
        logger.error(str(e))
        return JsonResponseStatus("error", message="Invalid file supplied.")

    return JsonResponseStatus("success", data=data)

    #request.session["upload_finished"] = True
    #return JsonResponse({"success": True})

@login_required()
def api_uploadMoleculesStatus(request):
    if "last_added_mols_index" in request.session.keys():
        try:
            new_added_mols = request.session["uploaded_mols"][request.session["last_added_mols_index"]:]
            request.session["last_added_mols_index"] = len(new_added_mols)
            response = {"success": True,
                        "mol_number": len(request.session["uploaded_mols"]),
                        "uploaded_mols": new_added_mols}

            if "upload_finished" in request.session.keys():
                if request.session["upload_finished"]:
                    response["upload_finished"] = True
                    request.session["upload_finished"] = False
                    request.session["uploaded_mols"] = []
                    request.session["last_added_mols_index"] = 0

            return JsonResponse(response)
        except IndexError:
            return JsonResponse({"success": False, "error": "IndexError"})

def api_downloadMolecules(request):
    folder = "static/temp/"

    if "download_all" in request.GET.keys():
        mols = Molecule.objects.all()
        filename = "molecules_all_{}_{}.sdf".format(time.strftime("%Y-%m-%d"), random.randint(1, 10000))
    else:
        mol_ids = [int(x) for x in request.GET.getlist("mol_ids[]")]
        mols = Molecule.objects.filter(id__in=mol_ids)
        filename = "molecules_{}_{}.sdf".format(time.strftime("%Y-%m-%d"), random.randint(1, 10000))

    path = folder + filename

    writer = AllChem.SDWriter(path)
    for mol in mols:
        writer.write(AllChem.MolFromSmiles(mol.smiles))
    writer.close()

    with open(path, mode="r", encoding="utf-8") as f:
        response = HttpResponse(FileWrapper(f), content_type='application/download')
        response['Content-Disposition'] = 'attachment; filename={}'.format(filename)
        return response

def api_searchMoleculesByStructure(request):
    if "search_paginators" not in request.session.keys():
        request.session["search_paginators"] = []

    if "page" in request.POST.keys():
        return search_pagination(request, page=request.POST["page"], paginator_number=request.POST["paginator_number"])
    else:
        try:
            mol = Chem.MolFromMolBlock(request.POST["molfile"])

            if mol:
                if request.POST["type"] == "exact":
                    mols_list = Molecule.objects.filter(inchi=Chem.MolToInchi(mol))
                elif request.POST["type"] == "substructure":
                    mols_list = Molecule.objects.filter(rdmol__hassubstruct=Chem.MolToSmiles(mol))

                paginator = SearchPaginator(Molecule,
                                            mols_list.values_list("id", flat=True),
                                            MOLECULES_PER_SEARCH_PAGE)
                request.session["search_paginators"].append({
                    "paginator": paginator,
                    "query": mols_list.query
                })

                return search_pagination(request,
                                         paginator=paginator,
                                         paginator_number=len(request.session["search_paginators"]) - 1)
            else:
                raise Molecule.MoleculeCreationError
        except Molecule.MoleculeCreationError as e:
            return JsonResponseStatus("error", message=e.message)

def api_searchMoleculesByFilter(request):
    if "search_paginators" not in request.session.keys():
        request.session["search_paginators"] = []

    if "page" in request.POST.keys():
        return search_pagination(request, page=request.POST["page"], paginator_number=request.POST["paginator_number"])
    else:
        search_fields = ["rdmol", "sum_formula", "inchi", "inchi_key", "smarts", "mw", "mw__lt", "mw__gt"]
        search_type = request.POST["search_type"]

        fields = request.POST.getlist("fields[]")
        values = request.POST.getlist("values[]")
        unknown_fields = []
        search_q = Q()

        if fields and values:
            if "smarts" in fields:
                smarts_index = fields.index("smarts")
                search_q = Q(rdmol__hassubstruct=QMOL(Value(values.pop(smarts_index))))
                del fields[smarts_index]

            if "name" in fields:
                name_index = fields.index("name")
                search_q = Q(name__contains=values.pop(name_index))
                del fields[name_index]

            for field in fields:
                if field not in search_fields:
                    unknown_fields.append(field)

            if unknown_fields:
                return JsonResponseStatus("error",
                                      message="Unknown search fields: {}. Allowed fields are {}".format(
                                          ["'{}'".format(x) for x in unknown_fields],
                                          ["'{}'".format(y) for y in search_fields]
                                      ))

            predicates = list(zip(fields, values))
            q_list = [Q(x) for x in predicates]

            if search_type == "or":
                if fields and values:
                    search_q = search_q | reduce(operator.or_, q_list)
            elif search_type == "and":
                if fields and values:
                    search_q = search_q & reduce(operator.and_, q_list)

            try:
                mols_list = Molecule.objects.filter(search_q)

                paginator = SearchPaginator(Molecule,
                                            mols_list.values_list("id", flat=True),
                                            MOLECULES_PER_SEARCH_PAGE)
                request.session["search_paginators"].append({
                    "paginator": paginator,
                    "query": mols_list.query
                })

                return search_pagination(request,
                                         paginator=paginator,
                                         paginator_number=len(request.session["search_paginators"]) - 1)
            except DataError_psycopg2 as e:
                return JsonResponseStatus("error", message=str(e))
            except DataError as e:
                return JsonResponseStatus("error", message=str(e))
        else:
            return JsonResponseStatus("error", message="No search fields supplied.")

@login_required()
def api_addOrder(request):
    try:
        amount = float(request.POST["amount"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'amount' value in POST")

    try:
        mol_id = int(request.POST["mol_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'mol_id' value in POST")

    message = addOrder(request, amount, mol_id)
    request.session["order_progress"] = True

    return JsonResponseStatus("success", message=message)

@login_required()
def api_addPurchase(request):
    try:
        amount = float(request.POST["amount"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'amount' value in POST")

    try:
        mol_id = int(request.POST["mol_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'mol_id' value in POST")

    purchase_only = request.POST["purchase_only"]

    if purchase_only == "true":
        purchase_only = True
    elif purchase_only == "false":
        purchase_only = False
    else:
        return JsonResponseStatus("error", message="Invalid 'purchase_only' value in POST")

    amount_db = Molecule.objects.get(id=mol_id).amount
    message_order = ""

    if amount > amount_db and not purchase_only:
        amount_order = amount - amount_db
        amount_purchase = amount_db
        amount_new = 0
        message_order = addOrder(request, amount_order, mol_id)
    else:
        amount_purchase = amount
        amount_new = amount_db - amount

    try:
        purchase_db = Purchase.objects.get(status__status="in progress", person=request.user.id)

        try:
            purchase_compound_db = purchase_db.purchasecompound_set.get(compound=mol_id)
            purchase_compound_db.amount += amount_purchase
            purchase_compound_db.save()
        except ObjectDoesNotExist:
            purchase_compound_db = purchase_db.purchasecompound_set.create(
                compound=Molecule.objects.get(id=mol_id),
                purchase=purchase_db,
                amount=amount_purchase
            )
        finally:
            mol = Molecule.objects.get(id=mol_id)
            mol.amount -= amount_purchase
            mol.save(update=True)

            request.session["purchase_progress"] = True

            return JsonResponseStatus("success",
                                      message="Added {}\u03BCl to current 'in progress' purchase (total {}\u03BCl). "\
                                        .format(amount_purchase, purchase_compound_db.amount) + message_order,
                                      data={"new_amount": amount_new})
    except ObjectDoesNotExist:
        purchase_db = Purchase(
            person=Person.objects.get(id=request.user.id),
            status=PurchaseStatus.objects.get(id=1),
        )

        mol = Molecule.objects.get(id=mol_id)

        purchase_db.save(update=False)
        purchase_db.purchasecompound_set.create(
            compound=mol,
            purchase=purchase_db,
            amount=amount
        )

        mol.amount -= amount_purchase
        mol.save(update=True)

        request.session["purchase_progress"] = True

        return JsonResponseStatus("success",
                                  message="Created new purchase and added {} \u03BCl to it. "
                                    .format(amount_purchase) + message_order,
                                  data={"new_amount": amount_new})

@login_required()
def api_rmMolFromInProgressPurchase(request):
    try:
        mol_id = int(request.POST["mol_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'mol_id' value in POST")

    try:
        purchase_id = int(request.POST["purchase_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'purchase_id' value in POST")

    try:
        purchase_progress = Purchase.objects.get(status__status="in progress", person=request.user.id, id=purchase_id)
    except ObjectDoesNotExist:
        return JsonResponseStatus("error", message="Cannot remove molecule from purchase: purchase not in progress.")

    try:
        purchase_compound = purchase_progress.purchasecompound_set.get(compound=mol_id)
        purchase_compound.compound.amount += purchase_compound.amount
        purchase_compound.compound.save(update=True)
        purchase_compound.delete()
        return JsonResponseStatus("success", message="Successfully removed from current in progress purchase.")
    except ObjectDoesNotExist:
        return JsonResponseStatus("error", message="Cannot remove molecule from purchase: molecule not in purchase.")

@login_required()
def api_rmInProgressPurchase(request):
    try:
        purchase_id = int(request.POST["purchase_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'purchase_id' value in POST")

    try:
        purchase_progress = Purchase.objects.get(status__status="in progress", person=request.user.id, id=purchase_id)

        for purchase_compound in purchase_progress.purchasecompound_set.all():
            purchase_compound.compound.amount += purchase_compound.amount
            purchase_compound.compound.save(update=True)
            purchase_compound.delete()

        purchase_progress.delete()
        request.session["purchase_progress"] = False
        return JsonResponseStatus("success", message="Successfully removed current in progress purchase.")
    except ObjectDoesNotExist:
        return JsonResponseStatus("error", message="Cannot remove molecule from purchase: purchase not in progress.")

@login_required()
def api_finishInProgressPurchase(request):
    try:
        purchase_id = int(request.POST["purchase_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'purchase_id' value in POST")

    try:
        purchase_progress = Purchase.objects.get(status__status="in progress", person=request.user.id, id=purchase_id)
        purchase_progress.status = PurchaseStatus.objects.get(id=2)
        purchase_progress.save()
        request.session["purchase_progress"] = False
        return JsonResponseStatus("success", message="Successfully finished current in progress purchase.")
    except ObjectDoesNotExist:
        return JsonResponseStatus("error", message="Cannot finish the purchase: purchase not in progress.")

@login_required()
def api_rmMolFromInProgressOrder(request):
    try:
        mol_id = int(request.POST["mol_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'mol_id' value in POST")

    try:
        order_id = int(request.POST["order_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'order_id' value in POST")

    try:
        order_progress = Order.objects.get(status__status="in progress", person=request.user.id, id=order_id)
    except ObjectDoesNotExist:
        return JsonResponseStatus("error", message="Cannot remove molecule from order: order not in progress.")

    try:
        order_compound = order_progress.ordercompound_set.get(compound=mol_id)
        order_compound.compound.amount += order_compound.amount
        order_compound.compound.save(update=True)
        order_compound.delete()
        return JsonResponseStatus("success", message="Successfully removed from current in progress order.")
    except ObjectDoesNotExist:
        return JsonResponseStatus("error", message="Cannot remove molecule from order: molecule not in order.")

@login_required()
def api_rmInProgressOrder(request):
    try:
        order_id = int(request.POST["order_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'order_id' value in POST")

    try:
        order_progress = Order.objects.get(status__status="in progress", person=request.user.id, id=order_id)

        for order_compound in order_progress.ordercompound_set.all():
            order_compound.compound.amount += order_compound.amount
            order_compound.compound.save(update=True)
            order_compound.delete()

        order_progress.delete()
        request.session["order_progress"] = False
        return JsonResponseStatus("success", message="Successfully removed current in progress order.")
    except ObjectDoesNotExist:
        return JsonResponseStatus("error", message="Cannot remove molecule from order: order not in progress.")

@login_required()
def api_finishInProgressOrder(request):
    try:
        order_id = int(request.POST["order_id"])
    except ValueError:
        return JsonResponseStatus("error", message="Invalid 'order_id' value in POST")

    try:
        order_progress = Order.objects.get(status__status="in progress", person=request.user.id, id=order_id)
        order_progress.status = OrderStatus.objects.get(id=2)
        order_progress.save()

        for order_compound in order_progress.ordercompound_set.all():
            order_compound.compound.amount += order_compound.amount
            order_compound.compound.save(update=True)

        request.session["order_progress"] = False
        return JsonResponseStatus("success", message="Successfully finished current in progress order.")
    except ObjectDoesNotExist:
        return JsonResponseStatus("error", message="Cannot finish the order: order not in progress.")

# helper functions

def saveMol(smiles=None, molfile=None, rdmol=None, name=None):
    if smiles or molfile or rdmol:
        mol = Molecule(amount=randint(1, 1000))

        try:
            mol.save(smiles=smiles, molfile=molfile, rdmol=rdmol, name=name)
        except Molecule.MoleculeExistsInDatabase as e:
            return addMoleculeDictSerialize(False, error=e.message, smiles=e.smiles)
        except Molecule.MoleculeCreationError as e:
            if smiles:
                return addMoleculeDictSerialize(False, error=e.message, smiles=smiles)
            else:
                return addMoleculeDictSerialize(False, error=e.message, smiles="-")

        return addMoleculeDictSerialize(True, internal_id=mol.internal_id, smiles=mol.smiles)
    else:
        return addMoleculeDictSerialize(False, error="Cannot add empty molecule.", smiles="-")

def addMoleculeDictSerialize(success, internal_id=None, error=None, smiles=None):
    if success:
        response = {"success": True, "internal_id": internal_id}
    else:
        response = {"success": False, "error": error}

    if smiles:
        response.update({"smiles": smiles})

    return response

def JsonResponseStatus(status, data=None, message=None):
    return JsonResponse({"status": status, "message": message, "data": data})

def search_pagination(request, page=1, paginator=None, paginator_number=0):
    paginator_number = int(paginator_number)

    if not paginator:
        paginator = request.session["search_paginators"][paginator_number]["paginator"]

    try:
        page = paginator.get_page(page)
        return JsonResponseStatus("success", data={
            "table": render_to_string(
                '_molecules_table.html',
                context={"mols": page.objects_list,
                         "type": "search",
                         "paginator": page,
                         "paginator_number": paginator_number},
                request=request),
            "mols_count": paginator.objects_count,
        })
    except SearchPaginator.PageDoesNotExist:
        return JsonResponseStatus("error", message="No molecules found in database.")

def addOrder(request, amount, mol_id):
    try:
        order_db = Order.objects.get(Q(status__status="in progress") & Q(person=request.user.id))

        try:
            order_compound_db = order_db.ordercompound_set.get(compound=mol_id)
            order_compound_db.amount += amount
            order_compound_db.save()
        except ObjectDoesNotExist:
            order_compound_db = order_db.ordercompound_set.create(
                compound=Molecule.objects.get(id=mol_id),
                order=order_db,
                amount=amount
            )
        finally:
            return "Added {}\u03BCl to current 'in progress' order (total {}\u03BCl)."\
                .format(amount, order_compound_db.amount)
    except ObjectDoesNotExist:
        order_db = Order(
            person=Person.objects.get(id=request.user.id),
            status=OrderStatus.objects.get(id=1),
        )

        order_db.save(update=False)
        order_db.ordercompound_set.create(
            compound=Molecule.objects.get(id=mol_id),
            order=order_db,
            amount=amount
        )

        return "Created new order and added {} \u03BCl to it.".format(amount)