from django.core.exceptions import ObjectDoesNotExist
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from django_rdkit import models
import datetime
from django.contrib.auth.models import User

class Molecule(models.Model):
    """
    Represents one molecule.
    """

    # fields, which can be calculated on save():
    rdmol = models.MolField()
    internal_id = models.CharField(max_length=32, db_index=True)
    image_svg = models.TextField(null=True)
    mw = models.FloatField(db_index=True)
    sum_formula = models.CharField(max_length=32, db_index=True)
    fingerprint = models.CharField(max_length=1024, db_index=True)
    inchi = models.TextField(db_index=True)
    inchi_key = models.CharField(max_length=27, db_index=True)
    name = models.TextField(db_index=True, null=True)
    smiles = models.TextField(db_index=True)
    amount = models.FloatField()
    created = models.DateTimeField(auto_now_add=True)

    # excluded molecules SMILES (they cause rdKit stuck)
    EXCLUDED_MOLECULES = ["C", "CH3", "CH4", "[CH3]", "[C]", "[CH4]"]

    def __str__(self):
        return "Molecule ({id}): '{name}', formula: '{formula}'".format(id=self.internal_id, name=self.name, formula=self.sum_formula)

    def save(self, smiles=None, molfile=None, rdmol=None, inchi=None, name=None, update=False, *args, **kwargs):
        if not update:
            if molfile:
                mol = AllChem.MolFromMolBlock(molfile)
            elif smiles:
                mol = AllChem.MolFromSmiles(smiles)
            elif rdmol:
                mol = rdmol
            elif inchi:
                mol = AllChem.MolFromInchi(inchi)

            if mol:
                inchi = AllChem.MolToInchi(mol)
                smiles = AllChem.MolToSmiles(mol)

                if inchi and Molecule.objects.filter(inchi=inchi).count() == 0 and len(inchi) > 1:
                    self.inchi = inchi

                    self.mw = float("{0:.2f}".format(AllChem.CalcExactMolWt(mol)))
                    self.sum_formula = AllChem.CalcMolFormula(mol)
                    self.fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=1024).ToBitString()
                    self.inchi_key = AllChem.InchiToInchiKey(self.inchi)
                    self.molfile = AllChem.MolToMolBlock(mol)
                    self.smiles = smiles
                    self.rdmol = mol

                    # generating SVG image
                    if self.smiles not in self.EXCLUDED_MOLECULES:
                        binMol = AllChem.Mol(self.rdmol.ToBinary())

                        if not binMol.GetNumConformers():
                            rdDepictor.Compute2DCoords(self.rdmol)

                        drawer = rdMolDraw2D.MolDraw2DSVG(100, 100)
                        drawer.DrawMolecule(self.rdmol)
                        drawer.FinishDrawing()
                        svg = drawer.GetDrawingText().replace('svg:', '')

                        # remove first line containg XML meta information
                        self.image_svg = "\n".join(svg.split("\n")[1:]).strip()
                    else:
                        self.image_svg = None

                    if name:
                        self.name = name
                    else:
                        try:
                            self.name = mol.GetProp("LONGNAME")
                        except KeyError:
                            self.name= None

                    if Molecule.objects.all().count() == 0:
                        self.internal_id = "MI-J-1"
                    else:
                        self.internal_id = "MI-J-{}".format(Molecule.objects.latest("id").id + 1)

                    super(Molecule, self).save(*args, **kwargs)
                else:
                    raise self.MoleculeExistsInDatabase(smiles)
            else:
                raise self.MoleculeCreationError
        else:
            super(Molecule, self).save(*args, **kwargs)

    class Meta:
        ordering = ['id']

    class MoleculeExistsInDatabase(Exception):
        def __init__(self, smiles):
            super(Exception, self).__init__(smiles)
            self.smiles = smiles
            self.message = "Cannot add the molecule: it already exists in database."

    class MoleculeCreationError(Exception):
        def __init__(self):
            super(Exception, self).__init__()
            self.message = "Cannot add the molecule: check your structure (valence etc.)."

class Building(models.Model):
    """
    Represents one building.
    """

    name = models.CharField(max_length=2, unique=True)

    def __str__(self):
        return "Building ({}): {}".format(self.id, self.name)

class Room(models.Model):
    """
    Represents one room.
    """

    code = models.CharField(max_length=8, unique=True)

    def __str__(self):
        return "Room ({}): {}".format(self.id, self.code)

class Place(models.Model):
    """
    Represents one place which consists of room and building.
    """

    building = models.ForeignKey(Building)
    room = models.ForeignKey(Room)

    class Meta:
        unique_together = ("building", "room")

    def __str__(self):
        return "Place ({}): {}/{}".format(self.id, self.building.name, self.room.code)

class Person(models.Model):
    """
    Represents one person.
    """

    user = models.OneToOneField(User, on_delete=models.CASCADE)
    name = models.CharField(max_length=64)
    surname = models.CharField(max_length=128)
    place = models.ForeignKey(Place)
    telephone = models.CharField(max_length=16)

    def __str__(self):
        return "Person ({}): {}, {} @ {}".format(self.id, self.surname, self.name, str(self.place))

    class Meta:
        ordering = ['id']

class OrderStatus(models.Model):
    """
    Represents the order status.
    """

    status = models.CharField(max_length=32, unique=True)

    def __str__(self):
        return "OrderStatus ({}): {}".format(self.id, self.status)

class Order(models.Model):
    """
    Represents one order.
    """

    compounds = models.ManyToManyField(Molecule, through="OrderCompound")
    person = models.ForeignKey(Person)
    status = models.ForeignKey(OrderStatus)
    order_internal_id = models.CharField(max_length=32, unique=True)
    created = models.DateTimeField(auto_now_add=True)

    def save(self, update=True, *args, **kwargs):
        if not update:
            today_date = datetime.date.today()
            today_min = datetime.datetime.combine(datetime.date.today(), datetime.time.min)
            today_max = datetime.datetime.combine(datetime.date.today(), datetime.time.max)

            try:
                today_last_order = Order.objects.filter(created__range=(today_min, today_max)).latest("id")
                self.order_internal_id = "{:%Y-%m-%d}/{}".format(datetime.date.today(), today_last_order.id + 1)
            except ObjectDoesNotExist:
                self.order_internal_id = "{:%Y-%m-%d}/1".format(today_date)

        super(Order, self).save(*args, **kwargs)

    def __str__(self):
        return "Order ({} / {}): for {} | compounds: {} | {}".format(self.id,
                                                                     self.order_internal_id,
                                                                     str(self.person),
                                                                     self.compounds.count(),
                                                                     str(self.status))

    class Meta:
        ordering = ['id']

class OrderCompound(models.Model):
    """
    Join table for order and its compounds.
    """

    compound = models.ForeignKey(Molecule)
    order = models.ForeignKey(Order)
    amount = models.FloatField()

    def __str__(self):
        return "OrderCompound ({}): compound {} in order {} ".format(self.id, self.compound.id, self.order.id)

class PurchaseStatus(models.Model):
    """
    Represents the purchase status.
    """

    status = models.CharField(max_length=32, unique=True)

    def __str__(self):
        return "PurchaseStatus ({}): {}".format(self.id, self.status)

class Purchase(models.Model):
    """
    Represents one purchase.
    """

    compounds = models.ManyToManyField(Molecule, through="PurchaseCompound")
    person = models.ForeignKey(Person)
    status = models.ForeignKey(PurchaseStatus)
    purchase_internal_id = models.CharField(max_length=32, unique=True)
    created = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return "Purchase ({} / {}): for {} | compounds: {} | {}".format(self.id,
                                                                     self.purchase_internal_id,
                                                                     str(self.person),
                                                                     self.compounds.count(),
                                                                     str(self.status))

    def save(self, update=True, *args, **kwargs):
        if not update:
            today_date = datetime.date.today()
            today_min = datetime.datetime.combine(datetime.date.today(), datetime.time.min)
            today_max = datetime.datetime.combine(datetime.date.today(), datetime.time.max)

            try:
                today_last_purchase = Purchase.objects.filter(created__range=(today_min, today_max)).latest("id")
                print(today_last_purchase.id)
                self.purchase_internal_id = "{:%Y-%m-%d}/{}".format(datetime.date.today(), today_last_purchase.id + 1)
            except ObjectDoesNotExist:
                self.purchase_internal_id = "{:%Y-%m-%d}/1".format(today_date)

        super(Purchase, self).save(*args, **kwargs)

    class Meta:
        ordering = ['id']

class PurchaseCompound(models.Model):
    """
    Join table for purchase and its compounds.
    """

    compound = models.ForeignKey(Molecule)
    purchase = models.ForeignKey(Purchase)
    amount = models.FloatField()

    def __str__(self):
        return "PurchaseCompound ({}): compound {} in purchase {} ".format(self.id, self.compound.id, self.purchase.id)