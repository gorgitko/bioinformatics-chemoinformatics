import os
import django
from django.core.exceptions import ObjectDoesNotExist
from random import sample

os.environ["DJANGO_SETTINGS_MODULE"] = "_settings"
django.setup()

from moldb.models import *

def add_places(places):
    places_added = []

    for place in places:
        try:
            building_db = Building.objects.get(name=place["b_name"])
        except ObjectDoesNotExist:
            building_db = Building(name=place["b_name"])
            building_db.save()

        try:
            room_db = Room.objects.get(code=place["r_code"])
        except ObjectDoesNotExist:
            room_db = Room(code=place["r_code"])
            room_db.save()

        try:
            place_db = Place.objects.get(building__name=place["b_name"], room__code=place["r_code"])
        except ObjectDoesNotExist:
            place_db = Place(building=building_db, room=room_db)
            place_db.save()

        places_added.append(place_db)

    return places_added

def add_persons(persons, places):
    persons_added = []

    for i, person in enumerate(persons):
        try:
            person_db = Person.objects.get(name=person["name"], surname=person["surname"])
        except ObjectDoesNotExist:
            user_db = User.objects.create_user(person["username"], email=person["email"], password=person["password"])

            person_db = Person(name=person["name"], surname=person["surname"],
                               telephone=person["telephone"], place=places[i], user=user_db)
            person_db.save()

        persons_added.append(person_db)

    return persons_added

def add_molecules(molecules):
    molecule_added_ids = []

    for mol in molecules:
        try:
            mol_db = Molecule.objects.get(smiles=mol["smiles"])
        except ObjectDoesNotExist:
            mol_db = Molecule(amount=mol["amount"])
            mol_db.save(smiles=mol["smiles"])

        molecule_added_ids.append(mol_db)

    return molecule_added_ids

def add_order_statuses(order_statuses):
    order_statuses_added_ids = []

    for status in order_statuses:
        try:
            status_db = OrderStatus.objects.get(status=status["status"])
        except ObjectDoesNotExist:
            status_db = OrderStatus(status=status["status"])
            status_db.save()

        order_statuses_added_ids.append(status_db)

    return order_statuses_added_ids

def add_purchase_statuses(purchase_statuses):
    purchase_statuses_added_ids = []

    for status in purchase_statuses:
        try:
            status_db = PurchaseStatus.objects.get(status=status["status"])
        except ObjectDoesNotExist:
            status_db = PurchaseStatus(status=status["status"])
            status_db.save()

        purchase_statuses_added_ids.append(status_db)

    return purchase_statuses_added_ids

def add_orders(orders):
    for order in orders:
        order_db = Order(person=order["person"], status=order["status"])
        order_db.save(update=False)

        for comp in order["compounds"]:
            order_compound_db = OrderCompound(compound=Molecule.objects.get(id=comp["mol_id"]),
                                              order=order_db,
                                              amount=comp["amount"])

            order_compound_db.save()

buildings = [
    {"name": "A"},
    {"name": "B"},
    {"name": "C"}
]

rooms = [
    {"code": "211"},
    {"code": "250"},
    {"code": "347"},
    {"code": "201"},
    {"code": "233"},
    {"code": "147"},
    {"code": "101"},
    {"code": "136"},
]

places = [
    {"b_name": "A", "r_code": "211"},
    {"b_name": "A", "r_code": "120"},
    {"b_name": "A", "r_code": "324"},
    {"b_name": "B", "r_code": "340"},
    {"b_name": "B", "r_code": "225"},
    {"b_name": "B", "r_code": "229"}
]

persons = [
    {"name": "Hana", "surname": "Dohnalova", "telephone": "775425746",
     "email": "hannicka12@seznam.cz", "username": "hanicka", "password": "hanicka"},
    {"name": "Eda", "surname": "Ehler", "telephone": "798625112",
     "email": "eda@seznam.cz", "username": "eda", "password": "eda"},
    {"name": "Slavek", "surname": "Tretjacenko", "telephone": "656874996",
     "email": "slavek@gmail.com", "username": "slavek", "password": "slavek"},
    {"name": "Daniel", "surname": "Svozil", "telephone": "778986332",
     "email": "dan.svozil@gmail.com", "username": "dan", "password": "dan"}
]

molecules = [
    {"smiles": "CCCC", "amount": 10},
    {"smiles": "CCc1ccccc1", "amount": 50},
    {"smiles": "Cl", "amount": 1000},
    {"smiles": "FCCc1ccccc1", "amount": 5},
    {"smiles": "ClCCCCCCF", "amount": 20},
    {"smiles": "FCC=O", "amount": 20},
    {"smiles": "FF", "amount": 20}
]

order_statuses = [
    {"status": "in progress"},
    {"status": "finished"}
]

purchase_statuses = [
    {"status": "in progress"},
    {"status": "finished"}
]

#building_ids = add_buildings(buildings)
#room_ids = add_rooms(rooms)
#print(building_ids, room_ids)

places_added = add_places(places)
persons_added = add_persons(persons, sample(places_added, len(persons)))
molecule_ids = add_molecules(molecules)
order_status_ids = add_order_statuses(order_statuses)
purchase_status_ids = add_purchase_statuses(purchase_statuses)

orders = [
    {
        "status": OrderStatus.objects.get(status="in progress"),
        "person": persons_added[0],
        "compounds": [
            {"mol_id": 1, "amount": 5},
            {"mol_id": 2, "amount": 50},
            {"mol_id": 3, "amount": 10}
        ]
    },
        {
        "status": OrderStatus.objects.get(status="finished"),
        "person": persons_added[1],
        "compounds": [
            {"mol_id": 3, "amount": 20},
            {"mol_id": 4, "amount": 30},
            {"mol_id": 5, "amount": 60}
        ]
    }
]

#print(orders)

add_orders(orders)