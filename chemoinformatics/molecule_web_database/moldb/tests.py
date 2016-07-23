from django.test import TestCase
from .models import Molecule
import json
from rdkit.Chem import AllChem

class AddMoleculeTestCase(TestCase):
    def setUp(self):
        self.propane = """Molecule from ChemDoodle Web Components

http://www.ichemlabs.com
  3  2  0  0  0  0            999 V2000
   -0.9330   -0.2500    0.0000 C   0  0  0  0  0  0
    0.0670   -0.2500    0.0000 C   0  0  0  0  0  0
    0.9330    0.2500    0.0000 C   0  0  0  0  0  0
  1  21  0     0  0
  2  31  0     0  0
M  END"""

    def test_add(self):
        mol = Molecule()
        mol.save(molfile=self.propane)

        self.assertEqual(mol.internal_id, "MI-J-1")
        self.assertEqual(mol.smiles, "CCC")
        self.assertEqual(mol.sum_formula, "C3H8")
        self.assertEqual(mol.mw, 44.06)
        self.assertEqual(mol.inchi, "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3")
        self.assertEqual(mol.inchi_key, "ATUOYWHBWRKTHZ-UHFFFAOYSA-N")

    def test_api_addMolecule(self):
        response = self.client.post(path="/api/addMolecule", data={"molfile": self.propane})
        self.assertEqual(response.status_code, 200)

        mol = AllChem.MolFromMolBlock(self.propane)
        mol_added = Molecule.objects.last()

        self.assertEqual(float("{0:.2f}".format(AllChem.CalcExactMolWt(mol))), mol_added.mw)
        self.assertEqual(AllChem.MolToSmiles(mol), mol_added.smiles)
        self.assertEqual(AllChem.CalcMolFormula(mol), mol_added.sum_formula)
        inchi = AllChem.MolToInchi(mol)
        self.assertEqual(inchi, mol_added.inchi)
        self.assertEqual(AllChem.InchiToInchiKey(inchi), mol_added.inchi_key)

    def test_api_molConverter(self):
        response = self.client.post(path="/api/molConverter", data={"data": self.propane,
                                                                "format_from": "molfile",
                                                               "format_to": "smiles"})
        self.assertEqual(response.status_code, 200)
        data = json.loads(str(response.content, encoding="utf-8"))
        self.assertEqual(data["data"], "CCC")

        response = self.client.post(path="/api/molConverter", data={"data": "NotValidMolfile",
                                                                "format_from": "molfile",
                                                                "format_to": "smiles"})
        self.assertEqual(response.status_code, 200)
        data = json.loads(str(response.content, encoding="utf-8"))
        self.assertTrue(data["error"])

    def test_api_uploadMolecules(self):
        with open('moldb/test_data/smiles_test.txt') as f:
            response = self.client.post('/api/uploadMolecules', {'file': f})

        self.assertEqual(response.status_code, 200)
        data = json.loads(str(response.content, encoding="utf-8"))

        first_mol = data[0]
        second_mol = data[1]

        self.assertEqual(first_mol["internal_id"], "MI-J-1")
        self.assertEqual(first_mol["smiles"], "FCCCl")
        self.assertTrue(first_mol["success"])

        self.assertTrue(second_mol["error"])
        self.assertEqual(second_mol["smiles"], "C#C=C")
        self.assertFalse(second_mol["success"])