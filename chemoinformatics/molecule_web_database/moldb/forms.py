from django import forms
from .models import Molecule

class AddMoleculeForm(forms.Form):
    name = forms.CharField(label="Name", )
    smiles = forms.CharField(label="SMILES")

    class Meta:
        model = Molecule