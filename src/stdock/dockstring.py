# Created by roy.gonzalez-aleman at 20/11/2023

import rdkit.Chem as Chem
from dockstring import load_target

smiles = 'CC1=C(C(=O)N2CCCCC2=N1)CCN3CCC(CC3)C4=NOC5=C4C=CC(=C5)F'
Chem.MolFromSmiles(smiles)
target = load_target('DRD2')
score, aux = target.dock(smiles)
