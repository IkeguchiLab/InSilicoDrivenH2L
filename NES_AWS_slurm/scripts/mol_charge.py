import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops

mol_pdb = sys.argv[1]
mol1 = Chem.MolFromPDBFile(mol_pdb, removeHs=False)
charge = rdkit.Chem.rdmolops.GetFormalCharge(mol1)
print(int(charge))
