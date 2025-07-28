from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import fegrow
from fegrow import RGroups, RLinkers
import pandas as pd

def read_sdf_file(sdf_file):
    """
    Read an SDF file and return a list of RDKit molecules.

    Parameters:
        sdf_file (str): Path to the input SDF file.

    Returns:
        list: List of RDKit molecules.
    """
    molecules = []
    suppl = Chem.SDMolSupplier(sdf_file,removeHs=False)
    for mol in suppl:
        if mol is not None:
            molecules.append(mol)
    return molecules

def make_init_mol(mol):
    init_mol = Chem.RWMol(mol)
    for i,atom in enumerate(init_mol.GetAtoms()):
        if atom.GetAtomicNum() ==0:
            attachment_index = i
            atom.SetAtomicNum(0)
    return init_mol, attachment_index

selected_rgroups = read_sdf_file('Rgroup_select_30504.sdf')
title_list = [mol.GetProp("_Name") for mol in selected_rgroups]

merge_smiles_list=[]
suppl = Chem.SDMolSupplier("core4_tol_R1.sdf")
print(suppl)
mol_core = Chem.RWMol(suppl[0])
init_mol, attachment_index = make_init_mol(mol_core)
attachment_index = [attachment_index]
rmols = fegrow.build_molecules(init_mol, selected_rgroups, attachment_index)
df_rmols = rmols.dataframe
smiles_list = df_rmols["Smiles"].to_list()
merge_smiles_list.extend(smiles_list)

save_mol_list = []
if len(merge_smiles_list) != len(title_list):
    raise ValueError("length error")
for smiles, title in zip(merge_smiles_list, title_list):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol.SetProp('_Name', title)
        save_mol_list.append(mol)

writer = SDWriter('core4_tol_R1_Rgroup_select_30477.sdf')
for mol in save_mol_list:
    writer.write(mol)
writer.close()


