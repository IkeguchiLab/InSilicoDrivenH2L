from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import fegrow
from fegrow import RGroups, RLinkers
import pandas as pd
import glob
import os
from concurrent.futures import ProcessPoolExecutor

def read_sdf_file(sdf_file):
    """
    Read an SDF file and return a list of RDKit molecules.

    Parameters:
        sdf_file (str): Path to the input SDF file.

    Returns:
        list: List of RDKit molecules.
    """
    molecules = []
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    for mol in suppl:
        if mol is not None:
            molecules.append(mol)
    return molecules

def make_init_mol(mol):
    """
    Prepare the initial molecule by identifying the attachment point.

    Parameters:
        mol (RDKit Mol): The input molecule.

    Returns:
        tuple: The modified molecule and the attachment index.
    """
    init_mol = Chem.RWMol(mol)
    attachment_index = None
    for i, atom in enumerate(init_mol.GetAtoms()):
        if atom.GetSymbol() == "R":
            attachment_index = i
            atom.SetAtomicNum(0)
    return init_mol, attachment_index

def process_file(file_path, selected_rgroups, title_list):
    """
    Process a single file to generate new molecules with attached R-groups.

    Parameters:
        file_path (str): Path to the input ligand file.
        selected_rgroups (list): List of selected R-groups.
        title_list (list): List of titles for the R-groups.

    Returns:
        str: Path to the output file.
    """
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    variable_part = base_name.split('_', 1)[1]
    output_file = f"core2_{variable_part}_R2_30477.sdf"

    if os.path.exists(output_file):
        print(f"The file '{output_file}' exists.")
        return output_file

    merge_smiles_list = []
    suppl = Chem.SDMolSupplier(file_path, removeHs=False)
    mol_core = Chem.RWMol(suppl[0])
    init_mol, attachment_index = make_init_mol(mol_core)
    attachment_index = [attachment_index]
    rmols = fegrow.build_molecules(init_mol, selected_rgroups, attachment_index)
    df_rmols = rmols.dataframe
    smiles_list = df_rmols["Smiles"].to_list()
    merge_smiles_list.extend(smiles_list)

    save_mol_list = []
    if len(merge_smiles_list) != len(title_list):
        raise ValueError("The length of the lists does not match.")
    for smiles, title in zip(merge_smiles_list, title_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            title = f"{variable_part}+{title}"
            mol.SetProp('_Name', title)
            mol.SetProp('core', "core2")
            mol.SetProp('SMILES_fegrow', smiles)
            save_mol_list.append(mol)

    writer = SDWriter(output_file)
    for mol in save_mol_list:
        writer.write(mol)
    writer.close()

    return output_file

def main(max_workers=None):
    # Read selected R-groups
    selected_rgroups = read_sdf_file('./Rgroup_30477.sdf')
    title_list = [mol.GetProp("_Name") for mol in selected_rgroups]

    # File pattern for ligand files
    file_pattern = "lig_*.sdf"

    # Search for files matching the pattern
    files = glob.glob(file_pattern)

    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_file, file_path, selected_rgroups, title_list) for file_path in files]
        for future in futures:
            try:
                result = future.result()
                print(f"Processed file: {result}")
            except Exception as e:
                print(f"Error processing file: {e}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process SDF files in parallel.")
    parser.add_argument('--max_workers', type=int, default=None, help='Number of worker processes to use.')
    args = parser.parse_args()

    main(max_workers=args.max_workers)
