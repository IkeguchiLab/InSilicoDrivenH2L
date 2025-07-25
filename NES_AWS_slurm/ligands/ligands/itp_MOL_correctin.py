import os

current_directory = os.getcwd()

for dir_name in os.listdir(current_directory):
    if dir_name.startswith("lig_") and os.path.isdir(dir_name):
        mol_itp_path = os.path.join(dir_name, "MOL.itp")
        if os.path.exists(mol_itp_path):
            with open(mol_itp_path, 'r') as file:
                lines = file.readlines()

            lines = lines[:-3]

            with open(mol_itp_path, 'w') as file:
                file.writelines(lines)

            print(f"Updated: {mol_itp_path}")
print("All done")
