import csv
import glob
import numpy as np
from collections import defaultdict

csv_files = glob.glob('fep_mapper*.csv')

fep_data = defaultdict(list)

for csv_file in csv_files:
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            ligand_i = int(row['Ligand1'])
            ligand_j = int(row['Ligand2'])
            fep_value = float(row['FEP'])
            
            fep_data[(ligand_i, ligand_j)].append(fep_value)

results = []

for (ligand_i, ligand_j), fep_values in fep_data.items():
    fep_mean = np.mean(fep_values)
    fep_std = np.std(fep_values)
    
    if fep_std == 0:
        fep_std = 0.01
    
    results.append((ligand_i, ligand_j, fep_mean, fep_std))

with open('final_results_core4_cycle_FEP.tsv', 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['ligand_i', 'ligand_j', 'DDG(i->j) (kcal/mol)', 'uncertainty (kcal/mol)'])
    
    for result in results:
        writer.writerow(result)
