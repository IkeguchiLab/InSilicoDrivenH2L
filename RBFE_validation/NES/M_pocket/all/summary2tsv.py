import pandas as pd

df = pd.read_csv('results_summary_core4_Mpocket_cycle.csv')

ligands = df['edge_run'].str.extract(r'(lig_\d+-?\d*)_(lig_\d+-?\d*)')
df['ligand_i'] = ligands[0]
df['ligand_j'] = ligands[1]

df_output = df[['ligand_i', 'ligand_j', 'ddG_mean', 'ddG_standard_dev']].copy()
df_output.columns = ['ligand_i', 'ligand_j', 'DDG(i->j) (kcal/mol)', 'uncertainty (kcal/mol)']

df_output["uncertainty (kcal/mol)"] = df_output["uncertainty (kcal/mol)"].replace(0, 0.01)

df_output.to_csv("final_results_core4_Mpocket_cycle.tsv", sep='\t', index=False)
