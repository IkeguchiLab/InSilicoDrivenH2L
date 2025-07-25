import sys
import os
import glob
import pandas as pd
import statistics

proj_home=sys.argv[1]
results_csv=sys.argv[2]
print(proj_home)

edge_file=f"{proj_home}workpath/transition/edgelist.txt"
df=pd.read_csv(results_csv)

edge_list=[]
with open(edge_file, 'r') as f:
    for line in f:
        line=line.rstrip()
        line=line.split(",")
        edge_line="edge_"+line[0]+"_"+line[1]
        edge_list.append(edge_line)

print(edge_list)

summary_list=[]
for edge in edge_list:
    filtered_df = df[df["edge_run"].str.contains(edge)]

    if filtered_df.empty:
        continue

    filtered_df_protein = filtered_df[filtered_df["edge_run"].str.contains('protein')]
    filtered_df_water = filtered_df[filtered_df["edge_run"].str.contains('water')]
    print(filtered_df_protein)
    print(filtered_df_water)

    ddG_list=[]
    for run in ['run1', 'run2', 'run3']:
        protein_run = filtered_df_protein[filtered_df_protein["edge_run"].str.contains(run)]
        water_run = filtered_df_water[filtered_df_water["edge_run"].str.contains(run)]
        
        # skip the calculation when no. num_forward < 10 or no. num_backward < 10
        if not protein_run.empty and not water_run.empty:
            if protein_run['num_forward'].values[0] < 10 or protein_run['num_backward'].values[0] < 10:
                continue
            if water_run['num_forward'].values[0] < 10 or water_run['num_backward'].values[0] < 10:
                continue
            
            ddG = protein_run['CGI_dG'].values[0] - water_run['CGI_dG'].values[0]
            ddG_list.append(ddG)

    if len(ddG_list) == 0:
        continue
    try:
        mean_value = statistics.mean(ddG_list)
    except:
        print(f'{edge} is not succeeded.')
        continue
    try:
        variance_value = statistics.variance(ddG_list)
    except:
        variance_value = 0
    try:
        standard_deviation = statistics.stdev(ddG_list)
    except:
        standard_deviation = 0
    print(ddG_list)
    print("mean: ", mean_value)
    print("variance: ", variance_value)
    print("standard_dev: ", standard_deviation)
    _summary_list=[edge, mean_value, variance_value, standard_deviation, len(ddG_list)]
    summary_list.append(_summary_list)

df_summary = pd.DataFrame(summary_list, columns=["edge_run", "ddG_mean", "ddG_variance", "ddG_standard_dev", "num_run"])

df_summary.to_csv(f"{proj_home}workpath/summary/results_summary.csv")
