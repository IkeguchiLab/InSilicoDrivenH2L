import sys
import os
import glob
import pandas as pd

proj_home=sys.argv[1]
print(proj_home)
results_file="results_re.txt"
results_file_paths=sorted(glob.glob(f'{proj_home}workpath/transition/edge_*/analyse_*/{results_file}', recursive=True))

run_dir = None
edge_dir = None

out=[]

for results_file_path in results_file_paths:
    path_comonents = results_file_path.split("/")
    for directory in path_comonents:
        if "analyse" in directory:
            run_dir = directory
            run_dir = run_dir.replace("analyse_", "")
        if "edge" in directory:
            edge_dir = directory
    print(edge_dir+"_"+run_dir)
    edge_run=edge_dir+"_"+run_dir

    with open(results_file_path, 'r') as f:
        CGI_dG=0
        CGI_Std_Err_parametric=0
        CGI_Std_Err_bootstrip=0
        BAR_dG=0
        BAR_Std_Err_bootstrip=0
        BAR_Std_Err_analytical=0
        num_forward=0
        num_backward=0
        for l in f:
            l = l.rstrip()
            foo = l.split()
            if 'CGI: dG' in l:
                CGI_dG=float(foo[-2])
            elif 'CGI: Std Err (bootstrap:parametric)' in l:
                CGI_Std_Err_parametric=float(foo[-2])
            elif 'CGI: Std Err (bootstrap)'in l:
                CGI_Std_Err_bootstrip=float(foo[-2])
            elif 'BAR: dG' in l:
                BAR_dG=float(foo[-2])
            elif 'BAR: Std Err (bootstrap)' in l:
                BAR_Std_Err_bootstrip=float(foo[-2])
            elif 'BAR: Std Err (analytical)' in l:
                BAR_Std_Err_analytical=float(foo[-2])
            elif '0->1' in l:
                num_forward=int(foo[-1])
            elif '1->0' in l:
                num_backward=int(foo[-1])
    result_list=[edge_dir+"_"+run_dir, CGI_dG, CGI_Std_Err_parametric, CGI_Std_Err_bootstrip, BAR_dG, BAR_Std_Err_bootstrip, BAR_Std_Err_analytical, num_forward, num_backward]
    out.append(result_list)
    print(result_list)

df = pd.DataFrame(out, columns=["edge_run", "CGI_dG", "CGI_Std_Err_parametric", "CGI_Std_Err_bootstrip", "BAR_dG", "BAR_Std_Err_bootstrip", "BAR_Std_Err_analytical", "num_forward", "num_backward"])

df.to_csv(f"{proj_home}workpath/summary/results.csv")
