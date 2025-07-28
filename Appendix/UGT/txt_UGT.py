import pandas as pd
import os

columns_to_check = ["UGT1A1","UGT1A3","UGT1A4","UGT1A6","UGT1A8","UGT1A9","UGT1A10","UGT2B7","UGT2B15"]
txt_files = [f for f in os.listdir('.') if f.endswith('PDE9A_core4_22_R2_AP.txt')]

results = []

for txt_file in txt_files:
    print(txt_file)
    df = pd.read_csv(txt_file, delimiter='\t')
    df['UGT_Yes'] = df[columns_to_check].apply(lambda row: sum('Yes' in str(val) for val in row), axis=1)
    output_df = df[['*molname', 'UGT_Yes']]
    results.append(output_df)

final_df = pd.concat(results, ignore_index=True)

output_file = 'merged_core4_22_R2_AP_UGT.csv'

final_df.to_csv(output_file, index=False, sep='\t')
