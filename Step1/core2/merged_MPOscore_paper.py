import pandas as pd
import numpy as np

# ファイルの読み込み
file_path = 'merged_core2_R2_30477_AP.txt'
df = pd.read_csv(file_path, delimiter='\t')

# 新しいカラムの作成
df['S+Sw_uM'] = df['S+Sw'] * 1000000 / df['MWt']
df['log10_S+Sw_uM'] = np.log10(df['S+Sw_uM'])
df['log10_CYP_HLM_CLint'] = np.log10(df['CYP_HLM_CLint'])

# point_S+Peffの計算
df['point_S+Peff'] = df['S+Peff'].apply(lambda x: 0.33 * x if x <= 3 else 1)

# point_log10_CYP_HLM_CLintの計算
def calculate_point_log10_CYP_HLM_CLint(value):
    if value <= 1.17:
        return 1
    elif 1.17 < value <= 1.9:
        return (-value / 0.73) + (1.9 / 0.73)
    else:
        return 0

df['point_log10_CYP_HLM_CLint'] = df['log10_CYP_HLM_CLint'].apply(calculate_point_log10_CYP_HLM_CLint)

# point_log10_S+Sw_uMの計算
def calculate_point_log10_S_Sw_uM(value):
    if value <= 1:
        return 0
    elif 1 < value <= 1.95:
        return (value / 0.95) + 1 - (1.95 / 0.95)
    else:
        return 1

df['point_log10_S+Sw_uM'] = df['log10_S+Sw_uM'].apply(calculate_point_log10_S_Sw_uM)

# MPO_scoreの計算
df['MPO_score'] = df['point_S+Peff'] + df['point_log10_CYP_HLM_CLint'] + df['point_log10_S+Sw_uM']

# ビン分けと件数・割合の計算
def bin_and_calculate(df, column, bins):
    bin_labels = range(len(bins) - 1)
    df[f'{column}_bin'] = pd.cut(df[column], bins=bins, labels=bin_labels, include_lowest=True)
    bin_counts = df[f'{column}_bin'].value_counts().sort_index()
    bin_percentages = bin_counts / len(df) * 100
    return bin_counts, bin_percentages

# ビンの設定
bins_logP = [-np.inf, 0, 1, 2, 3, 4, np.inf]
bins_Sw_uM = [-np.inf, 10, 50, 90, np.inf]
bins_CYP_HLM_CLint = [-np.inf, 15, 30, 80, np.inf]
bins_Peff = [0, 1, 2, 3, np.inf]
bins_sum_points = [-np.inf, 1.0, 1.5, 2, 2.5, 3]

# ビン分けと結果の表示
columns_to_bin = ['S+logP', 'S+Sw_uM', 'CYP_HLM_CLint', 'S+Peff', 'MPO_score']
bins = [bins_logP, bins_Sw_uM, bins_CYP_HLM_CLint, bins_Peff, bins_sum_points]

for column, bin_limits in zip(columns_to_bin, bins):
    counts, percentages = bin_and_calculate(df, column, bin_limits)
    print(f'Column: {column}')
    print('Counts:')
    print(counts)
    print('Percentages:')
    print(percentages)
    print()

# 結果の保存
#df.to_csv('processed_merged_core1_R2_30477_AP.txt', sep='\t', index=False)
