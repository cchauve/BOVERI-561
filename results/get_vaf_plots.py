"""
Generate figures of scatter plot of TP/FN indels expected VAF versus VAF
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys

EXT = sys.argv[1]

GRID = [
    'VAF_NG_0.0_0.5_0.0_4.0', 'VAF_NG_0.0_0.5_4.0_8.0',
    'VAF_NG_0.0_0.5_8.0_16.0', 'VAF_NG_0.0_0.5_16.0_32.0',
    'VAF_NG_0.0_0.5_32.0_100.0',
    'VAF_NG_0.5_1.0_0.0_4.0', 'VAF_NG_0.5_1.0_4.0_8.0',
    'VAF_NG_0.5_1.0_8.0_16.0', 'VAF_NG_0.5_1.0_16.0_32.0',
    'VAF_NG_0.5_1.0_32.0_100.0',
    'VAF_NG_1.0_5.0_0.0_4.0', 'VAF_NG_1.0_5.0_4.0_8.0',
    'VAF_NG_1.0_5.0_8.0_16.0', 'VAF_NG_1.0_5.0_16.0_32.0',
    'VAF_NG_1.0_5.0_32.0_100.0'
]

PREF = 'v4MiSeq_commercial_samples_expected_indels_ExpectedVAFs_MissingSeraSeq'

for setting in GRID:
    setting_split = setting.replace('VAF_NG_', '').split('_')
    vaf_min = setting_split[0]
    vaf_max = setting_split[1]
    ng_min = setting_split[2]
    ng_max = setting_split[3].replace('100.0', '64.0')
    title = f"VAF=({vaf_min},{vaf_max}] DNA=({ng_min},{ng_max}]"
    in_file = '_'.join([PREF, setting, EXT, 'vaf.tsv'])
    out_file = in_file.replace('.tsv', '.png')
    in_df = pd.read_csv(in_file, dtype=str, sep='\t')
    in_df['vaf'] = in_df['vaf'].astype(np.float64)
    in_df['exp_vaf'] = in_df['exp_vaf'].astype(np.float64)
    correlation = in_df['vaf'].corr(in_df['exp_vaf'], method='pearson')
    plt.figure()
    vaf_plot = sns.scatterplot(
        data=in_df, x='vaf', y='exp_vaf', hue='status',
        palette={'TP': 'blue', 'FN': 'orange'}
    )
    plt.title(f"{title}, Pearson correlation={round(correlation, 3)}")
    plt.xticks(rotation=90)
    # fig = vaf_plot.get_figure()
    plt.savefig(out_file)
