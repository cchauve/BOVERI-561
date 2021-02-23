"""
Extract the best results from the grid parameters exploration
"""
import collections
import sys
import numpy as np
import pandas as pd

EXT = sys.argv[1]
FACTOR = float(sys.argv[2])
COL1 = sys.argv[3]
COL2 = sys.argv[4]
MIN_FREQUENCY = int(sys.argv[5])

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

header = ['LLOD', 'score']
header += ['w_comp', 'TP', 'FP', 'TN', 'FN_u', 'FN_d']
header += ['sens.', 'spec.', 'acc.', 'prec.', 'recall', 'F1', 'FDR']

COMBINATIONS = []
for setting in GRID:
    setting_split = setting.replace('VAF_NG_', '').split('_')
    vaf_min = setting_split[0]
    vaf_max = setting_split[1]
    ng_min = setting_split[2]
    ng_max = setting_split[3].replace('100.0', '64.0')
    in_file = '_'.join([PREF, setting, EXT, 'out.tsv'])
    in_df = pd.read_csv(in_file, dtype=str, sep='\t')
    in_df_columns = list(in_df.columns)
    in_df[COL1] = in_df[COL1].astype(np.float64)
    in_df[COL2] = in_df[COL2].astype(np.float64)
    COL1_max = in_df[COL1].max()
    COL2_max = in_df[COL2].max()
    print(f"\n### VAF=({vaf_min},{vaf_max}] DNA=({ng_min},{ng_max}]\n#### FILE={in_file}")
    print('```')
    print('\t'.join(header))
    for _, row in in_df.iterrows():
        if (row[COL1] >= FACTOR * COL1_max) and (row[COL2] >= FACTOR * COL2_max):
            info = [row[c] for c in in_df_columns]
            print('\t'.join([str(x) for x in info]))
            COMBINATIONS.append((row['score'], row['w_comp']))
    print('```')

print('\n### Combinations frequencies\n```\nscore\tw_comp\tfrequency')
FREQUENCIES = collections.Counter(COMBINATIONS)
FREQUENCIES_KEYS = list(FREQUENCIES.keys())
FREQUENCIES_KEYS.sort(key=lambda c: FREQUENCIES[c])
for key in FREQUENCIES_KEYS:
    if FREQUENCIES[key] >= MIN_FREQUENCY:
        print(f"{key[0]}\t{key[1]}\t{FREQUENCIES[key]}")
print('```')
