"""
Analyze FP and FN for the optimal thresholds
"""
from collections import defaultdict
import sys
import numpy as np
import pandas as pd

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
SUFF = 'optimal_thresholds_with_blacklist'
FEATURES = ['vaf', 'comp.', 'supp.', 'overlap', 'ctrl']
ERROR_TYPES = ['FP', 'FN_u', 'FN_d']

header = '\t'.join(
    ['freq.', 'chr', 'pos', 'ref', 'alt'] +
    [f"avg_{x}\tmax_{x}" for x in FEATURES]
)
FEATURES_VAL = {
    error_type: {feature: defaultdict(list) for feature in FEATURES}
    for error_type in ERROR_TYPES
}
FEATURES_AVG = {
    error_type: {feature: defaultdict(float) for feature in FEATURES}
    for error_type in ERROR_TYPES
}
FEATURES_MAX = {
    error_type: {feature: defaultdict(float) for feature in FEATURES}
    for error_type in ERROR_TYPES
}
VARIANTS_COUNT = {error_type: defaultdict(int) for error_type in ERROR_TYPES}
OUT_FILES = {
    error_type: f"{PREF}_{SUFF}_{error_type}.txt" for error_type in ERROR_TYPES
}


for setting in GRID:
    in_file = '_'.join([PREF, setting, SUFF, 'errors.txt'])
    print(in_file)
    in_df = pd.read_csv(in_file, dtype=str, sep='\t')
    for _, row in in_df.iterrows():
        v_key = '_'.join([
            str(x) for x in [row['chr'], row['pos'], row['ref'], row['alt']]
        ])
        error_type = row['w_comp']
        print(error_type)
        # Bug in output: LLOD used as index, error type not in header
        VARIANTS_COUNT[error_type][v_key] += 1
        for feature in FEATURES:
            FEATURES_VAL[error_type][feature][v_key].append(float(row[feature]))

for error_type in ERROR_TYPES:
    out_file = open(OUT_FILES[error_type], 'w')
    out_file.write(header)
    VARIANTS_KEYS = list(VARIANTS_COUNT[error_type])
    VARIANTS_KEYS.sort(
        key=lambda x: VARIANTS_COUNT[error_type][x],
        reverse=True
    )
    for v_key in VARIANTS_KEYS:
        for feature in FEATURES:
            FEATURES_AVG[error_type][feature][v_key] = (
                round(np.mean(FEATURES_VAL[error_type][feature][v_key]), 3)
            )
            FEATURES_MAX[error_type][feature][v_key] = (
                round(np.max(FEATURES_VAL[error_type][feature][v_key]), 3)
            )
        out_str_frequency = [str(VARIANTS_COUNT[error_type][v_key])]
        out_str_variant = v_key.split('_')
        out_str_features = [
            f"{FEATURES_AVG[error_type][feature][v_key]}\t{FEATURES_MAX[error_type][feature][v_key]}"
            for feature in FEATURES
        ]
        out_str = '\t'.join(
            out_str_frequency + out_str_variant + out_str_features)
        out_file.write(f"\n{out_str}")
    out_file.close()
