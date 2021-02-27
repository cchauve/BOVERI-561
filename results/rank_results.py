"""
Aggregate results from the grid parameters exploration
"""
from collections import defaultdict
import sys
import numpy as np
import pandas as pd

IN_FILE = 'v51NextSeq_commercial_samples_expected_indels_NO_WT_report_with_blacklist_out.tsv'

GRID_NG = ['4.0_64.0', '8.0_64.0']
GRID_VAF = ['0.1_0.125', '0.25', '0.375', '0.5', '0.7_0.75', '1.0']
NB_VAF = len(GRID_VAF)
GRID = [(vaf, ng) for vaf in GRID_VAF for ng in GRID_NG]
GRID_SCORE = [round(0.1 * x, 2) for x in range(1, 21)]
GRID_W_COMP = [round(0.1 * x, 2) for x in range(0, 11)]
GRID_W_SCORE_VAF = [
    (ng, score, w_comp, i)
    for ng in GRID_NG
    for score in GRID_SCORE
    for w_comp in GRID_W_COMP
    for i in range(NB_VAF)
]

PREF = 'v51NextSeq_commercial_samples_expected_indels_NO_WT_report_with_blacklist'

SETTINGS_1 = ['vaf_values']
SETTINGS_2 = ['ng_range', 'score', 'w_comp', "{:<10}".format('vaf_range')]
STATS_1 = ['TP', 'FP', 'TN', 'FN', 'FN_u', 'FN_o']
AGG_STATS_1 = [f"agg_{stat}" for stat in STATS_1]
STATS_2 = ['sens.', 'spec.', 'acc.', 'youden', 'prec.', 'recall', 'F1', 'FDR']
AGG_STATS_2 = [f"agg_{stat}" for stat in STATS_2]
STATS = STATS_1 + STATS_2

IN_DF = pd.read_csv(IN_FILE, dtype=str, sep='\t')
for col in STATS_1:
    IN_DF[col] = IN_DF[col].astype(np.int64)
for col in STATS_2 + ['score', 'w_comp']:
    IN_DF[col] = IN_DF[col].astype(np.float64)

AGG_STATS_DICT_1 = {x: defaultdict(list) for x in GRID_W_SCORE_VAF}
AGG_STATS_DICT_2 = {x: defaultdict(list) for x in GRID_W_SCORE_VAF}
for _, row in IN_DF.iterrows():
    vaf_idx = GRID_VAF.index(row['vaf_values'])
    key_1 = (row['ng_range'], row['score'], row['w_comp'], vaf_idx)
    for stat in STATS:
        AGG_STATS_DICT_1[key_1][stat].append(row[stat])
    for i in range(vaf_idx + 1):
        key_2 = (row['ng_range'], row['score'], row['w_comp'], i)
        for stat in STATS:
            AGG_STATS_DICT_2[key_2][stat].append(row[stat])

def compute_ratio(x, y, precision=3):
    """
    Computes the ratio x/y with precision precision, but if y=0 in whih case
    returns 0.0.
    """
    return (0.0 if y == 0 else  round(x / y, precision))

for key, stats_dict in AGG_STATS_DICT_1.items():
    stats_dict['mean_youden'] = round(np.mean(stats_dict['youden']), 3)
    for stat in STATS_1:
        stats_dict[f"agg_{stat}"] = sum(stats_dict[stat])
    nb_tp = stats_dict['agg_TP']
    nb_fp = stats_dict['agg_FP']
    nb_tn = stats_dict['agg_TN']
    nb_fn = stats_dict['agg_FN']
    nb_all = nb_tp + nb_fp + nb_tn + nb_fn
    stats_dict['agg_sens.'] = compute_ratio(nb_tp, nb_tp + nb_fn)
    stats_dict['agg_spec.'] = compute_ratio(nb_tn, nb_tn + nb_fp)
    stats_dict['agg_youden'] = round(stats_dict['agg_sens.'] + stats_dict['agg_spec.'] - 1.0, 3)
    stats_dict['agg_acc.'] = compute_ratio(nb_tp + nb_tn, nb_all)
    stats_dict['agg_prec.'] = compute_ratio(nb_tp, nb_tp + nb_fp) # = PPV
    precision = stats_dict['agg_prec.']
    stats_dict['agg_recall'] = stats_dict['agg_sens.']
    recall = stats_dict['agg_recall']
    stats_dict['agg_F1'] = compute_ratio(2.0 * precision * recall, precision + recall)
    stats_dict['agg_FDR'] = round(1.0 - precision, 3)

for key, stats_dict in AGG_STATS_DICT_2.items():
    stats_dict['mean_youden'] = round(np.mean(stats_dict['youden']), 3)
    for stat in STATS_1:
        stats_dict[f"agg_{stat}"] = sum(stats_dict[stat])
    nb_tp = stats_dict['agg_TP']
    nb_fp = stats_dict['agg_FP']
    nb_tn = stats_dict['agg_TN']
    nb_fn = stats_dict['agg_FN']
    nb_all = nb_tp + nb_fp + nb_tn + nb_fn
    stats_dict['agg_sens.'] = compute_ratio(nb_tp, nb_tp + nb_fn)
    stats_dict['agg_spec.'] = compute_ratio(nb_tn, nb_tn + nb_fp)
    stats_dict['agg_youden'] = round(stats_dict['agg_sens.'] + stats_dict['agg_spec.'] - 1.0, 3)
    stats_dict['agg_acc.'] = compute_ratio(nb_tp + nb_tn, nb_all)
    stats_dict['agg_prec.'] = compute_ratio(nb_tp, nb_tp + nb_fp) # = PPV
    precision = stats_dict['agg_prec.']
    stats_dict['agg_recall'] = stats_dict['agg_sens.']
    recall = stats_dict['agg_recall']
    stats_dict['agg_F1'] = compute_ratio(2.0 * precision * recall, precision + recall)
    stats_dict['agg_FDR'] = round(1.0 - precision, 3)

SORTED_KEYS_1 = list(AGG_STATS_DICT_1.keys())
SORTED_KEYS_1.sort(key=lambda x: AGG_STATS_DICT_1[x]['agg_youden'], reverse=True)

SORTED_KEYS_2 = list(AGG_STATS_DICT_2.keys())
SORTED_KEYS_2.sort(key=lambda x: AGG_STATS_DICT_2[x]['agg_youden'], reverse=True)

HEADER = [x.replace('agg_', '') for x in SETTINGS_2 + AGG_STATS_1 + AGG_STATS_2 + ['mean_youden']]

print('\t'.join(HEADER))
for key in SORTED_KEYS_1:
    stats_dict = AGG_STATS_DICT_1[key]
    out_fields = [key[0], key[1], key[2], "{:<10}".format(f"{GRID_VAF[key[3]]}")]
    out_fields += [stats_dict[stat] for stat in AGG_STATS_1 + AGG_STATS_2]
    out_fields += [stats_dict['mean_youden']]
    print('\t'.join([str(x) for x in out_fields]))

for key in SORTED_KEYS_2:
    stats_dict = AGG_STATS_DICT_2[key]
    out_fields = [key[0], key[1], key[2], "{:<10}".format(f">={GRID_VAF[key[3]]}")]
    out_fields += [stats_dict[stat] for stat in AGG_STATS_1 + AGG_STATS_2]
    out_fields += [stats_dict['mean_youden']]
    print('\t'.join([str(x) for x in out_fields]))
