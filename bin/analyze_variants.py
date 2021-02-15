import csv
import os
import pandas as pd
import sys

RUN_CSV_FILE = sys.argv[1]
ALL_INDELS_DF = pd.read_csv('data/v4MiSeq_commercial_samples_expected_indels.csv')
ALL_INDELS_DF.rename(columns={'chromosome': 'chr'}, inplace=True)

CTRL_VAF_DIFF = float(sys.argv[2])

RUN_ID_LIST = []
with open(RUN_CSV_FILE) as csvfile:
        runs_data = csv.reader(csvfile, delimiter=',')
        for row in runs_data:
            RUN_ID_LIST.append(row[1])

def difference_df(df1, df2):
    df1_aux = df1[['sample', 'chr', 'pos', 'ref', 'alt']]
    df2_aux = df2[['sample', 'chr', 'pos', 'ref', 'alt']]
    df2_m_df1 = pd.concat([df2_aux, df1_aux, df1_aux]).drop_duplicates(keep=False)
    df1_m_df2 = pd.concat([df1_aux, df2_aux, df2_aux]).drop_duplicates(keep=False)
    return list(df1_m_df2.index), list(df2_m_df1.index)

def find_indel(sample, chr, pos, ref, alt, df):
    return list(df.loc[(df['sample']==sample) &
                       (df['chr']==chr) &
                       (df['pos']==pos) &
                       (df['ref']==ref) &
                       (df['alt']==alt)
                       ].index)

def compare_indels(true_df, run_df_in, score, verbose=False):
    run_df = run_df_in.loc[run_df_in['score']<=score]
    index_true = list(true_df.index)
    index_run = list(run_df.index)
    index_only_true, index_only_run = difference_df(true_df, run_df)
    index_tp = [x for x in index_true if x not in index_only_true]
    index_fn = index_only_true
    index_fp = [x for x in index_run if x in index_only_run]
    nb_tp, nb_fn, nb_fp = len(index_tp), len(index_fn), len(index_fp)
    if verbose:
        stats = f"TP:{nb_tp}\tFN:{nb_fn}\tFP:{nb_fp}\tSPEC:{specificity}"
        print(f"\tscore_min:{score}\t{stats}")
        for _, row in true_df[true_df.index.isin(index_fn)].iterrows():
            print(f"\tFN\t{row['sample']}\t{row['chr']}\t{row['pos']}\t{row['ref']}\t{row['alt']}")
        for _, row in run_df[run_df.index.isin(index_fp)].iterrows():
            out_str = (
                f"\tFP\t{row['sample']}\t{row['chr']}\t{row['pos']}\t{row['ref']}\t{row['alt']}\t"
                f"{row['vaf']}\t{row['score']}\t{row['complexity']}\t"
                f"{row['support']}\t{row['overlap']}\t{row['control']}\t"
                f"{row['repeats']}"
            )
            print(out_str)
    fp_df = run_df.loc[index_fp]
    fp_mean_vaf = round(fp_df['vaf'].mean(), 2)
    fp_mean_score = round(fp_df['score'].mean(), 2)
    nb_fp_1 = len(fp_df.loc[fp_df['vaf']<1.0].index)
    nb_fp_05 = len(fp_df.loc[fp_df['vaf']<0.5].index)
    index_tp = [x for x in index_run if x not in index_only_run]
    tp_df = run_df.loc[index_tp]
    tp_mean_vaf = round(tp_df['vaf'].mean(), 2)
    tp_mean_score = round(tp_df['score'].mean(), 2)
    nb_tp_1 = len(tp_df.loc[tp_df['vaf']<1.0].index)
    nb_tp_05 = len(tp_df.loc[tp_df['vaf']<0.5].index)
    return nb_tp, nb_fn, nb_fp, fp_mean_vaf, fp_mean_score, nb_fp_1, nb_fp_05, tp_mean_vaf, tp_mean_score, nb_tp_1, nb_tp_05


def assess_indels(true_df, run_df, run_id, out_file):
    result = {}
    NB_TP, NB_FN, NB_FP = 0, 0, 0
    for index, row in true_df.iterrows():
        indel_info = [run_id, row['sample'], row['chr'], row['pos'], row['ref'], row['alt']]
        indel_str = '\t'.join([str(x) for x in indel_info])
        index_run = find_indel(row['sample'], row['chr'], row['pos'], row['ref'], row['alt'], run_df)
        if len(index_run) == 0:
            result[index] = None
            stat_str = f"nan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan\tnan"
        else:
            indel_run = run_df.loc[index_run[0]]
            score = indel_run['score']
            vaf = indel_run['vaf']
            nb_tp, nb_fn, nb_fp, fp_vaf, fp_score, nb_fp_1, nb_fp_05, tp_vaf, tp_score, nb_tp_1, nb_tp_05 = compare_indels(true_df, run_df, score)
            stat_str = '\t'.join([str(x) for x in [vaf, score, nb_tp, tp_vaf, tp_score, nb_tp_1, nb_tp_05, nb_fn, nb_fp, fp_vaf, fp_score, nb_fp_1, nb_fp_05]])
            result[index] = (score, nb_tp, nb_fn, nb_fp)
            NB_TP += nb_tp
            NB_FN += nb_fn
            NB_FP += nb_fp
        out_file.write(f"\n{indel_str}\t{stat_str}")
    return NB_TP, NB_FN, NB_FP

# Analysis: expected indels by minimum score to call them
_, DATASET_NAME = os.path.split(RUN_CSV_FILE)
OUT_FILE_1 = os.path.join('results', DATASET_NAME.replace('.csv', '_out_1.tsv'))
OUT_1 = open(OUT_FILE_1, 'w')
header = 'run_id\tsample\tchr\tpos\tref\talt\tvaf\tscore\tTP:nb\tTP:mean_vaf\t'
header += 'TP:mean_score\tTP:nb:vaf<1\tTP:nb:vaf<0.5\tFN:nb\tFP:nb\tFP:mean_vaf\t'
header += 'FP:mean_score\tFP:nb:vaf<1\tFP:nb:vaf<0.5'
OUT_1.write(header)
for RUN_ID in RUN_ID_LIST:
    RUN_INDELS_FILE = os.path.join('results', RUN_ID, f"{RUN_ID}_indels.tsv")
    RUN_INDELS_DF_AUX = pd.read_csv(RUN_INDELS_FILE, sep='\t')
    RUN_INDELS_DF_AUX['sample'] = RUN_INDELS_DF_AUX.apply(
        lambda row: row['sample'].split('_S')[0], axis=1
    )
    RUN_INDELS_DF = RUN_INDELS_DF_AUX[~RUN_INDELS_DF_AUX['sample'].str.startswith('QMRS')].round(3)
    TRUE_INDELS_DF = ALL_INDELS_DF.loc[ALL_INDELS_DF['run_id']==RUN_ID]
    RUN_INDELS_DF_CTRL = RUN_INDELS_DF.loc[RUN_INDELS_DF['control']>CTRL_VAF_DIFF]

    _, _, _ = assess_indels(TRUE_INDELS_DF, RUN_INDELS_DF_CTRL, RUN_ID, OUT_1)
OUT_1.close()

# Analysis: statistics by minimum score
OUT_FILE_2 = os.path.join('results', DATASET_NAME.replace('.csv', '_out_2.tsv'))
OUT_2 = open(OUT_FILE_2, 'w')
header = 'run_id\tmax_score\tTP:nb\tTP:mean_vaf\tTP:mean_score\tTP:nb:vaf<1\tTP:nb:vaf<0.5\t'
header += 'FN:nb\tFP:nb\tFP:mean_vaf\tFP:mean_score\tFP:nb:vaf<1\tFP:nb:vaf<0.5'
OUT_2.write(header)
for RUN_ID in RUN_ID_LIST:
    RUN_INDELS_FILE = os.path.join('results', RUN_ID, f"{RUN_ID}_indels.tsv")
    RUN_INDELS_DF_AUX = pd.read_csv(RUN_INDELS_FILE, sep='\t')
    RUN_INDELS_DF_AUX['sample'] = RUN_INDELS_DF_AUX.apply(
        lambda row: row['sample'].split('_S')[0], axis=1
    )
    RUN_INDELS_DF = RUN_INDELS_DF_AUX[~RUN_INDELS_DF_AUX['sample'].str.startswith('QMRS')].round(3)
    TRUE_INDELS_DF = ALL_INDELS_DF.loc[ALL_INDELS_DF['run_id']==RUN_ID]
    RUN_INDELS_DF_CTRL = RUN_INDELS_DF.loc[RUN_INDELS_DF['control']>CTRL_VAF_DIFF]
    for score in [0.25, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]:
        nb_tp, nb_fn, nb_fp, fp_vaf, fp_score, nb_fp_1, nb_fp_05, tp_vaf, tp_score, nb_tp_1, nb_tp_05 = compare_indels(TRUE_INDELS_DF, RUN_INDELS_DF_CTRL, score)
        stat_str = '\n' + '\t'.join([str(x) for x in [RUN_ID, score, nb_tp, tp_vaf, tp_score, nb_tp_1, nb_tp_05, nb_fn, nb_fp, fp_vaf, fp_score, nb_fp_1, nb_fp_05]])
        OUT_2.write(stat_str)
OUT_2.close()
