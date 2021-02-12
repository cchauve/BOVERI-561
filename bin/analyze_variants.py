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
    df2_minus_df1 = pd.concat([df2_aux, df1_aux, df1_aux]).drop_duplicates(keep=False)
    df1_minus_df2 = pd.concat([df1_aux, df2_aux, df2_aux]).drop_duplicates(keep=False)
    return list(df1_minus_df2.index), list(df2_minus_df1.index)

def find_indel(sample, chr, pos, ref, alt, df):
    return list(df.loc[(df['sample']==sample) & (df['chr']==chr) & (df['pos']==pos) & (df['ref']==ref) & (df['alt']==alt)].index)

def compare_indels(true_df, run_df_in, vaf, score, verbose=False):
    run_df = run_df_in.loc[(run_df_in['vaf']>=vaf) & (run_df_in['score']<=score)]
    index_true = list(true_df.index)
    index_run = list(run_df.index)
    index_only_true, index_only_run = difference_df(true_df, run_df)
    index_tp = [x for x in index_true if x not in index_only_true]
    index_fn = index_only_true
    index_fp = [x for x in index_run if x in index_only_run]
    nb_tp = len(index_tp)
    nb_fn = len(index_fn)
    nb_fp = len(index_fp)
    specificity = round(nb_tp / (nb_tp + nb_fn), 2)
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
    return nb_tp, nb_fn, nb_fp, specificity


def assess_indels(true_df, run_df, run_id):
    result = {}
    NB_TP, NB_FN, NB_FP = 0, 0, 0
    for index, row in true_df.iterrows():
        indel_info = [run_id, row['sample'], row['chr'], row['pos'], row['ref'], row['alt']]
        indel_str = '\t'.join([str(x) for x in indel_info])
        index_run = find_indel(row['sample'], row['chr'], row['pos'], row['ref'], row['alt'], run_df)
        if len(index_run) == 0:
            result[index] = None
            stat_str = f"na\tna\tna\tna\tna"
        else:
            indel_run = run_df.loc[index_run[0]]
            score = indel_run['score']
            vaf = indel_run['vaf']
            run_df_aux = run_df.loc[run_df['score']>=score]
            nb_tp, nb_fn, nb_fp, _ = compare_indels(true_df, run_df_aux, vaf, score)
            stat_str = '\t'.join([str(x) for x in [vaf, score, nb_tp, nb_fn, nb_fp]])
            result[index] = (score, nb_tp, nb_fn, nb_fp)
            NB_TP += nb_tp
            NB_FN += nb_fn
            NB_FP += nb_fp
        print(f"{indel_str}\t{stat_str}")
    return NB_TP, NB_FN, NB_FP

print('run_id\tsample\tchr\tpos\tref\talt\vaf\tscore\tTP\tFN\tFP')
for RUN_ID in RUN_ID_LIST:
    RUN_INDELS_FILE = os.path.join('results', RUN_ID, f"{RUN_ID}_indels.tsv")
    RUN_INDELS_DF_AUX = pd.read_csv(RUN_INDELS_FILE, sep='\t')
    RUN_INDELS_DF_AUX['sample'] = RUN_INDELS_DF_AUX.apply(
        lambda row: row['sample'].split('_S')[0], axis=1
    )
    RUN_INDELS_DF = RUN_INDELS_DF_AUX[~RUN_INDELS_DF_AUX['sample'].str.startswith('QMRS')].round(3)

    TRUE_INDELS_DF = ALL_INDELS_DF.loc[ALL_INDELS_DF['run_id']==RUN_ID]

    RUN_INDELS_DF_CTRL = RUN_INDELS_DF.loc[RUN_INDELS_DF['control']>CTRL_VAF_DIFF]

    _, _, _ = assess_indels(TRUE_INDELS_DF, RUN_INDELS_DF_CTRL, RUN_ID)
