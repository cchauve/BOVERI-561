#!/usr/bin/env python3

"""
Analyzing pipeline calls compared to expected calls
"""
import argparse
import csv
import os
import pandas as pd

# ------------------------------------------------------------------------------
# Global input
# Input: Expected indels dataframe
ALL_EXPECTED_INDELS_DF = pd.read_csv(
    'data/v4MiSeq_commercial_samples_expected_indels_ExpectedVAFs.tsv',
    sep='\t'
)
ALL_EXPECTED_INDELS_DF.rename(columns={'chromosome': 'chr'}, inplace=True)
# Input: amplicons manifest
MANIFEST_DF = pd.read_csv(
    'data/CG001v4.0_Amplicon_Manifest_Panel4.0.3_20181101.tsv', sep='\t'
)
MANIFEST_DICT = {}
for _, row in MANIFEST_DF.iterrows():
    MANIFEST_DICT[row['Amplicon_ID']] = row['Amplicon']
# Analysis: Features defining an indel
INDEL_FEATURES_EXPVAF = ['sample', 'chr', 'pos', 'ref', 'alt', 'exp_vaf']
INDEL_FEATURES_NOVAF = ['sample', 'chr', 'pos', 'ref', 'alt']

# ------------------------------------------------------------------------------
# Auxiliary functions
def read_run_list(run_csv_file):
    """
    Extracts from a csv file a list of run IDs
    """
    run_id_list = []
    with open(run_csv_file) as csvfile:
        runs_data = csv.reader(csvfile, delimiter=',')
        for row in runs_data:
            run_id_list.append(row[1])
    return run_id_list

def find_indel(row, df):
    """
    Computes the index of all indels defined by sample, chr, pos, ref, alt
    """
    return list(df.loc[(df['sample']==row['sample']) &
                       (df['chr']==row['chr']) &
                       (df['pos']==row['pos']) &
                       (df['ref']==row['ref']) &
                       (df['alt']==row['alt'])
                       ].index)

def find_amplicon_coordinates(row, manifest_df):
    """
    Returns a dictionary indexed by amplicon ID of the coordinate in the
    amplicon of row['pos']
    """
    result = {}
    for amplicon_id, amplicon in manifest_df.iterrows():
        test_chr = amplicon['Chr'] == row['chr']
        test_pos = amplicon['End'] >= row['pos'] >= amplicon['Start']
        if test_chr and test_pos:
            result[amplicon['Amplicon_ID']] = row['pos'] - amplicon['Start']
    return result

def get_cluster_consensus_sequences(clusters_file):
    """
    Returns a dictionary cluster_id -> cluster sequence for a given clusters
    file
    """
    clusters_sequences = {}
    clusters = open(clusters_file, 'r').readlines()
    cluster_to_read = False
    for line in clusters:
        if line[0] == '>':
            cluster_id = line.split()[0][1:]
            cluster_to_read = True
        elif cluster_to_read:
            clusters_sequences[cluster_id] = line.rstrip()
            cluster_to_read = False
    return clusters_sequences

def check_indel(row, manifest_dict, amp_pos, run_id, out_file):
    """
    Returns True if the variant described in row matches with some cluster
    consensus sequence.
    """
    sample = f"{row['sample']}_S{row['sample'].split('-')[-1]}"
    indel_info = [run_id] + [row[x] for x in INDEL_FEATURES_EXPVAF]
    indel_str = '\t'.join([str(x) for x in indel_info])
    for amplicon_id, pos in amp_pos.items():
        clusters_file_name = f"{sample}_{amplicon_id}_clusters.txt"
        clusters_file_path = os.path.join('results', run_id, clusters_file_name)
        clusters_sequences = get_cluster_consensus_sequences(clusters_file_path)
        amp_seq = manifest_dict[amplicon_id]
        if len(row['ref']) == 1: # insertion
            v_seq = f"{amp_seq[:pos]}{row['alt']}{amp_seq[pos + 1:]}"
        else:
            v_seq = f"{amp_seq[:pos + 1]}{amp_seq[pos + len(row['ref']):]}"
        for cluster_id, c_seq in clusters_sequences.items():
            len_diff = abs(len(c_seq) - len(v_seq))
            test_1 = len_diff > 2 and len(c_seq) == len(v_seq)
            test_2 = len_diff <= 2 and c_seq == v_seq
            if test_1 or test_2:
                out_file.write(f"\n{amplicon_id}_{cluster_id}\t{indel_str}")
                return True
    out_file.write(f"\nnan\t{indel_str}")
    return False

def exclude_qmrs(df):
    """
    Returns a variants dataframe excluding all calls from QMRS
    """
    return df.loc[~df['sample'].str.startswith('QMRS')]

def exclude_control_calls(df, ctrl_vaf_diff):
    """
    Returns a variants dataframe with all calls present in control with a close
    VAF filtered out
    """
    return df.loc[df['control']>ctrl_vaf_diff]

def get_run_data(run_id, expected_indels_df, ctrl_vaf_diff):
    """
    Returns the dataframe of expected and called indels for run run_id
    """
    # Reading the pipeline results
    run_indels_file = os.path.join('results', run_id, f"{run_id}_indels.tsv")
    called_df_aux = pd.read_csv(run_indels_file, sep='\t')
    # Reformatting the sample column to match the true indels format
    called_df_aux['sample'] = called_df_aux.apply(
        lambda row: row['sample'].split('_S')[0], axis=1
    )
    # Excluding the QMRS sample and rounding all numerical values
    run_indels_df_1 = exclude_qmrs(called_df_aux).round(3)
    # Excluding indels seen in the control samples with a close VAF
    run_called_indels_df = exclude_control_calls(run_indels_df_1, ctrl_vaf_diff)
    # Reading the expected (true) indels
    df_aux = expected_indels_df
    run_expected_indels_df = df_aux.loc[df_aux['run_id']==run_id]
    return run_expected_indels_df, run_called_indels_df

def diff_df(df1, df2, columns=['sample', 'chr', 'pos', 'ref', 'alt']):
    """
    Computes df1-df2 and df2-df1 on a set of common columns
    """
    df1_aux = df1[columns]
    df2_aux = df2[columns]
    df2_m_df1 = pd.concat(
        [df2_aux, df1_aux, df1_aux]
    ).drop_duplicates(keep=False)
    df1_m_df2 = pd.concat(
        [df1_aux, df2_aux, df2_aux]
    ).drop_duplicates(keep=False)
    return list(df1_m_df2.index), list(df2_m_df1.index)

def compare_indels(expected_df, called_df_in, score, verbose=False):
    """
    Computes the index of true positive, false positive, false negative
    If verbose, prints stats about TP, FN, FP and prints all FN and FP indels.
    Return index of TP and FP indels in called_df_in and FN indels in
    expected_df
    """
    called_df = called_df_in.loc[called_df_in['score']<=score]
    expected_index = list(expected_df.index)
    called_index = list(called_df.index)
    expected_only_index, called_only_index = diff_df(expected_df, called_df)
    index_tp = [x for x in called_index if x not in called_only_index]
    index_fn = expected_only_index
    index_fp = [x for x in called_index if x in called_only_index]
    nb_tp, nb_fn, nb_fp = len(index_tp), len(index_fn), len(index_fp)
    if verbose:
        stats = f"TP:{nb_tp}\tFN:{nb_fn}\tFP:{nb_fp}"
        print(f"\tscore_min:{score}\t{stats}")
        for _, row in expected_df[expected_df.index.isin(index_fn)].iterrows():
            out_str = (
                f"\tFN\t{row['sample']}\t{row['chr']}\t{row['pos']}"
                f"\t{row['ref']}\t{row['alt']}\t{row['exp_vaf']}"
            )
            print(out_str)
        for _, row in called_df[called_df.index.isin(index_fp)].iterrows():
            out_str = (
                f"\tFP\t{row['sample']}\t{row['chr']}\t{row['pos']}"
                f"\t{row['ref']}\t{row['alt']}"
                f"\t{row['vaf']}\t{row['score']}\t{row['complexity']}"
                f"\t{row['support']}\t{row['overlap']}\t{row['control']}"
                f"\t{row['repeats']}"
            )
            print(out_str)
    return index_tp, index_fn, index_fp

def print_stats(
    called_indels_df, expected_indels_df, index_tp, index_fn, index_fp, prefix
):
    """
    Creates a string with all statistics for a set of TP, FN and FP
    """
    fp_df = called_indels_df.loc[index_fp]
    tp_df = called_indels_df.loc[index_tp]
    fn_df = expected_indels_df.loc[index_fn]
    nb_tp, nb_fn, nb_fp = len(index_tp), len(index_fn), len(index_fp)
    # FP statistics
    fp_mean_vaf = round(fp_df['vaf'].mean(), 2)
    fp_mean_score = round(fp_df['score'].mean(), 2)
    nb_fp_1 = len(fp_df.loc[fp_df['vaf']<1.0].index)
    nb_fp_05 = len(fp_df.loc[fp_df['vaf']<0.5].index)
    # TP statistics
    tp_mean_vaf = round(tp_df['vaf'].mean(), 2)
    tp_mean_score = round(tp_df['score'].mean(), 2)
    nb_tp_1 = len(tp_df.loc[tp_df['vaf']<1.0].index)
    nb_tp_05 = len(tp_df.loc[tp_df['vaf']<0.5].index)
    # FN statistics
    fn_mean_vaf = round(fn_df['exp_vaf'].mean(), 2)
    nb_fn_100 = len(fn_df.loc[fn_df['exp_vaf']>=1.0].index)
    nb_fn_1 = len(fn_df.loc[fn_df['exp_vaf']<1.0].index)
    nb_fn_05 = len(fn_df.loc[fn_df['exp_vaf']<0.5].index)
    # Printing
    stat_fields = prefix
    stat_fields += [nb_tp, tp_mean_vaf, tp_mean_score, nb_tp_1, nb_tp_05]
    stat_fields += [nb_fn, fn_mean_vaf, nb_fn_100, nb_fn_1, nb_fn_05]
    stat_fields += [nb_fp, fp_mean_vaf, fp_mean_score, nb_fp_1, nb_fp_05]
    stat_str = '\t'.join([str(x) for x in stat_fields])
    return stat_str

# ------------------------------------------------------------------------------
HEADER_STATS = (
    f"TP:nb\tTP:mean_vaf\tTP:mean_score\tTP:nb:vaf<1\tTP:nb:vaf<0.5\t"
    f"FN:nb\tFN:mean_vaf\tFN:nb:vaf>=1\tFN:nb:vaf<1\tFN:nb:vaf<0.5\t"
    f"FP:nb\tFP:mean_vaf\tFP:mean_score\tFP:nb:vaf<1\tFP:nb:vaf<0.5"
)

def run_cmd_a1(dataset_name, run_id_list, ctrl_vaf_diff):
    """
    Analysis 1: expected indels by minimum score to call them: for each
    expected indel, print statistics about indels called with a score at
    least as good.
    """
    out_file_name = dataset_name.replace('.csv', '_out_1.tsv')
    out_file = open(os.path.join('results', out_file_name), 'w')
    header_1 = 'run_id\tsample\tchr\tpos\tref\talt\t'
    header_2 = 'vaf\t'
    header_2_len = header_2.count('\t')
    header_3 = 'exp_vaf'
    header_4 = f"\tscore\tcomplexity\tsupport\toverlap\t{HEADER_STATS}"
    header_4_len = header_4.count('\t')
    out_file.write(header_1 + header_2 + header_3 + header_4)
    for run_id in run_id_list:
        expected_indels_df, called_indels_df = get_run_data(
            run_id, ALL_EXPECTED_INDELS_DF, ctrl_vaf_diff
        )
        # loop on all expected indels
        for index, row in expected_indels_df.iterrows():
            exp_vaf = row['exp_vaf']
            indel_info = [run_id] + [row[x] for x in INDEL_FEATURES_NOVAF]
            indel_str = '\t'.join([str(x) for x in indel_info])
            called_index = find_indel(row, called_indels_df)
            if len(called_index) == 0:
                stat_str = '\t'.join(['nan' for i in range(header_2_len)])
                stat_str += f"\t{str(exp_vaf)}\t"
                stat_str += '\t'.join(['nan' for i in range(header_4_len)])
            else:
                indel_called = called_indels_df.loc[called_index[0]]
                score = indel_called['score']
                vaf = indel_called['vaf']
                complexity = indel_called['complexity']
                support = indel_called['support']
                overlap = indel_called['overlap']
                index_tp, index_fn, index_fp = compare_indels(
                    expected_indels_df, called_indels_df, score
                )
                stat_str = print_stats(
                    called_indels_df, expected_indels_df,
                    index_tp, index_fn, index_fp,
                    [vaf, exp_vaf, score, complexity, support, overlap]
                )
            out_file.write(f"\n{indel_str}\t{stat_str}")
    out_file.close()

def run_cmd_a2(dataset_name, run_id_list, ctrl_vaf_diff):
    """
    Analysis 2: statistics by minimum score
    For a range of scores we print statistics per run about indels called with
    this score
    """
    out_file_name = dataset_name.replace('.csv', '_out_2.tsv')
    out_file = open(os.path.join('results', out_file_name), 'w')
    out_file.write(f"run_id\tmax_score\tf{HEADER_STATS}")
    for run_id in run_id_list:
        expected_indels_df, called_indels_df = get_run_data(
            run_id, ALL_EXPECTED_INDELS_DF, ctrl_vaf_diff
        )
        for score in [0.25, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]:
            index_tp, index_fn, index_fp = compare_indels(
                expected_indels_df, called_indels_df, score
            )
            stat_str = print_stats(
                called_indels_df, expected_indels_df,
                index_tp, index_fn, index_fp,
                [run_id, score]
            )
            out_file.write(f"\n{stat_str}")
    out_file.close()

def run_cmd_fn(dataset_name, run_id_list, ctrl_vaf_diff):
    """
    Checking if FN are present in cluster sequences
    """
    out_file_name = dataset_name.replace('.csv', '_out_3.tsv')
    out_file = open(os.path.join('results', out_file_name), 'w')
    header = '\t'.join(INDEL_FEATURES_EXPVAF)
    out_file.write(f"found\t{header}")
    for run_id in run_id_list:
        expected_indels_df, called_indels_df = get_run_data(
            run_id, ALL_EXPECTED_INDELS_DF, ctrl_vaf_diff
        )
        for index, row in expected_indels_df.iterrows():
            indel_info = [run_id] + [row[x] for x in INDEL_FEATURES_EXPVAF]
            indel_str = '\t'.join([str(x) for x in indel_info])
            called_index = find_indel(row, called_indels_df)
            #if len(called_index) > 0 and not indel_present:
            #    print(f"TP missing: {indel_str}")
            if len(called_index) == 0:
                amp_pos = find_amplicon_coordinates(row, MANIFEST_DF)
                indel_present = check_indel(
                    row,  MANIFEST_DICT, amp_pos, run_id, out_file
                )
    out_file.close()

# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # Analysis parameters
    # Input file
    ARGS_RUNS_FILE = ['runs_csv_file', None, 'Runs CSV file']
    # Command to execute
    ARGS_CMD = ['cmd', None, 'Command to execute']
    # Max difference of VAF with control to filter out an indel call
    ARGS_CTRL_VAF_DIFF = ['-c', '--ctrl_vaf_diff', 'control VAF difference']
    parser = argparse.ArgumentParser(description='Indels testing: analyze')
    parser.add_argument(ARGS_RUNS_FILE[0], type=str, help=ARGS_RUNS_FILE[2])
    parser.add_argument(ARGS_CMD[0], type=str, help=ARGS_CMD[2])
    parser.add_argument(ARGS_CTRL_VAF_DIFF[0],
                        ARGS_CTRL_VAF_DIFF[1],
                        type=float,
                        default=5.0,
                        help=ARGS_CTRL_VAF_DIFF[2])
    args = parser.parse_args()

    # Input: Run IDs list
    run_id_list = read_run_list(args.runs_csv_file)
    # Dataset name
    _, dataset_name = os.path.split(args.runs_csv_file)
    if args.cmd == 'a1':
        run_cmd_a1(dataset_name, run_id_list, args.ctrl_vaf_diff)
    elif args.cmd == 'a2':
        run_cmd_a2(dataset_name, run_id_list, args.ctrl_vaf_diff)
    elif args.cmd == 'fn':
        run_cmd_fn(dataset_name, run_id_list, args.ctrl_vaf_diff)
