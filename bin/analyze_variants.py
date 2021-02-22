#!/usr/bin/env python3

"""
Analyzing pipeline calls compared to expected calls to generate accuracy
statistics and identify False Positive (FP) and False Negative (FN) cases
"""
import argparse
import csv
import os
import pandas as pd
import yaml

# ------------------------------------------------------------------------------
# Global input

## Analysis: Features defining an indel
## Basic indel features
INDEL_FEATURES = ['sample', 'chr', 'pos', 'ref', 'alt']
# Expected indels features
INDEL_FEATURES_EXPVAF = INDEL_FEATURES + ['exp_vaf']
# Detected indels features
INDEL_FEATURES_VAF_1 = INDEL_FEATURES + ['vaf']
INDEL_FEATURES_VAF_2 = ['score', 'complexity', 'support', 'overlap', 'control']
INDEL_FEATURES_VAF = INDEL_FEATURES_VAF_1 + INDEL_FEATURES_VAF_2

# Paremeters keys for the YAML file
GRID_KEY = 'grid'
LLOD_KEY = 'llod'
SCORES_KEY = 'scores'
W_COMP_KEY = 'w_comp'
BLACKLIST_KEY = 'blacklist'
EXPECTED_INDELS_FILE = 'expected_indels_file'
OUT_SUFFIX_KEY = 'out_suffix'

# ------------------------------------------------------------------------------
# Auxiliary functions

def read_parameters(yaml_file_path):
    """
    Read a YAML parameters file and returns a dictionary of parameters values:
    - thresholds grid
    - llod parameters
    - blacklist
    """
    parameters = {BLACKLIST_KEY: [], GRID_KEY: None, OUT_SUFFIX_KEY: ''}
    with open(yaml_file_path) as c_file:
        parameters_dict = yaml.safe_load(c_file)
        for key, value in parameters_dict.items():
            if key == GRID_KEY:
                parameters[GRID_KEY] = [
                    (float(x.split('_')[0]), float(x.split('_')[1]))
                    for x in value.split()
                ]
            elif key == SCORES_KEY:
                scores_range = value.split()
                start, end = int(scores_range[0]), int(scores_range[1])
                step = float(scores_range[2])
                scores = [round(x * step, 2) for x in range(start, end)]
            elif key == W_COMP_KEY:
                w_comp_range = value.split()
                start, end = int(w_comp_range[0]), int(w_comp_range[1])
                step = float(w_comp_range[2])
                w_comp = [round(x * step, 2) for x in range(start, end)]
            elif key == BLACKLIST_KEY:
                parameters[BLACKLIST_KEY] = read_blacklist(value)
            elif key == LLOD_KEY:
                parameters[LLOD_KEY] = [float(x) for x in value.split()]
            elif key == OUT_SUFFIX_KEY:
                parameters[OUT_SUFFIX_KEY] = value
    if parameters[GRID_KEY] is None:
        parameters[GRID_KEY] = [(s, w) for s in scores for w in w_comp]
    return parameters


def read_blacklist(blacklist_csv_file):
    """
    Extracts from a CSV file a list of (chr, pos, ref, alt) describing indels
    :param: blacklist_csv_file (str): path to a CSV file describing indels in
    the format chr,pos,ref,alt
    :return: list(str, int, str, str): list of blacklisted indels in the format
    (chr, pos, ref, alt)
    """
    blacklist = []
    with open(blacklist_csv_file) as csvfile:
        blacklist_data = csv.reader(csvfile, delimiter=',')
        for row in blacklist_data:
            blacklist.append((row[0], int(row[1]), row[2], row[3]))
    return blacklist

def filter_blacklist(df, blacklist):
    """
    Returns a dataframe of variants from df by filtering all calls that match
    the features (chr, pos, ref, alt) in blacklist
    :param: df (DataFrame): dataframe of indels with a column 'sample'
    :param: blacklist (list(str, int, str, str)): list of blacklisted indels
    in the format (chr, pos, ref, alt)
    :return: DataFrame: input dataframe without any row where the entries
    'chr', 'pos', 'ref', 'alt' match a blacklisted indel from blacklist
    """
    index = []
    for (chrom, pos, ref, alt) in blacklist:
        index += list(
            df.loc[
                (df['chr']==chrom) &
                (df['pos']==pos) &
                (df['ref']==ref) &
                (df['alt']==alt)
            ].index
        )
    return  df.loc[~df.index.isin(index)]

def get_run_data(
    run_id, samples, expected_indels_df, blacklist=[]
):
    """
    Returns the dataframes of expected and detected indels for run run_id from
    data and results files
    :param: run_id (str): ID of the run
    :param: samples (list(str)): list of samples in the run for which there is
    at least one expected indel
    :param: expected_indels_df (DataFrame): dataframe of expected indels
    over a set of runs
    :param: blacklist (list(str, int, str, str)): list of blacklisted indels
    in the format (chr, pos, ref, alt)
    :return DataFRame, DataFrame: dataframe of expected indels for run run_ID,
    dataframe of detected indels for run run_id
    """
    # Reading the pipeline results
    run_indels_file = os.path.join('results', run_id, f"{run_id}_indels.tsv")
    detected_df_aux = pd.read_csv(run_indels_file, sep='\t')
    # Reformatting the sample column to match the true indels format
    detected_df_aux['sample'] = detected_df_aux.apply(
        lambda row: row['sample'].split('_S')[0], axis=1
    )
    run_indels_df = detected_df_aux.loc[
        detected_df_aux['sample'].isin(samples)
    ].round(3)
    # Excluding indels in the blasklist
    run_detected_indels_df = filter_blacklist(run_indels_df, blacklist)
    # Reading the expected (true) indels
    df_aux = expected_indels_df
    run_expected_indels_df = df_aux.loc[df_aux['run_id']==run_id]
    return run_expected_indels_df, run_detected_indels_df

def find_indel(row, df):
    """
    Computes the index of all indels defined by sample, chr, pos, ref, alt
    :param: row (dict()): dictionary with keys sample, chr, pos, ref, alt
    encoding an indel
    :param: df (DataFrame): dataframe of indels with columns matching the
    keys of row
    :return: list(int): list of index in df where the entries match the values
    of row
    """
    return list(df.loc[(df['sample']==row['sample']) &
                       (df['chr']==row['chr']) &
                       (df['pos']==row['pos']) &
                       (df['ref']==row['ref']) &
                       (df['alt']==row['alt'])
                       ].index)

def diff_df(df1, df2, columns=INDEL_FEATURES):
    """
    Computes df1-df2 and df2-df1 on a set of common columns, by default the
    basic features of an indel.
    :param: df1, df2 (DataFrame): dataframes with columns including columns
    :param: columns (list(str)): columns on which df1 and df2 are compared
    :return: DataFrame, DataFrame: df1-df2, df2-df1
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

def compare_indels(expected_df, detected_df, score):
    """
    Computes the index of true positive, false positive, false negative
    :param: expected_df (DataFrame): dataframe of all expected indels
    :param: detected_df (DataFrame): dataframe of all dected indels
    :param: score (float): max penalty to call a detected indel
    :return: list(int), list(int), list(int), list(int):
    - index of TP calls in detected_df
    - index of FN indels in expected_df
    - index of FP calls in detected_df
    - index of TN calls in detected _df
    """
    called_df = detected_df.loc[detected_df['score']<=score]
    uncalled_df = detected_df.loc[detected_df['score']>score]
    expected_index = list(expected_df.index)
    called_index = list(called_df.index)
    expected_only_index, called_only_index = diff_df(expected_df, called_df)
    _, uncalled_only_index = diff_df(expected_df, uncalled_df)
    index_tp = [x for x in called_index if x not in called_only_index]
    index_fn = expected_only_index
    index_fp = [x for x in called_index if x in called_only_index]
    index_tn = uncalled_only_index
    return index_tp, index_fn, index_fp, index_tn

# ------------------------------------------------------------------------------
STAT_PRECISION = 3

def compute_ratio(x, y, precision=STAT_PRECISION):
    """
    Computes the ratio x/y with precision precision, but if y=0 in whih case
    returns 0.0
    """
    return (0.0 if y == 0 else  round(x / y, precision))

def process_indels(parameters):
    """
    Computing sensitivity, specificity, precision, recall, F1, FDR, accuracy
    for a grid of thresholds and parameters
    :param: parameters (dict): parameters dictionary containing the LLOD,
    grid of thresholds, expected indels input file and blacklist list
    :return: None
    Write two files. If exp_indels_file is NAME.tsv:
    - NAME_out.tsv: statistics
    - NAME_errors.tsv: FP and FN indels
    """
    exp_indels_file = parameters[EXPECTED_INDELS_FILE]
    llod_min = parameters[LLOD_KEY][0]
    llod_lowering_factor = parameters[LLOD_KEY][1]
    grid = parameters[GRID_KEY]
    blacklist = parameters[BLACKLIST_KEY]
    out_suffix = parameters[OUT_SUFFIX_KEY]
    # Reading expected indels
    ALL_EXPECTED_INDELS_DF = pd.read_csv(exp_indels_file, sep='\t')
    ALL_EXPECTED_INDELS_DF.rename(columns={'chromosome': 'chr'}, inplace=True)
    # LLOD
    LLOD = ALL_EXPECTED_INDELS_DF['exp_vaf'].min()
    if pd.isnull(LLOD):
        LLOD = 0.0
    LLOD_LOWERED = max(llod_min, llod_lowering_factor * LLOD)
    # List of runs
    RUNS_ID_LIST = list(ALL_EXPECTED_INDELS_DF['run_id'].unique())
    # Output file
    _, dataset_name = os.path.split(exp_indels_file)
    out_file_name = dataset_name.replace('.tsv', f"{out_suffix}_out.tsv")
    out_file = open(os.path.join('results', out_file_name), 'w')
    header_1 = ['LLOD', 'score', 'w_comp']
    header_2 = ['TP', 'FP', 'TN', 'FN_u', 'FN_d']
    header_3 = ['sens.', 'spec.', 'acc.', 'prec.', 'recall', 'F1', 'FDR']
    out_file.write('\t'.join(header_1 + header_2 + header_3))
    out_file_errors_name = out_file_name.replace('_out.tsv', '_errors.txt')
    out_file_errors = open(os.path.join('results', out_file_errors_name), 'w')
    header_4 = ['run', 'sample', 'chr', 'pos', 'ref', 'alt', 'vaf']
    header_5 = ['score', 'comp.', 'supp.', 'overlap', 'ctrl']
    out_file_errors.write('\t'.join(header_1 + header_4 + header_5))
    # Going through pipeline results
    for (score, w) in grid:
        w1 = round(1.0 - w, 2)
        print(LLOD, LLOD_LOWERED, score, w)
        # Numbers of TP, FP and TN
        nb_tp, nb_fp, nb_tn = 0, 0, 0
        # Number of FN, FN detected but not called, FN undetected
        nb_fn, nb_fn_detected, nb_fn_undetected = 0, 0, 0
        # Looping on all runs for which at least one indel s expected
        for run_id in RUNS_ID_LIST:
            # Samples in which at least one indel is expected
            samples = ALL_EXPECTED_INDELS_DF.loc[
                ALL_EXPECTED_INDELS_DF['run_id']==run_id
            ]['sample'].unique()
            # Extracting all expected indels and all detected indels for the
            # selected samples of the run, but the black-listed ones
            expected_indels_df, run_indels_df = get_run_data(
                run_id, samples, ALL_EXPECTED_INDELS_DF, blacklist=blacklist)
            # Reducing the weight of complexity penalty by factor w
            run_indels_df['score'] = run_indels_df.apply(
                lambda row: row['score'] - (w1 * row['complexity']),
                axis=1
            )
            # Filtering out detected indels with a VAF below LLOD_LOWERED
            detected_indels_df = run_indels_df.loc[
                (run_indels_df['vaf']>=(LLOD_LOWERED))
            ]
            # Index in dataframes of TP (detected_indels_df),
            # FN (expected_indels_df), FP (detected_indels_df) and
            # TN (detected_indels_df)
            index_tp, index_fn, index_fp, index_tn = compare_indels(
                expected_indels_df, detected_indels_df, score
            )
            # Updating statistics on numbers of TP, FP, FN, TN
            nb_tp += len(index_tp)
            nb_fn += len(index_fn)
            nb_fp += len(index_fp)
            nb_tn += len(index_tn)
            # Dataframe of FP indels and FN indels
            fp_df = detected_indels_df.loc[index_fp]
            fn_df = expected_indels_df.loc[index_fn]
            # Output of FP indels
            for _, row in fp_df.iterrows():
                indel_info = (
                    [round(LLOD, 2), round(score, 2), round(w, 2)] +
                    ['FP', run_id] +
                    [row[x] for x in INDEL_FEATURES_VAF_1] +
                    [round(row[x], 3) for x in INDEL_FEATURES_VAF_2]
                )
                indel_str = '\t'.join([str(x) for x in indel_info])
                out_file_errors.write('\n' + indel_str)
            # Output of FN indels, separated in two cases, FN detected but not
            # called and FN undetected
            for _, row in fn_df.iterrows():
                index_indel_in_detected = find_indel(row, run_indels_df)
                if len(index_indel_in_detected) == 1:
                    # Detected but uncalled FN
                    index_found = index_indel_in_detected[0]
                    nb_fn_detected += 1
                    fn_status = 'FN_d'
                    features = [
                        round(run_indels_df.at[index_found, x], 3)
                        for x in INDEL_FEATURES_VAF_2
                    ]
                else:
                    # Undetected FN
                    nb_fn_undetected += 1
                    fn_status = 'FN_u'
                    features = ['nan' for x in INDEL_FEATURES_VAF_2]
                indel_info = [round(LLOD, 2), round(score, 2), round(w, 2)]
                if fn_status == 'FN_u':
                    indel_info += [row[x] for x in INDEL_FEATURES_EXPVAF]
                else:
                    indel_info += [
                        run_indels_df.at[index_found, x]
                        for x in INDEL_FEATURES_VAF_1
                    ]
                indel_info += features
                indel_str = '\t'.join([str(x) for x in indel_info])
                out_file_errors.write('\n' + indel_str)
        # Checking all FNs are either detected or undected
        assert(nb_fn == nb_fn_undetected + nb_fn_detected)
        # Computing statistics over all samples
        sensitivity = compute_ratio(nb_tp, nb_tp + nb_fn)
        specificity = compute_ratio(nb_tn, nb_tn + nb_fp)
        accuracy = compute_ratio(nb_tp + nb_tn, nb_tp + nb_tn + nb_fp + nb_fn)
        precision = compute_ratio(nb_tp, nb_tp + nb_fp)
        recall = sensitivity
        F1 = compute_ratio(2.0 * precision * recall, precision + recall)
        FDR = round(1.0 - precision, STAT_PRECISION)
        # Output of statistics
        param = [LLOD, score, round(w, 2)]
        stats = [nb_tp, nb_fp, nb_tn, nb_fn_undetected, nb_fn_detected]
        stats += [sensitivity, specificity, accuracy, precision, recall]
        stats += [F1, FDR]
        results = [str(x) for x in param + stats]
        out_file.write('\n' + '\t'.join(results))
    out_file.close()
    out_file_errors.close()

# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # Analysis parameters
    # Input file
    ARGS_EXPECTED_FILE = ['exp_indels_file', None, 'Expected indels file']
    # Black list file
    ARGS_PARAMETERS_FILE = ['parameters_file', None, 'Parameters YAML file']
    parser = argparse.ArgumentParser(description='Indels testing: explore scores')
    parser.add_argument(ARGS_EXPECTED_FILE[0],
                        type=str,
                        help=ARGS_EXPECTED_FILE[2])
    parser.add_argument(ARGS_PARAMETERS_FILE[0],
                        type=str,
                        help=ARGS_PARAMETERS_FILE[2])
    args = parser.parse_args()

    parameters = read_parameters(args.parameters_file)
    parameters[EXPECTED_INDELS_FILE] = args.exp_indels_file
    process_indels(parameters)
