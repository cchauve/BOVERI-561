"""
Functions to postprocess list of variants and features
"""
from collections import defaultdict
from itertools import groupby

from vcf_utils import (
    COMPLEXITY,
    SUPPORT,
    OVERLAP,
    CONTROL,
)
from bin.features_utils import (
    MAX_COV,
    SCORE,
    SOURCE,
)
from vcf_utils import (
    VAF,
    CHR_COL,
    POS_COL,
    REF_COL,
    ALT_COL,
    INFO_COL,
    ID_COL,
)


def get_control_samples_feature(parameters, variants_features, ctrl_variants_features):
    """
    :param: parameters (dict(str -> int/float)): filtering parameters
    :param: variants_features (list(Variant, VariantFeatures)): input list of
    variants and their features
    :param: ctrl_variants_features dict(str, (list(Variant, VariantFeatures)):
    dictionary of list of variants and their features in control samples, indexed
    by control samples ID

    :return:  list(Variant, VariantFeatures): input list where every variant
    occuring in a control sample with a VAF differing by at most
    parameters[FILTER_CTRL_VAF_DIFF] has been filtered out
    """
    # Control variants VAF indexed by variant string
    ctrl_variants_vaf = defaultdict(list)
    for ctrl_sample_id, ctrl_variants in ctrl_variants_features.items():
        # Recording variants for a control sample in a dict. indexed by variant string
        for (variant, features, _) in ctrl_variants_features[ctrl_sample_id]:
            ctrl_variants_vaf[variant.to_str(novaf=True)].append(variant.get_vaf())
    ctrl_variants_str = list(ctrl_variants_vaf.keys())
    # Filtering input variants
    out_list = []
    for (variant, features, score_dict) in variants_features:
        v_str, closest_vaf, closest_vaf_diff = variant.to_str(novaf=True), 0.0, 1.0
        if v_str in ctrl_variants_str:
            ctrl_vaf_list = ctrl_variants_vaf[v_str]
            in_vaf = variant.get_vaf()
            for ctrl_vaf in ctrl_vaf_list:
                vaf_diff = abs(in_vaf - ctrl_vaf)
                if vaf_diff < closest_vaf_diff:
                    closest_vaf = ctrl_vaf
                    closest_vaf_diff = vaf_diff
        score_dict[CONTROL] = closest_vaf
        out_list.append((variant, features, score_dict))
    return out_list

## Computing confidence score

def compute_complexity_seq(seq, kmin, kmax):
    def kmers(seq, k):
        i_max = len(seq) - k
        coords = [(i, i + k) for i in range(i_max + 1)]
        return set([seq[i:j] for (i, j) in coords])
    k_range = list(range(kmin, kmax + 1))
    nb_kmers = sum([len(kmers(seq, k)) for k in k_range])
    nb_kmers_max_1 = sum([len(seq) - k + 1 for k in k_range])
    nb_kmers_max_2 = sum([4 ** k for k in k_range])
    nb_kmers_max = min(nb_kmers_max_1, nb_kmers_max_2)
    return float(nb_kmers) / float(nb_kmers_max)


def compute_complexity_variant(ref,
                               alt,
                               pos,
                               amplicon_seq,
                               amplicon_start,
                               kmin,
                               kmax,
                               l):
    v_amp_start = pos - amplicon_start + (1 if ref[0] == alt[0] else 0)
    flanking_start = max(0, v_amp_start - l)
    v_amp_end = pos + len(ref) - 1 - amplicon_start
    flanking_end = min(len(amplicon_seq) - 1, v_amp_end + l)
    left_flanking_seq = amplicon_seq[flanking_start:v_amp_start]
    right_flanking_seq = amplicon_seq[v_amp_end+1:flanking_end + 1]
    ref_complexity = compute_complexity_seq(
        f"{left_flanking_seq}{ref}{right_flanking_seq}", kmin, kmax
    )
    alt_complexity = compute_complexity_seq(
        f"{left_flanking_seq}{alt}{right_flanking_seq}", kmin, kmax
    )
    return ref_complexity, alt_complexity


def compute_complexity_row(row, manifest, kmin, kmax, l):
    ref, alt, pos = row[REF_COL], row[ALT_COL], row[POS_COL]
    info_list = row[INFO_COL].split(';')
    for info in info_list:
        info_split = info.split('=')
        if info_split[0] == SOURCE:
            amplicon_id = info_split[1].split(',')[0]
            amplicon = manifest.get_amplicon(amplicon_id)
            amplicon_start = amplicon.get_start()
            amplicon_seq = amplicon.get_seq()
    return compute_complexity_variant(
        ref, alt, pos, amplicon_seq, amplicon_start, kmin, kmax, l
    )

COMP_REF_COL = 'comp_ref'
COMP_ALT_COL = 'comp_alt'

def compute_complexity_score(row, weight):
    score = weight * (1 - min(row[COMP_REF_COL], row[COMP_ALT_COL]))
    return score

def add_complexity_score(df, amplicons_data, kmin, kmax, l, weight=1.0):
    def comp_ref_alt(row, i):
        return compute_complexity_row(row, amplicons_data, kmin, kmax, l)[i]
    df[COMP_REF_COL] = df.apply(lambda row: comp_ref_alt(row, 0), axis=1)
    df[COMP_ALT_COL] = df.apply(lambda row: comp_ref_alt(row, 1), axis=1)
    df[COMPLEXITY] = df.apply(
        lambda row: compute_complexity_score(row, weight), axis=1
    )
    df.drop(columns=[COMP_REF_COL, COMP_ALT_COL], inplace=True)

def compute_support_score(row, support_min, weight):
    # Penalty due to size of largest supporting cluster
    info_list = row[INFO_COL].split(';')
    for info in info_list:
        info_split = info.split('=')
        if info_split[0] == MAX_COV:
            max_cluster = int(info_split[1])
    score = (support_min ** (1.0 / max_cluster) - 1) / (support_min - 1)
    return  weight * score

def add_support_score(df, support_min, weight=1.0):
    df[SUPPORT] = df.apply(
        lambda row: compute_support_score(row, support_min, weight), axis=1
    )


def compute_overlapping_calls(df):
    overlaps = {index: [] for index in df.index}
    opened_indels = []
    for current_index, row in df.iterrows():
        info_list = row[INFO_COL].split(';')
        for info in info_list:
            info_split = info.split('=')
            if info_split[0] == VAF:
                current_vaf = float(info_split[1])
        current_chr, current_pos = row[CHR_COL], row[POS_COL]
        current_ref = row[REF_COL]
        nb_opened_indels = len(opened_indels)
        current_end = current_pos + len(current_ref) - 1
        if nb_opened_indels == 0 or opened_indels[-1][3] != current_chr:
            opened_indels = [
                (current_index, current_end, current_vaf, current_chr)
            ]
        else:
            indels_to_remove = []
            for (index, end, vaf, chrom) in opened_indels:
                if end < current_pos:
                    indels_to_remove.append((index, end, vaf, chrom))
                else:
                    overlaps[current_index].append((index, vaf))
                    overlaps[index].append((current_index, current_vaf))
            for indel_to_remove in indels_to_remove:
                opened_indels.remove(indel_to_remove)
            current_indel = (
                current_index, current_end, current_vaf, current_chr
            )
            opened_indels.append(current_indel)
    return overlaps

def compute_overlapping_score(df, overlaps_dict, weight):
    for index, overlaps in overlaps_dict.items():
        if len(overlaps) == 0:
            df.at[index, OVERLAP] = 0
        else:
            info_list = df.at[index, INFO_COL].split(';')
            for info in info_list:
                info_split = info.split('=')
                if info_split[0] == VAF:
                    vaf = float(info_split[1])
            overlap_vaf = min(100.0, sum([x[1] for x in overlaps]))
            score =  weight * min(1.0, overlap_vaf / vaf)
            df.at[index, OVERLAP] = score

def add_overlap_score(df, weight=1.0):
    overlaps = compute_overlapping_calls(df)
    compute_overlapping_score(df, overlaps, weight)

def compute_confidence_score(row):
    return round(row[COMPLEXITY] + row[SUPPORT] + row[OVERLAP], 3)

def add_confidence_score(df):
    df[SCORE] = df.apply(
        lambda row: compute_confidence_score(row), axis=1
    )