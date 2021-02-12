"""
Functions to manipulate VCF files
"""
import copy
from datetime import date
import io
import os
import pandas as pd
import subprocess
from textwrap import dedent

import vcf

from bin.features_utils import (
    HP_LEFT_BASE,
    HP_LEFT_LEN,
    HP_RIGHT_BASE,
    HP_RIGHT_LEN,
    MAX_COV,
    SCORE,
    SOURCE,
    SOURCE_COV,
    TOTAL_COV,
    V_RU,
    V_RU_CNB,
    V_RU_LEFT_CNB,
    V_RU_RIGHT_CNB,
    WT_RU,
    WT_RU_CNB,
    WT_RU_LEFT_CNB,
    WT_RU_RIGHT_CNB,
    VariantFeatures,
)
from bin.parameters_utils import (
    SNPEFF_LOG,
    SNPEFF_PATH,
)
from bin.pipeline_file_names import vcf_sample_file
from bin.variants_utils import (
    DEL,
    INS,
    VAF_PRECISION,
    Variant,
)


# Sorting criterion for list of triples (variant_str, variant, features)
def sort_by_position(x):
    """
    :param: x (Variant, VariantFeatures)
    :return: (str, int): source, variant start position
    """
    return (x[0].get_source(), x[0].get_start())


def sort_by_vaf(x):
    """
    :param: x (Variant, VariantFeatures)
    :return: float: variant VAF
    """
    return x[0].get_vaf()


VCF_SORT = {}  # Dictionary of sorting function and boolean indicating the reverse value
VCF_SORT_POS, VCF_SORT_VAF, = 'pos', 'vaf'
VCF_SORT[VCF_SORT_POS] = (sort_by_position, False)
VCF_SORT[VCF_SORT_VAF] = (sort_by_vaf, True)


def vcf_sorting_criteria():
    return list(VCF_SORT.keys())


# VCF file writing

COMPLEXITY = 'COMPLEXITY_SCORE'
SUPPORT = 'SUPPORT_SCORE'
OVERLAP = 'OVERLAP_SCORE'
CONTROL = 'CONTROL'
VAF = 'VAF'
V_TYPE = 'TYPE'
ANNOTATION = 'ANN'
SAMPLE = 'SAMPLE'
RUN_ID = 'RUN_ID'
RUN_NAME = 'RUN_NAME'

CHR_COL = 'CHROM'
POS_COL = 'POS'
REF_COL = 'REF'
ALT_COL = 'ALT'
INFO_COL = 'INFO'
ID_COL = 'ID'

VCF_DESC = {}
VCF_DESC[V_TYPE] = '\"Variant type (DEL, INS, DELINS, MNV, SNP)\"'
VCF_DESC[VAF] = '\"Variant Allele Frequency\"'
VCF_DESC[SCORE] = '\"Confidence score\"'
VCF_DESC[COMPLEXITY] = '\"Complexity score\"'
VCF_DESC[SUPPORT] = '\"Support score\"'
VCF_DESC[OVERLAP] = '\"Overlap score\"'
VCF_DESC[CONTROL] = '\"Closest VAF in a control sample (0.0 means not present in control samples)\"'
VCF_DESC[WT_RU] = '\"Length of minimal repeat unit (RU) in reference\"'
VCF_DESC[WT_RU_CNB] = '\"Copy number of ref RU in reference\"'
VCF_DESC[WT_RU_LEFT_CNB] = '\"Copy number of ref RU left of reference\"'
VCF_DESC[WT_RU_RIGHT_CNB] = '\"Copy number of ref RU right of reference\"'
VCF_DESC[V_RU] = '\"Length of minimal RU in alternate\"'
VCF_DESC[V_RU_CNB] = '\"Copy number of alt RU in alternate\"'
VCF_DESC[V_RU_LEFT_CNB] = '\"Copy number of alt RU left of reference\"'
VCF_DESC[V_RU_RIGHT_CNB] = '\"Copy number of alt RU right of reference\"'
VCF_DESC[HP_LEFT_LEN] = '\"Length of homopolymer left of reference\"'
VCF_DESC[HP_LEFT_BASE] = '\"Base of homopolymer left of reference\"'
VCF_DESC[HP_RIGHT_LEN] = '\"Length of homopolymer right of reference\"'
VCF_DESC[HP_RIGHT_BASE] = '\"Base of homopolymer right of reference\"'
VCF_DESC[SOURCE_COV] = '\"Source coverage\"'
VCF_DESC[TOTAL_COV] = '\"Total variant coverage\"'
VCF_DESC[MAX_COV] = '\"Largest cluster coverage of variant\"'
VCF_DESC[SOURCE] = '\"Amplicons where the variant is observed\"'
VCF_DESC[SAMPLE] = '\"Sample where the variant is observed\"'
VCF_DESC[RUN_ID] = '\"ID of run where the variant is observed\"'
VCF_DESC[RUN_NAME] = '\"Name of run where the variant is observed\"'

def vcf_header():
    """
    String of the metadata of a VCF file
    """
    vcf_str = dedent(
        f"""
        ##fileformat=VCFv4.1
        ##fileDate={str(date.today().strftime("%Y%m%d"))}
        ##source=Canexia_indels
        ##INFO=<ID={V_TYPE},Number=1,Type=String,Description={VCF_DESC[V_TYPE]}>
        ##INFO=<ID={VAF},Number=1,Type=Float,Description={VCF_DESC[VAF]}>
        ##INFO=<ID={SCORE},Number=1,Type=Float,Description={VCF_DESC[SCORE]}>
        ##INFO=<ID={COMPLEXITY},Number=1,Type=Float,Description={VCF_DESC[COMPLEXITY]}>
        ##INFO=<ID={SUPPORT},Number=1,Type=Float,Description={VCF_DESC[SUPPORT]}>
        ##INFO=<ID={OVERLAP},Number=1,Type=Float,Description={VCF_DESC[OVERLAP]}>
        ##INFO=<ID={CONTROL},Number=1,Type=Float,Description={VCF_DESC[CONTROL]}>
        ##INFO=<ID={WT_RU},Number=1,Type=Integer,Description={VCF_DESC[WT_RU]}>
        ##INFO=<ID={WT_RU_CNB},Number=1,Type=Integer,Description={VCF_DESC[WT_RU_CNB]}>
        ##INFO=<ID={WT_RU_LEFT_CNB},Number=1,Type=Integer,Description={VCF_DESC[WT_RU_LEFT_CNB]}>
        ##INFO=<ID={WT_RU_RIGHT_CNB},Number=1,Type=Integer,Description={VCF_DESC[WT_RU_RIGHT_CNB]}>
        ##INFO=<ID={V_RU},Number=1,Type=Integer,Description={VCF_DESC[V_RU]}>
        ##INFO=<ID={V_RU_CNB},Number=1,Type=Integer,Description={VCF_DESC[V_RU_CNB]}>
        ##INFO=<ID={V_RU_LEFT_CNB},Number=1,Type=Integer,Description={VCF_DESC[V_RU_LEFT_CNB]}>
        ##INFO=<ID={V_RU_RIGHT_CNB},Number=1,Type=Integer,Description={VCF_DESC[V_RU_RIGHT_CNB]}>
        ##INFO=<ID={HP_LEFT_LEN},Number=1,Type=Integer,Description={VCF_DESC[HP_LEFT_LEN]}>
        ##INFO=<ID={HP_LEFT_BASE},Number=1,Type=String,Description={VCF_DESC[HP_LEFT_BASE]}>
        ##INFO=<ID={HP_RIGHT_LEN},Number=1,Type=Integer,Description={VCF_DESC[HP_RIGHT_LEN]}>
        ##INFO=<ID={HP_RIGHT_BASE},Number=1,Type=String,Description={VCF_DESC[HP_RIGHT_BASE]}>
        ##INFO=<ID={SOURCE_COV},Number=1,Type=Integer,Description={VCF_DESC[SOURCE_COV]}>
        ##INFO=<ID={TOTAL_COV},Number=1,Type=Integer,Description={VCF_DESC[TOTAL_COV]}>
        ##INFO=<ID={MAX_COV},Number=1,Type=Integer,Description={VCF_DESC[MAX_COV]}>
        ##INFO=<ID={SOURCE},Number=.,Type=String,Description={VCF_DESC[SOURCE]}>
        ##INFO=<ID={SAMPLE},Number=.,Type=String,Description={VCF_DESC[SAMPLE]}>
        ##INFO=<ID={RUN_ID},Number=.,Type=String,Description={VCF_DESC[RUN_ID]}>
        ##INFO=<ID={RUN_NAME},Number=.,Type=String,Description={VCF_DESC[RUN_NAME]}>
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""
    )[1:]
    return vcf_str

def tsv_header():
    return 'sample\trun_id\trun\tchr\tpos\tref\talt\tVAF\ttype\tscore\tcomplexity\tsupport\toverlap\tcontrol\tcov\talt_cov\tmax_cov\tannotation'


def vcf_variant(variant, features, extra_features, precision=VAF_PRECISION):
    """
    String of a VCF entry for one variant with fixed precision for VAF

    :param: variant (Variant): variant
    :param: features (VariantFeatures): features of the variant
    :param: precision (int): precision of floating numbers

    :return: str: VCF entry
    """
    hp_left = features.get_hp_left()
    hp_right = features.get_hp_right()

    vcf_str = (
        f"{variant.get_source()}"
        f"\t{variant.get_start()}"
        f"\t."
        f"'\t{variant.get_wt_seq()}"
        f"\t{variant.get_variant_seq()}"
        f"\t.\tFAIL"
        f"\tTYPE={variant.get_type()}"
        f";{VAF}={round(variant.get_vaf(), precision)}"
        f";{SCORE}={round(extra_features[SCORE], precision)}"
        f";{COMPLEXITY}={round(extra_features[COMPLEXITY], precision)}"
        f";{SUPPORT}={round(extra_features[SUPPORT], precision)}"
        f";{OVERLAP}={round(extra_features[OVERLAP], precision)}"
        f";{CONTROL}={round(extra_features[CONTROL], precision)}"
        f";{WT_RU}={len(features.get_wt_ru_seq())}"
        f";{WT_RU_CNB}={features.get_wt_ru_cnb()}"
        f";{WT_RU_LEFT_CNB}={features.get_wt_ru_cnb_left()}"
        f";{WT_RU_RIGHT_CNB}={features.get_wt_ru_cnb_right()}"
        f";{V_RU}={len(features.get_v_ru_seq())}"
        f";{V_RU_CNB}={features.get_v_ru_cnb()}"
        f";{V_RU_LEFT_CNB}={features.get_v_ru_cnb_left()}"
        f";{V_RU_RIGHT_CNB}={features.get_v_ru_cnb_right()}"
        f";{HP_LEFT_LEN}={hp_left[1]}"
        f";{HP_LEFT_BASE}={hp_left[0]}"
        f";{HP_RIGHT_LEN}={hp_right[1]}"
        f";{HP_RIGHT_BASE}={hp_right[0]}"
        f";{SOURCE_COV}={features.get_source_coverage()}"
        f";{TOTAL_COV}={features.get_total_support()}"
        f";{MAX_COV}={features.get_max_support()}"
        f";{SOURCE}={','.join(features.get_source())}"
        f";{SAMPLE}={extra_features[SAMPLE]}"
        f";{RUN_ID}={extra_features[RUN_ID]}"
        f";{RUN_NAME}={extra_features[RUN_NAME]}"
    )
    return vcf_str


def vcf_write_variants_features(out_file, variants_list, precision=VAF_PRECISION, sorting=None):
    """
    Write a list of variants and features into a vcf file

    :param: out_file (str): path to output file
    :param: variants_list (list(Variant, VariantFeatures)): input list
    :param: precision (int): precision of floating numbers
    :param: sorting (None or str in VCF_SORT.keys()): None implies no sorting
    """
    out_vcf = open(out_file, 'w')
    out_vcf.write(vcf_header())
    if sorting is not None:
        sort_arg = VCF_SORT[sorting]
        variants_list.sort(key=lambda x: sort_arg[0](x), reverse=sort_arg[1])
    for (variant, features, extra_features) in variants_list:
        out_vcf.write(f"\n{vcf_variant(variant, features, extra_features, precision=precision)}")
    out_vcf.close()

def vcf_write_df(out_file, variants_df, precision=VAF_PRECISION, sorting=None, append=False):
    """
    Write a list of variants and features into a vcf file

    :param: out_file (str): path to output file
    :param: variants_df (DataFrame): input dataframe
    :param: precision (int): precision of floating numbers
    :param: sorting (None or str in VCF_SORT.keys()): None implies no sorting
    """
    if not append:
        out_vcf = open(out_file, 'w')
        out_vcf.write(vcf_header())
    else:
        out_vcf = open(out_file, 'a')
    if sorting is not None:
        sort_arg = VCF_SORT[sorting]
        variants_df.sort_values(by=['CHROM', 'POS'])
    for _, row in variants_df.iterrows():
        row_str = '\t'.join([str(value[1]) for value in row.items()])
        out_vcf.write(f"\n{row_str}")
    out_vcf.close()

def vcf_import_df(vcf_file):
    """
    Import a VCF file into a dataframe
    :param: vcf_file (str): path to input VCF file
    :return: DataFrame
    """
    with open(vcf_file, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM', })



def dump_vcf_to_tsv(vcf_file, out_tsv_file, append=False):
    if not append:
        out_tsv = open(out_tsv_file, 'w')
        out_tsv.write(tsv_header())
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        chrom = record.CHROM
        pos = record.POS
        ref = str(record.REF)
        alt = str(record.ALT[0])
        sample = record.INFO[SAMPLE][0]
        run_id = record.INFO[RUN_ID][0]
        run_name = record.INFO[RUN_NAME][0]
        vaf = round(record.INFO[VAF], VAF_precision)
        v_type = record.INFO[V_TYPE]
        score = round(record.INFO[SCORE], VAF_precision)
        complexity = round(record.INFO[COMPLEXITY], VAF_precision)
        support = round(record.INFO[SUPPORT], VAF_precision)
        overlap = round(record.INFO[OVERLAP], VAF_precision)
        control = round(record.INFO[CONTROL], VAF_precision)
        cov = record.INFO[SOURCE_COV]
        alt_cov = record.INFO[TOTAL_COV]
        max_cov = record.INFO[MAX_COV]
        annotation = ','.join(record.INFO[ANNOTATION])
        fields = [sample, run_id, run_name, chrom, pos, ref, alt, vaf, v_type, score, complexity, support, overlap, control, cov, alt_cov, max_cov, annotation]
        out_tsv.write('\n' + '\t'.join([str(x) for x in fields]))
    return False


def run_snpeff(parameters, in_vcf_file):
    """
    Runs snpEff to annotate a VCF file

    :param: parameters (Parameters): pipeline parameters
    :param: in_vcf_file (str): path to input VCF file
    """
    snpEff_parameters = parameters.get_snpeff_parameters()
    snpEff_dir = snpEff_parameters[SNPEFF_PATH]
    if not os.path.isdir(os.path.join(snpEff_dir, 'data', 'hg19')):
        wd = os.getcwd()
        os.chdir(snpEff_dir)
        cmd = ['java', '-jar', 'snpEff.jar', 'download', 'hg19']
        subprocess.call(cmd)
        os.chdir(wd)
    out_vcf_file = in_vcf_file.replace('.vcf', '_snpeff.vcf')
    cmd = ['java', '-Xmx4g', '-jar', os.path.join(snpEff_dir, 'snpEff.jar')]
    cmd += ['eff', '-dataDir', os.path.join(snpEff_dir, 'data')]
    cmd += ['-nodownload', '-noInteraction', '-noMotif', '-noNextProt', '-noLog']
    cmd += ['-v', '-hgvs', 'hg19', in_vcf_file]
    with open(out_vcf_file, 'w') as f:
        subprocess.call(cmd, stdout=f)
    # Deleting snpEff log files
    for snpEff_log in snpEff_parameters[SNPEFF_LOG]:
        if os.path.isfile(snpEff_log):
            os.remove(snpEff_log)
