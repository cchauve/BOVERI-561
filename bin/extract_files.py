#!/usr/bin/env python3
"""
Functions to extract the results files from AWS
"""

# Standard imports
import argparse
import csv
import fileinput
import os
import re
import shutil
import subprocess
import tarfile

# Third-party imports
import boto3

# Manifests
MANIFESTS = {
    'CG001Qv4': 'CG001v4.0_Amplicon_Manifest_Panel4.0.3_20181101.tsv',
    'CG001Qv5': 'CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv'
}
MANIFEST_KEY_LG = len(list(MANIFESTS.keys())[0])

# Default AWS parameters
AWS_DIR = 'cchauve-orchestration-ch'
AWS_CP = ['aws', 's3', 'cp']
RESULTS_DIR = 'results'

DATA_DIR = 'data'

def get_files_in_s3(prefix, s3_bucket):
    """
    Get a list of files from the indels pipeline output of a run.
    :param: prefix (str): path to run data, e.g.
    input/201014_M03829_0366_000000000-JBV6Y
    :param: s3_bucket (str): S3 bucket containing files to read,
    e.g. 'cchauve-orchestration-ch', 'ch-testdata'

    :return: list(str): file paths of files in directory prefix;
    None if the directory is empty or does not exist
    """
    s3_client = boto3.client('s3')
    s3_objects = s3_client.list_objects_v2(Bucket=s3_bucket, Prefix=prefix)
    if s3_objects['KeyCount'] == 0:
        return None
    else:
        return [obj['Key'] for obj in s3_objects['Contents']]

def get_runs_manifests_list(runs_csv_file):
    """
    Get from a CSV file containing runs ID and runs name a list of pairs
    (run ID, manifest)
    :param: runs_csv_file (str): path to the input CSV file
    :assumption: column 0 is run name, column 1 is run ID

    :return: list((str, str)): list (run_id, manifest)
    """
    result = []
    with open(runs_csv_file) as csvfile:
        runs_data = csv.reader(csvfile, delimiter=',')
        for row in runs_data:
            run_id = row[1]
            run_name = row[0]
            manifest = MANIFESTS[run_name[0:MANIFEST_KEY_LG]]
            result.append((run_id, manifest, run_name))
    return result

PARAM_EXT = '.yaml'
MAIN_EXT = '_main.tar.gz'
VCF_EXT = '_indels_unfiltered.vcf'

def extract_files(files_list, s3_bucket, out_dir, manifest_file):
    for file_path in files_list:
        file_dir, file_name = os.path.split(file_path)
        if file_name.endswith(MAIN_EXT):
            sample_id = file_name.replace(MAIN_EXT, '')
            s3_file_path = f"s3://{s3_bucket}/{file_path}"
            aws_cp_cmd = AWS_CP + [s3_file_path, out_dir]
            subprocess.call(aws_cp_cmd)
            tgz_file_path = os.path.join(out_dir, file_name)
            tarfile.open(tgz_file_path, 'r:gz').extractall(path=out_dir)
            os.remove(tgz_file_path)
        elif file_name.endswith(PARAM_EXT):
            s3_file_path = f"s3://{s3_bucket}/{file_path}"
            parameters_file_path = os.path.join(out_dir, file_name)
            aws_cp_cmd = AWS_CP + [s3_file_path, out_dir]
            subprocess.call(aws_cp_cmd)
            for line in fileinput.input(parameters_file_path, inplace=True):
                if line[0] == '#':
                    print(line.rstrip())
                else:
                    line_split = line.rstrip().split()
                    if line_split[0] == 'path_input:':
                        print(f"path_input: {out_dir}")
                    elif line_split[0] == 'path_output:':
                        print(f"path_output: {out_dir}")
                    elif line_split[0] == 'snpeff_path:':
                        print('snpeff_path: /home/cchauve/bin/snpEff')
                    elif line_split[0] == 'path_manifest:':
                        manifest_path = os.path.join(out_dir, manifest_file)
                        print(f"manifest_path: {manifest_path}")
                    else:
                        print(line.rstrip())


if __name__ == "__main__":
    """
    Arguments:
    - runs_csv_file: CSV file with 2 fields <run_name>,<run_id>
      run_name is used to define the amplicon manifest
    - s3_output: S3 bucket where to store the results
      results for run_id go into <s3_output>/<run_id>
    """
    # Input file
    ARGS_RUNS_FILE = ['runs_csv_file', None, 'Runs CSV file']
    # S3 directory where to store the results
    ARGS_OUTPUT_BUCKET = ['-o', '--s3_output', 'output S3 bucket directory']


    parser = argparse.ArgumentParser(description='Indels testing: extract results')
    parser.add_argument(ARGS_RUNS_FILE[0], type=str, help=ARGS_RUNS_FILE[2])
    parser.add_argument(ARGS_OUTPUT_BUCKET[0],
                        ARGS_OUTPUT_BUCKET[1],
                        type=str,
                        default=AWS_DIR,
                        help=ARGS_OUTPUT_BUCKET[2])
    args = parser.parse_args()
    s3_bucket = args.s3_output

    # Creating a log file located in the same directory than the YAML
    # configuration file and with the same name with .yaml replaced by .log
    _, run_file_name = os.path.split(args.runs_csv_file)
    log_file_path = os.path.join('log',
                                 run_file_name.replace('.csv', '_input.log'))
    log_file = open(log_file_path, 'w')

    runs_manifests_list = get_runs_manifests_list(args.runs_csv_file)
    for (run_id, manifest, run_name) in runs_manifests_list:
        out_dir = os.path.join(RESULTS_DIR, run_id)
        os.makedirs(out_dir, exist_ok=True)
        log_file.write(f"RUN:{run_id}.{run_name}\n")
        run_files = get_files_in_s3(run_id, s3_bucket)
        extract_files(run_files, s3_bucket, out_dir, manifest)
        manifest_file_path1 = os.path.join(DATA_DIR, manifest)
        manifest_file_path2 = os.path.join(out_dir, manifest)
        shutil.copyfile(manifest_file_path1, manifest_file_path2)
    log_file.close()
