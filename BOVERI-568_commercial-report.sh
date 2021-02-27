#!/usr/bin/env bash

# python bin/extract_files.py data/BOVERI-568_commercial.csv
# python bin/dump_variants.py data/BOVERI-568_commercial.csv -m CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv
python bin/analyze_variants.py data/v51NextSeq_commercial_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial-report.yaml

