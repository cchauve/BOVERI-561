#!/usr/bin/env bash

# python bin/extract_files.py data/BOVERI-532.csv
# python bin/dump_variants.py data/BOVERI-532.csv

INPUT=`ls data/v4MiSeq_commercial_samples_expected_indels_ExpectedVAFs_MissingSeraSeq_VAF_NG*`

for I in ${INPUT}
do
    python bin/analyze_variants.py ${I} BOVERI-532_parameters-1.yaml
    python bin/analyze_variants.py ${I} BOVERI-532_parameters-2.yaml
    python bin/analyze_variants.py ${I} BOVERI-532_parameters-3.yaml
done
