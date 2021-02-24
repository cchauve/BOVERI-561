"""
Separate expected indels into wildtype indels and not wildtype indels
"""

import sys
import csv

INPUT_FILE = sys.argv[1]
HEADER = open(INPUT_FILE, 'r').readlines()[0].rstrip()

NO_WILDTYPE_FILE_NAME = INPUT_FILE.replace('.tsv', '_NO_WT.tsv')
NO_WILDTYPE_FILE = open(NO_WILDTYPE_FILE_NAME, 'w')
NO_WILDTYPE_FILE.write(HEADER)
WILDTYPE_FILE_NAME = INPUT_FILE.replace('.tsv', '_WT.tsv')
WILDTYPE_FILE = open(WILDTYPE_FILE_NAME, 'w')
WILDTYPE_FILE.write(HEADER)


with open(INPUT_FILE) as csvfile:
    expected_indels_data = csv.DictReader(csvfile, delimiter='\t')
    for row in expected_indels_data:
        vaf = float(row['exp_vaf'])
        row_str = '\t'.join([str(x) for x in row.values()])
        if vaf > 0:
            NO_WILDTYPE_FILE.write('\n' + row_str)
        else:
            WILDTYPE_FILE.write('\n' + row_str)

NO_WILDTYPE_FILE.close()
WILDTYPE_FILE.close()
