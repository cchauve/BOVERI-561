import sys
import csv

VAF_RANGES = [(0.0, 0.5), (0.5, 1.0), (1.0, 5.0), (5.0, 100.0)]
NG_RANGES = [(0.0, 4.0), (4.0, 8.0), (8.0, 16.0), (16.0, 32.0), (32.0, 100.0)]

INPUT_FILE = sys.argv[1]
HEADER = open(INPUT_FILE, 'r').readlines()[0].rstrip()

GRID = [(vaf, ng) for vaf in VAF_RANGES for ng in NG_RANGES]

def format_setting(vaf_range, ng_range):
    (v1, v2) = vaf_range
    (n1, n2) = ng_range
    return '_'.join([str(x) for x in [v1, v2, n1, n2]])

def value_to_range(val, range_list):
    for (x1, x2) in range_list:
        if x1 < val <= x2:
            return (x1, x2)
    return None

SETTING_2_FILE_NAME = {}
for (vaf_range, ng_range) in GRID:
    key = format_setting(vaf_range, ng_range)
    SETTING_2_FILE_NAME[key] = INPUT_FILE.replace('.tsv', f"_VAF_NG_{key}.tsv")

SETTING_2_FILE = {}
for key, file_name in SETTING_2_FILE_NAME.items():
    SETTING_2_FILE[key] = open(SETTING_2_FILE_NAME[key], 'w')
    SETTING_2_FILE[key].write(HEADER)

with open(INPUT_FILE) as csvfile:
    expected_indels_data = csv.DictReader(csvfile, delimiter='\t')
    for row in expected_indels_data:
        vaf = float(row['exp_vaf'])
        if vaf > 0:
            ng = float(row['ng_est'])
            vaf_range = value_to_range(vaf, VAF_RANGES)
            ng_range = value_to_range(ng, NG_RANGES)
            key = format_setting(vaf_range, ng_range)
            SETTING_2_FILE[key].write('\n' + '\t'.join(row.values()))

for key, file_name in SETTING_2_FILE.items():
    SETTING_2_FILE[key].close()
