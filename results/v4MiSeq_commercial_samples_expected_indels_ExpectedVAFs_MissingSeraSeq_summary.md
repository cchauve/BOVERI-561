# Analysis of the new indel caller on Horizon and SeraSeq commercial dilution series
Cedric Chauve
February 20, 2021

This document describes the results of the new indel caller on a set of runs
obtained from Horizon and SeraSeq commercial samples, diluted at various levels
of DNA input and expected VAF.

## Methods and data

### Data
There are 15 different set of samples that were analyzed in this experiment,
defined by a range of DNA input and a range of expected VAF:
- DNA input (in ng): (0, 4], (4, 8], (8, 16], (16, 32], (32, 64]
- expected VAF: (0.0, 0.5], (0.5, 1.0], (1.0, 5.0]
Note: there were no expected indel with expected VAF above 5%.

For each setting defined by a DNA input range and an expected VAF, all samples
sequenced with that setting were considered as a dataset. The LLOD was defined
as max(0.25, 0.5 * min_VAF) where min_VAF is the lowest expected VAF among all
expected indels in these samples.

For these samples, the indels detected by the new indel caller were extracted
from the pipeline result files. An indel that was dectected by the indel caller
is said to be a "detected indel", independent of its VAF (that might make it be
filtered out if it is below the LLOD threshold) or score.

### Methods
To each detected indel that has a VAF at least the LLOD threshold, a score
(that should  actually be called a penalty) is associated defined as the sum
of 4 terms:
- weighted sequence complexity penalty, defined as the sequence complexity
  penalty weighted by a multiplicative factor W in [0.0, 1.0]
- support penalty (nb reads in largest supporting cluster / number of supporting
  reads),
- overlap penalty (min(1.0, cumulative VAF of overlapping indels / VAF)),
- control penalty (min(1.0, closest VAF observed in a control sample / VAF)),

For a given score threshold S, every detected with a VAF >= LLOD and a
score <= X is said ot be a "called indel".

The performance of the indel caller was explored, for each setting, within a
grid of values for S and W:
- S ranges from 0.0 to 2.0 by increments of 0.1,
- W ranges from 0.0 to 1.0 by increments of 0.1.

The following classes of indels are then recorder:
- TP: an expected indel that is called,
- FP: a called indel which is not an expected indel,
- TN: a detected indel with VAF >= LLOD which is not expected and not called,
- FN_u: an expected indel which is not dectected (and thus not called),
- FN_d: an expected indel which is detected but not called (score > S or VAF < LLOD).
- FN = FN_u + FN_d

Remark. Three indels observed as high frequency FP indels in a previous
experiment were also filtered out from detected indels:
- chr19 4117563 GCTT  G
- chr15 90631917  TC T
- chr12 25398297  C CCG

From these statistics, for each setting the following quantities were computed:
- sensitivity: TP / (TP + FN),
- specificity: TN / (TN + FP),
- accuracy: (TP + TN) / (TP + TN + FP + FN),
- precision / Positive Predictive Value (PPV): TP / (TP + FP),
- recall: TP / (TP + FN)
- F1: 2 * precision * recall / (precision + recall),
- FDR: FP / (TP + FP) = 1.0 - precision = 1.0 - PPV,

Note that these statistics are somewhat biased by the expected indels that are
undetected FNs, as this implies these indels are very likely not present in the
data and so lower artificially the statistics.

## Results

For the various combinations of DNA input and expected VAF I show below the best
results obtained with the new indel caller. To do so, I recorder for each
setting the best F1 score and list all combinations of (score, sequence
complexity weight) that lead to an F1 score at least 95% of the best one.
I use the F1 score as it does not depend on the number of TN that is generally
high and flattens the specificity.

In the tables below I show the results for all such combinations (score, w_comp)
split per setting (expected VAF, DNA input amount). It is followed by the list
of combinations that appear in at least 13 settings, suggesting these are
choices for the score threshold and sequence complexity weight that work well
over a large number of settings (there are 15 settings).

The resulting table is provided in the file
```
v4MiSeq_commercial_samples_expected_indels_ExpectedVAFs_MissingSeraSeq_VAF_NG_F1_095.txt
```

The top-ranked combinations, that are called *optimal thresholds*  are
```
score   w_comp  frequency
0.6     0.5     13
0.6     0.4     13
0.7     0.7     13
0.6     0.0     13
0.6     0.1     13
0.6     0.2     13
0.7     0.3     13
0.7     0.4     13
0.7     0.5     13
0.8     0.8     13
0.9     1.0     13
0.5     0.2     14
0.5     0.3     14
0.5     0.0     14
0.5     0.1     14
0.6     0.3     14
0.8     0.7     14
0.4     0.0     15
```

The score threshold ranges from 0.4 (with sequence complexity not accounted for
with a weight of 0) to 0.9 with sequence complexity fully accounted for.
Generally it shows that the score threshold needs to be stringent (the maximum
possible score is 4.0). For low settings (expected VAF and DNA input), higher
scores and lower sequence complexity weight lead to a much higher rate of FPs.

As expected as the score threshold increases, the sequence complexity weight
decreases, showing a trade-off between both parameters. The range of scores
in [0.5, 0.6] dominates this list, with a sequence complexity weight below 0.5.

Last I looked at all three types of errors, FP, FN_u and FN_d, over all optimal
thresholds, and recorded for each error its frequency (how many times it does
appear in the list of errors), average and max VAF, complexity score, support
score, overlap score and control score. The results are in the files
```
v4MiSeq_commercial_samples_expected_indels_ExpectedVAFs_MissingSeraSeq_optimal_thresholds_with_blacklist_FP.txt
v4MiSeq_commercial_samples_expected_indels_ExpectedVAFs_MissingSeraSeq_optimal_thresholds_with_blacklist_FN_u.txt
v4MiSeq_commercial_samples_expected_indels_ExpectedVAFs_MissingSeraSeq_optimal_thresholds_with_blacklist_FN_d.txt
```

For FP, there are 107 different indel that are called at least once as a FP over
all runs with optimal thresholds, for a total of 2,782 FP calls. Among them, a
set of 4 indel calls have a much higher frequency than the other ones:
```
freq.   chr     pos     ref     alt     avg_vaf max_vaf avg_comp.       max_comp.       avg_supp.       max_supp.       avg_overlap     max_overlap   avg_ctrl        max_ctrl
181     chr22   33559506        C       CCA     0.313   0.413   0.348   0.348   0.117   0.25    0.215   0.454   0.0     0.0
147     chr7    55228049        G       GC      0.541   1.128   0.143   0.143   0.069   0.286   0.038   0.167   0.009   0.485
113     chr7    116411729       AT      A       0.375   0.545   0.286   0.286   0.111   0.333   0.0     0.0     0.268   0.552
109     chr7    116411671       GT      G       0.385   0.548   0.429   0.429   0.006   0.6     0.0     0.0     0.243   0.358
96      chr11   534294  CCCA    C       0.438   0.78    0.44    0.44    0.137   0.286   0.057   0.251   0.058   0.379
```

Generally FP calls have a low VAF, with only 4 indels having a maximum VAF of
at least 1.0, accounting for 213 FPs. If the maximum VAF threshold is lowered
to 0.5, there are 31 indels accounting for 1,169 of the FPs.

The number of FN is naturally smaller and they are shown below: first the
undetected ones.
```
freq.   chr     pos     ref     alt     avg_vaf max_vaf avg_comp.       max_comp.       avg_supp.       max_supp.       avg_overlap     max_overlap     avg_ctrl        max_ctrl
270     chr7    55242464        AGGAATTAAGAGAAGC        A       0.547   2.0     nan     nan     nan     nan     nan     nan     nan     nan
252     chr7    55248998        A       ATGGCCAGCG      0.357   0.75    nan     nan     nan     nan     nan     nan     nan     nan
108     chr7    55249012        C       CGGT    0.312   1.0     nan     nan     nan     nan     nan     nan     nan     nan
90      chr7    55242465        GGAATTAAGAGAAGCA        G       0.225   0.5     nan     nan     nan     nan     nan     nan     nan     nan
90      chr17   7577557 AG      A       0.575   1.0     nan     nan     nan     nan     nan     nan     nan     nan
90      chr17   37880981        A       AGCATACGTGATG   0.525   1.0     nan     nan     nan     nan     nan     nan     nan     nan
54      chr4    55141048        T       TA      0.458   0.75    nan     nan     nan     nan     nan     nan     nan     nan
18      chr17   7579419 AG      A       1.0     1.0     nan     nan     nan     nan     nan     nan     nan     nan
```
Next the detected ones (note: there is a mistake in my scripts and the shown VAF
is the expected VAF not the one from the pipeline):
```
freq.   chr     pos     ref     alt     avg_vaf max_vaf avg_comp.       max_comp.       avg_supp.       max_supp.       avg_overlap     max_overlap     avg_ctrl        max_ctrl
655     chr7    55248998        A       ATGGCCAGCG      0.659   3.0     0.351   0.351   0.151   0.7     0.041   1.0     0.0     0.0
407     chr17   7579419 AG      A       0.491   2.0     0.381   0.381   0.302   0.778   0.361   1.0     0.29    1.0
361     chr7    55242464        AGGAATTAAGAGAAGC        A       1.001   5.0     0.388   0.388   0.141   0.375   0.407   1.0     0.0     0.0
306     chr17   37880981        A       AGCATACGTGATG   0.478   1.0     0.302   0.302   0.16    0.5     0.078   0.502   0.0     0.0
216     chr4    55141048        T       TA      0.417   1.0     0.158   0.158   0.083   0.267   0.015   0.066   0.0     0.0
181     chr17   7577557 AG      A       0.499   1.0     0.053   0.053   0.11    0.429   0.0     0.0     0.0     0.0
154     chr7    55242465        GGAATTAAGAGAAGCA        G       0.591   1.0     0.388   0.388   0.094   0.333   0.341   0.996   0.0     0.0
144     chr7    55249012        C       CGGT    0.516   1.0     0.211   0.211   0.167   0.333   0.082   0.659   0.0     0.0
```
