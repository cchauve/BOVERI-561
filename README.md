# BOVERI-561
Extract indel calls from commercial samples.

The list of commercial samples run IDs is in BOVERI-532. The corresponding input
file is in data/BOVERI-532.csv.

The code in this repo generates one TSV file per run, containing the list of
all (almost) unfiltered indels. The only filters that have been applied are
- discarding alignments with at least 5 gaps from variants support,
- filtering out variants supported by an alignment in a low-quality region
  defined as a region of 3 bases pre and post breakpoint with an average quality
  within a read cluster that is at most 20 for all bases of the considered
  region.

Each entry in the TSV file for a run contains the usual fields of a variant
(chromosome, position, reference sequence, alternate sequence, VAF) plus the
following fields
- sample: sample ID where the variant was observed
- run_id: run ID containing this sample
- run_name: run name
- v_type: DEL, INS, DELINS, MNV
- score: penalty associated to the variant (high penalty implies low confidence)
- complexity: penalty associated to the complexity of the sequence around the
              breakpoint (5 bases pre and post)
- support: support penalty associated to the size of the largest cluster
           supporting the variant
- overlap: overlap penalty associated to the VAF of overlapping variants
- control: closest VAF of the same variant in a control sample (0.0 if not
           present in a control sample)
- cov: coverage of the amplicon(s) by merged reads (a merged read corresponds
       to two reads)
- alt_cov: coverage of the variant by merged reads
- max_cov: size (in number of merged reads) of the largest supporting cluster
- repeats: WT1:WT2:WT3:WT4,V1:V2:V3:V4
    WT1 = length of the repeat unit of the reference sequence
    WT2 = number of copies of the repeat unit in the reference sequence
    WT3 = number of copies of the repeat unit left of variant breakpoint
    WT4 = number of copies of repeat unit right of variant breakpoint
    V1, V2, V3,  V4: same for alternate sequence
- annotation: snpEff annotation

Moreover, the file results/BOVERI-532_out_1.tsv contains an analysis
of the expected indels and if they are found by the indel pipeline.
Note: the pipeline results were filtered to discard any indel call that was
found in a non-control sample and in a control sample (Blank or normal female)
with a VAF differeing by at most 5%.

Each line of this file corresponds to one expected indel, as indicated by the
vendor.
It contains the information about the sample expected to have this indel and the
indel itself, complemented by information about the indel as found by the
indels calling pipeline (na means the indel was not found):
- vaf: VAF of the indel as given by he pipeline
- score: penalty P associated to the indel call (0 means no penalty, very high
  confidence call)
- TP/FN/FP: number of true positive (expected indels found), false negative
  (missed expected indels) and false positive if one would call all indels
  with a score at most P.
- for the TP and FP categories, the file shows the average VAF and score of the
  corresponding indels and the number of indels with a VAF < 1% and the number
  of indels with a VAF < 0.5%.

Altogether, out of 726 expected indels, 156 are not found. A further analysis
would be to discard the ones from the wild-type samples and to see if they are
correlated with low input or low VAF samples.

A second interesting point is that the wide majority of FP are indels with a
VAF <0.5%.

The file results/BOVERI-532_out_2.tsv shows the same statistics for each run
as a function of the maximum score to accept a call. Here again we can see that
the number of FPs is large but FPs are mostly indels of very low VAF. Again the
large number of FNs raise the question to know if they are actually present in
the sequencing data.
