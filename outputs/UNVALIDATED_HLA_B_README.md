# HLA-B tumor-vs-normal — NOT FOR PRESENTATION

Quarantined per Phase 1 (`outputs/PHASE1_BC_strand_check.md`): the exon5_skip/intron5_retention classifier could not be validated for HLA-B -- only 2 candidate AS events were found in the entire 935-tumor/106-normal cohort at the standard 10%-local-rate threshold (vs 60 for HLA-C, ~200+ for HLA-A), far too few to tell whether the naive or strand-flipped geometry is correct, or whether the exon5_skip/intron5_retention labels mean anything at all for this gene.

These files are the same computations Phase 2 runs for validated genes, kept here for audit/reference ONLY. **Do not fold any number from these files into RESULTS_SUMMARY.md, PRESENTATION_NOTES.md, or figures_overview.pdf** unless a manual strand/geometry review specifically validates HLA-B.

Pending manual review: (1) confirm with the alignment/reference-construction pipeline (outside this repo) whether personalized HLA-B contigs are built the same way as HLA-C's (Phase 1 found HLA-C's coordinates are already transcript-oriented, contradicting the initial reverse-strand-flip hypothesis) and (2) investigate why HLA-B's maximum observed Local_rate in this cohort is only 12.6% (barely above the 10% detection floor) vs 52-60% for A/C -- real biology or a technical/reference-quality issue specific to HLA-B.
