# Results index

[Project overview](../README.md) · [Results narrative](../docs/RESULTS_SUMMARY.md) · [Methods](../docs/METHODS.md)

## Main figures and tables

| Analysis | Figure | Table |
|---|---|---|
| HLA-A headline events | [View](fig_classes.png) | [Carrier rates and tests](table_headline_classes.csv) |
| HLA-C headline events | [View](fig_classes_HLA-C.png) | [Carrier rates and tests](table_headline_classes_HLA-C.csv) |
| HLA-A allele comparisons | [View](fig_per_allele.png) | [Per-allele results](table_per_allele.csv) |
| HLA-C allele comparisons | [View](fig_per_allele_HLA-C.png) | [Per-allele results](table_per_allele_HLA-C.csv) |
| HLA-A recurrent patterns | [Volcano](fig_volcano.png) | [Detailed tests and FDR](results_BRCA_tumor_vs_normal_HLA-A.csv) |
| HLA-C recurrent patterns | [Volcano](fig_volcano_HLA-C.png) | [Detailed tests and FDR](results_BRCA_tumor_vs_normal_HLA-C.csv) |
| Threshold sensitivity | [View](fig_threshold_sensitivity.png) | [Sweep results](table_threshold_sensitivity.csv) |
| Cohort-matched pan-cancer comparisons | [View](fig_phase5_matched_normal.png) | [Summary, including flagged HLA-B rows](table_phase5_summary.csv) |


## Supporting reports

- [Gene-level review](FINAL_REVIEW.md): findings, limitations, and excluded results.
- [HLA-B/C classifier checks](PHASE1_BC_strand_check.md): strand and classification diagnostics.
- [Pan-cancer descriptive analysis](PHASE3_pancancer_descriptive.md): cohort summaries and HLA-A*11:01 recurrence.
- [Cohort-matched comparisons](PHASE5_pancancer_matched_normal.md): normal-sample mapping and results.
- [Annotation artifact investigation](../docs/ANNOTATION_ARTIFACT_NOTE.md): rationale for the HLA-A*11:01 exclusion.
- [Extended write-up](MINI_PAPER.md) and [figure audit](FIGURE_AUDIT.md): supporting interpretation and figure checks.
- [Exploratory diagnostics](scratch/OPEN_QUESTIONS_DIAGNOSTIC.md): unresolved questions; not a standalone reproduction script.

## Filename conventions

| Name | Meaning |
|---|---|
| `fig_*` | Figures, typically PDF for export and PNG for browsing |
| `table_*` | Summaries, counts, diagnostics, or statistical comparisons |
| `results_BRCA_*` | Detailed per-pattern BRCA tests |
| `_HLA-C` | HLA-C results; original unsuffixed headline figures/tables refer to HLA-A |
| `phase1`, `phase3`, `phase5` | Classifier checks, descriptive pan-cancer analysis, and cohort-matched comparisons |
| `*_all_alleles` | Unfiltered audit copy, before headline allele exclusions |
| `UNVALIDATED_HLA_B_*` | HLA-B audit outputs; classifier unvalidated and results excluded from presented conclusions |
| `archive/` | Retained presentation/export copies; use the indexed figures above |

The [HLA-B status note](UNVALIDATED_HLA_B_README.md) explains the exclusion. Tables that combine genes may still contain HLA-B rows; respect their quarantine flags when selecting results.
