# Rabadan HLA splicing project — canonical folder

ALL work lives here: `~/Downloads/rabadan-splicing-project/`
- `data/`      input TSVs (do not edit) — `deprecated/` holds superseded files, kept for reference only
- `outputs/`   every generated figure/CSV from the scripts below
- root         scripts + the presentation-ready write-up (.md/.pdf)

## Data files (13 cols, tab-separated, NO header)
`Cohort, Type, Barcode, Allele, Contig, First_base, Last_base, Strand,
Motif, Annotated, Uniquely_mapping, Multimapping, Overhang`

- `data/export_BRCA_tumor_HLA-A.tsv`   935 BRCA tumor samples  (Type=tumor)
- `data/export_BRCA_normal_HLA-A.tsv`  106 BRCA normal samples (Type=normal)
- `data/pan_HLA-A_SJ_combined.tsv`     pan-cancer, 6 cohorts (NORMAL pool +
  TCGA-CHOL/HNSC/KIRC/LGG/UCEC — no BRCA) — used only as a "does the pipeline
  run unmodified on a second dataset" check, not a BRCA result
- `data/deprecated/brca_HLA-A_SJ_combined.tsv` — IGNORE, no Type column, all-tumor

## Scripts (run in this order to reproduce everything)
1. `analyze_brca.py` — main analysis: cohort overview, headline tumor-vs-normal
   comparison (CI, odds ratio, Fisher p) for exon5_skip/intron5_retention,
   per-allele view, volcano plot, per-sample event burden, FDR-corrected
   per-pattern table, heatmaps. All thresholds are CLI flags
   (`--local-rate-threshold`, `--exon-skip-min-jump`,
   `--intron-retention-max-shift`, `--min-allele-samples`, `--min-carriers`,
   `--exclude-alleles`) — see its module docstring for the full list and for
   how to point it at HLA-B/C or the pan-cancer file later.
2. `diagnostic_annotation_check.py` — investigates any allele whose top event
   looks suspiciously universal (currently HLA-A\*11:01): raw SJ rows,
   variant-vs-annotated read share, per-allele penetrance, plain verdict.
   Diagnostic only, does not modify `outputs/`.
3. `additional_analyses.py` — power analysis (achieved power + minimum
   detectable effect for the headline test), threshold-sensitivity sweep
   (does the conclusion hold across reasonable cutoff choices), and the
   HLA-A allele frequency spectrum (cohort background context).
4. `make_figures_overview.py` — composites the headline/per-allele/volcano
   results into `figures_overview.pdf` (panels A/B/C) for slides.
5. `make_background_figures.py` — draws `background_method.pdf`, a pure
   schematic (not data-driven) explaining the pipeline steps and the
   exon5_skip/intron5_retention coordinate geometry, for audiences who need
   the method explained before the results.

## Write-up (project root)
- `RESULTS_SUMMARY.md` — the full narrative: cohort → pipeline generalization
  → headline test (with power analysis + threshold sensitivity) → HLA-A\*11:01
  diagnosis → bottom line.
- `PRESENTATION_NOTES.md` — 7 claim/stat/caveat bullets + a next-steps slide,
  framed against the MHC Hammer paper.
- `ANNOTATION_ARTIFACT_NOTE.md` — the HLA-A\*11:01 writeup on its own.
- `figures_overview.pdf` — results slide panel (A: headline classes, B:
  per-allele contrast, C: volcano).
- `background_method.pdf` — background/method schematic (pipeline flowchart
  + splicing-geometry cartoon).

## Key outputs/ files
- `table_headline_classes.csv` (presented, A\*11:01 excluded) /
  `table_headline_classes_all_alleles.csv` (audit copy, unfiltered)
- `table_power_analysis.csv`, `table_threshold_sensitivity.csv` +
  `fig_threshold_sensitivity.pdf`
- `table_per_allele.csv` + `fig_per_allele.pdf`
- `results_BRCA_tumor_vs_normal_HLA-A.csv` (full per-pattern table + FDR q)
- `table_volcano.csv` + `fig_volcano.pdf`
- `table_allele_frequency_spectrum.csv` + `fig_allele_frequency_spectrum.pdf`
- `annotation_artifact_check.csv` + `annotation_artifact_raw_rows_hla_a_11_01.csv`
- `fig_cohort_overview.pdf`, `fig_event_burden.pdf`, `fig_heatmap_tumor.pdf`,
  `fig_heatmap_normal.pdf`

Always write new generated files into `outputs/` (or the project root for
slide-ready composites, matching the pattern above) — never into `data/`.
