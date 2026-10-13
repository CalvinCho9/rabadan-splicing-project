# Analysis guide

See [README.md](../README.md) for the project overview and quick start.
Paths below are relative to the repository root.
- `data/`      input TSVs (do not edit) — `deprecated/` holds superseded files, kept for reference only
- `outputs/`   every generated figure/CSV from the scripts below
- `scripts/`   analysis and figure-generation code
- `docs/`      methods, findings, and analysis records
- `presentation/` presentation notes and optional figure-generation instructions

## Data files (13 cols, tab-separated, NO header)
`Cohort, Type, Barcode, Allele, Contig, First_base, Last_base, Strand,
Motif, Annotated, Uniquely_mapping, Multimapping, Overhang`

- `data/export_BRCA_tumor_HLA-A.tsv`   935 BRCA tumor samples  (Type=tumor)
- `data/export_BRCA_normal_HLA-A.tsv`  106 BRCA normal samples (Type=normal)
- `data/export_BRCA_{tumor,normal}_HLA-{B,C}.tsv` — same schema, HLA-B/C
  (935/106 tumor/normal for B; 933/106 for C). See
  `outputs/PHASE1_BC_strand_check.md` before trusting Event_class labels on
  these — **HLA-B is quarantined** (classifier unvalidated), HLA-C is
  validated.
- `data/pan_HLA-A_SJ_combined.tsv`     pan-cancer, 6 cohorts (NORMAL pool +
  TCGA-CHOL/HNSC/KIRC/LGG/UCEC — no BRCA) — used only as a "does the pipeline
  run unmodified on a second dataset" check, not a BRCA result
- `data/pan_HLA-{B,C}_SJ_combined.tsv` — pan-cancer, 8 cohorts each,
  including a `Cohort=="BRCA_NORMAL"` block. **That block's `Type` column is
  mislabeled "tumor"** (row count matches the real BRCA normal export
  exactly) — see `docs/ANALYSIS_LOG.md`. `scripts/phase3_pancancer_descriptive.py`
  excludes it; do the same in any new script against these files.
- `data/pancancer_status.txt` — per-cohort sample-count reference table,
  not machine-parsed by any script here.
- `data/deprecated/brca_HLA-A_SJ_combined.tsv` — IGNORE, no Type column, all-tumor
- `data/gdc_sample_sheet_normal_with_hla.tsv` — GDC sample sheet (17 cols,
  WITH header) for 740 normal samples across 24 TCGA projects, with
  HLA-A/B/C typing per sample. Used to resolve the pan-cancer
  files' pooled `Cohort=="NORMAL"` samples (see above) back to their true
  TCGA project of origin via `Case ID` -> `Project ID` (145/145 barcodes
  resolved for the HLA-A pan file, 318/318 for HLA-B/C). Used by
  `scripts/phase5_pancancer_matched_normal.py` — see that script's docstring and
  `outputs/PHASE5_pancancer_matched_normal.md`.

## Scripts (run in this order to reproduce everything)
1. `scripts/analyze_brca.py` — main analysis: cohort overview, headline tumor-vs-normal
   comparison (CI, odds ratio, Fisher p) for exon5_skip/intron5_retention,
   per-allele view, volcano plot, per-sample event burden, FDR-corrected
   per-pattern table, heatmaps. All thresholds are CLI flags
   (`--local-rate-threshold`, `--exon-skip-min-jump`,
   `--intron-retention-max-shift`, `--min-allele-samples`, `--min-carriers`,
   `--exclude-alleles`) — see its module docstring for the full list and for
   how to point it at HLA-B/C or the pan-cancer file later.
2. `scripts/diagnostic_annotation_check.py` — investigates any allele whose top event
   looks suspiciously universal (currently HLA-A\*11:01): raw SJ rows,
   variant-vs-annotated read share, per-allele penetrance, plain verdict.
   Writes diagnostic CSVs to `outputs/`; does not change the main result tables.
3. `scripts/additional_analyses.py` — power analysis (achieved power + minimum
   detectable effect for the headline test), threshold-sensitivity sweep
   (does the conclusion hold across reasonable cutoff choices), and the
   HLA-A allele frequency spectrum (cohort background context).
4. `scripts/make_figures_overview.py` — composites the headline/per-allele/volcano
   results into `presentation/figures_overview.pdf` (panels A/B/C) for slides.
5. `scripts/make_background_figures.py` — draws `presentation/background_method.pdf`, a pure
   schematic (not data-driven) explaining the pipeline steps and the
   exon5_skip/intron5_retention coordinate geometry, for audiences who need
   the method explained before the results.

### HLA-B/C extension + pan-cancer pass (added, see `docs/ANALYSIS_LOG.md`)

6. `scripts/hla_common.py` — shared helpers (not a standalone script): generic
   tumor+normal loader, a strand-aware variant of the classifier
   (`classify_events_strand_aware`, tested and NOT needed for either B or
   C — see #7), and the generalized (allele-parametrized) version of
   `scripts/diagnostic_annotation_check.py`'s constitutive-vs-variable check.
   Imports `analyze_brca` as a library; never modifies it.
7. `scripts/phase1_bc_strand_check.py` — validates the exon5_skip/intron5_retention
   classifier for HLA-B/C by testing naive vs. a strand-mirrored
   classification against both genes and checking which produces a
   coherent signal. Writes the per-gene verdict to
   `outputs/PHASE1_BC_strand_check.md`. **Run this, and read its verdict,
   before trusting Event_class labels on any new HLA-B/C data** — the
   result was gene-specific and counter to the naive strand hypothesis.
8. `scripts/phase2_bc_headline.py` — HLA-B/C tumor-vs-normal (headline + per-allele
   + volcano + FDR), same statistics as `scripts/analyze_brca.py`, gated by Phase
   1's verdict: HLA-C's outputs are presented normally
   (`outputs/*_HLA-C.*`); HLA-B's are quarantined
   (`outputs/UNVALIDATED_HLA_B_*`, not for presentation).
9. `scripts/phase3_pancancer_descriptive.py` — per-cohort, per-allele penetrance
   across the pan-cancer files, descriptive only (no tumor-vs-normal
   claims). Includes the A\*11:01 cross-cohort recurrence check. Excludes
   the mislabeled `BRCA_NORMAL` block from the B/C pan files. Same
   gene-level quarantine as Phase 2.
10. `scripts/phase5_pancancer_matched_normal.py` — uses
    `data/gdc_sample_sheet_normal_with_hla.tsv` to resolve Phase 3's pooled
    `NORMAL` cohort to each sample's true TCGA project of origin, turning
    the pan-cancer pass from purely descriptive into a real matched-cohort
    tumor-vs-normal test (same Fisher/Wilson/OR machinery as Phase 2) for
    every cohort with 3+ samples on both sides after resolution: TCGA-CHOL
    and TCGA-HNSC for HLA-A, plus TCGA-KIRC for HLA-B/C. Re-applies the
    A\*11:01 exclusion unconditionally (catches it even in cohorts too
    small for the per-cohort re-derivation to flag it on its own — see the
    script's docstring for a near-miss this caught) and re-runs the
    per-cohort constitutive check on top for any new artifact. Same
    gene-level quarantine as Phase 2/3.

Running order for a from-scratch reproduction of everything:
`scripts/analyze_brca.py` → `scripts/diagnostic_annotation_check.py` →
`scripts/additional_analyses.py` → `scripts/phase1_bc_strand_check.py` →
`scripts/phase2_bc_headline.py` → `scripts/phase3_pancancer_descriptive.py` →
`scripts/phase5_pancancer_matched_normal.py` →
`scripts/make_figures_overview.py` → `scripts/make_background_figures.py`.

## Write-up and presentation

The two presentation PDFs listed below are optional generated outputs, not
committed files. Historical reports retain references to their slide panels.
The research figures remain available in `outputs/`.
- `docs/RESULTS_SUMMARY.md` — the full narrative: cohort → pipeline generalization
  → headline test (with power analysis + threshold sensitivity) → HLA-A\*11:01
  diagnosis → HLA-C headline test → A\*11:01 pan-cancer generalization →
  matched-normal pan-cancer extension (Finding 6) → bottom line. HLA-B is
  deliberately absent (quarantined).
- `presentation/PRESENTATION_NOTES.md` — 12 claim/stat/caveat bullets + a next-steps
  slide, framed against the MHC Hammer paper (bullets 1-7 original HLA-A,
  8-10 the B/C extension, 11-12 the Phase 5 matched-normal extension).
- `docs/ANNOTATION_ARTIFACT_NOTE.md` — the HLA-A\*11:01 writeup on its own
  (BRCA-only; see `outputs/PHASE3_pancancer_descriptive.md` for the
  pan-cancer follow-up).
- `docs/ANALYSIS_LOG.md` — chronological run log for the B/C extension +
  Phase 5: every assumption, judgment call, and thing to double-check, plus
  a "top things to review" list.
- `outputs/FINAL_REVIEW.md` — one paragraph per gene: what's real / null /
  artifact / underpowered / quarantined.
- `presentation/figures_overview.pdf` — results slide panel, HLA-A + HLA-C (A1/A2:
  headline classes per gene, B: HLA-A per-allele contrast, C: combined
  volcano, D: A\*11:01 pan-cancer recurrence). HLA-B intentionally excluded.
  **Not updated for Phase 5** (all-null results there don't earn a 6th panel
  in an already-tuned figure) — see `outputs/fig_phase5_matched_normal.pdf`
  and `outputs/PHASE5_pancancer_matched_normal.md` instead.
- `presentation/background_method.pdf` — background/method schematic (pipeline flowchart
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
- `table_phase5_summary.csv` (every Phase 5 gene x cohort x class headline
  test), `table_phase5_cohort_meta.csv` (N + excluded alleles per pair),
  `table_phase5_hla_c_candidate_cross_cohort.csv` (Finding-4 candidate
  pooled across BRCA/CHOL/HNSC/KIRC) + `fig_phase5_matched_normal.pdf`

Always write new generated files into `outputs/` (or `presentation/` for
slide-ready composites and `docs/` for the artifact note) — never into `data/`.
