# Figure audit — HLA splicing talk

Audit performed by reading the plotting code in `scripts/analyze_brca.py`,
`scripts/additional_analyses.py`, `scripts/phase2_bc_headline.py`,
`scripts/phase5_pancancer_matched_normal.py`, `scripts/make_figures_overview.py`, and by
rendering page 1 of every candidate figure to PNG and reading the title text
directly (not just tracing code).

Legend: **A** = figure is broken (mislabeled / overcrowded / wrong chart
type) → fix. **B** = result is genuinely null (few/zero events, nothing
crosses significance) → correct emptiness, improve legibility only, do not
fill or drop.

| filename | data gene | what it shows | why it looks empty / suspect | category | proposed fix | off-limits change it would need |
|---|---|---|---|---|---|---|
| fig_cohort_overview | HLA-A | total samples vs AS-event carriers, tumor/normal | Doesn't actually look empty — bars + %-carrier annotations render correctly, title correctly says "HLA-A". | none | no fix needed | — |
| fig_classes | HLA-A | headline exon5_skip/intron5_retention carrier % tumor vs normal, Wilson CI, Fisher p | Renders correctly, title correct. Small bars + wide CIs are honest (small effect, n=935 tumor/106 normal, few carriers) — not a bug. | none | no fix needed | — |
| fig_per_allele | HLA-A | 18-panel grid, one panel per allele (n≥20), exon5_skip/intron5_retention Wilson CI | Title correct. Most of 18 panels show 0 carriers on one/both sides (thin CI whisker only, no colored bar) — that's real (most alleles have literally 0 events in one arm), but 18 panels is unreadable projected on a slide. | A (readability; underlying panels are correct nulls) | Keep full grid as-is (appendix). Add `fig_per_allele_talk.pdf`: subset to largest-n alleles (hla_a_02_01, hla_a_03_01, hla_a_01_01, hla_a_24_02) ∪ alleles with ≥1 nonzero carrier bar (hla_a_01_01, hla_a_03_01, hla_a_23_01, hla_a_24_02, hla_a_30_01) = 6 alleles, 2×3 grid. | — |
| fig_volcano | HLA-A | log2FC vs -log10(Fisher p) for every AS_pattern with ≥3 carriers | Title correct. 6 patterns tested; **0 pass BH-FDR q<0.05** (`table_volcano.csv`: Tested=6, FDR_q<0.05 count=0), even though one point sits above the nominal (uncorrected) p=0.05 dashed line. The plot currently draws that line but never labels it as *nominal, not FDR* and never states the tested/significant counts — so a reader can misread the one point above the dashed line as "significant." | B (genuinely null after correction) | Label the dashed line explicitly as "nominal p=0.05 (uncorrected)"; annotate "N patterns tested (≥3 carriers): 6, 0 significant after BH-FDR q<0.05" pulled from the already-computed `Tested`/`FDR_q` columns in the same dataframe passed to the plotting function — no new computation. | — |
| fig_event_burden | HLA-A | boxplot of AS events/sample, tumor vs normal, Mann-Whitney U | **Wrong chart type.** >82% of samples in both arms have exactly 0 events (max observed = 2), so a boxplot collapses to a flat line at 0 with fliers — real distribution, unreadable presentation. | A (chart-type bug, not a data problem) | Replace with a grouped bar chart of % of samples at each observed event count (0 / 1 / 2) by tumor/normal — same `N_events` numbers, same Mann-Whitney p kept in the title, y-axis starts at 0. | — |
| fig_heatmap_tumor | HLA-A | samples × alleles, tumor, colored by AS_pattern | Title correct, legible (118 tumor samples with ≥1 event × alleles). Renders fine. | none | no fix needed | — |
| fig_heatmap_normal | HLA-A | samples × alleles, normal, colored by AS_pattern | Title correct, legible (19 normal samples). Renders fine. | none | no fix needed | — |
| fig_threshold_sensitivity | HLA-A | intron5_retention OR vs 3 swept thresholds | Renders correctly, conclusion visibly stable (all CIs cross OR=1, direction consistent) — this is itself an honest, informative null/stability result, not a bug. | none | no fix needed | — |
| fig_allele_frequency_spectrum | HLA-A | top-20 allele carrier counts, tumor vs normal | Renders correctly, legible. | none | no fix needed | — |
| fig_classes_HLA-C | HLA-C | headline classes, HLA-C data | **Title says "HLA-A alternative splicing: tumor vs normal"** while plotting HLA-C's headline table. Cause: `analyze_brca.plot_classes()` hardcodes the string `"HLA-A"` in `ax.set_title(...)`; `scripts/phase2_bc_headline.py` calls `ab.plot_classes(headline, fname(...))` without any gene parameter. | A (mislabel) | Parameterize `gene_label` in `plot_classes()`; pass `"HLA-C"` from `scripts/phase2_bc_headline.py`. | — |
| fig_per_allele_HLA-C | HLA-C | 19-panel grid, per-allele HLA-C | **Title says "HLA-A exon5_skip / intron5_retention by allele"** over HLA-C panels (verified: allele labels are `hla_c_*`). Same hardcoded-string cause in `plot_per_allele()`'s `fig.suptitle(...)`. Also 19 panels, unreadable on a slide (same readability issue as the HLA-A grid). | A (mislabel + readability) | Parameterize `gene_label`; pass `"HLA-C"`. Add `fig_per_allele_HLA-C_talk.pdf` restricted to `hla_c_03_03`, `hla_c_03_04` (the Finding-4 real-variable-splicing candidate — the only story-carrying HLA-C panels; everything else is a 0-carrier null, per `table_per_allele_HLA-C.csv`). | — |
| fig_volcano_HLA-C | HLA-C | volcano, HLA-C patterns ≥3 carriers | **Title says "HLA-A AS_pattern recurrence, tumor vs normal"** over HLA-C's 4 tested patterns. Same hardcoded string in `plot_volcano()`. Additionally: 4 patterns tested, **0 pass BH-FDR q<0.05** (`table_volcano_HLA-C.csv`) — genuinely null at the corrected level, same as HLA-A's volcano. | A (mislabel) + B (genuinely null after correction) | Parameterize `gene_label` → `"HLA-C"`. Same nominal-vs-FDR line label + tested/significant count annotation as fig_volcano. | — |
| fig_phase5_matched_normal | pan-cancer, HLA-A + HLA-C (HLA-B computed but excluded from the plot) | grouped bars, % samples carrying event, per gene×cohort×class, matched tumor-vs-normal | Not mislabeled (built from per-gene/per-cohort labels, no hardcoded gene string) and not sparse in the sense of missing data — but **0 of the 10 presented (gene, cohort, class) tests have an OR 95% CI excluding 1** (`table_phase5_summary.csv`, presented rows only: HLA-A×{CHOL,HNSC}×2 classes + HLA-C×{CHOL,HNSC,KIRC}×2 classes = 4+6=10; HLA-B's 6 rows are quarantined out), i.e. it is a genuinely null pan-cancer result, and the figure gave no explicit tested/significant count to tell a viewer that at a glance. | B (genuinely null) | Keep the figure exactly as computed. Add a summary annotation ("10 gene×cohort×class tests shown, 0 with OR 95% CI excluding 1 (nominal, not FDR-corrected)") pulled from the already-computed `Quarantined`/`OR_CI_low`/`OR_CI_high` columns of the `summary` dataframe already in memory in `plot_summary()` — identical arithmetic to what `write_report()` already does elsewhere in the same script, not a new computation. | — |
| UNVALIDATED_HLA_B_fig_classes | HLA-B (quarantined) | headline classes, HLA-B data | **Title says "HLA-A alternative splicing: tumor vs normal"** over HLA-B's headline bars. Same `plot_classes()` hardcoded-string cause. Bars are near-zero because HLA-B has only 2 candidate AS events in the entire 935-tumor/106-normal cohort (Phase 1 finding) — that emptiness is correct and is exactly why the gene is quarantined, not a plotting bug. | A (mislabel only; the near-zero content itself is category B and correctly quarantined) | Parameterize `gene_label` → `"HLA-B"`. Do not un-quarantine; do not touch `UNVALIDATED_HLA_B_README.md`'s quarantine language. | — |
| UNVALIDATED_HLA_B_fig_per_allele | HLA-B (quarantined) | 24-panel grid, per-allele HLA-B | **Title says "HLA-A exon5_skip / intron5_retention by allele"** over `hla_b_*` panels. Same cause. 24 panels, unreadable; only 2 alleles (`hla_b_15_01`, `hla_b_27_05`) have any nonzero carrier bar at all (`UNVALIDATED_HLA_B_table_per_allele.csv`). | A (mislabel + readability) | Parameterize `gene_label` → `"HLA-B"`. Add `UNVALIDATED_HLA_B_fig_per_allele_talk.pdf` (still quarantine-prefixed, never referenced from any presented document) restricted to the 4 largest-n alleles ∪ the 2 alleles with a visible bar = 6 panels, for audit convenience only. | — |

## Phase 1 summary

- **7 figures** need a title fix (hardcoded `"HLA-A"` string bleeding into
  HLA-B/HLA-C plots): `fig_classes_HLA-C`, `fig_per_allele_HLA-C`,
  `fig_volcano_HLA-C`, `UNVALIDATED_HLA_B_fig_classes`,
  `UNVALIDATED_HLA_B_fig_per_allele`. (`fig_volcano` for HLA-B is never
  generated — 0 patterns reach the ≥3-carrier floor — so there is no
  mislabeled HLA-B volcano to fix.)
- **1 figure** needs a chart-type fix: `fig_event_burden` (boxplot → grouped
  bar of event-count distribution).
- **3 per-allele grids** get an additional `*_talk.pdf` slide-sized subset
  (full grid kept untouched as appendix): `fig_per_allele`,
  `fig_per_allele_HLA-C`, `UNVALIDATED_HLA_B_fig_per_allele`.
- **3 figures** are genuinely null results that get legibility-only
  annotations (no numbers changed): `fig_volcano`, `fig_volcano_HLA-C`,
  `fig_phase5_matched_normal`.
- **8 figures** need no change: `fig_cohort_overview`, `fig_classes`,
  `fig_heatmap_tumor`, `fig_heatmap_normal`, `fig_threshold_sensitivity`,
  `fig_allele_frequency_spectrum` (all HLA-A, all already correct), plus the
  full-grid appendix versions of the three per-allele figures once retitled.
- **0 figures** require an off-limits change (no threshold, parameter, or
  statistic needs to move to fix any of the above).

---

## Phase 2/3 — what changed, verified before/after

All fixes were presentation-only: title text, chart type, legend placement,
figure width, and added annotations pulled from already-computed columns.
**No threshold, parameter, statistic, or plotted value was changed.** Every
pipeline (`scripts/analyze_brca.py`, `scripts/phase2_bc_headline.py`,
`scripts/phase5_pancancer_matched_normal.py`, `scripts/make_figures_overview.py`) was
re-run end to end after the edits; all regenerated cleanly with the same
printed numbers as before (verified by diffing console `Fisher p` / `OR` /
carrier-count lines against the pre-fix run logs). Every corrected PDF's
page 1 was then rasterized (QuickLook `qlmanage -t`, not code tracing) and
the title text read directly off the image.

### Title fixes (`gene_label` parameterized in `analyze_brca.plot_classes()` /
`plot_per_allele()` / `plot_volcano()`; every caller in `scripts/analyze_brca.py`'s
`main()` and `scripts/phase2_bc_headline.py`'s `run_gene()` now passes its own gene)

| figure | before | after |
|---|---|---|
| `fig_classes_HLA-C.pdf` | "HLA-A alternative splicing: tumor vs normal" | "**HLA-C** alternative splicing: tumor vs normal" |
| `fig_per_allele_HLA-C.pdf` | "HLA-A exon5_skip / intron5_retention by allele (95% Wilson CI)" | "**HLA-C** exon5_skip / intron5_retention by allele (95% Wilson CI)" |
| `fig_volcano_HLA-C.pdf` | "HLA-A AS_pattern recurrence, tumor vs normal" | "**HLA-C** AS_pattern recurrence, tumor vs normal" |
| `UNVALIDATED_HLA_B_fig_classes.pdf` | "HLA-A alternative splicing: tumor vs normal" | "**HLA-B** alternative splicing: tumor vs normal" |
| `UNVALIDATED_HLA_B_fig_per_allele.pdf` | "HLA-A exon5_skip / intron5_retention by allele (95% Wilson CI)" | "**HLA-B** exon5_skip / intron5_retention by allele (95% Wilson CI)" |

`fig_classes.pdf` / `fig_per_allele.pdf` / `fig_volcano.pdf` (the actual
HLA-A figures) were already correctly titled and are unchanged by this fix
(gene_label="HLA-A" passed explicitly, same string as before).

### Chart-type fix

`fig_event_burden.pdf`: boxplot (collapsed to a flat line at 0, >80% of
samples zero-inflated) → grouped bar of % of samples at each observed event
count (0/1/2), same `N_events` numbers, Mann-Whitney U p=0.159 kept
unchanged in the title, y-axis 0–100%. Cross-checked: the new figure's
">=1 event: tumor 13%, normal 18%" annotation matches `fig_cohort_overview`'s
independently-computed `Pct_carriers` exactly.

### New "talk" figures (full grids untouched, kept as appendix)

- `fig_per_allele_talk.pdf` (HLA-A, 6 panels: hla_a_02_01, hla_a_03_01,
  hla_a_01_01, hla_a_24_02, hla_a_23_01, hla_a_30_01)
- `fig_per_allele_talk_HLA-C.pdf` (HLA-C, 2 panels: hla_c_03_03, hla_c_03_04
  — the Finding-4 real-variable-splicing candidate)
- `UNVALIDATED_HLA_B_fig_per_allele_talk.pdf` (HLA-B, 6 panels — audit
  convenience only, still quarantine-prefixed, not referenced from any
  presented document)

While building these, a fig-level legend collided with the wrapped 2-line
title on narrow (≤2-column) grids specifically — fixed generally in
`plot_per_allele()` by anchoring the legend inside the last subplot's own
corner instead of at the figure level whenever `ncols < 3`; the ≥3-column
full grids keep the original figure-level legend placement unchanged.

### Null-legibility annotations (no numbers changed, only surfaced)

- `fig_volcano.pdf` / `fig_volcano_HLA-C.pdf`: the existing dashed
  significance line is now labeled "nominal p=0.05 (uncorrected)", and the
  title gains a line reading "N patterns tested, M significant after BH-FDR
  q<0.05" — pulled from the `Tested`/`FDR_q` columns already present in the
  `volcano` dataframe (computed by `detailed_table()`'s existing BH-FDR
  correction; not recomputed). HLA-A: 6 tested, 0 significant. HLA-C: 4
  tested, 0 significant.
- `fig_phase5_matched_normal.pdf`: suptitle gains "10 gene×cohort×class
  tests shown, 0 with OR 95% CI excluding 1 (nominal, not FDR-corrected)" —
  same `OR_CI_low`/`OR_CI_high` arithmetic `write_report()` already runs on
  the same `summary` dataframe elsewhere in the script, just surfaced on the
  figure too.

### Deliberately NOT changed

- `fig_cohort_overview`, `fig_classes`, `fig_heatmap_tumor`,
  `fig_heatmap_normal`, `fig_threshold_sensitivity`,
  `fig_allele_frequency_spectrum` — already correct, left untouched.
- The full per-allele grids (`fig_per_allele.pdf`, `fig_per_allele_HLA-C.pdf`,
  `UNVALIDATED_HLA_B_fig_per_allele.pdf`) — kept complete (no panel
  dropped), only retitled; the new `*_talk.pdf` files are additions, not
  replacements.
- HLA-B remains fully quarantined: not un-quarantined, no
  `UNVALIDATED_HLA_B_*` file was moved, renamed out of the `UNVALIDATED_`
  prefix, or referenced from `presentation/figures_overview.pdf` / `docs/RESULTS_SUMMARY.md` /
  `presentation/PRESENTATION_NOTES.md`. `UNVALIDATED_HLA_B_README.md`'s quarantine
  language is untouched.
- No p-value, odds ratio, confidence interval, carrier count, threshold
  (`local_rate_threshold`, `exon_skip_min_jump`,
  `intron_retention_max_shift`, `min_carriers`, `min_allele_samples`,
  `MIN_SIDE_N`, `LOW_N_FLAG`, `HIGH_PENETRANCE`, `LOW_PENETRANCE`,
  `DOMINANT_SHARE`), exclusion list (`exclude_alleles`,
  `KNOWN_ARTIFACT_ALLELES`), or the `wilson_ci` / `odds_ratio_ci` / `bh_fdr`
  functions were touched. All null results (0 events crossing FDR
  significance anywhere in the project) remain null in every regenerated
  figure — verified by diffing the pre-fix and post-fix console logs'
  `Fisher p` / `Odds ratio` / carrier-count lines, which are identical.

### Figures that would have needed an off-limits change

None. Every figure in the audit table was fixable with presentation-only
changes; no figure's emptiness required moving a threshold or filter to
"fix."
