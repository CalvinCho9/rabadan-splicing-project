# Overnight run log — HLA-B/C extension + pan-cancer descriptive pass

Run started 2026-09-11, unattended, per instructions: run Phases 1-4 start to
finish, document everything, don't stop to ask. This file is written
chronologically as each phase runs; the **TOP THINGS TO REVIEW** section at
the bottom is the fastest way to catch anything that might be wrong.

Context carried in from the existing HLA-A work (not redone, see
`RESULTS_SUMMARY.md` / `PRESENTATION_NOTES.md` / `ANNOTATION_ARTIFACT_NOTE.md`):
pipeline generalizes; BRCA HLA-A tumor-vs-normal is a null, underpowered
result; HLA-A\*11:01's recurrent event is a constitutive artifact, excluded
from headline numbers.

---

## Pre-flight: data inventory

New files present in `data/` beyond what HLA-A used:
- `export_BRCA_{tumor,normal}_HLA-B.tsv`, same for HLA-C — same 13-col
  schema, no header, confirmed identical column semantics to the HLA-A files.
- `pan_HLA-B_SJ_combined.tsv`, `pan_HLA-C_SJ_combined.tsv` — pan-cancer,
  multi-cohort, like the existing `pan_HLA-A_SJ_combined.tsv`.
- `data/pancancer_status.txt` — a per-cohort sample-count table (2 numeric
  columns per TCGA cohort, presumably tumor/normal counts feeding the
  pan-cancer export). Not machine-parsed by any script below; included here
  as background only, not verified against the pan files row-for-row.

Row counts (`wc -l`, headerless so this equals sample-row count):
| file | rows |
|---|---|
| export_BRCA_tumor_HLA-A.tsv | 28,434 |
| export_BRCA_tumor_HLA-B.tsv | 31,484 |
| export_BRCA_tumor_HLA-C.tsv | 27,261 |
| export_BRCA_normal_HLA-A.tsv | 2,868 |
| export_BRCA_normal_HLA-B.tsv | 3,101 |
| export_BRCA_normal_HLA-C.tsv | 2,618 |
| pan_HLA-A_SJ_combined.tsv | 9,738 |
| pan_HLA-B_SJ_combined.tsv | 24,808 |
| pan_HLA-C_SJ_combined.tsv | 21,263 |

Distinct barcodes: HLA-B tumor 935 / normal 106 (matches HLA-A's 935/106
exactly — same cohort, as expected). **HLA-C tumor is 933, not 935** — 2
tumor barcodes present for A/B are absent from the HLA-C tumor export.
Normal is 106/106, matching. **Assumption:** treated as a real (if minor)
per-gene sample-availability difference in the source export, not
investigated further (2/935 samples, does not change any conclusion) — noted
here rather than silently reconciled.

**Judgment call — flagged, not acted on:** `pan_HLA-B_SJ_combined.tsv` and
`pan_HLA-C_SJ_combined.tsv` both contain a `Cohort == "BRCA_NORMAL"` block
whose `Type` column reads **"tumor"** — 3,101 rows for HLA-B (exactly equal
to `wc -l export_BRCA_normal_HLA-B.tsv`) and 2,618 rows for HLA-C (exactly
equal to `wc -l export_BRCA_normal_HLA-C.tsv`). The row counts match the
*normal* exports exactly, which strongly suggests this is the BRCA normal
cohort with `Type` mislabeled as "tumor" in the pan-cancer combined file —
an upstream export bug, not something in these scripts. (`pan_HLA-A_SJ_combined.tsv`
has no BRCA rows at all, so this is new to the B/C pan files and not
something the original HLA-A analysis had to deal with.) Per the "never
fabricate, note and move on" rule, **`Type` was NOT silently corrected**.
Phase 3 (pan-cancer, descriptive-only, no tumor-vs-normal claims) **excludes
the `BRCA_NORMAL` cohort rows entirely** rather than either trusting a
label that looks wrong or guessing a fix — see Phase 3 section below. This
does not affect Phase 2 (which uses the clean, correctly-labeled
`export_BRCA_{tumor,normal}_HLA-{B,C}.tsv` files, not the pan files).

---

## Phase 1 — HLA-B/C strand validation

Script: `phase1_bc_strand_check.py` (new). Shared helpers factored into
`hla_common.py` (new) — imports `analyze_brca` as a library, does not modify
it, so all existing HLA-A outputs/behavior are untouched.

**What was tested.** The task brief's working hypothesis was that HLA-B/C
being reverse-strand genes might mean the A classifier's geometry (donor =
First_base, acceptor = Last_base, "downstream" = increasing coordinate — only
strictly correct for plus-strand junctions) mislabels B/C events, and would
need a strand-aware coordinate flip. Verified first that this is even a live
question: `First_base <= Last_base` holds for every row in A/B/C (0
violations — coordinates are always stored in left-to-right order in this
export, not pre-oriented per transcript), and the `Strand` column confirms
HLA-A is ~97.8% `Strand==1`, HLA-B ~97.8% and HLA-C ~94.7% `Strand==2` — so
if genomic order and transcript order disagree for minus-strand genes, B/C
are overwhelmingly the case that would matter.

Built `hla_common.classify_events_strand_aware(flip_minus=True)`: for
`Strand==2` rows only, mirror coordinates (`W1'=-W2, W2'=-W1`, same for
V1/V2) before applying the *unmodified* A geometry — this swaps donor/
acceptor roles and reverses "downstream," which is the correct transform
**if** First_base/Last_base in this export are genomic-order coordinates.
Ran classification both ways (naive/strand-blind vs strand-aware/flipped)
on both genes and compared which produced a geometrically coherent,
recurrent exon5_skip/intron5_retention signal.

**Result — the flip hypothesis was empirically wrong for both genes.**
Flipping did not help; it destroyed what signal existed:

| Gene | naive exon5/intron5 carriers | flipped carriers | verdict |
|---|---|---|---|
| HLA-B | 1 | 0 | **untrustworthy** (insufficient data either way) |
| HLA-C | 30 | 0 | **valid_as_is** (naive is correct, flip is wrong) |

For HLA-C, the naive classifier finds a coherent, A-like signal — e.g.
`hla_c_03_04`'s top pattern `1273_1711->1291_1711` cleanly fits the
intron5_retention geometry (acceptor fixed at 1711, donor shifts +18bp) and
shows real, non-constitutive penetrance (19.3% of tumor carriers, 11.1% of
normal carriers — well inside the 5-80% "real variable splicing" band, not
the ~100%-in-both-groups A\*11:01 pattern). Under the flip, 59/60 candidate
event rows relabel to "other" and the entire signal disappears. **Working
explanation** (documented as an inference, not confirmed against the
upstream alignment/reference-construction code, which is outside this
repo): the personalized per-patient reference contigs for HLA-B/C were
likely already reverse-complemented at reference-construction time (as the
task brief describes), which would make First_base/Last_base in this
*export* already transcript-oriented (same convention as A) — and the
`Strand` column plausibly just carries the gene's known genomic strand as
metadata (HLA-A sits on the + strand of chr6, HLA-B/C on the -, which is
correct human-genome biology) rather than instructing a per-row coordinate
transform. This is inferred from the data behavior, not verified against
pipeline source — flagged as a top review item below.

For HLA-B, this couldn't be tested with any confidence either way: only
**2 candidate AS events total** were found across all 935+106 samples (vs
60 for HLA-C, and roughly 200+ for HLA-A at the same 10%-local-rate
threshold). This traces to a real, data-driven ceiling, not a script bug —
checked the *maximum* `Local_rate` among any novel-junction-downstream-of-
annotated row in HLA-B: **12.6%**, barely above the 10% detection
threshold at all, vs 52.5% for A and 60.2% for C. HLA-B in this cohort
essentially doesn't produce alternative-splicing signal at the 10%
local-read-share bar the pipeline uses. Whether that reflects real biology
(HLA-B genuinely less alternatively spliced / more strongly selected here)
or a technical issue upstream (allele diversity at B spreading reads too
thin per bucket, personalized-reference quality, etc.) is **not
determined** — noted, not invented an answer.

**Constitutive-artifact check (A\*11:01-style), both genes:** ran the same
penetrance-vs-tumor/normal check used to catch A\*11:01 (this check is
AS_pattern-based, not Event_class-based, so it's independent of the
strand-classifier verdict above and trustworthy regardless). **No B or C
allele shows an A\*11:01-style near-universal (>=80% in both tumor and
normal) constitutive pattern.** HLA-C's alleles are either
`real_variable_splicing_candidate` (6 alleles, all well under 20% peak
penetrance) or `rare_private_event`/too-few-samples. HLA-B has too few
events to call anything beyond "too rare." So the quarantine rule's
concern about hidden constitutive artifacts is a non-issue for whichever
genes proceed to Phase 2 — but note this check could not run meaningfully
for HLA-B either, for the same low-N reason.

**Outputs:** `outputs/PHASE1_BC_strand_check.md` (full verdict + evidence
table), `outputs/table_phase1_top_patterns_HLA_{B,C}.csv`,
`outputs/table_phase1_penetrance_HLA_{B,C}.csv`,
`outputs/table_phase1_classifier_comparison_HLA_{B,C}.csv`.

**Quarantine rule applied from here on:**
- **HLA-C**: verdict `valid_as_is` → Phase 2 proceeds normally with the
  standard (unflipped) classifier, results go into the presented outputs
  (`RESULTS_SUMMARY.md`, `PRESENTATION_NOTES.md`, `figures_overview.pdf`)
  once Phase 4 folds them in.
- **HLA-B**: verdict `untrustworthy` (insufficient data to validate either
  scheme) → per the non-negotiable quarantine rule, Phase 2's HLA-B
  tumor-vs-normal output goes to `outputs/UNVALIDATED_HLA_B_*` only, marked
  **NOT FOR PRESENTATION — pending manual strand review**, and is excluded
  from `RESULTS_SUMMARY.md`, `PRESENTATION_NOTES.md`, and
  `figures_overview.pdf`.

Committed to git after this phase.

---

## Phase 2 — HLA-B/C tumor-vs-normal

Script: `phase2_bc_headline.py` (new). Runs the exact same headline
(Fisher + Wilson CI + OR + Haldane-Anscombe correction), per-allele
(min_allele_samples=20), and detailed/volcano (min_carriers=3, BH-FDR)
comparisons `analyze_brca.py` runs for A, pointed at HLA-B/C, with the
standard (unflipped) classifier for both genes per Phase 1's finding. Also
re-ran the constitutive-artifact (A\*11:01-style) check independently within
this script (not just trusting Phase 1's numbers) as the `--exclude-alleles`
equivalent — **0 alleles flagged for either gene**, so no allele exclusions
were applied to either gene's headline/per-allele comparison.

**Bug caught and fixed before results were trusted:** `analyze_brca.py`'s
`plot_per_allele()`/`plot_volcano()` always write to the fixed filenames
`<outdir>/fig_per_allele.pdf` / `<outdir>/fig_volcano.pdf`. First draft of
this script called them with `OUTDIR` directly (outputs/) and renamed
afterwards — which clobbered/mislabeled HLA-A's own canonical
`outputs/fig_per_allele.pdf` and `outputs/fig_volcano.pdf` (deleted them, and
in one case relabeled HLA-A's original volcano plot as
`UNVALIDATED_HLA_B_fig_volcano.pdf`, since HLA-B produces no volcano data of
its own and the rename silently picked up the stale HLA-A file sitting at
the default path). Caught by diffing outputs/ against the Phase 1 git commit
before trusting anything downstream. Fixed by routing both plot calls
through an isolated `outputs/_phase2_scratch/` directory and moving the
result into the gene-tagged filename, then verified byte-for-byte
(`md5`) that `outputs/fig_per_allele.pdf` / `outputs/fig_volcano.pdf` match
the Phase 1 commit exactly after the fix. **Flagging this prominently below
in TOP THINGS TO REVIEW** since it's exactly the kind of silent-clobber bug
that's easy to miss in an unattended run — worth a spot-check that HLA-A's
figures still look right.

### HLA-C (presented)

933 tumor / 106 normal, 60 candidate AS events (2 exon5_skip / 2
intron5_retention normal; 1 exon5_skip / 54 intron5_retention / 1 other
tumor, by row — not distinct carriers).

| Event | Tumor | Normal | Fisher p | OR (95% CI) | Direction |
|---|---|---|---|---|---|
| exon5_skip | 0.1% (1/933, CI 0.0–0.6%) | 0.9% (1/106, CI 0.2–5.2%) | 0.194 | 0.11 (0.01–1.81) | depleted in tumor |
| intron5_retention | 3.1% (29/933, CI 2.2–4.4%) | 0.9% (1/106, CI 0.2–5.2%) | 0.354 | 3.37 (0.45–24.98) | enriched in tumor |

Both CIs cross 1 — **neither result is statistically confident**, same
honesty rule as A. Note for the room: intron5_retention's *direction* here
(enriched in tumor) is the same direction the MHC Hammer paper predicts,
unlike A's BRCA result (which trended the opposite way) — but with only 1
normal carrier total, this is an extremely fragile point estimate (CI runs
0.45 to nearly 25x), not a replication. 0/4 tested per-pattern combinations
survive BH-FDR q<0.05 in the detailed table. No allele shows a
constitutive-artifact pattern, so nothing was excluded from these numbers.

**Outputs:** `outputs/table_headline_classes_HLA-C.csv`,
`outputs/fig_classes_HLA-C.pdf`, `outputs/table_per_allele_HLA-C.csv`,
`outputs/fig_per_allele_HLA-C.pdf`,
`outputs/results_BRCA_tumor_vs_normal_HLA-C.csv`,
`outputs/table_volcano_HLA-C.csv`, `outputs/fig_volcano_HLA-C.pdf`,
`outputs/table_penetrance_recheck_HLA-C.csv`.

### HLA-B (quarantined — NOT FOR PRESENTATION)

Computed for completeness/audit but **excluded from every presented
document** per the quarantine rule (Phase 1 verdict: untrustworthy,
insufficient data — only 2 candidate AS events in the whole cohort). Both
headline classes show 0-1 total carriers combined tumor+normal — there is
essentially nothing here to report either way, which is itself consistent
with "not enough signal to say anything," not a null result worth stating
confidently. All outputs prefixed `outputs/UNVALIDATED_HLA_B_*`, with
`outputs/UNVALIDATED_HLA_B_README.md` explaining the quarantine and what
would need to happen (manual strand/reference review) before these numbers
could ever be presented.

Committed to git after this phase.

---

## Phase 3 — pan-cancer descriptive pass

Script: `phase3_pancancer_descriptive.py` (new). Descriptive only, as
instructed — no tumor-vs-normal test anywhere in this phase. The `NORMAL`
cohort in these files is a separate cross-cohort normal-tissue pool, not
matched normals per cancer type, so even cohorts with both `tumor` and
`normal` rows present are not treated as a paired comparison.

**BRCA_NORMAL handling:** excluded the mislabeled `Cohort=="BRCA_NORMAL"`
block from both pan_HLA-B/C files entirely (3,101 / 2,618 rows respectively)
rather than trust or silently fix the `Type=="tumor"` label — logged in
Pre-flight above. No information is lost: BRCA is already covered correctly
by Phase 2's dedicated export files.

**Quarantine applied at the gene level**, per the "when in doubt, quarantine"
rule — HLA-B's entire Phase 3 output goes to
`outputs/UNVALIDATED_HLA_B_pancancer_*`, even though the per-cohort
penetrance check itself doesn't use the exon5_skip/intron5_retention
classifier (it's AS_pattern-based, same as Phase 1's). Decided this way for
simplicity/safety rather than carving out a classifier-independent
exception — HLA-B's data is sparse pan-cancer too (620 samples, only 2
candidate AS events total, one allele reaches the 10-sample bar and shows
no real signal), so there's nothing being lost by quarantining it wholesale.

**Headline finding — A\*11:01's constitutive isoform generalizes beyond
BRCA.** Tested whether the specific junction established as constitutive in
BRCA (`2365_2806->2383_2806`, `ANNOTATION_ARTIFACT_NOTE.md`) recurs in other
cohorts. It does, strongly: **NORMAL pool 20/20 (100%), TCGA-CHOL 3/3
(100%), TCGA-LGG 6/7 (85.7%)** — only TCGA-HNSC is an outlier (1/3, but n=3
is too small to weigh against the other three). Pooled: 30/33 (91%) of
pan-cancer A\*11:01 carriers show the pattern. This is a genuinely useful
strengthening of the existing A\*11:01 finding — it's a fixed property of
the allele's reference transcript across tissue/cancer types, not a
BRCA-cohort quirk. TCGA-KIRC and TCGA-UCEC had 0 A\*11:01 carriers to test
(this cohort's allele frequency spectrum apparently doesn't include A*11:01
much, or samples are too few — not investigated further).

**HLA-A general cross-cohort scan:** besides A\*11:01, the broader per-allele
penetrance-by-cohort table (alleles with >=10 pan-cancer samples) found 3
(allele, cohort) combinations at >=80% penetrance — all of them are, in
fact, A\*11:01 itself in different cohorts (NORMAL/CHOL/LGG, the same 3
cohorts as above) — i.e. this scan independently rediscovers the same
signal rather than turning up a second constitutive allele. No other A
allele showed a similar pattern.

**HLA-C:** no (allele, cohort) combination reaches 80% penetrance pan-cancer
(max observed 50%, on n=2 — noise). Consistent with Phase 1/2's BRCA-only
finding that no C allele looks constitutive. `hla_c_03_04`'s intron5_retention
pattern (the same one flagged as `real_variable_splicing_candidate` in
Phase 1/2) shows modest, variable penetrance across cohorts (0-21%) — looks
like real, biologically variable splicing wherever it's seen, not an
artifact, consistent with the BRCA-only read.

**HLA-B (quarantined):** essentially no signal pan-cancer either (2
candidate events total across 620 samples) — consistent with the Phase 1
finding that this gene just doesn't produce detectable AS events at the
10%-local-rate threshold in this data, not specific to BRCA.

**Outputs:** `outputs/PHASE3_pancancer_descriptive.md`,
`outputs/table_phase3_cohort_overview_HLA-{A,C}.csv`,
`outputs/table_phase3_allele_penetrance_by_cohort_HLA-{A,C}.csv`,
`outputs/table_phase3_A11_01_cross_cohort.csv`,
`outputs/UNVALIDATED_HLA_B_pancancer_*.csv`.

Verified HLA-A's canonical outputs/ files remain untouched (git status shows
zero deletions) before committing.

Committed to git after this phase.

---

## Phase 4 — consolidation

Updated `RESULTS_SUMMARY.md` (added Finding 4: HLA-C headline test, same
honest-null treatment as HLA-A; Finding 5: A\*11:01 pan-cancer
generalization; updated Bottom line; added a gene-coverage note up top
making the HLA-B quarantine explicit before anyone reads further) and
`PRESENTATION_NOTES.md` (added bullets 8-10 for HLA-C's headline result,
the A\*11:01 pan-cancer finding, and an explicit "HLA-B: validation
attempted, inconclusive" bullet so it's addressed proactively rather than
left as a gap if asked live; updated Next Steps). Only HLA-C and the
pan-cancer A\*11:01/HLA-C descriptive results were folded in, per the
quarantine rule — no HLA-B number appears anywhere in either file.

Rewrote `make_figures_overview.py` to a two-gene, five-panel layout (A1/A2:
HLA-A and HLA-C headline classes; B: HLA-A per-allele constitutive-artifact
contrast, unchanged; C: combined HLA-A+HLA-C volcano; D: new panel, the
A*11:01 pan-cancer recurrence bar chart) and regenerated
`figures_overview.pdf`. Titled explicitly as an HLA-A/HLA-C figure with
HLA-B's absence called out in the suptitle and panel C's caption, rather
than silently omitting it. Rendered a PNG preview before finalizing to
visually confirm the layout wasn't broken (all 4 subplot groups legible,
labels not overlapping) — see the preview check before the final
`figures_overview.pdf` write in this session's tool history.

Wrote `outputs/FINAL_REVIEW.md`: one paragraph per gene, each claim tagged
real/null/artifact/underpowered/quarantined, plus its own "TOP THINGS TO
REVIEW" list (6 items, focused on what a reader of that specific file needs
to double-check).

Lightly updated `README_FOR_CLAUDE.md` (data file list, script list, running
order, write-up section) so it stays an accurate entry point for reproducing
everything — not explicitly required by the brief but it's the project's
canonical map and leaving it silently stale after adding 4 new scripts and
~2 dozen new output files seemed like the wrong kind of shortcut.

Final verification before this commit: re-diffed `outputs/fig_per_allele.pdf`
and `outputs/fig_volcano.pdf` (HLA-A's originals) against the Phase 1 git
commit — still byte-identical (`md5` match). Confirmed no `git status`
deletions anywhere under `outputs/`, `RESULTS_SUMMARY.md`, or
`PRESENTATION_NOTES.md` (only additions/new files).

Committed to git after this phase.

---

## Phase 5 — pan-cancer NORMAL pool resolved to true cohort of origin

Not overnight/unattended like Phases 1-4 — triggered by Anqi's actual reply
to the "can we pull more BRCA normals?" question asked at the end of Phase
4. Two things in her message: (1) 106 BRCA normals is already enough for a
proportion comparison, no need to pull more; (2) she attached
`gdc_sample_sheet_normal_with_hla.tsv`, a GDC sample sheet + HLA typing for
740 normal samples across 24 TCGA projects. Read literally, that file looks
like it's *about* BRCA (her message was in response to a BRCA-normals
question) — it isn't, mostly. Checked: it maps `Case ID` -> `Project ID` for
every normal sample in it, and cross-referencing against the pan-cancer SJ
files' Barcode column shows every one of the pooled `Cohort=="NORMAL"`
barcodes (145 in the HLA-A pan file, 318 in HLA-B/C's) resolves to a real
TCGA project -- 0 unresolved in both cases. That pool is exactly what
`phase3_pancancer_descriptive.py`'s docstring calls "a separate cross-cohort
normal-tissue pool ... not per-cancer-type matched normals." Resolving it
turns that into real matched normals for whichever TCGA project also has
tumor rows in the same pan file.

Checked which cohorts that actually gives a usable test for (>=3 samples
both sides after resolving): TCGA-CHOL (33 tumor/8 normal) and TCGA-HNSC
(52 tumor/5 normal) for the HLA-A pan file; TCGA-CHOL (33/8), TCGA-HNSC
(67/35), TCGA-KIRC (92/37) for HLA-B/C's (same 7-cohort file, more tumor
cohorts). Cohorts that DON'T get a new test and why: TCGA-LGG and TCGA-ACC
have tumor rows but no TCGA project of that name anywhere in the 24-project
sheet at all (no normal exists to resolve to); TCGA-THCA has 47 resolved
normals but only 1 tumor sample; TCGA-UCEC has up to 79 tumor samples but
only 1 resolved normal. Wrote `phase5_pancancer_matched_normal.py` to run
the exact same headline (Fisher/Wilson/OR) machinery Phase 2 uses for BRCA,
restricted to each of these cohorts one at a time, gated by the same Phase
1 quarantine (HLA-B -> `outputs/UNVALIDATED_HLA_B_phase5_*`).

**Near-miss, caught before writing up results:** the first version excluded
constitutive-artifact alleles by re-deriving the exclusion per cohort, same
as Phase 2 does for BRCA (`hla_common.top_event_penetrance` +
`add_verdict`, which won't flag an allele with fewer than 10-20 carriers).
TCGA-CHOL and TCGA-HNSC only have 2-3 tumor samples carrying A\*11:01 each
-- nowhere near that bar -- so the re-derivation said "nothing to exclude,"
and the resulting headline numbers for HLA-A in both cohorts were "enriched
in tumor." Before writing any of that into `RESULTS_SUMMARY.md`, opened the
underlying per-pattern detail CSV to see which allele/pattern was actually
driving each headline number (worth doing as a habit, not just when
something looks suspicious) and found: 100% of both cohorts' original
intron5_retention carriers were `hla_a_11_01`'s already-known
`2365_2806->2383_2806` artifact junction (Finding 3/5), not new biology.
Fixed by excluding A\*11:01 unconditionally in every cohort (it's already
established project-wide, not something that should need to clear a
per-cohort sample-size bar again), keeping the per-cohort re-derivation
running on top in case a *different* cohort-specific artifact shows up
(none did). Re-ran; HLA-A's TCGA-HNSC comparison now has zero qualifying
events on either side, and TCGA-CHOL's intron5_retention count dropped from
3 to the 1 real (non-artifact) carrier.

Also ran a bonus check in the same script: does HLA-C's own Finding-4
candidate (`hla_c_03_03`/`hla_c_03_04`, a real, non-artifact
intron5_retention pattern that in BRCA rested on a single fragile
normal-side carrier) show up again in the newly-resolvable cohorts? It
does, in all three (TCGA-CHOL 1/0, TCGA-HNSC 4/0, TCGA-KIRC 3/0
tumor/normal carriers) -- pooled with BRCA, 42 tumor vs 2 normal across 4
independent cohorts, directionally consistent every time. Written up
carefully as "worth a dedicated follow-up," explicitly NOT as a replication
-- no individual test reaches significance, this wasn't pre-registered, and
raw-count pooling across cohorts isn't a proper stratified test.

All 16 gene x cohort x event-class headline tests from this phase come back
null (every OR 95% CI spans 1) -- consistent with, not contradicting, the
BRCA-only story. Didn't touch `figures_overview.pdf` given an all-null
result set doesn't earn a 6th panel in an already-tuned 5-panel figure;
`outputs/fig_phase5_matched_normal.pdf` stands alone instead. Updated
`RESULTS_SUMMARY.md` (Finding 6, Bottom line), `PRESENTATION_NOTES.md`
(bullets 11-12, corrected the now-stale "more normals" next-step to reflect
Anqi's actual answer), and `README_FOR_CLAUDE.md` (data file, script #10,
key outputs). Moved `gdc_sample_sheet_normal_with_hla.tsv` into `data/`
(it arrived at the project root) to match the project's existing
data/root/outputs layout convention.

Also fixed a latent bug surfaced by this phase, not introduced by it:
`hla_common.top_event_penetrance()` crashed with a `KeyError` on any
input where zero alleles reach `min_allele_samples` with a qualifying
event (an empty `rows` list produces a columnless `pd.DataFrame`, and
`.sort_values()` on that raises rather than returning empty) -- never
triggered before because Phase 1/2 always call it on the full 935+106 BRCA
cohort, where this never happens, but several of Phase 5's small
single-cohort slices (as few as 33+8 samples) hit it immediately. Fixed by
returning an empty, correctly-columned frame in that case; Phase 1/2's
BRCA-scale behavior is unaffected (never exercises that branch).

Committed to git after this phase.

---

## TOP THINGS TO REVIEW

Ranked most-likely-to-be-wrong first. (A shorter, gene-organized version of
this same list is also in `outputs/FINAL_REVIEW.md`.)

1. **The whole HLA-C validation rests on an inference about the upstream
   alignment pipeline that this repo can't confirm.** Phase 1 empirically
   rejected the strand-flip hypothesis for HLA-C (flipping destroyed a
   coherent signal; naive classification produced one), then explained
   *why* with an inference — "the personalized reference contigs were
   probably already reverse-complemented at construction time" — that was
   never checked against the actual reference-construction code (outside
   this repo). The empirical result (naive > flipped) is solid; the
   causal story for *why* is a plausible guess. If that guess is wrong,
   HLA-C's classifier could still be behaving correctly for the wrong
   reason, or there could be a third explanation not considered. **Worth a
   5-minute check with whoever built the personalized HLA-B/C references
   before this goes beyond a lab meeting.**
2. **HLA-B's near-total absence of detectable splicing signal is
   unexplained.** Confirmed it's not a script bug (checked the raw
   Local_rate distribution directly — the ceiling really is ~12.6% vs
   52-60% for A/C), but whether that reflects real biology or a technical/
   reference-quality problem specific to HLA-B was explicitly left open.
   This is the reason HLA-B is quarantined at all — resolving it either
   validates HLA-B for a future pass or turns into its own finding.
3. **A mid-run bug (Phase 2) would have silently corrupted HLA-A's own
   canonical figures** (`outputs/fig_per_allele.pdf`, `outputs/fig_volcano.pdf`)
   by deleting or mislabeling them with B/C data, due to
   `analyze_brca.py`'s plotting functions writing to fixed filenames. Caught
   by diffing against the Phase 1 git commit rather than assuming success,
   fixed with an isolated scratch directory, and re-verified byte-identical
   after every subsequent phase. Flagging because it's exactly the kind of
   thing that's easy to miss when running unattended — worth a human glance
   at those two files to confirm they still look like the familiar
   single-gene HLA-A figures.
4. **The pan-cancer `BRCA_NORMAL`/`Type=tumor` mislabeling was worked
   around, not fixed at the source.** If this export is reused elsewhere
   without reading this project's notes, someone will silently double-count
   106-108 real BRCA normal samples as tumor. Worth a note to whoever
   maintains that export pipeline.
5. **Every tumor-vs-normal number in this entire project (HLA-A and HLA-C
   alike) is underpowered and none survive multiple-testing correction.**
   Both `RESULTS_SUMMARY.md` and `PRESENTATION_NOTES.md` were written to
   make this hard to drop accidentally in a live presentation, but it's
   worth a deliberate reminder before presenting — especially HLA-C's
   intron5_retention number, whose direction happens to match the MHC
   Hammer paper's prediction and is therefore the single most tempting
   number in the whole deck to overstate. It rests on n=1 on the normal
   side.
6. **The exclusion thresholds reused from HLA-A** (`min_allele_samples=20`,
   `min_carriers=3`, `local_rate_threshold=0.1`, the 80%/5% penetrance
   gates for the constitutive-artifact check) were applied to B/C
   unchanged, without their own gene-specific sensitivity sweep (Phase 1/2
   scope didn't include re-running `additional_analyses.py`'s
   threshold-sensitivity sweep for B/C — only HLA-A's original sweep from
   the prior session covers this). If HLA-C's headline numbers need to
   survive real scrutiny (not just a lab meeting), that sweep should be
   extended to HLA-C before citing it further.
7. **Phase 3's "3 (allele, cohort) combinations at >=80% penetrance"
   headline for HLA-A** turned out, on inspection, to just be A\*11:01
   itself in 3 different cohorts, not 3 independent constitutive alleles —
   correctly caught and written up that way in
   `outputs/PHASE3_pancancer_descriptive.md`/`RESULTS_SUMMARY.md`, but
   flagging because a less careful read of just the summary count ("3
   combinations flagged") could be misread as 3 separate findings if
   someone skims only the number.
8. **Phase 5's HLA-C bonus finding (bullet 12/Finding 6) pools raw carrier
   counts across 4 cohorts (42 tumor vs 2 normal) without a proper
   stratified test (e.g. Cochran-Mantel-Haenszel).** Deliberately presented
   as descriptive corroboration, not a pooled significance claim — same
   restraint Finding 5 used for A\*11:01's cross-cohort recurrence — but
   worth a second look before anyone is tempted to compute a naive pooled
   Fisher test on those numbers and call it significant; the four cohorts
   have different allele frequencies and ascertainment, so a naive pool
   would overstate confidence.
9. **The exclusion-derivation near-miss documented in Phase 5 above** (the
   per-cohort constitutive-artifact re-check missing A\*11:01 in small
   cohorts) is a structural gap in `hla_common.top_event_penetrance`'s
   design, not just a one-off Phase 5 bug: any *future* script that calls
   it on a small single-cohort slice without also applying the
   project-wide known-artifact list will make the same mistake. Worth a
   comment or guard in `hla_common.py` itself if this pattern gets reused
   again, rather than re-deriving the fix ad hoc each time.
