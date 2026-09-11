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
