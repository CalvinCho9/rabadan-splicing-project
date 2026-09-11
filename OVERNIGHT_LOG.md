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
