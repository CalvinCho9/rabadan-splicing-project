# Final review — HLA-A/B/C, per gene

Written at the end of the overnight B/C run (Phases 1-4). One paragraph per
gene, explicitly tagged with what kind of result each claim is: **real**
(a finding worth stating with confidence), **null** (a proper test that
found no effect), **artifact** (a signal that turned out not to be
biological), **underpowered** (can't tell either way — not the same as
null), or **quarantined** (not trusted enough to present at all). Full
detail and every judgment call behind these calls: `OVERNIGHT_LOG.md`.
Presentation framing: `PRESENTATION_NOTES.md`. Full narrative:
`RESULTS_SUMMARY.md`.

---

## HLA-A

- **Pipeline generalization: REAL.** The A pipeline ran unmodified on a
  second, non-BRCA multi-cohort file (`data/pan_HLA-A_SJ_combined.tsv`) and
  correctly loaded/bucketed/classified events — a genuine capability
  demonstration, not a biological claim.
- **BRCA tumor-vs-normal (exon5_skip, intron5_retention): UNDERPOWERED
  NULL, not a contradiction of the paper.** exon5_skip essentially
  unobserved (1 tumor case, 9.8% achieved power). intron5_retention trends
  opposite the paper's predicted enrichment (4.2% tumor vs 7.5% normal,
  OR 0.53, CI 0.24–1.17) but at only 29.4% achieved power — the CI crosses
  1, so this is honestly "can't tell," not "paper is wrong." Stable across
  a 12-setting threshold sensitivity sweep (direction and rough OR range
  don't move). 0/6 patterns survive FDR.
- **HLA-A\*11:01's recurrent event: ARTIFACT (constitutive, not
  tumor-associated) — and this round confirmed it's not even
  BRCA-specific.** 96.7% of BRCA tumor carriers / 100% of normal carriers
  show the same junction, annotated junction stays dominant in reads
  (rules out a simple coordinate swap). New this round (Phase 3,
  descriptive): the same junction recurs at 85–100% penetrance in 3 of 4
  other pan-cancer cohorts with A*11:01 carriers (NORMAL pool, TCGA-CHOL,
  TCGA-LGG; TCGA-HNSC is a small-n outlier). Excluded from all headline/
  per-allele comparisons; the strongest single non-null finding in this
  entire project, but it's a reference-transcript/annotation finding, not
  a cancer-splicing finding.

## HLA-B

- **Classifier validation: FAILED — QUARANTINED, all results.** Phase 1
  tested both the naive (as-used-for-A) and a strand-aware
  (reverse-strand-mirrored) classifier; neither could be validated because
  the gene produces almost no detectable alternative splicing in this
  cohort at all — only 2 candidate AS events across 935 tumor + 106 normal
  samples (max local-junction read-share 12.6%, barely above the 10%
  detection floor, vs 52–60% for A/C). This is a genuine, checked
  data-driven ceiling (not a script bug — verified against the raw
  Local_rate distribution), but its cause (real biology vs. a
  reference/alignment quality issue specific to B) is **not determined**.
- **Tumor-vs-normal, per-allele, pan-cancer: all computed, all
  QUARANTINED.** Every number exists under `outputs/UNVALIDATED_HLA_B_*`
  for audit, and none of it should be read as a result — headline classes
  show essentially 0-1 total carriers either way, which is itself a symptom
  of "not enough data," not evidence of anything biological.
- **What would un-quarantine this gene:** (1) confirm with the
  alignment/reference-construction pipeline (outside this repo) how HLA-B's
  personalized references are built, specifically whether it matches
  HLA-C's (which Phase 1 found behaves as already transcript-oriented); (2)
  either find more HLA-B signal at a lower/adjusted detection threshold
  with its own justification, or confirm the low yield is real biology.

## HLA-C

- **Classifier validation: PASSED, evidence-based.** Phase 1 directly
  tested the strand-flip hypothesis (motivated by HLA-C being ~94.7%
  Strand==2/minus) and found it was WRONG: flipping destroyed a coherent,
  geometrically sensible signal that the naive (unflipped, same-as-A)
  classifier produces cleanly — including a real, non-constitutive
  recurrent allele (`hla_c_03_04`). Naive classifier validated as-is.
- **BRCA tumor-vs-normal: UNDERPOWERED NULL, same honest treatment as A.**
  exon5_skip 0.1% tumor vs 0.9% normal (OR 0.11, CI 0.01–1.81);
  intron5_retention 3.1% tumor vs 0.9% normal (OR 3.37, CI 0.45–25.0). Both
  CIs cross 1. intron5_retention's *direction* happens to match the paper's
  prediction (unlike A), but rests on a single normal-side carrier — a
  fragile point estimate, explicitly flagged as not a replication. 0/4
  tested patterns survive FDR.
- **Constitutive-artifact check: NULL (no artifact found) — this is
  itself informative.** No HLA-C allele shows an A\*11:01-style
  near-universal pattern, in BRCA or pan-cancer (max cross-cohort
  penetrance 50%, on n=2 — noise). `hla_c_03_04`'s top pattern looks like
  genuine, biologically variable splicing (5–20% penetrance range,
  cohort-dependent) rather than an artifact.
- **Pan-cancer: DESCRIPTIVE only, as instructed.** No tumor-vs-normal claim
  made; per-cohort penetrance reported as context.

---

## TOP THINGS TO REVIEW

Ranked by how likely each is to be wrong or need a second look, most
important first:

1. **The core Phase 1 empirical result is counter-intuitive and rests on
   inference, not ground truth.** The strand-flip hypothesis was rejected
   for HLA-C because flipping destroyed signal and naive classification
   produced a coherent one — but *why* the naive (A-style) geometry is
   correct for a majority-minus-strand gene was inferred ("the
   personalized reference contigs were probably already
   reverse-complemented at construction time, so First_base/Last_base are
   already transcript-oriented") rather than confirmed against the actual
   alignment/reference-construction pipeline, which is outside this repo.
   This inference is the single largest assumption the whole B/C extension
   rests on. **Recommend:** a 10-minute conversation with whoever built the
   personalized HLA-B/C references to confirm or refute this before citing
   HLA-C's numbers anywhere more permanent than a lab meeting.
2. **HLA-B's low event yield is unexplained.** Real biology vs. a technical
   issue specific to HLA-B was explicitly left undetermined (Phase 1,
   Phase 2 quarantine note). If it's technical, HLA-B might have real,
   currently-invisible signal; if it's biological, that's itself a
   result worth investigating and presenting properly. Don't let this sit
   unresolved past this round.
3. **The `BRCA_NORMAL`-mislabeled-as-tumor issue in the pan-cancer B/C
   files was worked around (excluded), not fixed.** If anyone else reuses
   `pan_HLA-B_SJ_combined.tsv` / `pan_HLA-C_SJ_combined.tsv` without
   reading this project's notes, they will silently double-count BRCA
   normals as tumor samples. Worth flagging to whoever maintains that
   export.
4. **A bug was caught and fixed mid-run (Phase 2)** where the multi-gene
   figure-generation loop would have silently deleted/mislabeled HLA-A's
   own canonical `outputs/fig_per_allele.pdf` / `outputs/fig_volcano.pdf`.
   Verified byte-identical to the pre-bug committed version after the fix,
   and re-verified after every subsequent phase — but worth a human glance
   at `outputs/fig_per_allele.pdf` and `outputs/fig_volcano.pdf` (the
   plain, HLA-A-only ones, still referenced by `RESULTS_SUMMARY.md`
   Findings 1-3) to confirm they look like the familiar single-gene A
   figures and not anything B/C-shaped.
5. **Every tumor-vs-normal number in this project (A and C alike) is
   underpowered, none significant after correction.** Easy to accidentally
   overstate in a live presentation by dropping the "not significant"
   caveat under time pressure — Findings 2 and 4 in `RESULTS_SUMMARY.md`
   and bullets 3/4/8 in `PRESENTATION_NOTES.md` are written to make the
   caveat hard to drop, but worth a deliberate reminder before Friday.
6. **HLA-C intron5_retention's paper-matching direction rests on n=1**
   (a single normal-side carrier, 1/106). It would take only one different
   sample to flip this point estimate's direction entirely — presented
   with that caveat, but it's the kind of number that's tempting to
   overstate live because it "agrees with the paper."
