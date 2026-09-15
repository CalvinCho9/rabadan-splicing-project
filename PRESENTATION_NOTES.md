# Presentation notes — HLA-A/C splicing, TCGA-BRCA + pan-cancer

Talking points for the lab meeting. Each bullet: **claim** → supporting stat
→ honest caveat. Full detail in `RESULTS_SUMMARY.md`; slide panel in
`figures_overview.pdf`; method/background reference sheet in
`background_method.pdf` if the room needs the mechanics explained first.

**Gene coverage: A and C only.** HLA-B's classifier could not be validated
(Phase 1: only 2 candidate AS events in the whole BRCA cohort — nowhere near
enough to test) and is quarantined out of every number below. If asked
about HLA-B, the honest answer is "not enough signal in this data yet to
say anything, still under review" — see `outputs/FINAL_REVIEW.md` /
`OVERNIGHT_LOG.md` for the full reasoning, do not improvise a number.

1. **The analysis generalizes past BRCA, by design, not by accident.**
   The same code (unmodified) ran against a separate multi-cohort file and
   correctly loaded/bucketed/classified 46 AS events across 301 samples, 6
   cohorts.
   *Caveat:* that file doesn't include BRCA tumor data, so this demonstrates
   pipeline reuse, not a BRCA pan-cancer result — no biological claim is
   being made about those cohorts.

2. **BRCA has a modest AS-event carrier rate overall, across a diverse
   allele pool.** 13.0% of tumor samples (122/935) and 17.9% of normal
   samples (19/106) carry at least one qualifying novel-junction event; 39
   distinct HLA-A alleles are present, none overwhelmingly dominant (top
   allele `hla_a_02_01` is ~45% of samples).
   *Caveat:* "carrier of any event" mixes the two target classes with
   everything else that passes the >10%-local-reads filter — cohort context,
   not itself a tumor-vs-normal test.

3. **exon5_skip is essentially unobserved here.** 0.1% tumor (1/935) vs 0.0%
   normal (0/106), Fisher p=1, OR 0.34 (95% CI 0.01–8.45); achieved power to
   detect this was 9.8%.
   *Caveat:* one tumor case total. Not evidence for or against the paper's
   predicted depletion — there's not enough data to test it, and the power
   number quantifies exactly how little.

4. **intron5_retention trends the opposite direction from the paper's
   prediction, but the test was underpowered to say much either way.** 4.2%
   tumor (39/935) vs 7.5% normal (8/106) after excluding the A\*11:01
   artifact, Fisher p=0.13, OR 0.53 (95% CI 0.24–1.17), achieved power 29.4%.
   At 80% power we'd have needed a tumor rate outside roughly [1.8%, 16.8%]
   to detect anything at all.
   *Caveat:* don't present this as "we didn't reproduce the paper" — present
   it as "underpowered to test it here," with a number attached.

5. **The conclusion doesn't hinge on where we drew the two geometry
   thresholds or the 10% local-read cutoff.** A 12-setting one-at-a-time
   sweep (local-rate 5–20%, exon-skip jump 50–150bp, intron-retention shift
   20–70bp) keeps the direction "depleted in tumor" throughout, OR staying
   in 0.41–0.60.
   *Caveat:* one of the 12 settings crosses nominal p<0.05 — expected noise
   given 12 tests on a small effect, not a robust independent finding, and
   nothing here survives multiple-testing correction (0/6 tested patterns
   pass FDR q<0.05 in the main per-pattern table).

6. **HLA-A\*11:01 carries a near-universal minor splice isoform, unrelated
   to tumor status — caught because we had matched normals.** 96.7% of
   A\*11:01 tumor carriers (118/122) and 100% of normal carriers (19/19)
   show the same junction; the read-count share (annotated junction stays
   dominant, ~72% of reads, in both groups) rules out a simple reference
   coordinate swap and points to a real but constitutive minor isoform.
   *Caveat:* this is a data-quality finding about one allele's reference
   transcript, not a splicing-biology finding — it's excluded from the
   headline numbers above (`--exclude-alleles`), not folded into them.

7. **Against MHC Hammer specifically:** we reproduce their *framework*
   (exon5_skip / intron5_retention as the two events worth testing, and the
   value of testing per-allele rather than only pooled) but our BRCA cohort
   is too small on the normal side to reproduce their reported effect sizes
   in either direction — quantified above as 10–29% achieved power. Where
   our numbers most visibly diverge from a naive read of the paper is
   intron5_retention's direction (depleted vs. predicted enriched), flagged
   honestly as noise-consistent given the power ceiling, not as a
   disagreement with their conclusion.

8. **HLA-C extends to a second gene, with the same honest-null result.**
   933 tumor / 106 normal, classifier independently validated (see below).
   exon5_skip 0.1% tumor vs 0.9% normal (Fisher p=0.19, OR 0.11, CI
   0.01–1.81); intron5_retention 3.1% tumor vs 0.9% normal (Fisher p=0.35,
   OR 3.37, CI 0.45–25.0). *Caveat:* both CIs cross 1 — not significant.
   intron5_retention's direction here happens to match the paper's
   prediction, but it rests on a single normal-side carrier — flag as
   fragile, not a replication, if asked.

9. **A\*11:01's constitutive isoform isn't a BRCA quirk — it recurs
   pan-cancer.** Same junction, tested in every other cohort with A*11:01
   carriers: NORMAL pool 20/20 (100%), TCGA-CHOL 3/3 (100%), TCGA-LGG 6/7
   (85.7%); TCGA-HNSC 1/3 is the outlier but n=3 is too small to weigh
   against the rest. Pooled: 30/33 (91%). *Caveat:* descriptive only — the
   pan-cancer `NORMAL` cohort is a separate cross-cohort pool, not matched
   normals per cancer type, so this strengthens the "constitutive, not
   tumor-associated" read but isn't itself a new significance test.

10. **HLA-B: validation attempted, inconclusive — quarantined, not
    presented.** Ran the same strand/geometry validation used for C; HLA-B
    simply doesn't produce enough alternative-splicing signal in this
    cohort to test it (2 candidate events total, vs 60 for C and ~200+ for
    A at the same threshold — max observed junction read-share was 12.6%,
    barely above the 10% detection floor, vs 52–60% for A/C). *Caveat:*
    this is presented as an open question, not a finding — could be real
    biology (HLA-B genuinely less alternatively spliced here) or a
    technical/reference-quality issue specific to B; not determined. If the
    room asks for HLA-B numbers, the honest answer is "not yet — see
    `outputs/UNVALIDATED_HLA_B_*` for what exists, but don't present it as
    a result."

11. **Anqi's fix for "we need more normals" wasn't more normals — it was
    knowing where the ones we have come from.** She confirmed 106 BRCA
    normals is already enough for a proportion comparison, and sent
    `gdc_sample_sheet_normal_with_hla.tsv` (740 normal samples, 24 TCGA
    projects) specifically so the pan-cancer files' pooled `NORMAL` samples
    could be traced back to their real cancer type. Every one resolves
    (145/145 for HLA-A, 318/318 for HLA-B/C) — turning three more cohorts
    (TCGA-CHOL, TCGA-HNSC, and for HLA-C also TCGA-KIRC) into genuine
    matched tumor-vs-normal tests, the same statistics as the BRCA headline
    test. *Caveat:* all 16 tests come back null (every OR CI spans 1) —
    this is a methods win (a real question answered: no, we don't need more
    raw BRCA normals) more than a new biological result.

12. **HLA-C's one real candidate allele (not A\*11:01's artifact) shows the
    same tumor-only pattern in every cohort tested so far.** The
    `hla_c_03_03`/`hla_c_03_04` intron5_retention junction from bullet 8's
    candidate, which rested on a single BRCA normal-side carrier, recurs in
    TCGA-CHOL, TCGA-HNSC, and TCGA-KIRC — 42 tumor carriers vs 2 normal
    carriers pooled across all 4 cohorts, never seen in 80 non-BRCA normal
    samples carrying either allele. *Caveat:* directionally consistent
    ≠ significant — every individual cohort test still has p>0.5, this
    wasn't a pre-registered comparison, and pooling raw counts across
    cohorts isn't a proper stratified test. Presented as "worth a dedicated
    follow-up," not as confirmed variable splicing. **Also worth
    mentioning if asked:** a first pass at this analysis almost let
    A\*11:01's already-known artifact back in as a fake "cohort-specific"
    HLA-A finding in TCGA-CHOL/TCGA-HNSC — the per-cohort re-derivation
    used elsewhere in this project requires more carriers than those small
    cohorts have. Caught before presenting anything, fixed by excluding
    A\*11:01 unconditionally everywhere, same as BRCA.

---

## Next steps / what would strengthen this

- ~~More normals.~~ **Answered by Anqi directly: no.** 106 is enough for
  the proportion comparison this analysis is built around; what actually
  moved things forward was resolving pan-cancer normal samples to their
  true cancer type (bullets 11-12), not adding raw BRCA count. The CI/power
  numbers above are still wide, and a larger BRCA cohort would still
  tighten them — but that's a "would be nice," not the blocking gap it was
  framed as before this conversation.
- **More normals specifically carrying `hla_c_03_03`/`hla_c_03_04`.**
  Unlike the point above, this one's still open: bullet 12's cross-cohort
  candidate has only 2 normal-side carriers total across 4 cohorts (all in
  BRCA) — a targeted pull of more normals known to carry either allele
  (typeable in advance from a sheet like `gdc_sample_sheet_normal_with_hla.tsv`)
  would directly test whether the tumor-only pattern holds up.
- **Confirm the A\*11:01 isoform against the actual reference transcript**
  (not just inferred from read-share) — worth a manual look at that allele's
  GTF/contig to settle whether it's a known, characterized minor isoform;
  now doubly motivated given it recurs pan-cancer (bullet 9).
- **Resolve HLA-B.** Two independent paths: (a) manually confirm with
  whoever built the personalized HLA-B references whether the
  reference-construction step for B matches C's (Phase 1 found C's
  coordinates behave as already transcript-oriented, contradicting the
  initial reverse-strand-flip hypothesis — worth confirming this applies to
  B too, not just inferring it from C), and (b) investigate why HLA-B's
  detectable-splicing ceiling is so much lower than A/C's in this cohort —
  could reflect a real, interesting biological difference (stronger
  purifying selection against alternative splicing at B?) worth its own
  slide once explained, or a fixable data issue.
- **The pan-cancer files' `BRCA_NORMAL` mislabeling** (see
  `OVERNIGHT_LOG.md`) should be fixed at the export/upstream step, not
  worked around per-analysis indefinitely.
