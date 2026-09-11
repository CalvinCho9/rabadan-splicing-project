# Presentation notes — HLA-A splicing, TCGA-BRCA

Talking points for the lab meeting. Each bullet: **claim** → supporting stat
→ honest caveat. Full detail in `RESULTS_SUMMARY.md`; slide panel in
`figures_overview.pdf`; method/background reference sheet in
`background_method.pdf` if the room needs the mechanics explained first.

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

---

## Next steps / what would strengthen this

- **More normals.** 106 is the binding constraint on every CI and power
  number above — even a doubling would meaningfully tighten the
  intron5_retention estimate and raise achieved power well past the current
  ~30%.
- **Confirm the A\*11:01 isoform against the actual reference transcript**
  (not just inferred from read-share) — worth a manual look at that allele's
  GTF/contig to settle whether it's a known, characterized minor isoform.
- **Extend to HLA-B/C and the full pan-cancer file** once alignment finishes
  for those — the pipeline (this repo, `analyze_brca.py` +
  `additional_analyses.py`) needs no code changes, only new
  `--tumor`/`--normal` file paths and a gene-appropriate re-check of the
  exon5_skip/intron5_retention coordinate thresholds (already demonstrated
  working unmodified on a second, non-BRCA multi-cohort file — Finding 1).
