# A focused reproduction of the MHC Hammer exon-5 splicing axis in TCGA-BRCA HLA-A/B/C: an underpowered null and a characterized constitutive-splicing confounder

*Internal report. Sources: `METHODS.md`, `RESULTS_SUMMARY.md`, `PRESENTATION_NOTES.md`,
`outputs/FINAL_REVIEW.md`, `OVERNIGHT_LOG.md`, `outputs/PHASE1_BC_strand_check.md`,
`outputs/PHASE3_pancancer_descriptive.md`, `outputs/PHASE5_pancancer_matched_normal.md`,
`outputs/FIGURE_AUDIT.md`, `ANNOTATION_ARTIFACT_NOTE.md`,
`outputs/scratch/OPEN_QUESTIONS_DIAGNOSTIC.md`, and the CSVs under `outputs/`.
Every number below traces to one of these files; where a number is not available
from them, it is marked `[verify: …]` rather than estimated.*

---

## Abstract

We built a schema-first, statistically explicit pipeline to test two specific
alternative-splicing (AS) events at the HLA class I exon-5/intron-5 boundary —
complete exon-5 skipping and partial intron-5 retention — reported by MHC
Hammer (Puttick et al., *Nat Genet* 2024) as pervasive in normal tissue and,
within-patient, enriched in a subset of tumors. Using TCGA-BRCA STAR
splice-junction calls for HLA-A, HLA-B, and HLA-C (935 tumor / 106 normal
samples) and a deliberately strict, high-confidence detection threshold, we
find a rigorously power-characterized **null** result for both target events
in both testable genes (HLA-A, HLA-C): no tumor-vs-normal difference survives
Benjamini-Hochberg FDR correction anywhere in the project, and achieved
statistical power for the two headline comparisons ranges from 9.8% to 29.4%
— this is an honestly-quantified "cannot tell," not "no effect." Alongside
the null, the project's strongest finding is methodological: a
constitutive-annotation-artifact check caught HLA-A\*11:01 carrying a
near-universal (96.7–100% of carriers), tumor-independent minor splice
isoform that would have looked like a promising cancer signal in a
tumor-only design, and — on closer inspection prompted by this report — very
likely a second such constitutive event in HLA-C (`C*03:04`/`C*03:03`,
truncated by our detection threshold rather than genuinely rare). This is a
reproduction that returns an honestly power-characterized null plus the
identification of allele-specific constitutive splicing as a confounder for
this kind of analysis — **not** a positive replication of the paper's
tumor-enrichment finding, and it should not be read or presented as one.

---

## Background

**Biology.** HLA class I genes (HLA-A/B/C) undergo alternative splicing at
their exon-5/intron-5 boundary. Complete exon-5 skipping removes the
transmembrane-anchoring exon, producing a soluble HLA isoform; partial
intron-5 retention is a partial retention of intron 5 near the same
boundary. MHC Hammer (Puttick et al., *Nat Genet* 2024) profiled HLA
splicing across large normal and tumor cohorts and reported: (a) HLA
splicing is pervasive in **normal**
tissue, with exon-5 skipping the single most common event (76% of normal
lung samples, 69% of normal breast samples show it); (b) the large majority
of individual splicing events (~86–90%) occur at low abundance — below 10%
of novel-transcript proportion; (c) HLA-C\*03:04 shows an allele-specific,
essentially constitutive (~100% of carriers, 70/70 lung and 57/57 breast)
in-frame partial intron-5 retention in normal tissue; and (d) tumor
enrichment of AS, when it occurs, is a **within-patient** phenomenon —
compared against that same patient's matched adjacent-normal sample, present
in an estimated 15% of ER+ breast tumors — and the paper states matched
normals are required to make that call.

**What this project set out to test.** A focused reproduction of two of the
paper's four reported AS classes — complete exon-5 skip (`exon5_skip`) and
partial intron-5 retention (`intron5_retention`) — in TCGA-BRCA, extended to
all three classical class I genes (HLA-A/B/C) and, where the data allow, to
several pan-cancer cohorts. The paper's other two classes (partial exon
skip, complete intron retention) and other splice axes it reports (e.g.
exon-3 skipping) are **out of scope here**; complete intron retention was
not observed at all in the paper's own data and is not tested here either.

**Design difference, stated up front.** The paper's tumor-enrichment claim
is **within-patient**, tumor vs that same patient's matched adjacent-normal.
This project's BRCA comparison is an **unpaired cohort proportion** — tumor
samples vs a normal cohort that is not per-patient matched to the tumor
cohort (935 tumor samples vs 106 normal samples, not 106 matched pairs).
This is a materially different, weaker design for detecting a within-patient
effect, and every tumor-vs-normal number in this report should be read with
that in mind. This project does not attempt and does not claim a replication
of the paper's within-patient result.

**The contribution.** Beyond re-running the paper's target events, this
project's addition is a formal normal-tissue comparator with full statistical
accounting (Wilson CIs, odds ratios, Fisher exact tests, BH-FDR across
patterns, and an explicit power/minimum-detectable-effect analysis) plus a
constitutive-annotation-artifact check that is not part of the paper's
published method as reproduced here. That check is what surfaced this
project's two most important findings (below).

---

## Methods

Full detail, thresholds, and line-level code references: `METHODS.md`. Summary,
in pipeline order:

1. **Upstream, outside this repo (not verified here).** Raw TCGA-BRCA BAMs
   sliced to the MHC region, personalized per-patient HLA-A/B/C reference
   contigs built from each patient's own genotype, STAR-aligned to produce
   one `SJ.out.tab`-style junction file per sample × allele, then exported to
   the flat 13-column TSVs this repo reads. **Two steps here are explicitly
   unconfirmed inferences, not verified facts**, and are labeled as such
   throughout this report:
   - *HLA typing tool.* The per-patient two-field HLA-A/B/C genotype was
     most plausibly generated by a tool such as arcasHLA, inferred only from
     the column shape of `data/gdc_sample_sheet_normal_with_hla.tsv` and
     common practice — **no typing-tool name, version, or invocation appears
     anywhere in this codebase.**
   - *HLA-B/C reference orientation.* HLA-B/C are minus-strand genes, yet
     the exon5_skip/intron5_retention geometry (written in transcript-relative,
     effectively plus-strand terms) produces a coherent signal for HLA-C
     **without** a strand-flip transform (Methods §3d below). The working
     explanation — that personalized HLA-B/C reference contigs were likely
     already reverse-complemented at construction time — is an **inference**,
     never confirmed against the actual reference-construction code, which is
     not present in this repo.

2. **Per-sample derived columns** (`analyze_brca.compute_derived_columns()`):
   cumulative-sum bucketing of each row under its nearest upstream annotated
   junction (`Annotated_cumsum`), and `Local_rate` — a novel junction's share
   of `Uniquely_mapping` reads among every row competing in its own bucket.

3. **AS-event detection** (`find_as_events()`): a row qualifies as a
   candidate AS event if `Annotated == 0` (novel junction), `Annotated_cumsum
   != 0` (occurs downstream of a real annotated junction), and **`Local_rate >
   0.10`** — the novel junction must carry more than 10% of its local
   bucket's reads. **This 10% floor is a deliberate, stricter-than-the-paper
   choice.** The paper finds ~86–90% of its splicing events occur *below* a
   0.1 novel-transcript-proportion cutoff; this project's threshold sits
   exactly where the paper says most of the signal is *excluded*. The
   tradeoff is explicit: high confidence that anything detected is a real,
   non-noise event, at the cost of low sensitivity to the low-abundance
   majority of events the paper actually reports. This tradeoff is a
   throughline of the Results below — it has two demonstrated downstream
   consequences (Results §5 and §4).

4. **Geometry classification** (`classify_events()`): parses each event as
   `W1_W2 -> V1_V2` (wild-type → variant junction) and calls `exon5_skip`
   when the donor is unchanged and the acceptor jumps ≥100 bp downstream
   (`exon_skip_min_jump = 100`), or `intron5_retention` when the acceptor is
   unchanged and the donor shifts 1–50 bp upstream (`intron_retention_max_shift
   = 50`). Both are coordinate-offset judgment calls, not independently
   re-derived per gene against measured exon-5 length — no reference/GTF file
   with exon boundaries exists in this repo to do so
   (`outputs/scratch/OPEN_QUESTIONS_DIAGNOSTIC.md`).

5. **Strand-orientation validation (Phase 1).** Because HLA-B (~97.7–97.8%
   minus-strand) and HLA-C (~94.4–94.7% minus-strand) junctions might need a
   strand-aware coordinate mirror before the plus-strand-derived geometry
   above applies correctly, both genes were classified two ways — naive and
   strand-mirrored — and compared empirically. **HLA-C: naive classifier
   validated (`valid_as_is`)** — it produces a coherent, plausible signal (30
   carriers vs 0 for the flipped version, which destroys nearly all of it).
   **HLA-B: neither scheme could be validated (`untrustworthy`)** — only 2
   candidate AS events exist in the entire 935-tumor/106-normal cohort at
   this threshold (max observed `Local_rate` 12.6%, barely above the 10%
   floor, vs 52–60% for A/C), too few to discriminate the two schemes.
   Consequently **HLA-B is quarantined**: every HLA-B output is written under
   an `outputs/UNVALIDATED_HLA_B_*` prefix, excluded from every presented
   figure/table in this report, and is not un-quarantined here.

6. **Constitutive-annotation-artifact check.** For alleles with ≥20
   tumor+normal carriers, find each allele's single most-recurrent AS
   pattern and its penetrance (% of that allele's own carriers showing it)
   in tumor and normal separately. An allele is flagged
   `constitutive_annotation_artifact` when penetrance is ≥80% in every group
   with usable data — i.e., essentially everyone who carries the allele
   shows the event, tumor or normal alike, which carries no tumor-vs-normal
   information regardless of statistical significance.

7. **Tumor-vs-normal statistics**, identical across every gene/cohort:
   Wilson 95% CIs on carrier proportions, risk difference and odds ratio
   (Haldane-Anscombe corrected when a 2×2 cell is 0) with CIs, Fisher exact
   tests, and Benjamini-Hochberg FDR across every tested per-pattern
   comparison (`min_carriers ≥ 3`).

8. **Power analysis.** Two-proportion z-test achieved power (Cohen's *h*) for
   the observed effect, plus the minimum detectable effect at 80% power given
   the actual n=935/106 — reported as an honest "how underpowered," not a
   pass/fail label.

### Every threshold used, and whether it was re-tuned per gene

| Parameter | Value | Re-tuned for B/C or pan-cancer? |
|---|---|---|
| `local_rate_threshold` | 0.10 | No — identical everywhere |
| `exon_skip_min_jump` | 100 bp | No — inherited from HLA-A |
| `intron_retention_max_shift` | 50 bp | No — inherited from HLA-A |
| `min_allele_samples` (BRCA scale) | 20 | No |
| `min_allele_samples` (pan-cancer per-cohort) | 10 | Lowered deliberately — smaller cohorts |
| `min_carriers` (per-pattern Fisher/FDR) | 3 | No |
| `exclude_alleles` (HLA-A) | `hla_a_11_01` | Re-derived per gene for B/C; none flagged |
| Constitutive-artifact gates | `HIGH_PENETRANCE=80%`, `LOW_PENETRANCE=5%`, `min_n=10` | No |

Full table with rationale: `METHODS.md` §9.

---

## Results

### 1. The pipeline generalizes across genes and cohorts

The HLA-A pipeline, unmodified, was pointed at a second multi-cohort file
(`data/pan_HLA-A_SJ_combined.tsv`, 9,738 rows — direct file count — across 6
cohorts, no BRCA) with zero code changes and correctly loaded, bucketed, and
classified 46 candidate AS events (`outputs/PHASE3_pancancer_descriptive.md`).
Distinct-sample count: 301 distinct barcodes by direct count of the file;
`outputs/PHASE3_pancancer_descriptive.md`'s own summary states "311 total
pan-cancer samples," which double-counts 10 barcodes that the file labels
under more than one `Cohort` value — both numbers trace to the same
underlying file, they are just counted two different ways, and the
discrepancy is noted here rather than silently resolved one way. This is a
capability demonstration, not a biological claim; it established the
pipeline was ready to be pointed at HLA-B/C and pan-cancer cohorts without
gene-specific rewrites — **an engineering generalization only**: this check
confirms the code runs on a new file, not that the exon5_skip/intron5_retention
analysis is valid for every gene it can run on. Per-gene analysis validity
was checked separately (Methods §5), and HLA-B failed that check (§6) even
though the same code executed without error on HLA-B's data too.
*Figures:* `figures_overview.pdf`; no dedicated figure for this check.

### 2. BRCA tumor-vs-normal: a null, led by its power, not by its p-value

**Lead with the power, not the p-value.** Both headline comparisons in BRCA
are underpowered to detect an effect of the size the paper reports, and we
can say precisely how underpowered (`outputs/table_power_analysis.csv`):

| Event class | Tumor | Normal | Fisher p | OR (95% CI) | Achieved power | MDE at 80% power |
|---|---|---|---|---|---|---|
| exon5_skip (A\*11:01 excluded) | 0.1% (1/935) | 0.0% (0/106) | 1.0 | 0.34 (0.01–8.45)† | **9.8%** | tumor rate would need to fall outside [0.0%, 2.0%] |
| intron5_retention (A\*11:01 excluded) | 4.2% (39/935) | 7.5% (8/106) | 0.13 | 0.53 (0.24–1.17) | **29.4%** | tumor rate would need to fall outside [1.8%, 16.8%] |

† Haldane-Anscombe corrected (a 2×2 cell was 0). Note: at 1 tumor carrier vs
0 normal carriers, exon5_skip is effectively uninterpretable as a direction
either way — the raw proportions (0.1% > 0.0%) point "enriched in tumor"
(as labeled on `outputs/fig_classes.pdf`), while the Haldane-corrected OR
(0.34, <1) points "depleted." Both are correct arithmetic on the same 1-vs-0
cell counts; neither should be read as a real direction. See §5 for why this
event class is near-absent at our threshold in the first place.

At 9.8% and 29.4% achieved power, this cohort's normal side (n=106) cannot
reliably detect an effect the size the paper reports in **either** direction.
`intron5_retention`'s point estimate trends *depleted* in tumor — the
opposite of the paper's predicted enrichment — but with a 95% CI that
crosses 1 (0.24–1.17) and 29.4% power, this is honestly "cannot tell," not a
contradiction of the paper. A 12-setting one-at-a-time sweep of all three
judgment-call thresholds keeps the direction "depleted in tumor" in all 12
settings (OR range 0.41–0.60), so the null is not an artifact of one
specific cutoff choice. **0 of 6 tested per-pattern comparisons survive
BH-FDR q<0.05** in the full per-pattern table.
*Figures:* `outputs/fig_classes.pdf`, `outputs/fig_threshold_sensitivity.pdf`,
`outputs/fig_volcano.pdf` (nominal-vs-FDR significance line labeled per
`outputs/FIGURE_AUDIT.md`), `figures_overview.pdf` panels A1/C.

### 3. HLA-A\*11:01: a constitutive splice isoform, not tumor biology — the strongest finding in the project

The `hla_a_11_01` allele's top recurrent event (`2365_2806->2383_2806`,
classified `intron5_retention`) is present in **96.7% of tumor carriers
(118/122)** and **100% of normal carriers (19/19)** — no other allele in the
cohort comes close, and the rate is identical whether the sample is tumor or
normal (`ANNOTATION_ARTIFACT_NOTE.md`). A read-share check rules out the
simplest explanation (a reference-coordinate swap): the **annotated**
junction stays the majority isoform at ~72.4% of local reads on average in
both groups, so this looks like a real, fixed, ~25–30%-penetrant minor
splice isoform of A\*11:01's reference transcript — a property of the
allele, not of cancer status. In a tumor-only design (no matched normal
comparator), this allele's high recurrence alone would have looked like a
promising cancer-associated signal; having normals is what caught it.
Consequently `hla_a_11_01` is excluded from the headline/per-allele
comparisons throughout (never from descriptive outputs). A pan-cancer
descriptive check (not a new significance test — no per-cancer-type matched
normal for this comparison) found the same junction recurs at 100% (20/20)
in the pooled normal-tissue reference, 100% (3/3) in TCGA-CHOL, and 85.7%
(6/7) in TCGA-LGG (TCGA-HNSC 1/3 is a small-n outlier); pooled across all
pan-cancer cohorts, 30/33 (91%) of A\*11:01 carriers show the pattern — this
is not a BRCA-specific quirk.
*Figures:* `outputs/fig_per_allele.pdf` / `outputs/fig_per_allele_talk.pdf`,
`figures_overview.pdf` panels B and D.

### 4. HLA-C: very likely a *second* constitutive artifact, truncated by our threshold — not tumor-enrichment [RESOLVED]

HLA-C's own recurrent candidate — `hla_c_03_04`/`hla_c_03_03`, pattern
`1273_1711->1291_1711`, classified `intron5_retention` — was previously
reported in this project's other documents (`RESULTS_SUMMARY.md` Finding 4/6)
as a `real_variable_splicing_candidate` with a tumor-enriched odds ratio
(OR 3.37) and cross-cohort directional corroboration. **That framing is
retired here.** MHC Hammer reports HLA-C\*03:04 partial intron-5 retention
as an essentially constitutive event — present in 100% of carriers (70/70
lung, 57/57 breast) in normal tissue (Puttick et al., *Nat Genet* 2024). Our
measured penetrance for the *same allele, same intron-5-retention geometry*
is far below that: 19.33% tumor (23/119) / 11.11% normal (1/9) for
`hla_c_03_04`, and 13.75% (11/80) / 7.14% (1/14) for `hla_c_03_03`
(`outputs/scratch/OPEN_QUESTIONS_DIAGNOSTIC.md`, cross-checked against
`outputs/table_per_allele_HLA-C.csv`). The gap between our ~11–19%
penetrance and the paper's ~100% is exactly in the direction our 10%
local-rate floor would produce if the true event is constitutive but
frequently expressed below that floor (Methods §3, the same tradeoff
already flagged there) — i.e., **very likely the same constitutive event
Finding 3 already found for A\*11:01, a second instance of the same
detection-threshold confound, not real tumor-enriched splicing.**

We state this as "very likely," not "is." The one remaining confirmation —
mapping coordinates 1291/1711 on the personalized HLA-C\*03:04 contig to
intron 5 to verify geometric identity with the paper's event — **could not
be done**: no reference contig, GTF, or feature-table file for the
personalized HLA-C\*03:04 sequence exists anywhere in this repository
(`outputs/scratch/OPEN_QUESTIONS_DIAGNOSTIC.md`). This is the single
open follow-up for this finding (see Follow-ups).

**Consequence for prior framing:** the previously-reported OR of 3.37 for
HLA-C `intron5_retention` (Results §2-style headline test, all C alleles
pooled) rested on a single normal-side carrier (1/106) and was already
non-significant (95% CI 0.45–25.0). Per-allele, the same pattern is
non-significant for both alleles individually (Fisher p=0.69 for
`hla_c_03_03`, p=1.0 for `hla_c_03_04`). None of this was ever a positive
finding by this project's own statistical bar; it is now better understood
as a **cautionary near-miss** — the constitutive-artifact check
(Results §3/Methods §6) would have flagged this allele the same way it
flagged A\*11:01, had the 10%-threshold-truncated penetrance not sat below
the check's 80% gate.
*Figures:* `outputs/fig_classes_HLA-C.pdf`, `outputs/fig_per_allele_HLA-C.pdf`
/ `outputs/fig_per_allele_talk_HLA-C.pdf`, `outputs/fig_volcano_HLA-C.pdf`.

### 5. exon5_skip near-absence: confirmed threshold-driven, not biological [RESOLVED]

HLA-A's `exon5_skip` is nearly unobserved at the presented 10% threshold (1
distinct tumor carrier, 0 normal). Because the paper reports exon-5 skipping
as its *single most common* event genome-wide (76% normal lung, 69% normal
breast) and ~86–90% of events occur below a 0.1 novel-transcript-proportion
cutoff, this near-absence was an open question: biological, or an artifact
of sitting at exactly the threshold the paper says excludes most signal?

An in-memory sensitivity diagnostic — the pipeline's own unmodified
`load_data` → `compute_derived_columns` → `find_as_events` →
`annotate_events` → `classify_events` functions, re-run with
`local_rate_threshold=0.02` passed as an explicit argument (no config file,
default, or output of record changed) — raised `exon5_skip` from **1 tumor /
0 normal distinct carriers at 0.10** (2 candidate rows collapsing to 1
distinct carrier, consistent with §2) to **44 tumor / 7 normal distinct
carriers at 0.02** (`outputs/scratch/OPEN_QUESTIONS_DIAGNOSTIC.md`). The
0.10-threshold counts in this diagnostic match the presented pipeline's own
output exactly, confirming the diagnostic used the identical code path, not
a divergent reimplementation.

**Conclusion: exon5_skip's near-absence at our 10% threshold is a
threshold/sensitivity effect, consistent with the paper's own low-abundance
regime for this event class — not biological absence.** We present this as
**confirmed by diagnostic**, not as a validated alternative analysis: no
0.02-threshold result set exists of record, no output under `outputs/` was
regenerated at 0.02, and no figure/table in this report uses the 0.02 run.
It answers "is the null real or a threshold artifact" without proposing 0.02
as a replacement threshold. The exon-5-length check that would have been an
independent sanity check on this remains unavailable (no GTF/reference file
in this repo) but is not needed given this result.

### 6. HLA-B: quarantined, insufficient signal at our threshold

HLA-B produces only 2 candidate AS events in the entire 935-tumor/106-normal
BRCA cohort at the 10% local-rate threshold (vs 60 for HLA-C, ~200+ for
HLA-A) — its maximum observed `Local_rate` in this cohort is 12.6%, barely
above the detection floor, vs 52–60% for A/C. Neither the naive nor the
strand-mirrored classifier could be validated against this little data
(Methods §5). Whether this ceiling reflects real biology (less alternative
splicing at HLA-B) or a technical/reference-quality issue specific to B is
**not determined** by this project. All HLA-B output exists for audit under
`outputs/UNVALIDATED_HLA_B_*` and is excluded from every result and figure
above; it is not un-quarantined in this report and should not be cited as a
finding, positive or negative.

### 7. Phase 5: a reusable matched-normal capability, not a solved power problem

A separate GDC sample sheet (740 normal samples, 24 TCGA projects, with
HLA-A/B/C typing) let this project resolve the pan-cancer files' pooled,
cohort-undifferentiated `NORMAL` block back to each sample's true TCGA
project of origin (145/145 barcodes resolved for HLA-A, 318/318 for HLA-B/C)
— turning a purely descriptive pan-cancer pass into genuine matched (same
cancer type) tumor-vs-normal tests, using the identical
Fisher/Wilson-CI/odds-ratio machinery as the BRCA headline test, wherever a
resolved cohort has ≥3 samples on both sides: TCGA-CHOL (33 tumor/8 normal)
and TCGA-HNSC (52 tumor/5 normal) for HLA-A; TCGA-CHOL (33/8), TCGA-HNSC
(67/35), and TCGA-KIRC (92/37) for HLA-C (HLA-B computed identically but
quarantined per §6).

**This should be read as a demonstrated *capability* — the pipeline can run
a real matched-normal test on any cohort a typing sheet resolves — not as
evidence the power problem is solved.** The newly-matched normal sides are
themselves small (as few as 5 samples, TCGA-HNSC/HLA-A) and, like BRCA, are
unpaired cohort proportions, not within-patient pairs. All 10 presented
(gene × cohort × class) tests come back null — every OR 95% CI spans 1, 0
have an OR CI excluding 1 at the nominal level (`outputs/table_phase5_summary.csv`,
`outputs/FIGURE_AUDIT.md`) — consistent with, not independent confirmation
beyond, the BRCA-scale power ceiling in §2. One methodological catch is
worth noting: an early version of this pass re-derived the A\*11:01
exclusion per cohort (requiring ≥10–20 carriers to flag anything) and missed
it in TCGA-CHOL/TCGA-HNSC, where only 2–3 tumor samples carry the allele —
100% of the *original* HLA-A TCGA-CHOL/TCGA-HNSC "enrichment" signal turned
out to be A\*11:01's already-known artifact reappearing. Fixed by excluding
A\*11:01 unconditionally, regardless of per-cohort carrier count.
*Figures:* `outputs/fig_phase5_matched_normal.pdf` (tested/significant count
annotated per `outputs/FIGURE_AUDIT.md`).

---

## Discussion

Three of the results above share one cause. The project's 10% local-rate
threshold was chosen deliberately for high confidence (anything detected is
almost certainly not alignment noise) at the acknowledged cost of low
sensitivity — the paper itself reports ~86–90% of real events sit below
that cutoff. That single choice explains **both** of this report's two
resolved open questions: exon5_skip's near-total absence at 10% (§5), which
a 0.02 diagnostic shows is a threshold effect, not biology; and HLA-C's
`C*03:04`/`C*03:03` candidate's measured penetrance (11–19%) sitting far
below the paper's reported ~100%, in exactly the direction a
threshold-truncated constitutive event would produce (§4). Put together with
the A\*11:01 finding (§3), the throughline is: **this project's real
methodological contribution is not a tumor-vs-normal effect size — it is a
constitutive-annotation-artifact check, now demonstrated on two independent
alleles across two genes, that would otherwise have been mistaken for
biology in a design without matched normals.** A tumor-only version of this
analysis would very plausibly have reported both A\*11:01 and
`C*03:04`/`C*03:03` as promising tumor-associated splicing hits; having
normals, and building a formal penetrance-based artifact check rather than
trusting recurrence alone, is what caught both.

None of this should be read as evidence against the paper's tumor-enrichment
finding. The paper's claim is within-patient, matched-adjacent-normal,
reported in an estimated 15% of ER+ breast tumors; this project's unpaired,
cohort-level BRCA comparison — even before the threshold and constitutive-
artifact issues above — was never powered to detect an effect of that shape
(§2). The honest summary is: this reproduction neither confirms nor refutes
the paper's tumor-enrichment claim, and identifies two mechanisms (threshold
truncation, allele-specific constitutive splicing) that a future, better-
powered, within-patient-paired attempt at this question would need to
account for.

---

## Limitations

- **Unpaired vs. within-patient design.** This project's BRCA (and Phase 5
  pan-cancer) comparisons are unpaired cohort proportions; the paper's
  tumor-enrichment claim is within-patient, matched-adjacent-normal. No
  result here is a replication of that claim, positive or negative.
- **Detection-threshold tradeoff.** The 10% local-rate floor is stricter
  than the regime the paper reports most events occur in (~86–90% below
  0.1). This is a deliberate high-confidence/low-sensitivity choice, and it
  is now demonstrated to materially suppress signal for at least two known
  events (exon5_skip generally, §5; HLA-C's `C*03:04`/`C*03:03` candidate,
  §4). No result-of-record in this project used a threshold other than 0.10.
- **Geometry thresholds not re-tuned per gene.** `exon_skip_min_jump`
  (100 bp) and `intron_retention_max_shift` (50 bp) were validated for HLA-A
  and carried over unchanged to HLA-B/C rather than independently re-derived
  against each gene's measured exon-5 length — no reference/GTF file with
  exon boundaries exists in this repo to do that re-derivation.
- **HLA-B is fully quarantined.** Insufficient signal at our threshold (2
  candidate events total, max `Local_rate` 12.6%) to validate a classifier
  or test anything; real biology vs. a technical/reference-quality issue
  specific to HLA-B is undetermined. No HLA-B number in this report should
  be cited as a finding.
- **A known export bug was worked around, not fixed.** The pan-cancer
  HLA-B/C files contain a `Cohort == "BRCA_NORMAL"` block with `Type`
  mislabeled `"tumor"`; every downstream script excludes this block
  explicitly rather than relying on a corrected upstream export. Anyone
  reusing `data/pan_HLA-{B,C}_SJ_combined.tsv` without this project's notes
  would silently double-count BRCA normals as tumor samples.
- **Two upstream inferences remain unconfirmed.** The HLA typing tool
  (plausibly arcasHLA, inferred only from a sample-sheet column layout) and
  the explanation for why HLA-B/C's naive (unflipped) coordinate geometry
  behaves correctly despite both genes being majority minus-strand (the
  working theory: personalized reference contigs were already
  reverse-complemented at construction time) are both inferences about
  steps outside this repo, never confirmed against the actual upstream
  reference-construction or typing code.

---

## Follow-ups

1. **Paired McNemar test on the 106 BRCA matched tumor/normal pairs.** The
   current headline test is an unpaired cohort proportion; a McNemar test on
   actually-paired samples would be a materially better-matched design
   against the paper's within-patient claim and does not require new data.
2. **Map HLA-C\*03:04 coordinates 1291/1711 to confirm RESOLVED finding
   §4.** Requires the personalized HLA-C\*03:04 reference contig and its
   exon/intron feature annotation from the upstream reference-construction
   step — not present in this repo; this is the single remaining
   confirmation for the "very likely a second constitutive artifact" call.
3. **Exon-3 skip and the paper's remaining AS classes** (partial exon skip;
   complete intron retention, not observed by the paper itself) — out of
   scope for this reproduction, open for a future pass.
4. **Confirm upstream reference construction and HLA typing directly** with
   whoever built the personalized per-patient HLA-A/B/C contigs: (a) whether
   HLA-B/C contigs were reverse-complemented at construction time (the
   inference this project's HLA-C classifier validation rests on), and (b)
   the actual HLA-typing tool/version used, before citing either upstream
   the way a publication would require.
5. **Resolve HLA-B's low event yield** — real biology vs. reference/alignment
   quality — as a prerequisite to un-quarantining any HLA-B result.
