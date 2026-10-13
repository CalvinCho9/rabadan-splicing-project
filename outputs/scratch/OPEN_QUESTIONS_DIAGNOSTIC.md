# Open-questions diagnostic (read-only)

Scratch file only. Nothing under `outputs/` referenced by a presented
document was modified. No pipeline default (`local_rate_threshold`,
`exon_skip_min_jump`, `intron_retention_max_shift`, or any other parameter)
was changed in place — the Q2 re-classification below passes
`local_rate_threshold=0.02` as an explicit function argument to
`analyze_brca.py`'s existing, unmodified `find_as_events()`, in a standalone
in-memory Python session; no source file was edited and no `outputs/*.csv`
was overwritten.

---

## Q1 — HLA-C candidate penetrance

Source files (both exist, both present in `outputs/`, and both contain
identical rows for the alleles in question):
- `outputs/table_phase1_penetrance_HLA_C.csv`
- `outputs/table_penetrance_recheck_HLA-C.csv`

Both list `1273_1711->1291_1711` as the **Top_AS_pattern** (i.e. the single
most-recurrent AS pattern) for both `hla_c_03_04` and `hla_c_03_03`, so the
Tumor/Normal carrier counts in these rows are specifically this pattern's
penetrance among each allele's carriers.

| Allele | Group | Carriers | Allele carriers (n) | Penetrance (%) |
|---|---|---|---|---|
| hla_c_03_04 | Tumor | 23 | 119 | 19.33% |
| hla_c_03_04 | Normal | 1 | 9 | 11.11% |
| hla_c_03_03 | Tumor | 11 | 80 | 13.75% |
| hla_c_03_03 | Normal | 1 | 14 | 7.14% |

From `outputs/table_per_allele_HLA-C.csv` (intron5_retention rows — this
pattern classifies as `intron5_retention` per both penetrance tables' `Top_event_class` column):

| Allele | Tumor carrier % (n) | Normal carrier % (n) | Fisher p |
|---|---|---|---|
| hla_c_03_04 | 19.33% (23/119) | 11.11% (1/9) | 1.0 |
| hla_c_03_03 | 13.75% (11/80) | 7.14% (1/14) | 0.686 |

These two sources agree exactly (same numerator/denominator/pct in both
tables), as expected since `intron5_retention` is each allele's top pattern.

**Placeholder for your confirmation:**
> MHC Hammer (Puttick 2024) reports C\*03:04 partial intron 5 retention at
> ~100% penetrance in carriers — [I will confirm this figure against the
> paper myself].

**Hypothesis to be confirmed (not resolved here):**
> IF our penetrance (19.33% tumor / 11.11% normal for C\*03:04; 13.75% /
> 7.14% for C\*03:03) is far below ~100%, this is consistent with the same
> constitutive event truncated by the 10% local-rate floor (i.e. a second
> A\*11:01-like artifact), NOT tumor-enrichment — REQUIRES my confirmation.

What I can say from the repo alone, without asserting the biological call:
the repo's own `add_verdict()` logic (`hla_common.py`) labeled both alleles
`real_variable_splicing_candidate`, not `constitutive_annotation_artifact` —
but that verdict is a `min_n=10`-gated, `HIGH_PENETRANCE=80`/`LOW_PENETRANCE=5`
threshold rule (see `hla_common.py:176-190`), not evidence against your
hypothesis; it only says the observed tumor/normal percentages in this data,
as measured at `local_rate_threshold=0.1`, didn't clear the 80% bar the repo
already uses elsewhere to call something constitutive-like. Whether that
threshold choice itself is masking a truly ~100%-penetrant event (as your
placeholder suggests, citing the paper) is exactly the open question — not
something this repo's code can adjudicate on its own.

---

## Q1-rigorous — coordinate mapping

Searched the entire repo (`find` on `*.gtf`, `*.gff*`, `*.fa`, `*.fasta`,
`*.bed`, and any file/dir name containing "annotation", "reference", or
"contig"; grepped every `.py`/`.md` file for `GTF`/`.gff`/"reference
contig"/"annotation file").

Result: **no personalized HLA-C\*03:04 reference contig, exon/intron
feature table, GTF, or FASTA exists anywhere in this repository.**
`data/` contains only the 13-column headerless STAR `SJ.out.tab`-style
export TSVs (junction coordinates + read counts) and the GDC sample sheet —
no sequence or annotation file. The repo's own documentation confirms this
directly and repeatedly:

- `METHODS.md:14`: "Per-patient personalized reference construction + STAR
  alignment... **Not verified in this repo** — the reference-construction
  and alignment code is not present here."
- `METHODS.md:135`: "the personalized per-patient reference contigs for
  HLA-B/C were likely already reverse-complemented at reference-construction
  time... **This is an unverified assumption about a step outside this
  repo**."
- `OVERNIGHT_LOG.md:110-120, 457-468`: documents this as an *inference*,
  never confirmed against the actual upstream reference-construction code,
  and lists confirming it with whoever built the personalized HLA-B/C
  references as an open follow-up item.

**CANNOT VERIFY — upstream C\*03:04 reference contig + annotation not
present in repo.** Stopping here per instructions; not substituting a
generic HLA-C annotation, and not guessing whether 1291/1711 fall at the
exon-5 donor / intron-5 region.

---

## Q2 — is exon5_skip absence threshold-driven

### Exon-5 length

Searched for any GTF/reference/feature-table file in the repo (same search
as Q1-rigorous, gene-agnostic). None exists — the repo contains junction
coordinate data only, never a feature/annotation table that states exon
boundaries or lengths.

**CANNOT VERIFY — exon-5 length source not in repo.** Not stating a length
from general knowledge.

### Diagnostic re-classification (HLA-A, in-memory only)

Ran `analyze_brca.py`'s existing, unmodified `load_data()` →
`compute_derived_columns()` → `find_as_events()` → `annotate_events()` →
`classify_events()` functions in a standalone Python session (not `main()`,
no CLI, no file writes), first at the real default
`local_rate_threshold=0.1`, then again at `local_rate_threshold=0.02`
(passed as an explicit function argument; `exon_skip_min_jump=100` and
`intron_retention_max_shift=50` left at their existing default values both
times — unchanged). Cohort: `data/export_BRCA_tumor_HLA-A.tsv` +
`data/export_BRCA_normal_HLA-A.tsv` (935 tumor / 106 normal samples, same
denominators the presented HLA-A headline uses).

| threshold | exon5_skip rows — tumor | exon5_skip rows — normal | exon5_skip distinct carriers — tumor | exon5_skip distinct carriers — normal |
|---|---|---|---|---|
| 0.10 (real/presented) | 2 | 0 | 1 | 0 |
| 0.02 (diagnostic only) | 69 | 8 | 44 | 7 |

The 0.10-threshold row/carrier counts here (2 rows / 1 distinct carrier,
tumor; 0/0, normal) match exactly what `analyze_brca.py`'s presented
pipeline already reports (`table_headline_classes_all_alleles.csv`: tumor
1/935 = 0.1%, normal 0/106 = 0.0%), confirming this diagnostic re-run used
the identical code path as the presented output and isn't a divergent
reimplementation.

**Interpretation (not resolved — yours to make):**
> IF exon5_skip events appear at 0.02 but not 0.10, the near-absence at
> 0.10 is a sensitivity/threshold effect (consistent with the paper
> reporting exon-5 skip as common but low-abundance), NOT biological
> absence. IF they remain absent even at 0.02, the cause is geometry or a
> coordinate/strand issue and needs the exon-5-length check above.

Since exon5_skip events go from ~0 (2 rows / 1 carrier at 0.10) to
substantially present (69 rows / 44 carriers at 0.02), the observed pattern
matches the **first** branch of your stated interpretation — but per your
instructions this diagnostic stops at reporting the counts; it does not
assert that the geometry/coordinate hypothesis is therefore ruled out, and
it does not make the biological call for you.

---

## What still needs your confirmation

1. **Q1** — the paper's actual reported penetrance figure for C\*03:04
   partial intron-5 retention (placeholder says "~100%"; not checked against
   Puttick 2024 in this diagnostic).
2. **Q1-rigorous** — the coordinate identity of 1291/1711 relative to the
   true exon-5 donor / intron-5 boundary on the personalized C\*03:04
   contig. Cannot be checked at all from this repo; requires the upstream
   reference-construction output (contig sequence + annotation), which is
   not present here.
3. **Q2** — whether the 0.10→0.02 threshold sensitivity shown above is
   itself evidence for "common but low-abundance" (paper-consistent) versus
   a coordinate/geometry artifact — this diagnostic reports the counts only
   and does not adjudicate between those two explanations. The exon-5
   length needed to sanity-check the geometry hypothesis independently is
   also not available in this repo.
