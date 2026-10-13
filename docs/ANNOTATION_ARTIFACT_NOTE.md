# Annotation-artifact check: is HLA-A*11:01 intron5_retention real?

Diagnostic only -- does not change the main analysis in `outputs/`.

## Verdict

**hla_a_11_01's recurrent event (2365_2806->2383_2806, classified as `intron5_retention`) looks constitutive, not tumor-associated.** It is present in 97% of tumor carriers (118/122) and 100% of normal carriers (19/19) of this allele -- almost everyone who carries A\*11:01 shows it, tumor or normal alike, which is the signature of an allele-intrinsic feature, not variable splicing that cancer turns on.

**But the read-share check argues against the simplest explanation.** If the reference/GTF had the exon-5 donor site at the wrong coordinate, we'd expect the 'novel' junction to dominate the read counts. It doesn't: the ANNOTATED junction remains the majority isoform at 72% of local reads, on average, vs only 27% for the 'novel' junction, consistently across every tumor and normal sample inspected in step 1 (annotated junction won in all 8/8 examples pulled, tumor and normal alike). So this isn't a coordinate swap -- it looks more like **A\*11:01's reference transcript has a real, fixed, ~25-30% minor splice isoform that essentially every carrier produces**, independent of tumor/normal status. Either way (coordinate error or genuine-but-universal minor isoform), the practical conclusion is the same: this allele's 'event' carries no tumor-vs-normal information and shouldn't be read as cancer-associated alternative splicing.

Across all 6 alleles checked (>= 20 total samples): **1 look like annotation artifacts**, **3 look like real variable splicing** (present in some carriers, absent in others), and 2 are too rare to call.

## How to read the table (`outputs/annotation_artifact_check.csv`)

For each allele, `Top_AS_pattern` is its single most-recurrent event. `Tumor_pct`/`Normal_pct` are the fraction of that allele's carriers (per Type) showing it; `Min_group_n` is the smaller of the two group sizes, so you can spot calls resting on a handful of samples. `Avg_variant_local_rate_pct` / `Read_share_pattern` describe *which kind* of constitutive signature it is, but do not gate the verdict -- the verdict is penetrance-only, per the question actually being asked (variable across carriers vs universal).

`Verdict` logic (penetrance-gated, matching 'present in some carriers' vs 'present in ~all carriers'):
- **constitutive_annotation_artifact**: penetrance >= 80% in every group with usable data -- basically everyone who carries this allele shows the event, tumor or normal alike, so it carries no tumor-vs-normal information regardless of how the reads split between the two junctions.
- **constitutive_annotation_artifact (one group underpowered)**: penetrance is near-universal in the one group with >= 10 samples, but the other group has too few carriers of this allele to confirm or refute the same -- treat as the same call, with that caveat attached.
- **real_variable_splicing_candidate**: penetrance is between roughly 5% and 80% -- present in some carriers, absent in others. This is what real, biologically variable alternative splicing should look like, and is the only category safe to present as a tumor-vs-normal splicing difference.
- **rare_private_event**: penetrance under 5% in every group with data -- too rare to characterize either way.

Separately, `Read_share_pattern` tells you *why* a constitutive call is constitutive: 'novel junction dominates' suggests a reference/GTF coordinate error (the annotated junction is actually the spurious one); 'annotated junction stays dominant' (A\*11:01's case) suggests a real but universal minor splice isoform baked into that allele's transcript, unrelated to any coordinate mistake.

## Practical implication

`hla_a_11_01` carriers should be **excluded before re-running the headline tumor-vs-normal comparison**, or at minimum footnoted -- as things stand, its near-universal 'event' inflates both groups' apparent intron5_retention carrier rate roughly equally, which *dilutes* rather than creates a tumor-vs-normal difference for that class, but it is not a real biological readout for this allele and shouldn't be presented as one. Any allele flagged `real_variable_splicing_candidate` is a legitimate candidate for the per-allele slide; alleles flagged as constitutive artifacts should either be dropped from that figure or footnoted, and (for A\*11:01 specifically) the reference contig is worth a manual look since the effect is consistent and large even though it isn't a simple coordinate swap.
