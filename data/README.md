# Input data

[Project overview](../README.md) · [Detailed methods](../docs/METHODS.md)

Inputs are preserved as supplied. Analysis scripts read these files and write results to `outputs/`.

| Files | Contents / use |
|---|---|
| `export_BRCA_{tumor,normal}_HLA-{A,B,C}.tsv` | BRCA splice-junction exports, one file per sample type and gene |
| `pan_HLA-{A,B,C}_SJ_combined.tsv` | Multi-cohort junction exports for descriptive and cohort-matched analyses |
| `gdc_sample_sheet_normal_with_hla.tsv` | Headered sample metadata; maps pooled normals to TCGA projects through `Case ID` and `Project ID` |
| `pancancer_status.txt` | Sample-count reference; not parsed by the scripts |
| `deprecated/brca_HLA-A_SJ_combined.tsv` | Superseded tumor-only export without a `Type` column; not used |

## Junction export schema

The current junction exports are tab-separated, with **no header**, and these 13 columns in order:

```text
Cohort  Type  Barcode  Allele  Contig  First_base  Last_base  Strand
Motif  Annotated  Uniquely_mapping  Multimapping  Overhang
```

`Type` distinguishes tumor from normal. Sample groupings use both `Type` and `Barcode`, since a barcode can occur in both groups. Preserve the input row order: junction bucketing depends on it. Definitions and row counts are in [Methods](../docs/METHODS.md#1-input-data).

## Known issues

- The HLA-B/C pan-cancer files contain a `BRCA_NORMAL` cohort labeled `tumor` in `Type`. The pan-cancer scripts exclude this block; BRCA normal comparisons use the dedicated exports.
- The pooled `NORMAL` cohort must be resolved to its originating TCGA project before cohort-specific comparisons. Phase 5 performs this mapping.
- HLA-B event classification remains unvalidated. See the [validation report](../outputs/PHASE1_BC_strand_check.md).
- Reference construction, alignment, and typing were performed upstream; those workflows are not supplied here.
