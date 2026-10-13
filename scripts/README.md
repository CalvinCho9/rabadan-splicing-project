# Analysis scripts

[Project overview and setup](../README.md#reproduce-the-analysis) · [Detailed analysis guide](../docs/ANALYSIS_GUIDE.md)

Run commands from the **repository root**, for example `python scripts/analyze_brca.py`. Inputs and outputs are resolved relative to that directory. The scripts share helpers through `analyze_brca.py` and `hla_common.py`.

| Order | Script | Purpose |
|---|---|---|
| 1 | [analyze_brca.py](analyze_brca.py) | Main HLA-A BRCA comparisons |
| 2 | [diagnostic_annotation_check.py](diagnostic_annotation_check.py) | HLA-A*11:01 annotation investigation |
| 3 | [additional_analyses.py](additional_analyses.py) | Power, threshold sensitivity, and allele frequencies |
| 4 | [phase1_bc_strand_check.py](phase1_bc_strand_check.py) | HLA-B/C classifier validation |
| 5 | [phase2_bc_headline.py](phase2_bc_headline.py) | HLA-B/C BRCA comparisons |
| 6 | [phase3_pancancer_descriptive.py](phase3_pancancer_descriptive.py) | Descriptive pan-cancer analysis |
| 7 | [phase5_pancancer_matched_normal.py](phase5_pancancer_matched_normal.py) | Cohort-matched pan-cancer comparisons |
| Optional | [make_figures_overview.py](make_figures_overview.py) | Combined presentation results |
| Optional | [make_background_figures.py](make_background_figures.py) | Method schematic |
| Library | [hla_common.py](hla_common.py) | Shared loading and classification helpers |

Phase numbers retain the original analysis sequence. See [outputs](../outputs/README.md) for saved results without running the code.
