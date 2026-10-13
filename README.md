# HLA alternative splicing in breast cancer and across TCGA cohorts

Analysis of allele-specific HLA-A, HLA-B, and HLA-C splice-junction exports, comparing tumor and normal samples in TCGA-BRCA and selected additional cancer cohorts. The project examines exon-5 skipping, partial intron-5 retention, and recurrent junction patterns, with checks for allele-specific annotation artifacts.

## Start here

| To explore… | Open |
|---|---|
| Findings and interpretation | [Results summary](docs/RESULTS_SUMMARY.md) |
| Analysis definitions and statistics | [Methods](docs/METHODS.md) |
| Figures, tables, and validation reports | [Results index](outputs/README.md) |
| Input files and known data issues | [Data guide](data/README.md) |
| Scripts and analysis order | [Analysis guide](docs/ANALYSIS_GUIDE.md) |
| Assumptions and follow-up questions | [Analysis log](docs/ANALYSIS_LOG.md) |

## Main findings

- **HLA-A:** neither headline event shows a statistically significant tumor–normal difference in BRCA. Sparse events and wide confidence intervals limit interpretation.
- **HLA-A*11:01:** a recurrent junction is consistent with a constitutive annotation mismatch. Its events are excluded from the headline comparison; unfiltered tables are retained for audit.
- **HLA-C:** headline comparisons also remain inconclusive. The classifier checks support its use in these exports, with upstream reference orientation still an unresolved assumption.
- **HLA-B:** the event classifier could not be validated with the available signal. Files prefixed `UNVALIDATED_HLA_B_` are exploratory audit outputs, excluded from the presented conclusions.
- **Pan-cancer:** comparisons using normals resolved to their cancer project remain inconclusive. “Matched” in these filenames means matched by cohort, not a patient-paired statistical analysis.

BRCA includes 935 tumor and 106 normal samples for HLA-A/B, and 933 tumor and 106 normal samples for HLA-C. See the [results summary](docs/RESULTS_SUMMARY.md) for effect sizes, uncertainty, and exclusions.

## Reproduce the analysis

The saved results can be browsed without running code. To regenerate them, use Python 3 and run the following from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python scripts/analyze_brca.py
python scripts/diagnostic_annotation_check.py
python scripts/additional_analyses.py
python scripts/phase1_bc_strand_check.py
python scripts/phase2_bc_headline.py
python scripts/phase3_pancancer_descriptive.py
python scripts/phase5_pancancer_matched_normal.py
```

These commands overwrite generated files in `outputs/` and the artifact note in `docs/`. Review the strand-check report before interpreting HLA-B/C classifications. `python scripts/analyze_brca.py --help` lists input paths, thresholds, and allele exclusions.

The dependency list reflects imports in the scripts; original package versions were not recorded, so it is not a locked reproduction environment. Upstream HLA typing, personalized reference construction, alignment, and export are outside this repository. Their limitations are documented in [Methods](docs/METHODS.md).

## Repository layout

```text
README.md                 Project overview and navigation
requirements.txt          Python dependencies
scripts/                  Analysis and figure-generation code
docs/                     Methods, findings, and analysis records
presentation/             Optional slide-figure regeneration instructions
data/                     Input exports and sample metadata
outputs/                  Generated figures, tables, and validation reports
  archive/                Previous presentation/export copies
  scratch/                Exploratory diagnostics
```

Script phase numbers retain the original analysis sequence; overview figure generation covers the presentation stage between phases 3 and 5. Detailed reports and scratch diagnostics are supporting records; begin with the results summary and indexed outputs above.
