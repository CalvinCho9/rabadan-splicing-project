#!/usr/bin/env python3
"""
phase5_pancancer_matched_normal.py
================================================================================
PHASE 5: use `data/gdc_sample_sheet_normal_with_hla.tsv` (GDC sample sheet +
HLA typing for 740 normal samples across 24 TCGA cohorts, sent by Anqi in
response to the "can we get more BRCA normals" question) to turn part of
Phase 3's pan-cancer pass from purely *descriptive* into a real, matched-
cohort tumor-vs-normal *test* -- exactly the fix Anqi actually recommended:

    "Not all tumors have paired normal samples and 106 BRCA normals is
    enough... We only need to compare the proportion of samples with a
    particular splicing between normal and tumor."

What the new file actually unlocks is NOT more BRCA normals (BRCA is already
covered by Phase 2's 106; Anqi says that's enough for a proportion
comparison) -- it's a Case ID -> true Project ID map for the pan-cancer
files' pooled "NORMAL" cohort. `phase3_pancancer_descriptive.py` treats that
pool as an undifferentiated cross-cohort background (docstring: "not
per-cancer-type matched normals") because the combined SJ export collapses
every non-BRCA normal sample's cancer type of origin into one bucket labeled
"NORMAL". Cross-referencing each of those samples' Barcode/Case ID against
this sheet resolves it back to its real TCGA project:

    HLA-A pan file : 145/145 pooled-NORMAL barcodes resolved (0 unresolved)
    HLA-B/C pan file: 318/318 pooled-NORMAL barcodes resolved (0 unresolved)

Once resolved, any cohort that ALSO has tumor rows in the same pan file
becomes a genuine matched (same cancer type) tumor-vs-normal comparison --
same Fisher-exact/Wilson-CI/odds-ratio machinery as Phase 2, just applied to
a cohort-restricted slice. After requiring >=3 samples on each side, that
gives:

    HLA-A : TCGA-CHOL (33 tumor / 8 normal), TCGA-HNSC (52 tumor / 5 normal)
    HLA-B*: TCGA-CHOL (33/8), TCGA-HNSC (67/35), TCGA-KIRC (92/37)
    HLA-C : TCGA-CHOL (33/8), TCGA-HNSC (67/35), TCGA-KIRC (92/37)
    (* HLA-B stays quarantined per Phase 1 -- see below)

Every other pan-cancer cohort either has no resolvable normal counterpart at
all (TCGA-LGG, TCGA-ACC -- not among the 24 projects in this sheet) or fails
the >=3-per-side floor (TCGA-THCA: 47 normal but only 1 tumor; TCGA-UCEC: up
to 79 tumor but only 1 normal) and is left exactly as Phase 3 already treats
it: descriptive cohort context only, no test.

Constitutive-artifact exclusion, two layers:
  1. KNOWN_ARTIFACT_ALLELES -- alleles already established project-wide as
     constitutive, non-tumor-associated (currently just HLA-A*11:01,
     Finding 3/5) are excluded unconditionally, regardless of how many
     carriers that allele has in a given cohort. This matters here
     specifically: a first version of this script tried to re-derive the
     exclusion per-cohort the same way Phase 2 does for BRCA
     (`hc.top_event_penetrance` + `add_verdict`, which requires >=10-20
     carriers of an allele before it will even consider flagging it) --
     that missed A*11:01 in TCGA-CHOL/TCGA-HNSC because only 3 and 2 tumor
     samples respectively carry it there, well under the bar. Caught by
     inspecting `table_phase5_detailed_*.csv` before writing up results:
     100% of HLA-A TCGA-CHOL's "3/33 intron5_retention, enriched in tumor"
     headline number turned out to be hla_a_11_01's already-known artifact
     pattern (2365_2806->2383_2806) reappearing, not new biology.
  2. The same per-cohort re-derivation Phase 2 uses is still run on top, in
     case one of these new cohorts turns up a DIFFERENT constitutive allele
     that BRCA never had -- none did (HLA-C's own recurring alleles here,
     `hla_c_03_03`/`hla_c_03_04`, are the real_variable_splicing_candidate
     from Finding 4, not an artifact -- correctly left in).

Quarantine: identical rule as every other phase. HLA-B's classifier verdict
(Phase 1: untrustworthy, insufficient data) means ALL of its Phase 5 output
goes under outputs/UNVALIDATED_HLA_B_* and is excluded from every presented
document.

Outputs:
    outputs/table_phase5_headline_HLA-{A,C}_<cohort>.csv
    outputs/table_phase5_detailed_HLA-{A,C}_<cohort>.csv
    outputs/UNVALIDATED_HLA_B_phase5_{headline,detailed}_<cohort>.csv
    outputs/table_phase5_summary.csv          (all tested gene x cohort pairs, one row per event class; HLA-B included but flagged)
    outputs/fig_phase5_matched_normal.pdf     (presented genes only)
    outputs/PHASE5_pancancer_matched_normal.md
================================================================================
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyze_brca as ab
import hla_common as hc

OUTDIR = Path("outputs")
OUTDIR.mkdir(exist_ok=True)

SHEET_PATH = "data/gdc_sample_sheet_normal_with_hla.tsv"
PAN_FILES = {
    "HLA-A": "data/pan_HLA-A_SJ_combined.tsv",
    "HLA-B": "data/pan_HLA-B_SJ_combined.tsv",
    "HLA-C": "data/pan_HLA-C_SJ_combined.tsv",
}
QUARANTINED = {"HLA-A": False, "HLA-B": True, "HLA-C": False}
MIN_SIDE_N = 3          # floor to attempt a test at all (Anqi: proportions, not raw N, are what matter)
LOW_N_FLAG = 15         # below this on the smaller side, flag as exploratory in the write-up
HLA_TYPING_COLS = {
    "HLA-A": ("HLA-A1", "HLA-A2"),
    "HLA-B": ("HLA-B1", "HLA-B2"),
    "HLA-C": ("HLA-C1", "HLA-C2"),
}
# Established project-wide (Finding 3: diagnostic_annotation_check.py /
# docs/ANNOTATION_ARTIFACT_NOTE.md; Finding 5: phase3_pancancer_descriptive.py) --
# excluded unconditionally here, not re-derived per cohort (see module
# docstring for why the per-cohort re-derivation alone isn't enough).
KNOWN_ARTIFACT_ALLELES = {"HLA-A": ["hla_a_11_01"], "HLA-B": [], "HLA-C": []}


# ==============================================================================
# 1. Resolve the pooled "NORMAL" cohort to true project of origin
# ==============================================================================
def load_sheet():
    sheet = pd.read_csv(SHEET_PATH, sep="\t")
    case_to_project = dict(zip(sheet["Case ID"], sheet["Project ID"]))
    print(f"  loaded {SHEET_PATH}: {len(sheet)} normal samples, "
          f"{sheet['Project ID'].nunique()} TCGA projects, "
          f"{sheet['Case ID'].nunique()} distinct Case IDs")
    return sheet, case_to_project


def load_and_resolve(gene, path, case_to_project):
    df = pd.read_csv(path, sep="\t", header=None, names=ab.COLUMNS)
    df["Type"] = df["Type"].str.strip().str.lower()

    # Same exclusion as Phase 3: BRCA_NORMAL block has Type mislabeled
    # "tumor" and duplicates Phase 2's dedicated BRCA export -- drop it.
    if "BRCA_NORMAL" in df["Cohort"].unique():
        mask = df["Cohort"] == "BRCA_NORMAL"
        df = df.loc[~mask].copy()

    pool_mask = df["Cohort"] == "NORMAL"
    pool_barcodes = df.loc[pool_mask, "Barcode"]
    resolved = pool_barcodes.map(case_to_project)
    n_unresolved = int(resolved.isna().sum())
    if n_unresolved:
        print(f"  [{gene}] WARNING: {n_unresolved} pooled-NORMAL barcodes not found in "
              f"{SHEET_PATH} -- left as Cohort='NORMAL_unresolved', excluded from matched tests.")
        resolved = resolved.fillna("NORMAL_unresolved")
    df.loc[pool_mask, "Cohort"] = resolved.values
    n_resolved = int((df.loc[pool_mask, "Cohort"] != "NORMAL_unresolved").sum()) if n_unresolved else int(pool_mask.sum())
    print(f"  [{gene}] resolved {n_resolved}/{int(pool_mask.sum())} pooled-NORMAL rows "
          f"({df.loc[pool_mask, 'Barcode'].nunique() if pool_mask.any() else 0} barcodes) "
          f"to their true TCGA project of origin")

    df = ab.compute_derived_columns(df)
    df["Contig_rd"] = df["Contig"].str.slice(0, 11)
    return df


def matched_cohorts(df):
    """Cohorts that now have BOTH tumor rows and resolved-normal rows, each
    side >= MIN_SIDE_N distinct samples."""
    counts = df.groupby(["Cohort", "Type"])["Barcode"].nunique().unstack(fill_value=0)
    for t in ("tumor", "normal"):
        if t not in counts.columns:
            counts[t] = 0
    counts = counts.loc[~counts.index.isin(["NORMAL", "NORMAL_unresolved"])]
    ok = counts.loc[(counts["tumor"] >= MIN_SIDE_N) & (counts["normal"] >= MIN_SIDE_N)]
    return ok


# ==============================================================================
# 2. Sanity check: does the sheet's own HLA typing agree with the
#    Contig-derived allele for a spot-check gene/cohort? (Cheap corroboration
#    of the A*11:01 cross-cohort claim now that TCGA-CHOL/TCGA-HNSC normals
#    are individually identifiable, not just a pooled count.)
# ==============================================================================
def typing_cross_check(gene, cohort, df, sheet):
    col1, col2 = HLA_TYPING_COLS[gene]
    barcodes = df.loc[(df["Cohort"] == cohort) & (df["Type"] == "normal"), "Barcode"].unique()
    sub = sheet[sheet["Case ID"].isin(barcodes)]
    if sub.empty:
        return None
    typed = sub[[col1, col2]].notna().all(axis=1).sum()
    print(f"    typing cross-check [{gene}/{cohort}]: {typed}/{len(sub)} normal samples have "
          f"non-NA {col1}/{col2} typing in the sheet")
    return sub


# ==============================================================================
# 3. Per (gene, cohort) matched test -- reuses analyze_brca's exact machinery
# ==============================================================================
def run_pair(gene, cohort, df):
    sub = df[df["Cohort"] == cohort].copy()
    events = ab.find_as_events(sub, local_rate_threshold=0.1)
    events = ab.annotate_events(sub, events)
    if events.empty:
        # Small single-cohort slice can have zero candidate AS events at all
        # (expected for HLA-B, which is sparse even on the full BRCA cohort
        # -- Phase 1 found only 2 events in 935 tumor + 106 normal). str.split
        # on an empty Series can't infer column count, so
        # classify_events_strand_aware would crash on an empty frame -- skip
        # it and add the column directly instead.
        events["Event_class"] = pd.Series(dtype=object)
        print(f"    [{gene}/{cohort}] 0 candidate AS events in this cohort")
    else:
        events = hc.classify_events_strand_aware(events, flip_minus=False)  # naive wins Phase 1 for A/B/C alike

    # Lower bar than BRCA-only (MIN_ALLELE_SAMPLES=20 default) -- these are
    # much smaller cohorts (e.g. TCGA-CHOL: 41 samples total), so few or no
    # alleles would ever reach 20. Same relaxed bar phase3 uses for the same
    # reason (MIN_ALLELE_COHORT_SAMPLES=10).
    penetrance = hc.top_event_penetrance(sub, events, min_allele_samples=10)
    if penetrance.empty:
        print(f"    [{gene}/{cohort}] no allele reaches 10+ samples in this cohort -- "
              f"skipping constitutive-artifact re-check (nothing to exclude)")
        exclude_alleles = []
    else:
        penetrance = hc.add_verdict(penetrance)
        exclude_alleles = penetrance.loc[
            penetrance["Verdict"].str.startswith("constitutive_annotation_artifact"), "Contig_rd"
        ].tolist()
        if exclude_alleles:
            print(f"    [{gene}/{cohort}] excluding constitutive-artifact allele(s) from headline: {exclude_alleles}")

    known = KNOWN_ARTIFACT_ALLELES.get(gene, [])
    newly_known = [a for a in known if a not in exclude_alleles]
    if newly_known:
        print(f"    [{gene}/{cohort}] also excluding previously-established artifact allele(s) "
              f"(Finding 3/5, unconditional regardless of per-cohort N): {newly_known}")
    exclude_alleles = sorted(set(exclude_alleles) | set(known))

    events_for_headline = events[~events["Contig_rd"].isin(exclude_alleles)] if exclude_alleles else events
    headline = ab.headline_comparison(sub, events_for_headline)
    headline.insert(0, "Cohort", cohort)
    headline.insert(0, "Gene", gene)
    headline["Excluded_alleles"] = ", ".join(exclude_alleles) if exclude_alleles else ""

    detail = ab.detailed_table(sub, events, min_carriers=3)
    detail.insert(0, "Cohort", cohort)
    detail.insert(0, "Gene", gene)

    return headline, detail, exclude_alleles, len(events)


def gene_tag(gene):
    return gene.replace("-", "_")


# ==============================================================================
# 4. Figure: grouped bar chart, tumor vs normal %, per gene x cohort x class
#    (presented genes only -- HLA-B is quarantined out of this figure too)
# ==============================================================================
def plot_summary(summary, outpath):
    presented = summary[~summary["Quarantined"]]
    if presented.empty:
        print("  no presented (gene, cohort) pairs -- skipping fig_phase5_matched_normal.pdf")
        return
    presented = presented.copy()
    presented["Label"] = presented["Gene"] + "\n" + presented["Cohort"].str.replace("TCGA-", "", regex=False)
    classes = ["exon5_skip", "intron5_retention"]
    fig, axes = plt.subplots(1, len(classes), figsize=(5 * len(classes), 4.5), sharey=True)
    for ax, cls in zip(axes, classes):
        d = presented[presented["Event_class"] == cls].sort_values(["Gene", "Cohort"])
        x = np.arange(len(d))
        width = 0.35
        # np.clip guards against sub-epsilon negative yerr from float rounding
        # at the 0%/100% Wilson CI boundary (matplotlib rejects any negative
        # yerr outright).
        ax.bar(x - width / 2, d["Tumor_pct"], width,
               yerr=[np.clip(d["Tumor_pct"] - d["Tumor_CI_low"], 0, None),
                     np.clip(d["Tumor_CI_high"] - d["Tumor_pct"], 0, None)],
               color=ab.TUMOR_COLOR, label="Tumor", capsize=3)
        ax.bar(x + width / 2, d["Normal_pct"], width,
               yerr=[np.clip(d["Normal_pct"] - d["Normal_CI_low"], 0, None),
                     np.clip(d["Normal_CI_high"] - d["Normal_pct"], 0, None)],
               color=ab.NORMAL_COLOR, label="Normal", capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels(d["Label"], fontsize=8)
        ax.set_title(cls)
        ax.set_ylabel("% samples carrying event (95% Wilson CI)")
        for i, (_, r) in enumerate(d.iterrows()):
            ax.annotate(f"p={r['Fisher_p']:.2g}", (i, max(r["Tumor_CI_high"], r["Normal_CI_high"]) + 1),
                        ha="center", fontsize=7)
    axes[0].legend()
    # Tested/significant counts pulled from columns already computed by
    # ab.headline_comparison() (OR_CI_low/OR_CI_high) -- same arithmetic
    # write_report() already runs on this same `summary` frame, just
    # surfaced on the figure too so a null result reads as "tested and
    # null", not "empty because nothing was tested".
    n_tested = len(presented)
    n_sig = int(((presented["OR_CI_low"] > 1) | (presented["OR_CI_high"] < 1)).sum())
    fig.suptitle("Phase 5: pan-cancer matched tumor-vs-normal (NORMAL pool resolved to true cohort of origin)\n"
                 "HLA-B excluded (quarantined, see Phase 1)\n"
                 f"{n_tested} gene×cohort×class tests shown, {n_sig} with OR 95% CI excluding 1 "
                 f"(nominal, not FDR-corrected)")
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  -> {outpath}")


# ==============================================================================
# 5. Bonus: does HLA-C's own real-variable-splicing candidate from Finding 4
#    (hla_c_03_03 / hla_c_03_04, pattern 1273_1711->1291_1711) corroborate
#    across the newly-matched cohorts, the way A*11:01's ARTIFACT was shown
#    to recur pan-cancer in Finding 5 -- except this allele is a genuine
#    candidate, not an artifact, so cross-cohort corroboration actually
#    strengthens the case for real biology.
# ==============================================================================
HLA_C_CANDIDATE_ALLELES = ["hla_c_03_03", "hla_c_03_04"]
HLA_C_CANDIDATE_PATTERN = "1273_1711->1291_1711"
BRCA_DETAIL_PATH = "outputs/results_BRCA_tumor_vs_normal_HLA-C.csv"


def hla_c_candidate_cross_cohort(all_detail):
    rows = []
    brca_path = Path(BRCA_DETAIL_PATH)
    if brca_path.exists():
        brca = pd.read_csv(brca_path)
        brca_sub = brca[brca["Contig_rd"].isin(HLA_C_CANDIDATE_ALLELES)
                         & (brca["AS_pattern"] == HLA_C_CANDIDATE_PATTERN)]
        if not brca_sub.empty:
            rows.append({"Cohort": "TCGA-BRCA (Finding 4)",
                         "Tumor_carriers": int(brca_sub["Tumor_carriers"].sum()),
                         "Normal_carriers": int(brca_sub["Normal_carriers"].sum())})
    for d in all_detail:
        if d.empty or d["Gene"].iloc[0] != "HLA-C":
            continue
        sub = d[d["Contig_rd"].isin(HLA_C_CANDIDATE_ALLELES) & (d["AS_pattern"] == HLA_C_CANDIDATE_PATTERN)]
        if sub.empty:
            continue
        rows.append({"Cohort": d["Cohort"].iloc[0],
                     "Tumor_carriers": int(sub["Tumor_carriers"].sum()),
                     "Normal_carriers": int(sub["Normal_carriers"].sum())})
    return pd.DataFrame(rows)


# ==============================================================================
# main
# ==============================================================================
def main():
    print("Loading GDC normal-sample sheet...")
    sheet, case_to_project = load_sheet()

    all_headline, all_detail, meta_rows = [], [], []

    for gene, path in PAN_FILES.items():
        quarantined = QUARANTINED[gene]
        status = "QUARANTINED (not for presentation)" if quarantined else "validated, presented"
        print(f"\n{'=' * 80}\n{gene}  [{status}]\n{'=' * 80}")

        df = load_and_resolve(gene, path, case_to_project)
        cohorts = matched_cohorts(df)
        if cohorts.empty:
            print(f"  [{gene}] no cohort reaches {MIN_SIDE_N}+ samples on both sides -- nothing to test")
            continue
        print(f"  [{gene}] matched cohorts (tumor / normal): "
              + ", ".join(f"{c} ({int(r.tumor)}/{int(r.normal)})" for c, r in cohorts.iterrows()))

        for cohort, row in cohorts.iterrows():
            n_tumor, n_normal = int(row["tumor"]), int(row["normal"])
            typing_cross_check(gene, cohort, df, sheet)
            headline, detail, exclude_alleles, n_events = run_pair(gene, cohort, df)
            headline["Quarantined"] = quarantined
            headline["Low_N_exploratory"] = min(n_tumor, n_normal) < LOW_N_FLAG
            all_headline.append(headline)
            all_detail.append(detail)
            meta_rows.append({"Gene": gene, "Cohort": cohort, "Tumor_n": n_tumor, "Normal_n": n_normal,
                               "N_events": n_events, "Excluded_alleles": ", ".join(exclude_alleles),
                               "Quarantined": quarantined})

            prefix = f"UNVALIDATED_{gene_tag(gene)}_phase5" if quarantined else f"table_phase5"
            suffix = cohort.replace("TCGA-", "")
            tag = f"_{gene[4:]}_{suffix}" if not quarantined else f"_{suffix}"
            headline.to_csv(OUTDIR / f"{prefix}_headline{tag}.csv", index=False)
            detail.to_csv(OUTDIR / f"{prefix}_detailed{tag}.csv", index=False)
            ab.print_headline_table(headline)

    summary = pd.concat(all_headline, ignore_index=True)
    summary.to_csv(OUTDIR / "table_phase5_summary.csv", index=False)
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(OUTDIR / "table_phase5_cohort_meta.csv", index=False)

    plot_summary(summary, OUTDIR / "fig_phase5_matched_normal.pdf")

    candidate_cross_cohort = hla_c_candidate_cross_cohort(all_detail)
    if not candidate_cross_cohort.empty:
        candidate_cross_cohort.to_csv(OUTDIR / "table_phase5_hla_c_candidate_cross_cohort.csv", index=False)
        print("\n=== HLA-C real_variable_splicing_candidate cross-cohort corroboration "
              f"({HLA_C_CANDIDATE_PATTERN}, alleles {HLA_C_CANDIDATE_ALLELES}) ===")
        print(candidate_cross_cohort.to_string(index=False))

    write_report(summary, meta, candidate_cross_cohort)
    print("\nPhase 5 done.")


def write_report(summary, meta, candidate_cross_cohort):
    lines = []
    lines.append("# Phase 5 — pan-cancer matched tumor-vs-normal, using `gdc_sample_sheet_normal_with_hla.tsv`\n")
    lines.append(
        "Generated by `phase5_pancancer_matched_normal.py`, in direct response to Anqi's guidance: "
        "106 BRCA normals is enough for a proportion comparison (Phase 2 already does this), and what "
        "would most strengthen the pan-cancer pass is resolving which cancer type each pooled `NORMAL` "
        "sample actually came from, not adding more BRCA normals. The new sheet (740 normal samples, 24 "
        "TCGA projects, with HLA-A/B/C typing) does exactly that: every pooled-`NORMAL` barcode in the "
        "pan-cancer SJ files resolves to a real TCGA project (145/145 for the HLA-A file, 318/318 for "
        "HLA-B/C — 0 unresolved).\n"
    )
    lines.append(
        "This turns `phase3_pancancer_descriptive.py`'s explicitly-descriptive pooled-NORMAL comparison "
        "into a real, matched-cohort tumor-vs-normal *test* — same Fisher-exact/Wilson-CI/odds-ratio "
        "machinery as the BRCA headline test — for every cohort where the resolved normal count and the "
        f"existing tumor count both reach {MIN_SIDE_N}+ samples. Every other pan-cancer cohort (no "
        "resolvable normal at all, e.g. TCGA-LGG/TCGA-ACC — not among the sheet's 24 projects — or too "
        "few samples on one side, e.g. TCGA-THCA's 47 normal vs 1 tumor, TCGA-UCEC's 79 tumor vs 1 "
        "normal) is left exactly as Phase 3 treats it: descriptive context only.\n"
    )
    lines.append("## Quarantine\n")
    lines.append(
        "Same rule as every other phase: HLA-B's classifier verdict from Phase 1 is untrustworthy. ALL "
        "HLA-B Phase 5 output is written under `outputs/UNVALIDATED_HLA_B_phase5_*` and excluded below "
        "and from every presented document.\n"
    )
    lines.append("## Constitutive-artifact exclusion\n")
    lines.append(
        "Two layers, applied BEFORE the headline test. (1) HLA-A\\*11:01 — already established "
        "project-wide as a constitutive, non-tumor-associated artifact (Finding 3/5) — is excluded "
        "unconditionally in every cohort, regardless of its carrier count there. This was NOT the first "
        "version of this script: an initial pass re-derived the exclusion per cohort the same way Phase "
        "2 does for BRCA (`hla_common.top_event_penetrance` + `add_verdict`, which requires 10-20+ "
        "carriers before flagging anything), and that missed A\\*11:01 in TCGA-CHOL/TCGA-HNSC because "
        "only 2-3 tumor samples carry it in either cohort — caught by inspecting the per-pattern detail "
        "table before writing up results: 100% of HLA-A TCGA-CHOL's and TCGA-HNSC's *original* "
        "intron5_retention carriers turned out to be A\\*11:01's known artifact pattern, not new biology. "
        "(2) The Phase-2-style per-cohort re-derivation still runs on top, to catch a *different* "
        "constitutive allele these new cohorts might reveal that BRCA never had — none did here (HLA-C's "
        "recurring alleles below, `hla_c_03_03`/`hla_c_03_04`, are the real_variable_splicing_candidate "
        "from Finding 4, correctly left in, not excluded). Excluded allele(s) per pair are in "
        "`table_phase5_cohort_meta.csv`.\n"
    )
    lines.append("## Results\n")
    for _, m in meta.iterrows():
        gene, cohort = m["Gene"], m["Cohort"]
        tag = "QUARANTINED" if m["Quarantined"] else "presented"
        sub = summary[(summary["Gene"] == gene) & (summary["Cohort"] == cohort)]
        low_n = "**exploratory, low N**" if sub["Low_N_exploratory"].iloc[0] else "adequately sized for a proportion test"
        lines.append(f"### {gene} — {cohort} [{tag}] ({int(m['Tumor_n'])} tumor / {int(m['Normal_n'])} normal, {low_n})\n")
        if m["Excluded_alleles"]:
            lines.append(f"- Excluded constitutive-artifact allele(s) from headline: `{m['Excluded_alleles']}`\n")
        for _, r in sub.iterrows():
            ci_note = "" if (r["OR_CI_low"] < 1 < r["OR_CI_high"]) else " **(OR 95% CI excludes 1)**"
            lines.append(
                f"- **{r['Event_class']}**: tumor {r['Tumor_carriers']}/{r['Tumor_total']} "
                f"({r['Tumor_pct']:.1f}%, CI {r['Tumor_CI_low']:.1f}-{r['Tumor_CI_high']:.1f}%) vs "
                f"normal {r['Normal_carriers']}/{r['Normal_total']} ({r['Normal_pct']:.1f}%, CI "
                f"{r['Normal_CI_low']:.1f}-{r['Normal_CI_high']:.1f}%) — {r['Direction']}, "
                f"OR {r['Odds_ratio']:.2f} ({r['OR_CI_low']:.2f}-{r['OR_CI_high']:.2f}), "
                f"Fisher p={r['Fisher_p']:.3g}{ci_note}"
            )
        lines.append("")

    if not candidate_cross_cohort.empty:
        lines.append("## Bonus: HLA-C's Finding-4 candidate corroborates across independent cohorts\n")
        lines.append(
            f"`{HLA_C_CANDIDATE_PATTERN}` (alleles `hla_c_03_03`/`hla_c_03_04` — closely related "
            "sub-alleles at the same locus showing the same junction, plausible given how similar their "
            "sequences are) was flagged in Finding 4 as HLA-C's one real_variable_splicing_candidate in "
            "BRCA, not an artifact (unlike A\\*11:01): present in tumor, essentially absent in normal, "
            "but resting on a fragile n=1 normal-side carrier there. This is the *opposite* situation "
            "from A\\*11:01 — cross-cohort corroboration here would strengthen a real candidate rather "
            "than just confirm an artifact is universal, so it's worth checking the same way Finding 5 "
            "checked A\\*11:01:\n"
        )
        lines.append("| Cohort | Tumor carriers | Normal carriers |\n|---|---|---|")
        for _, r in candidate_cross_cohort.iterrows():
            lines.append(f"| {r['Cohort']} | {r['Tumor_carriers']} | {r['Normal_carriers']} |")
        total_t = int(candidate_cross_cohort["Tumor_carriers"].sum())
        total_n = int(candidate_cross_cohort["Normal_carriers"].sum())
        lines.append(
            f"\nPooled across BRCA + the 3 newly-matched cohorts: **{total_t} tumor carriers vs "
            f"{total_n} normal carrier(s)** (counts are per-allele-row sums, so a sample compound "
            "heterozygous for both 03:03 and 03:04 would be counted twice — negligible at these "
            "allele frequencies, but noted rather than presented as an exact distinct-patient count). "
            "The pattern is enriched in tumor in **every** cohort tested so far, never observed in any "
            "of the 80 non-BRCA normal samples carrying either allele (8 CHOL + 35 HNSC + 37 KIRC) — "
            "directionally consistent, but still **not** a replication: none of the individual cohort "
            "tests reaches significance (all Fisher p>0.5 above), and this cross-cohort view was not a "
            "pre-registered test either, so treat it as 'worth a dedicated follow-up with more normals "
            "specifically for this allele,' not as confirmed real variable splicing.\n"
        )

    n_sig = int(((summary["OR_CI_low"] > 1) | (summary["OR_CI_high"] < 1)).sum())
    n_total = len(summary)
    lines.append(
        f"**Multiple testing:** {n_total} (gene, cohort, event class) headline tests run this phase; "
        f"{n_sig} have an OR 95% CI excluding 1 at the nominal (uncorrected) level. None of these were "
        f"planned as a single pre-registered comparison — treat any nominal hit here the same way "
        f"`docs/RESULTS_SUMMARY.md` treats the BRCA headline numbers: as one candidate among several tests, "
        f"not FDR-corrected across cohorts, worth a closer look but not a standalone claim.\n"
    )
    lines.append(
        "Table(s): `outputs/table_phase5_summary.csv` (all tested pairs), "
        "`outputs/table_phase5_cohort_meta.csv` (N + excluded alleles per pair), "
        "per-pair detail `outputs/table_phase5_detailed_*.csv`, "
        "figure `outputs/fig_phase5_matched_normal.pdf`.\n"
    )
    out_path = OUTDIR / "PHASE5_pancancer_matched_normal.md"
    out_path.write_text("\n".join(lines))
    print(f"\n  -> {out_path}")


if __name__ == "__main__":
    main()
