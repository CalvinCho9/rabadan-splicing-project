#!/usr/bin/env python3
"""
make_figures_overview.py
================================================================================
Composites results into a single multi-panel page (presentation/figures_overview.pdf,
project root) for dropping straight into slides.

UPDATED (Phase 4, HLA-B/C analysis): now a two-gene overview -- HLA-A (the
original, fully-worked analysis) plus HLA-C, the only other gene Phase 1
(outputs/PHASE1_BC_strand_check.md) validated the exon5_skip/intron5_retention
classifier for. HLA-B is NOT included anywhere in this figure -- its
classifier verdict was "untrustworthy" (insufficient data, only 2 candidate
AS events total), and per the project's quarantine rule, unvalidated B/C
results never enter this figure. This is "the validated three-gene
framework, currently showing the two genes that passed validation" rather
than a true three-gene panel -- see outputs/FINAL_REVIEW.md for the
per-gene status and docs/ANALYSIS_LOG.md for why.

Five panels:
    A1. HLA-A headline classes (tumor vs normal, A*11:01 excluded) --
        unchanged from the original single-gene figure.
    A2. HLA-C headline classes (tumor vs normal, no exclusions -- Phase 1/2
        found no constitutive-artifact allele for C).
    B.  Per-allele contrast: hla_a_11_01 (constitutive artifact) vs
        hla_a_01_01 / hla_a_03_01 (real variable splicing) -- unchanged,
        HLA-A only (this is specifically the A*11:01 diagnostic story; the
        Phase 3 pan-cancer panel below is where C's real-variable-splicing
        allele and A*11:01's cross-cohort recurrence are shown instead).
    C.  Volcano, HLA-A and HLA-C combined (colored by gene) -- all recurrent
        AS_patterns from both genes' detailed tables, unfiltered.
    D.  Pan-cancer generalization of hla_a_11_01's constitutive junction
        (Phase 3) -- new panel, the same 91% pooled / 100%-in-3-of-4-cohorts
        recurrence documented in outputs/PHASE3_pancancer_descriptive.md.

Does not recompute anything the main pipelines already produced (loads their
CSVs) except the panel B slice, which deliberately needs A*11:01 visible
without the standard exclusion applied.
================================================================================
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import analyze_brca as ab

OUTDIR = Path("outputs")
TARGET_ALLELE = "hla_a_11_01"
CONTRAST_ALLELES = [TARGET_ALLELE, "hla_a_01_01", "hla_a_03_01"]
PANEL_LABEL_KW = dict(fontsize=15, fontweight="bold", ha="left", va="top")
GENE_COLOR = {"HLA-A": "#2c3e50", "HLA-C": "#d35400"}


def build_frame_a():
    df = ab.load_data([Path("data/export_BRCA_tumor_HLA-A.tsv"), Path("data/export_BRCA_normal_HLA-A.tsv")])
    df = ab.compute_derived_columns(df)
    df["Contig_rd"] = df["Contig"].str.slice(0, 11)
    events = ab.find_as_events(df)
    events = ab.annotate_events(df, events)
    events = ab.classify_events(events)
    return df, events


def panel_headline(ax, headline, gene_label, note, panel_letter):
    x = np.arange(len(headline))
    width = 0.35
    t_err = [np.clip(headline["Tumor_pct"] - headline["Tumor_CI_low"], 0, None),
             np.clip(headline["Tumor_CI_high"] - headline["Tumor_pct"], 0, None)]
    n_err = [np.clip(headline["Normal_pct"] - headline["Normal_CI_low"], 0, None),
             np.clip(headline["Normal_CI_high"] - headline["Normal_pct"], 0, None)]
    ax.bar(x - width / 2, headline["Tumor_pct"], width, yerr=t_err, capsize=4,
           label="Tumor", color=ab.TUMOR_COLOR)
    ax.bar(x + width / 2, headline["Normal_pct"], width, yerr=n_err, capsize=4,
           label="Normal", color=ab.NORMAL_COLOR)
    ax.set_xticks(x)
    ax.set_xticklabels(headline["Event_class"])
    ax.set_ylabel("Carrier samples (%, 95% Wilson CI)")
    ax.set_title(f"{gene_label} headline classes\n{note}", fontsize=10)
    ymax = max(headline["Tumor_CI_high"].max(), headline["Normal_CI_high"].max(), 1)
    ax.set_ylim(0, ymax * 1.4)
    for i, row in headline.iterrows():
        top = max(row["Tumor_CI_high"], row["Normal_CI_high"])
        ax.text(i, top + 0.05 * ymax, f"p={row['Fisher_p']:.2g}, OR={row['Odds_ratio']:.2f}\n{row['Direction']}",
                ha="center", va="bottom", fontsize=7.5)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.12, 1.18, panel_letter, transform=ax.transAxes, **PANEL_LABEL_KW)


def panel_b(fig, gs_row, per_allele_full):
    axes = [fig.add_subplot(gs_row[0, i]) for i in range(len(CONTRAST_ALLELES))]
    classes = ("exon5_skip", "intron5_retention")
    for ax, allele in zip(axes, CONTRAST_ALLELES):
        sub = per_allele_full[per_allele_full["Contig_rd"] == allele].set_index("Event_class").loc[list(classes)]
        x = np.arange(len(classes))
        width = 0.35
        t_err = [np.clip(sub["Tumor_pct"] - sub["Tumor_CI_low"], 0, None),
                 np.clip(sub["Tumor_CI_high"] - sub["Tumor_pct"], 0, None)]
        n_err = [np.clip(sub["Normal_pct"] - sub["Normal_CI_low"], 0, None),
                 np.clip(sub["Normal_CI_high"] - sub["Normal_pct"], 0, None)]
        ax.bar(x - width / 2, sub["Tumor_pct"], width, yerr=t_err, capsize=3, color=ab.TUMOR_COLOR)
        ax.bar(x + width / 2, sub["Normal_pct"], width, yerr=n_err, capsize=3, color=ab.NORMAL_COLOR)
        ax.set_xticks(x)
        ax.set_xticklabels(classes, rotation=15, ha="right", fontsize=7.5)
        tag = " (constitutive artifact)" if allele == TARGET_ALLELE else " (real variable splicing)"
        ax.set_title(f"{allele}{tag}\nn_tumor={int(sub['Tumor_n'].iloc[0])}, n_normal={int(sub['Normal_n'].iloc[0])}",
                     fontsize=8.5)
        ax.set_ylabel("Carriers (%)", fontsize=8)
        ymax = max(sub["Tumor_CI_high"].max(), sub["Normal_CI_high"].max(), 1)
        ax.set_ylim(0, ymax * 1.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles = [Patch(facecolor=ab.TUMOR_COLOR, label="Tumor"), Patch(facecolor=ab.NORMAL_COLOR, label="Normal")]
    axes[-1].legend(handles=handles, frameon=False, fontsize=7.5, loc="upper right")
    axes[0].text(-0.35, 1.32, "B", transform=axes[0].transAxes, **PANEL_LABEL_KW)
    top_y = max(ax.get_position().y1 for ax in axes)
    fig.text(0.5, top_y + 0.05,
              "HLA-A per-allele contrast: constitutive artifact vs real variable splicing",
              ha="center", fontsize=10)


def panel_volcano_multigene(ax, volcano_a, volcano_c, label_n=4):
    volcano_a = volcano_a.copy(); volcano_a["Gene"] = "HLA-A"
    volcano_c = volcano_c.copy(); volcano_c["Gene"] = "HLA-C"
    v = pd.concat([volcano_a, volcano_c], ignore_index=True)
    marker = {"HLA-A": "o", "HLA-C": "^"}
    for gene, sub in v.groupby("Gene"):
        ax.scatter(sub["log2FC_tumor_vs_normal"], sub["neglog10_Fisher_p"],
                   c=GENE_COLOR[gene], marker=marker[gene], s=45, edgecolor="black",
                   linewidth=0.3, alpha=0.85, zorder=3, label=gene)
    top = v.nlargest(label_n, "neglog10_Fisher_p").sort_values("log2FC_tumor_vs_normal")
    for i, (_, r) in enumerate(top.iterrows()):
        dy = 12 + 14 * (i % 3)
        dx = -18 if i % 2 == 0 else 18
        ax.annotate(f"{r['Gene']} {r['Contig_rd']}", (r["log2FC_tumor_vs_normal"], r["neglog10_Fisher_p"]),
                    fontsize=6, xytext=(dx, dy), textcoords="offset points", ha="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="gray"), zorder=4)
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="gray", linestyle="-", linewidth=0.6)
    ax.set_xlabel("log2(tumor freq / normal freq)", fontsize=9)
    ax.set_ylabel("-log10(Fisher p)", fontsize=9)
    ax.set_title("All recurrent AS_patterns, HLA-A + HLA-C\n(unfiltered, >=3 carriers; HLA-B not shown -- quarantined)",
                 fontsize=9.5)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ymax = v["neglog10_Fisher_p"].max()
    ax.set_ylim(-0.05 * max(ymax, 1), max(ymax, 1) * 1.35)
    ax.text(-0.15, 1.15, "C", transform=ax.transAxes, **PANEL_LABEL_KW)


def panel_pancancer_a11(ax, a11):
    a11 = a11.sort_values("hla_a_11_01_carriers", ascending=False)
    x = np.arange(len(a11))
    ax.bar(x, a11["Pct"], color="#8e44ad", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(a11["Cohort"], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("% of A*11:01 carriers with\nthe BRCA-constitutive junction", fontsize=8.5)
    ax.set_title("Pan-cancer: A*11:01 constitutive junction\nrecurs outside BRCA (descriptive)", fontsize=9.5)
    ax.set_ylim(0, 115)
    ax.axhline(100, color="gray", linestyle=":", linewidth=0.7)
    for xi, (_, r) in zip(x, a11.iterrows()):
        ax.text(xi, r["Pct"] + 3, f"{int(r['Carriers_with_BRCA_pattern'])}/{int(r['hla_a_11_01_carriers'])}",
                ha="center", fontsize=7.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.2, 1.2, "D", transform=ax.transAxes, **PANEL_LABEL_KW)


def main():
    headline_a = pd.read_csv(OUTDIR / "table_headline_classes.csv")
    headline_c = pd.read_csv(OUTDIR / "table_headline_classes_HLA-C.csv")
    volcano_a = pd.read_csv(OUTDIR / "table_volcano.csv")
    volcano_c = pd.read_csv(OUTDIR / "table_volcano_HLA-C.csv")
    a11 = pd.read_csv(OUTDIR / "table_phase3_A11_01_cross_cohort.csv")

    print("Rebuilding HLA-A per-allele data WITHOUT exclusion for panel B contrast...")
    df, events = build_frame_a()
    per_allele_full, _ = ab.per_allele_view(df, events, min_allele_samples=20, exclude_alleles=())

    fig = plt.figure(figsize=(13, 15))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.85, 1.0], hspace=0.65)
    gs_top = gs[0].subgridspec(1, 2, wspace=0.35)
    gs_mid = gs[1].subgridspec(1, 3, wspace=0.5)
    gs_bot = gs[2].subgridspec(1, 2, wspace=0.35, width_ratios=[1.15, 1])

    ax_a1 = fig.add_subplot(gs_top[0, 0])
    ax_a2 = fig.add_subplot(gs_top[0, 1])
    panel_headline(ax_a1, headline_a, "HLA-A", "(A*11:01 excluded -- see panel B)", "A1")
    panel_headline(ax_a2, headline_c, "HLA-C", "(no exclusions -- no constitutive allele found)", "A2")

    panel_b(fig, gs_mid, per_allele_full)

    ax_c = fig.add_subplot(gs_bot[0, 0])
    ax_d = fig.add_subplot(gs_bot[0, 1])
    panel_volcano_multigene(ax_c, volcano_a, volcano_c)
    panel_pancancer_a11(ax_d, a11)

    fig.suptitle("HLA-A / HLA-C alternative splicing, TCGA-BRCA + pan-cancer\n"
                 "(HLA-B not shown: classifier unvalidated, quarantined — see outputs/FINAL_REVIEW.md)",
                 fontsize=14, y=0.995)
    outpath = Path("presentation/figures_overview.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")


if __name__ == "__main__":
    main()
