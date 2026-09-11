#!/usr/bin/env python3
"""
make_figures_overview.py
================================================================================
Composites three already-computed results into a single multi-panel page
(figures_overview.pdf, project root) for dropping straight into slides:

    A. Headline classes (tumor vs normal, --exclude-alleles applied, from
       outputs/table_headline_classes.csv -- this is the presented number)
    B. Per-allele contrast: hla_a_11_01 (constitutive artifact, excluded from
       the headline) side by side with the two alleles that show genuine
       variable splicing (hla_a_01_01, hla_a_03_01) -- recomputed WITHOUT the
       exclusion so A*11:01's near-100%-in-both-groups pattern is visible for
       the diagnostic story
    C. Volcano across all recurrent AS_patterns (outputs/table_volcano.csv,
       unfiltered -- descriptive, not the headline)

Does not recompute anything the main pipeline already produced (loads its
CSVs) except the panel B slice, which deliberately needs A*11:01 visible.
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


def build_frame():
    df = ab.load_data([Path("data/export_BRCA_tumor_HLA-A.tsv"), Path("data/export_BRCA_normal_HLA-A.tsv")])
    df = ab.compute_derived_columns(df)
    df["Contig_rd"] = df["Contig"].str.slice(0, 11)
    events = ab.find_as_events(df)
    events = ab.annotate_events(df, events)
    events = ab.classify_events(events)
    return df, events


def panel_a(ax, headline):
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
    ax.set_title("Tumor vs normal, headline classes\n(A*11:01 excluded -- see panel B)", fontsize=10)
    ymax = max(headline["Tumor_CI_high"].max(), headline["Normal_CI_high"].max(), 1)
    ax.set_ylim(0, ymax * 1.4)
    for i, row in headline.iterrows():
        top = max(row["Tumor_CI_high"], row["Normal_CI_high"])
        ax.text(i, top + 0.05 * ymax, f"p={row['Fisher_p']:.2g}, OR={row['Odds_ratio']:.2f}\n{row['Direction']}",
                ha="center", va="bottom", fontsize=7.5)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.12, 1.12, "A", transform=ax.transAxes, **PANEL_LABEL_KW)


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
    axes[0].text(-0.35, 1.28, "B", transform=axes[0].transAxes, **PANEL_LABEL_KW)
    top_y = max(ax.get_position().y1 for ax in axes)
    fig.text(0.5, top_y + 0.045,
              "Per-allele contrast: constitutive artifact vs real variable splicing",
              ha="center", fontsize=10)


def panel_c(ax, volcano, label_n=6):
    class_color = {"exon5_skip": "#8e44ad", "intron5_retention": "#16a085", "other": "#999999"}
    colors = volcano["Event_class"].map(class_color).fillna("#999999")
    ax.scatter(volcano["log2FC_tumor_vs_normal"], volcano["neglog10_Fisher_p"],
               c=colors, s=40, edgecolor="black", linewidth=0.3, alpha=0.85, zorder=3)
    top = volcano.nlargest(label_n, "neglog10_Fisher_p").sort_values("log2FC_tumor_vs_normal")
    for i, (_, r) in enumerate(top.iterrows()):
        dy = 12 + 14 * (i % 3)
        dx = -18 if i % 2 == 0 else 18
        ax.annotate(f"{r['Contig_rd']}", (r["log2FC_tumor_vs_normal"], r["neglog10_Fisher_p"]),
                    fontsize=6, xytext=(dx, dy), textcoords="offset points", ha="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="gray"), zorder=4)
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="gray", linestyle="-", linewidth=0.6)
    ax.set_xlabel("log2(tumor freq / normal freq)", fontsize=9)
    ax.set_ylabel("-log10(Fisher p)", fontsize=9)
    ax.set_title("All recurrent AS_patterns\n(unfiltered, >=3 carriers)", fontsize=10)
    handles = [Patch(facecolor=c, label=k) for k, c in class_color.items()]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ymax = volcano["neglog10_Fisher_p"].max()
    ax.set_ylim(-0.05 * max(ymax, 1), max(ymax, 1) * 1.35)
    ax.text(-0.15, 1.12, "C", transform=ax.transAxes, **PANEL_LABEL_KW)


def main():
    headline = pd.read_csv(OUTDIR / "table_headline_classes.csv")
    volcano = pd.read_csv(OUTDIR / "table_volcano.csv")

    print("Rebuilding per-allele data WITHOUT exclusion for panel B contrast...")
    df, events = build_frame()
    per_allele_full, _ = ab.per_allele_view(df, events, min_allele_samples=20, exclude_alleles=())

    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.85], hspace=0.55)
    gs_top = gs[0].subgridspec(1, 2, wspace=0.35, width_ratios=[1, 1.1])
    gs_bottom = gs[1].subgridspec(1, 3, wspace=0.5)

    ax_a = fig.add_subplot(gs_top[0, 0])
    ax_c = fig.add_subplot(gs_top[0, 1])
    panel_a(ax_a, headline)
    panel_c(ax_c, volcano)
    panel_b(fig, gs_bottom, per_allele_full)

    fig.suptitle("HLA-A alternative splicing, TCGA-BRCA: tumor vs normal", fontsize=15, y=0.985)
    outpath = Path("figures_overview.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")


if __name__ == "__main__":
    main()
