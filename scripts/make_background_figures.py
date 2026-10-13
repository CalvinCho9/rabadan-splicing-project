#!/usr/bin/env python3
"""
make_background_figures.py
================================================================================
Explanatory (not data-driven) graphics for the background/method section of
a talk: what the pipeline actually does (panel A) and what "exon5_skip" /
"intron5_retention" geometrically mean (panel B). Pure schematic -- no real
numbers are plotted here; see presentation/figures_overview.pdf / outputs/ for results.

Output: presentation/background_method.pdf (presentation folder)
================================================================================
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

TUMOR_COLOR = "#c0392b"
NORMAL_COLOR = "#2980b9"
EXON_SKIP_COLOR = "#8e44ad"
INTRON_RET_COLOR = "#16a085"
EXON_FACE = "#d9d9d9"
EXON_EDGE = "#333333"
PANEL_LABEL_KW = dict(fontsize=17, fontweight="bold", ha="left", va="top")


# ==============================================================================
# Panel A: pipeline workflow
# ==============================================================================
STEPS = [
    "STAR SJ.out.tab per sample x allele\n(2 files/sample: allele 1, allele 2)",
    "Concatenate + tag every row with\n[Cohort, Type, Barcode, Allele]",
    "Bucket rows by Annotated_cumsum\n(group each novel junction under the\nannotated junction it sits downstream of)",
    "Local_rate = junction's reads /\ntotal reads in its bucket",
    "Flag AS event:\nAnnotated==0  &  Local_rate > 10%  &  bucket != 0",
    "Classify by coordinate geometry:\nexon5_skip / intron5_retention / other",
    "Compare carrier frequency, tumor vs normal:\nWilson CI, Fisher exact, BH-FDR",
]


def panel_a(ax):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, len(STEPS) + 0.5)
    ax.axis("off")
    box_w, box_h = 9, 0.82
    ys = [len(STEPS) - i for i in range(len(STEPS))]
    for step, y in zip(STEPS, ys):
        box = FancyBboxPatch((0.5, y - box_h / 2), box_w, box_h,
                              boxstyle="round,pad=0.08,rounding_size=0.12",
                              facecolor="#eef2f7", edgecolor="#2c3e50", linewidth=1.2)
        ax.add_patch(box)
        ax.text(5, y, step, ha="center", va="center", fontsize=8.3)
        if y != ys[-1]:
            ax.add_patch(FancyArrowPatch((5, y - box_h / 2 - 0.02), (5, y - 1 + box_h / 2 + 0.05),
                                          arrowstyle="-|>", mutation_scale=14, color="#2c3e50", linewidth=1.2))
    ax.text(-0.02, 1.02, "A", transform=ax.transAxes, **PANEL_LABEL_KW)
    ax.set_title("Pipeline: from raw STAR output to a tested carrier-frequency\ndifference",
                 fontsize=11, pad=14)


# ==============================================================================
# Panel B: splicing-event geometry cartoon
# ==============================================================================
def exon(ax, x0, x1, y, face=EXON_FACE, edge=EXON_EDGE, hatch=None, label=None, height=0.5):
    ax.add_patch(Rectangle((x0, y - height / 2), x1 - x0, height,
                            facecolor=face, edgecolor=edge, linewidth=1.3, hatch=hatch, zorder=3))
    if label:
        ax.text((x0 + x1) / 2, y, label, ha="center", va="center", fontsize=8, zorder=4)


def intron_line(ax, x0, x1, y):
    ax.plot([x0, x1], [y, y], color="#888888", linewidth=1.2, zorder=1)


def splice_arc(ax, x0, x1, y, rad, color, lw=2.2, style="-"):
    arc = FancyArrowPatch((x0, y + 0.28), (x1, y + 0.28),
                           connectionstyle=f"arc3,rad={rad}", arrowstyle="-",
                           linestyle=style, linewidth=lw, color=color, zorder=5)
    ax.add_patch(arc)


def panel_b(ax):
    ax.set_xlim(-0.3, 13.3)
    ax.set_ylim(0, 10.5)
    ax.axis("off")
    ax.text(-0.02, 1.02, "B", transform=ax.transAxes, **PANEL_LABEL_KW)
    ax.set_title("What the two named events look like\n(coordinates: W = annotated junction, V = variant junction)",
                 fontsize=11, pad=14)

    # ---- Row 1: annotated / reference splicing --------------------------
    y1 = 8.1
    ax.text(0, y1 + 1.75, "Annotated (reference) splicing", fontsize=9.5, fontweight="bold")
    exon(ax, 0, 3, y1, label="Exon N")
    intron_line(ax, 3, 6, y1)
    exon(ax, 6, 9, y1, label="Exon N+1")
    exon(ax, 10.5, 13, y1, face="#f2f2f2", edge="#bbbbbb", label="Exon N+2")
    intron_line(ax, 9, 10.5, y1)
    splice_arc(ax, 3, 6, y1, rad=-0.6, color="#333333", style="--")
    ax.text(4.5, y1 + 1.1, "W1 → W2", ha="center", fontsize=8.5, color="#333333")
    ax.text(3, y1 - 0.55, "W1", ha="center", fontsize=7.5)
    ax.text(6, y1 - 0.55, "W2", ha="center", fontsize=7.5)

    # ---- Row 2: exon5_skip -----------------------------------------------
    y2 = 4.9
    ax.text(0, y2 + 1.75, "exon5_skip: same donor, acceptor jumps far downstream "
                            "(V1=W1, V2−W2≥100bp) → Exon N+1 skipped entirely",
            fontsize=9.5, fontweight="bold", color=EXON_SKIP_COLOR)
    exon(ax, 0, 3, y2, label="Exon N")
    intron_line(ax, 3, 6, y2)
    exon(ax, 6, 9, y2, face="#f2f2f2", edge="#bbbbbb", hatch="////", label="Exon N+1\n(skipped)")
    intron_line(ax, 9, 10.5, y2)
    exon(ax, 10.5, 13, y2, label="Exon N+2")
    splice_arc(ax, 3, 10.5, y2, rad=-0.45, color=EXON_SKIP_COLOR)
    ax.text(6.75, y2 + 1.1, "V1 → V2  (=W1 → far downstream)", ha="center", fontsize=8.5,
            color=EXON_SKIP_COLOR)
    ax.text(3, y2 - 0.55, "V1=W1", ha="center", fontsize=7.5)
    ax.text(10.5, y2 - 0.55, "V2", ha="center", fontsize=7.5)

    # ---- Row 3: intron5_retention -----------------------------------------
    y3 = 1.7
    ax.text(0, y3 + 1.75, "intron5_retention: same acceptor, donor shifts a few bp downstream "
                            "(V2=W2, 0<V1−W1≤50bp) → partial intron retained",
            fontsize=9.5, fontweight="bold", color=INTRON_RET_COLOR)
    exon(ax, 0, 3, y3, label="Exon N")
    exon(ax, 3, 3.6, y3, face="#cde9e2", edge=INTRON_RET_COLOR, label=None, height=0.5)
    ax.text(3.3, y3 - 0.85, "retained\n(≤50bp)", ha="center", fontsize=6.8, color=INTRON_RET_COLOR)
    intron_line(ax, 3.6, 6, y3)
    exon(ax, 6, 9, y3, label="Exon N+1")
    splice_arc(ax, 3.6, 6, y3, rad=-0.6, color=INTRON_RET_COLOR)
    ax.text(4.8, y3 + 1.1, "V1 → W2", ha="center", fontsize=8.5, color=INTRON_RET_COLOR)
    ax.text(3.6, y3 - 0.55, "V1", ha="center", fontsize=7.5)
    ax.text(6, y3 - 0.55, "W2=V2", ha="center", fontsize=7.5)

    legend_elems = [
        Line2D([0], [0], color="#333333", lw=2, linestyle="--", label="annotated (reference) junction"),
        Line2D([0], [0], color=EXON_SKIP_COLOR, lw=2, label="exon5_skip variant junction"),
        Line2D([0], [0], color=INTRON_RET_COLOR, lw=2, label="intron5_retention variant junction"),
    ]
    ax.legend(handles=legend_elems, loc="lower center", bbox_to_anchor=(0.5, -0.06),
              ncol=1, frameon=False, fontsize=8)


def main():
    fig = plt.figure(figsize=(15, 8.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.3], wspace=0.15)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    panel_a(ax_a)
    panel_b(ax_b)
    fig.suptitle("Background & method: how an HLA-A splice-junction file becomes a tumor-vs-normal test",
                 fontsize=14, y=1.01)
    outpath = Path("presentation/background_method.pdf")
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outpath}")


if __name__ == "__main__":
    main()
