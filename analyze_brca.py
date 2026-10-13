#!/usr/bin/env python3
"""
analyze_brca.py
================================================================================
Alternative-splicing comparison, tumor vs normal, for HLA-A (TCGA-BRCA).

See README_FOR_CLAUDE.md for the canonical folder layout:
    data/       input TSVs (do not edit)
    outputs/    everything this script generates
    analyze_brca.py   this file

Reproduces (vectorized, in pandas) the per-sample splicing math from
plot_HLA-A_BRCA.R / Anqi's original analysis, then tests two specific events
called out in the MHC Hammer paper (Nature Genetics):

    - exon5_skip          : complete exon-5 skipping   -> expected DEPLETED in tumor
    - intron5_retention   : partial intron-5 retention  -> expected ENRICHED in tumor

...and extends that headline comparison with a full presentation-ready
analysis: cohort overview, confidence intervals + effect sizes, an
allele-stratified view, a volcano plot across all recurrent AS_patterns,
per-sample event burden, and Benjamini-Hochberg FDR correction for the
many-pattern comparison.

Input
-----
Two headerless, tab-separated STAR SJ.out.tab-derived files, 13 columns each:

    Cohort  Type  Barcode  Allele  Contig  First_base  Last_base  Strand
    Motif  Annotated  Uniquely_mapping  Multimapping  Overhang

`Type` is "tumor" or "normal". A patient's tumor and normal sample can share
the same Barcode, so every per-sample grouping in this script keys on
[Type, Barcode, Allele] -- never on [Barcode, Allele] alone.

Output (all written into --outdir, default "outputs/")
----------------------------------------------------------------------
    table_cohort_overview.csv        + fig_cohort_overview.pdf
    table_events_per_sample.csv      (raw per-sample event counts; feeds burden plot)
    table_headline_classes.csv       + fig_classes.pdf              (presented: --exclude-alleles applied)
    table_headline_classes_all_alleles.csv                          (audit copy: no exclusions)
    table_per_allele.csv             + fig_per_allele.pdf           (allele-stratified view, excluded alleles dropped)
    table_volcano.csv                + fig_volcano.pdf              (all recurrent AS_patterns, no exclusions)
    fig_event_burden.pdf                                            (Mann-Whitney U)
    results_BRCA_tumor_vs_normal_HLA-A.csv                          (detailed table + FDR q-values, no exclusions)
    fig_heatmap_tumor.pdf / fig_heatmap_normal.pdf                  (samples x alleles, no exclusions)

--exclude-alleles (default "hla_a_11_01") drops an allele's events from the
tumor-vs-normal HEADLINE and PER-ALLELE comparisons only -- not from any
descriptive/overview output above. Use this once diagnostic_annotation_check.py
has shown an allele's top event is ~universal in both tumor and normal
carriers (constitutive), rather than tumor-associated. Both the excluded and
unfiltered headline numbers print to console and save to CSV so the effect of
exclusion is auditable, not hidden.

Extending to other inputs
--------------------------------
    - Pan-cancer file (data/pan_HLA-A_SJ_combined.tsv): same 13-column
      schema and already spans multiple Cohorts. Point --tumor/--normal at
      it (or split by Type first) -- everything below groups by
      Type/Barcode/Allele, not by Cohort, so multi-cohort input works
      unchanged. To break tumor-vs-normal out per Cohort too, add "Cohort"
      to the groupby keys in headline_comparison() / per_allele_view().
    - HLA-B / HLA-C: point --tumor/--normal at the new files. The
      exon5_skip / intron5_retention geometry thresholds
      (--exon-skip-min-jump / --intron-retention-max-shift) are coordinate
      offsets, not gene-specific constants -- re-check them against that
      gene's exon 5 length before trusting the labels. --min-allele-samples
      and --min-carriers are also CLI flags for the same reason: a smaller
      cohort (e.g. a rarer cancer type) will need lower cutoffs.
================================================================================
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless: never try to pop a GUI window
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.stats import fisher_exact, mannwhitneyu

# ==============================================================================
# Schema shared by every export_*.tsv file (no header row in the files).
# ==============================================================================
COLUMNS = [
    "Cohort", "Type", "Barcode", "Allele", "Contig",
    "First_base", "Last_base", "Strand", "Motif", "Annotated",
    "Uniquely_mapping", "Multimapping", "Overhang",
]

# A patient's tumor and normal sample can share a Barcode -> Type MUST be in
# every per-sample key, or their junctions get merged into one bucket.
GROUP_KEYS = ["Type", "Barcode", "Allele"]

TUMOR_COLOR = "#c0392b"   # red   (domain convention: tumor = red)
NORMAL_COLOR = "#2980b9"  # blue  (domain convention: normal = blue)
Z95 = 1.959963984540054   # z* for a 95% CI


# ==============================================================================
# 1. Load
# ==============================================================================
def load_data(paths):
    """Read one or more headerless 13-column SJ tsvs and concatenate.

    IMPORTANT: this preserves each file's internal row order. STAR already
    coordinate-sorts SJ.out.tab per sample/allele, and Annotated_cumsum below
    depends on that order -- never call sort_values() on this dataframe
    before compute_derived_columns().
    """
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep="\t", header=None, names=COLUMNS)
        df["Type"] = df["Type"].str.strip().str.lower()
        frames.append(df)
        print(f"  loaded {p}: {len(df):,} rows | "
              f"{df['Barcode'].nunique():,} distinct barcodes | "
              f"Type values = {sorted(df['Type'].unique())}")
    combined = pd.concat(frames, ignore_index=True)
    return combined


# ==============================================================================
# 2. Per-sample derived columns (cumsum bucketing + local splice rate)
# ==============================================================================
def compute_derived_columns(df):
    """Add Annotated_cumsum, Splice_coord, Local_rate.

    Annotated_cumsum: cumulative sum of Annotated within [Type, Barcode,
    Allele], in existing row order. This buckets each novel junction with the
    annotated junction it sits under (bucket N = "downstream of the Nth
    annotated junction").

    Local_rate: a junction's share of Uniquely_mapping reads among all rows
    in its own bucket -- i.e. how big this novel splice event is relative to
    the canonical junction it's competing with.
    """
    df = df.copy()
    g = df.groupby(GROUP_KEYS, sort=False)
    df["Annotated_cumsum"] = g["Annotated"].cumsum()
    df["Splice_coord"] = (
        df["First_base"].astype(str) + "_" + df["Last_base"].astype(str)
    )

    bucket_keys = GROUP_KEYS + ["Annotated_cumsum"]
    local_sum = df.groupby(bucket_keys, sort=False)["Uniquely_mapping"].transform("sum")
    df["Local_rate"] = df["Uniquely_mapping"] / local_sum
    return df


# ==============================================================================
# 3. Alternative-splicing event filter
# ==============================================================================
def find_as_events(df, local_rate_threshold=0.1):
    """Novel junction, carrying > local_rate_threshold of its local bucket's
    reads, downstream of at least one annotated junction (Annotated_cumsum
    != 0 excludes reads that fall before the first annotated junction in the
    transcript). local_rate_threshold is exposed (CLI: --local-rate-threshold,
    default 0.1 i.e. >10%) so the choice can be sensitivity-tested rather
    than trusted blindly -- see additional_analyses.py's threshold sweep."""
    events = df[
        (df["Annotated"] == 0)
        & (df["Local_rate"] > local_rate_threshold)
        & (df["Annotated_cumsum"] != 0)
    ].copy()
    return events


# ==============================================================================
# 4. Attach each event's annotated (wild-type) junction -> AS_pattern
# ==============================================================================
def annotate_events(df, events):
    """For each event, look up the Annotated==1 junction in the same
    [Type, Barcode, Allele, Annotated_cumsum] bucket; its Splice_coord
    becomes Annotated_coord. Builds AS_pattern = "{Annotated_coord}->{Splice_coord}"
    and Contig_rd = first 11 chars of Contig (2-field allele, e.g. hla_a_11_01)."""
    bucket_keys = GROUP_KEYS + ["Annotated_cumsum"]
    annotated_lookup = (
        df[df["Annotated"] == 1]
        .drop_duplicates(subset=bucket_keys, keep="first")  # 1 annotated junction per bucket, by construction
        [bucket_keys + ["Splice_coord"]]
        .rename(columns={"Splice_coord": "Annotated_coord"})
    )

    events = events.merge(annotated_lookup, on=bucket_keys, how="left")
    n_missing = events["Annotated_coord"].isna().sum()
    if n_missing:
        print(f"  WARNING: {n_missing} events have no matching annotated junction "
              f"in their bucket -- dropping them")
        events = events.dropna(subset=["Annotated_coord"])

    events["AS_pattern"] = events["Annotated_coord"] + "->" + events["Splice_coord"]
    events["Contig_rd"] = events["Contig"].str.slice(0, 11)
    return events


# ==============================================================================
# 5. Classify by coordinate geometry (generalizes across alleles/contigs,
#    since it never looks at absolute genomic position -- only at the
#    relative shift between the wild-type and variant junction).
# ==============================================================================
def classify_events(events, exon_skip_min_jump=100, intron_retention_max_shift=50):
    """Parse AS_pattern as "W1_W2->V1_V2" and classify:

        exon5_skip:        V1 == W1   and  V2 - W2 >= exon_skip_min_jump
                            (donor unchanged, acceptor jumps far downstream
                             -> the whole exon in between is skipped)

        intron5_retention: V2 == W2   and  0 < V1 - W1 <= intron_retention_max_shift
                            (acceptor unchanged, donor shifts a few bp
                             upstream into the intron -> partial retention)

        else: "other"
    """
    w = events["Annotated_coord"].str.split("_", expand=True).astype(int)
    v = events["Splice_coord"].str.split("_", expand=True).astype(int)
    w.columns = ["W1", "W2"]
    v.columns = ["V1", "V2"]
    events = pd.concat([events, w, v], axis=1)

    is_exon5_skip = (events["V1"] == events["W1"]) & (
        (events["V2"] - events["W2"]) >= exon_skip_min_jump
    )
    shift = events["V1"] - events["W1"]
    is_intron5_ret = (events["V2"] == events["W2"]) & (shift > 0) & (shift <= intron_retention_max_shift)

    events["Event_class"] = np.select(
        [is_exon5_skip, is_intron5_ret],
        ["exon5_skip", "intron5_retention"],
        default="other",
    )
    return events


# ==============================================================================
# Statistics helpers -- shared by the headline, per-allele, and per-pattern
# comparisons so every proportion/effect-size in this script is computed the
# same honest way.
# ==============================================================================
def wilson_ci(k, n, z=Z95):
    """Wilson score 95% CI for a binomial proportion k/n. Behaves sensibly
    (bounded in [0,1], sane width) even at small n or k=0/k=n, unlike the
    naive normal-approximation interval."""
    if n == 0:
        return float("nan"), float("nan")
    phat = k / n
    denom = 1 + z * z / n
    center = (phat + z * z / (2 * n)) / denom
    margin = (z * np.sqrt((phat * (1 - phat) + z * z / (4 * n)) / n)) / denom
    return max(0.0, center - margin), min(1.0, center + margin)


def risk_difference_ci(a, n1, c, n2, z=Z95):
    """Wald CI on the risk difference p1 - p2 (tumor - normal)."""
    if n1 == 0 or n2 == 0:
        return float("nan"), float("nan"), float("nan")
    p1, p2 = a / n1, c / n2
    diff = p1 - p2
    se = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    return diff, diff - z * se, diff + z * se


def odds_ratio_ci(a, b, c, d, z=Z95):
    """Odds ratio (tumor odds / normal odds) + Woolf-method log-OR CI.
    Applies a Haldane-Anscombe +0.5 continuity correction to all four cells
    whenever any cell is 0 (otherwise OR/CI would be 0 or infinite) and
    reports that a correction was applied so it can be flagged in tables."""
    corrected = min(a, b, c, d) == 0
    cc = 0.5 if corrected else 0.0
    a2, b2, c2, d2 = a + cc, b + cc, c + cc, d + cc
    or_ = (a2 * d2) / (b2 * c2)
    se_log = np.sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
    log_or = np.log(or_)
    return or_, float(np.exp(log_or - z * se_log)), float(np.exp(log_or + z * se_log)), corrected


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR q-values for an array of p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]  # enforce monotonicity
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out


# ==============================================================================
# 6. Headline comparison: class-level carrier fractions, CIs, effect sizes,
#    Fisher exact test
# ==============================================================================
def headline_comparison(df, events, classes=("exon5_skip", "intron5_retention")):
    """For each event class: fraction of tumor samples vs normal samples
    carrying >=1 such event (denominator = total distinct Barcodes of that
    Type, from the full dataset -- not just samples that happen to have an
    AS event). Reports Wilson 95% CIs on each proportion, the risk
    difference and odds ratio (both with 95% CIs), and a 2x2 Fisher exact
    test p-value."""
    totals = df.groupby("Type")["Barcode"].nunique()
    tumor_total = int(totals.get("tumor", 0))
    normal_total = int(totals.get("normal", 0))

    rows = []
    for cls in classes:
        carriers = (
            events.loc[events["Event_class"] == cls]
            .groupby("Type")["Barcode"].nunique()
        )
        a = int(carriers.get("tumor", 0))   # tumor carriers
        c = int(carriers.get("normal", 0))  # normal carriers
        b = tumor_total - a
        d = normal_total - c

        _, p_value = fisher_exact([[a, b], [c, d]])
        or_, or_lo, or_hi, or_corrected = odds_ratio_ci(a, b, c, d)
        rd, rd_lo, rd_hi = risk_difference_ci(a, tumor_total, c, normal_total)

        tumor_frac = a / tumor_total if tumor_total else float("nan")
        normal_frac = c / normal_total if normal_total else float("nan")
        t_lo, t_hi = wilson_ci(a, tumor_total)
        n_lo, n_hi = wilson_ci(c, normal_total)

        if tumor_frac > normal_frac:
            direction = "enriched in tumor"
        elif tumor_frac < normal_frac:
            direction = "depleted in tumor"
        else:
            direction = "no difference"

        rows.append({
            "Event_class": cls,
            "Tumor_carriers": a, "Tumor_total": tumor_total, "Tumor_pct": 100 * tumor_frac,
            "Tumor_CI_low": 100 * t_lo, "Tumor_CI_high": 100 * t_hi,
            "Normal_carriers": c, "Normal_total": normal_total, "Normal_pct": 100 * normal_frac,
            "Normal_CI_low": 100 * n_lo, "Normal_CI_high": 100 * n_hi,
            "Risk_diff_pct": 100 * rd, "Risk_diff_CI_low": 100 * rd_lo, "Risk_diff_CI_high": 100 * rd_hi,
            "Odds_ratio": or_, "OR_CI_low": or_lo, "OR_CI_high": or_hi, "OR_haldane_corrected": or_corrected,
            "Fisher_p": p_value,
            "Direction": direction,
        })
    return pd.DataFrame(rows)


def print_headline_table(headline):
    print("\n=== Headline comparison (class level, 95% CI, effect sizes) ===")
    for _, r in headline.iterrows():
        print(f"\n  {r['Event_class']}:")
        print(f"    Tumor : {r['Tumor_carriers']}/{r['Tumor_total']} = {r['Tumor_pct']:.1f}% "
              f"(95% CI {r['Tumor_CI_low']:.1f}-{r['Tumor_CI_high']:.1f}%)")
        print(f"    Normal: {r['Normal_carriers']}/{r['Normal_total']} = {r['Normal_pct']:.1f}% "
              f"(95% CI {r['Normal_CI_low']:.1f}-{r['Normal_CI_high']:.1f}%)")
        print(f"    Risk difference (tumor-normal): {r['Risk_diff_pct']:+.1f} pct pts "
              f"(95% CI {r['Risk_diff_CI_low']:+.1f} to {r['Risk_diff_CI_high']:+.1f})")
        or_note = " [Haldane-Anscombe corrected, a cell was 0]" if r["OR_haldane_corrected"] else ""
        print(f"    Odds ratio: {r['Odds_ratio']:.2f} (95% CI {r['OR_CI_low']:.2f}-{r['OR_CI_high']:.2f}){or_note}")
        print(f"    Fisher exact p = {r['Fisher_p']:.3g}  ->  {r['Direction']}")
        if r["Tumor_CI_low"] < r["Tumor_CI_high"] and (r["OR_CI_low"] < 1 < r["OR_CI_high"]):
            print(f"    NOTE: OR 95% CI spans 1 -- not a statistically confident effect at this sample size.")


# ==============================================================================
# 7. Detailed table: per [Contig_rd, AS_pattern, Event_class], with FDR
# ==============================================================================
def detailed_table(df, events, min_carriers=3):
    """Count distinct carriers (Barcodes) in tumor vs normal for every
    (Contig_rd, AS_pattern, Event_class) combo, with a per-allele denominator
    (how many samples of that Type carry that Contig_rd genotype at all --
    not every patient has every HLA-A allele). Runs a per-row Fisher exact
    test (only for rows with >= min_carriers total carriers -- anything
    rarer is too underpowered to say anything about) and applies
    Benjamini-Hochberg FDR correction across all tested rows, since we're
    running many tests here."""
    contig_denom = (
        df.groupby(["Contig_rd", "Type"])["Barcode"].nunique().unstack(fill_value=0)
    )
    for t in ("tumor", "normal"):
        if t not in contig_denom.columns:
            contig_denom[t] = 0

    carriers = (
        events.groupby(["Contig_rd", "AS_pattern", "Event_class", "Type"])["Barcode"]
        .nunique().unstack(fill_value=0)
    )
    for t in ("tumor", "normal"):
        if t not in carriers.columns:
            carriers[t] = 0
    carriers = carriers.rename(columns={"tumor": "Tumor_carriers", "normal": "Normal_carriers"}).reset_index()

    carriers["Tumor_denom"] = carriers["Contig_rd"].map(contig_denom["tumor"]).fillna(0).astype(int)
    carriers["Normal_denom"] = carriers["Contig_rd"].map(contig_denom["normal"]).fillna(0).astype(int)

    carriers["Tumor_freq"] = carriers["Tumor_carriers"] / carriers["Tumor_denom"].replace(0, np.nan)
    carriers["Normal_freq"] = carriers["Normal_carriers"] / carriers["Normal_denom"].replace(0, np.nan)
    carriers["Total_carriers"] = carriers["Tumor_carriers"] + carriers["Normal_carriers"]
    carriers["Tested"] = carriers["Total_carriers"] >= min_carriers

    # per-row Fisher exact test + odds ratio + Wilson CIs (only for tested rows)
    p_vals, ors, or_los, or_his = [], [], [], []
    t_los, t_his, n_los, n_his = [], [], [], []
    for _, r in carriers.iterrows():
        a, nt = r["Tumor_carriers"], r["Tumor_denom"]
        c, nn = r["Normal_carriers"], r["Normal_denom"]
        t_lo, t_hi = wilson_ci(a, nt)
        n_lo, n_hi = wilson_ci(c, nn)
        t_los.append(100 * t_lo); t_his.append(100 * t_hi)
        n_los.append(100 * n_lo); n_his.append(100 * n_hi)
        if r["Tested"]:
            _, p = fisher_exact([[a, nt - a], [c, nn - c]])
            orv, or_lo, or_hi, _ = odds_ratio_ci(a, nt - a, c, nn - c)
        else:
            p, orv, or_lo, or_hi = np.nan, np.nan, np.nan, np.nan
        p_vals.append(p); ors.append(orv); or_los.append(or_lo); or_his.append(or_hi)

    carriers["Fisher_p"] = p_vals
    carriers["Odds_ratio"] = ors
    carriers["OR_CI_low"] = or_los
    carriers["OR_CI_high"] = or_his
    carriers["Tumor_freq_CI_low"] = t_los
    carriers["Tumor_freq_CI_high"] = t_his
    carriers["Normal_freq_CI_low"] = n_los
    carriers["Normal_freq_CI_high"] = n_his

    q = np.full(len(carriers), np.nan)
    tested_mask = carriers["Tested"].to_numpy()
    if tested_mask.sum() > 0:
        q[tested_mask] = bh_fdr(carriers.loc[tested_mask, "Fisher_p"].to_numpy())
    carriers["FDR_q"] = q

    carriers = carriers.sort_values("Total_carriers", ascending=False).reset_index(drop=True)
    n_tested = int(carriers["Tested"].sum())
    n_sig = int((carriers["FDR_q"] < 0.05).sum())
    print(f"  detailed table: {len(carriers):,} distinct (Contig_rd, AS_pattern) combos; "
          f"{n_tested} tested (>= {min_carriers} total carriers); {n_sig} pass FDR q<0.05")
    return carriers


# ==============================================================================
# 8. Cohort overview: sample counts, carrier rate, event burden distribution,
#    distinct AS_pattern count
# ==============================================================================
def cohort_overview(df, events, outdir):
    """One row per Type: total samples, samples carrying >=1 AS event of any
    class, and the distribution (median/mean/max) of events-per-sample.
    Also returns the raw per-sample event counts table, reused by
    plot_event_burden()."""
    all_samples = df[["Type", "Barcode"]].drop_duplicates().set_index(["Type", "Barcode"])
    events_per_sample = events.groupby(["Type", "Barcode"]).size()
    counts = events_per_sample.reindex(all_samples.index, fill_value=0).reset_index(name="N_events")
    counts.to_csv(outdir / "table_events_per_sample.csv", index=False)

    rows = []
    for t in ("tumor", "normal"):
        sub = counts.loc[counts["Type"] == t, "N_events"]
        n_samples = len(sub)
        n_carriers = int((sub > 0).sum())
        rows.append({
            "Type": t,
            "N_samples": n_samples,
            "N_carriers_any_AS_event": n_carriers,
            "Pct_carriers": 100 * n_carriers / n_samples if n_samples else float("nan"),
            "Median_events_per_sample": sub.median() if n_samples else float("nan"),
            "Mean_events_per_sample": sub.mean() if n_samples else float("nan"),
            "Max_events_per_sample": int(sub.max()) if n_samples else 0,
            "N_distinct_AS_patterns": events.loc[events["Type"] == t, "AS_pattern"].nunique(),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "table_cohort_overview.csv", index=False)

    print("\n=== Cohort overview ===")
    print(summary.to_string(index=False))
    print(f"  {events['AS_pattern'].nunique()} distinct AS_patterns total across both groups")

    # -- figure: total samples vs carriers, tumor/normal side by side --
    fig, ax = plt.subplots(figsize=(5, 4.5))
    x = np.arange(len(summary))
    width = 0.35
    ax.bar(x - width / 2, summary["N_samples"], width, label="Total samples", color="#bdbdbd")
    bar_colors = [TUMOR_COLOR if t == "tumor" else NORMAL_COLOR for t in summary["Type"]]
    ax.bar(x + width / 2, summary["N_carriers_any_AS_event"], width,
           label="Carriers (>=1 AS event)", color=bar_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(summary["Type"])
    ax.set_ylabel("Samples")
    ax.set_title("HLA-A cohort overview")
    for i, r in summary.iterrows():
        ax.text(i + width / 2, r["N_carriers_any_AS_event"], f"{r['Pct_carriers']:.0f}%",
                ha="center", va="bottom", fontsize=8)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    outpath = outdir / "fig_cohort_overview.pdf"
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  wrote {outpath}")

    return summary, counts


# ==============================================================================
# 9. Per-sample event burden: total AS-event count per sample, tumor vs normal
# ==============================================================================
def plot_event_burden(counts, outdir):
    """AS events per sample (summed over both alleles), tumor vs normal, with
    a Mann-Whitney U test (non-parametric -- event counts are small
    non-negative integers, not remotely normal).

    Grouped bar of % of samples at each observed event count, not a boxplot:
    counts are so zero-inflated (>80% of samples in both arms have exactly 0
    events) that a boxplot collapses to a flat line at 0 with fliers -- an
    honest number, unreadable chart. The bar chart shows the exact same
    per-sample N_events distribution, just legibly."""
    tumor_vals = counts.loc[counts["Type"] == "tumor", "N_events"].to_numpy()
    normal_vals = counts.loc[counts["Type"] == "normal", "N_events"].to_numpy()

    u_stat, p_value = mannwhitneyu(tumor_vals, normal_vals, alternative="two-sided")

    max_events = int(max(tumor_vals.max(), normal_vals.max())) if len(tumor_vals) and len(normal_vals) else 0
    bins = np.arange(max_events + 1)
    tumor_pct = np.array([100 * (tumor_vals == b).sum() / len(tumor_vals) for b in bins])
    normal_pct = np.array([100 * (normal_vals == b).sum() / len(normal_vals) for b in bins])

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    width = 0.35
    ax.bar(bins - width / 2, tumor_pct, width, label="Tumor", color=TUMOR_COLOR)
    ax.bar(bins + width / 2, normal_pct, width, label="Normal", color=NORMAL_COLOR)
    ax.set_xticks(bins)
    ax.set_xlabel("AS events per sample (both alleles)")
    ax.set_ylabel("% of samples")
    ax.set_ylim(0, 100)
    tumor_carrier_pct = 100 * (tumor_vals > 0).sum() / len(tumor_vals) if len(tumor_vals) else float("nan")
    normal_carrier_pct = 100 * (normal_vals > 0).sum() / len(normal_vals) if len(normal_vals) else float("nan")
    ax.set_title(
        f"Per-sample AS-event burden (n_tumor={len(tumor_vals)}, n_normal={len(normal_vals)})\n"
        f"Mann-Whitney U p = {p_value:.3g}\n"
        f">=1 event: tumor {tumor_carrier_pct:.0f}%, normal {normal_carrier_pct:.0f}%",
        fontsize=11,
    )
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    outpath = outdir / "fig_event_burden.pdf"
    fig.savefig(outpath)
    plt.close(fig)

    result = {
        "U_statistic": u_stat, "p_value": p_value,
        "tumor_n": len(tumor_vals), "normal_n": len(normal_vals),
        "tumor_median": float(np.median(tumor_vals)) if len(tumor_vals) else float("nan"),
        "normal_median": float(np.median(normal_vals)) if len(normal_vals) else float("nan"),
        "tumor_mean": float(np.mean(tumor_vals)) if len(tumor_vals) else float("nan"),
        "normal_mean": float(np.mean(normal_vals)) if len(normal_vals) else float("nan"),
    }
    print("\n=== Per-sample event burden ===")
    print(f"  tumor  median={result['tumor_median']:.1f} mean={result['tumor_mean']:.2f} (n={result['tumor_n']})")
    print(f"  normal median={result['normal_median']:.1f} mean={result['normal_mean']:.2f} (n={result['normal_n']})")
    print(f"  Mann-Whitney U = {u_stat:.1f}, p = {p_value:.3g}")
    print(f"  wrote {outpath}")
    return result


# ==============================================================================
# 10. Per-allele view: is the effect allele-specific?
# ==============================================================================
def per_allele_view(df, events, min_allele_samples=20, classes=("exon5_skip", "intron5_retention"),
                     exclude_alleles=()):
    """For HLA-A alleles (Contig_rd) carried by >= min_allele_samples samples
    total (tumor+normal), compute exon5_skip / intron5_retention carrier
    frequency in tumor vs normal separately, with Wilson CIs and a per-allele
    Fisher exact test. Lets you see whether an effect that looks real at the
    whole-cohort level is actually concentrated in one or two alleles (the
    MHC Hammer paper reports this for HLA-A*11:01).

    exclude_alleles: Contig_rd values to leave out of this comparison (e.g.
    an allele diagnosed as a constitutive annotation artifact rather than
    real variable splicing -- see diagnostic_annotation_check.py). They are
    NOT removed from the input data, only from this particular table/figure;
    the descriptive/overview outputs (cohort overview, heatmaps, volcano,
    detailed table) still show them."""
    contig_denom = df.groupby(["Contig_rd", "Type"])["Barcode"].nunique().unstack(fill_value=0)
    for t in ("tumor", "normal"):
        if t not in contig_denom.columns:
            contig_denom[t] = 0
    contig_denom["Total"] = contig_denom["tumor"] + contig_denom["normal"]
    common_alleles = contig_denom.loc[contig_denom["Total"] >= min_allele_samples].index.tolist()
    excluded_present = [a for a in common_alleles if a in exclude_alleles]
    common_alleles = [a for a in common_alleles if a not in exclude_alleles]
    print(f"\n=== Per-allele view ===")
    print(f"  {len(common_alleles)} alleles with >= {min_allele_samples} total samples: {common_alleles}")
    if excluded_present:
        print(f"  excluded from this comparison (--exclude-alleles): {excluded_present}")

    rows = []
    for allele in common_alleles:
        nt = int(contig_denom.loc[allele, "tumor"])
        nn = int(contig_denom.loc[allele, "normal"])
        for cls in classes:
            sub = events[(events["Contig_rd"] == allele) & (events["Event_class"] == cls)]
            carriers = sub.groupby("Type")["Barcode"].nunique()
            a = int(carriers.get("tumor", 0))
            c = int(carriers.get("normal", 0))
            _, p = fisher_exact([[a, nt - a], [c, nn - c]])
            t_lo, t_hi = wilson_ci(a, nt)
            n_lo, n_hi = wilson_ci(c, nn)
            rows.append({
                "Contig_rd": allele, "Event_class": cls,
                "Tumor_carriers": a, "Tumor_n": nt,
                "Tumor_pct": 100 * a / nt if nt else float("nan"),
                "Tumor_CI_low": 100 * t_lo, "Tumor_CI_high": 100 * t_hi,
                "Normal_carriers": c, "Normal_n": nn,
                "Normal_pct": 100 * c / nn if nn else float("nan"),
                "Normal_CI_low": 100 * n_lo, "Normal_CI_high": 100 * n_hi,
                "Fisher_p": p,
            })
    result = pd.DataFrame(rows)
    return result, common_alleles


def plot_per_allele(per_allele, outdir, gene_label, classes=("exon5_skip", "intron5_retention"),
                     filename="fig_per_allele.pdf", title_suffix=""):
    outpath = outdir / filename
    if per_allele.empty:
        print(f"  no alleles meet --min-allele-samples -- skipping {outpath}")
        return
    alleles = sorted(per_allele["Contig_rd"].unique())
    ncols = min(3, len(alleles))
    nrows = int(np.ceil(len(alleles) / ncols))
    # extra headroom per subplot (4.6 vs 3.8) + constrained_layout so the
    # title of one row never collides with the p-value labels of the row
    # above, and the suptitle/legend don't leave a dead gap above row 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 4.6 * nrows),
                              squeeze=False, constrained_layout=True)

    for i, allele in enumerate(alleles):
        ax = axes[i // ncols][i % ncols]
        sub = per_allele[per_allele["Contig_rd"] == allele].set_index("Event_class").loc[list(classes)]
        x = np.arange(len(classes))
        width = 0.35
        # clip to >=0: Wilson bounds can be within float epsilon of the point
        # estimate (e.g. exactly 0% carriers), which produces a -1e-16
        # "negative" error bar that matplotlib rejects outright.
        t_err = [np.clip(sub["Tumor_pct"] - sub["Tumor_CI_low"], 0, None),
                 np.clip(sub["Tumor_CI_high"] - sub["Tumor_pct"], 0, None)]
        n_err = [np.clip(sub["Normal_pct"] - sub["Normal_CI_low"], 0, None),
                 np.clip(sub["Normal_CI_high"] - sub["Normal_pct"], 0, None)]
        ax.bar(x - width / 2, sub["Tumor_pct"], width, yerr=t_err, capsize=3,
               color=TUMOR_COLOR, label="Tumor")
        ax.bar(x + width / 2, sub["Normal_pct"], width, yerr=n_err, capsize=3,
               color=NORMAL_COLOR, label="Normal")
        ax.set_xticks(x)
        ax.set_xticklabels(classes, rotation=20, ha="right", fontsize=8)
        ax.set_title(f"{allele}\n(n_tumor={int(sub['Tumor_n'].iloc[0])}, n_normal={int(sub['Normal_n'].iloc[0])})",
                     fontsize=9, pad=16)
        ax.set_ylabel("Carriers (%)", fontsize=8)

        # headroom sized to THIS subplot's own tallest CI whisker, so the
        # p-value label always lands inside the axes instead of drifting
        # into the title above it
        ymax = max(sub["Tumor_CI_high"].max(), sub["Normal_CI_high"].max(), 1.0)
        ax.set_ylim(0, ymax * 1.28)
        for j, cls in enumerate(classes):
            nt_j, nn_j = sub.loc[cls, "Tumor_n"], sub.loc[cls, "Normal_n"]
            p_label = "n/a" if (nt_j == 0 or nn_j == 0) else f"p={sub.loc[cls, 'Fisher_p']:.2g}"
            ax.text(j, max(sub.loc[cls, "Tumor_CI_high"], sub.loc[cls, "Normal_CI_high"]) + 0.05 * ymax,
                    p_label, ha="center", va="bottom", fontsize=6.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # hide unused axes
    for j in range(len(alleles), nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")

    handles = [Patch(facecolor=TUMOR_COLOR, label="Tumor"), Patch(facecolor=NORMAL_COLOR, label="Normal")]
    if ncols >= 3:
        # wide enough that a fig-level legend "outside upper right" clears
        # the suptitle horizontally, even when title_suffix wraps to 2 lines
        fig.legend(handles=handles, loc="outside upper right", frameon=False)
    else:
        # narrow grid (talk-view subsets can be as few as 1-2 columns): a
        # fig-level legend at this width sits directly under a wrapped
        # 2-line suptitle and overlaps it -- anchor inside the last axes'
        # own corner instead, which stays clear of the suptitle regardless
        # of how many lines title_suffix adds
        axes[0][-1].legend(handles=handles, loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(f"{gene_label} exon5_skip / intron5_retention by allele (95% Wilson CI){title_suffix}", fontsize=13)
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  wrote {outpath}")


# ==============================================================================
# 11. Volcano plot: every recurrent AS_pattern, effect size vs significance
# ==============================================================================
def volcano_data(detailed, min_carriers=3):
    """Build the plotting frame from the detailed table: log2 fold-change of
    tumor vs normal carrier frequency (Haldane-Anscombe corrected so it's
    finite even when a pattern is private to one group) vs -log10(Fisher p).
    Only includes patterns with >= min_carriers total carriers -- anything
    rarer is a single-patient anecdote, not a recurrent pattern worth
    plotting on a significance axis."""
    d = detailed.loc[detailed["Total_carriers"] >= min_carriers].copy()
    tumor_freq_c = (d["Tumor_carriers"] + 0.5) / (d["Tumor_denom"] + 1)
    normal_freq_c = (d["Normal_carriers"] + 0.5) / (d["Normal_denom"] + 1)
    d["log2FC_tumor_vs_normal"] = np.log2(tumor_freq_c / normal_freq_c)
    p_floor = np.nextafter(0, 1)
    d["neglog10_Fisher_p"] = -np.log10(d["Fisher_p"].clip(lower=p_floor))
    d = d.sort_values("neglog10_Fisher_p", ascending=False).reset_index(drop=True)
    return d


def plot_volcano(volcano, outdir, gene_label, label_n=6, filename="fig_volcano.pdf"):
    outpath = outdir / filename
    if volcano.empty:
        print(f"  no patterns meet --min-carriers -- skipping {outpath}")
        return
    class_color = {"exon5_skip": "#8e44ad", "intron5_retention": "#16a085", "other": "#999999"}
    colors = volcano["Event_class"].map(class_color).fillna("#999999")

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(volcano["log2FC_tumor_vs_normal"], volcano["neglog10_Fisher_p"],
               c=colors, s=45, edgecolor="black", linewidth=0.3, alpha=0.85, zorder=3)

    # Label the top points, staggered into 3 vertical tiers (cycling) with a
    # thin connector line back to the marker -- with only a handful of
    # labeled points this is enough to keep them from stacking on top of
    # each other when several sit at a similar y (a common volcano-plot
    # failure mode when many patterns share ~the same p-value).
    top = volcano.nlargest(label_n, "neglog10_Fisher_p").sort_values("log2FC_tumor_vs_normal")
    for i, (_, r) in enumerate(top.iterrows()):
        dy = 14 + 16 * (i % 3)
        dx = -20 if i % 2 == 0 else 20
        ax.annotate(f"{r['Contig_rd']}\n{r['AS_pattern']}",
                    (r["log2FC_tumor_vs_normal"], r["neglog10_Fisher_p"]),
                    fontsize=6, xytext=(dx, dy), textcoords="offset points",
                    ha="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="gray"), zorder=4)

    nominal_line = -np.log10(0.05)
    ax.axhline(nominal_line, color="gray", linestyle="--", linewidth=0.8)
    ax.text(ax.get_xlim()[1] if ax.get_xlim()[1] > 0 else 0, nominal_line, "  nominal p=0.05 (uncorrected)",
            color="gray", fontsize=7, va="bottom", ha="left")
    ax.axvline(0, color="gray", linestyle="-", linewidth=0.6)
    ax.set_xlabel("log2( tumor carrier freq / normal carrier freq )")
    ax.set_ylabel("-log10( Fisher p )")

    # tested-vs-significant counts pulled straight from the already-computed
    # Tested/FDR_q columns (detailed_table()'s BH-FDR correction) -- no
    # re-thresholding or recomputation here, just reporting what's already
    # in the dataframe so a null result reads as "tested and null", not
    # "empty because nothing was tested".
    n_tested = int(volcano["Tested"].sum()) if "Tested" in volcano.columns else len(volcano)
    n_sig = int((volcano["FDR_q"] < 0.05).sum()) if "FDR_q" in volcano.columns else 0
    ax.set_title(f"{gene_label} AS_pattern recurrence, tumor vs normal\n"
                 f"(patterns with >= {int(volcano['Total_carriers'].min())} carriers shown)\n"
                 f"{n_tested} patterns tested, {n_sig} significant after BH-FDR q<0.05")
    # legend sits outside the plotting area entirely -- with this few points
    # spanning the full axes, any in-plot corner risks covering a marker/label
    handles = [Patch(facecolor=c, label=k) for k, c in class_color.items()]
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ymax = volcano["neglog10_Fisher_p"].max()
    ax.set_ylim(-0.05 * max(ymax, 1), max(ymax, 1) * 1.3)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outpath}")


# ==============================================================================
# 12a. Figure: grouped bar chart, tumor vs normal carrier %, with 95% CI
#      whiskers and Fisher p annotation
# ==============================================================================
def plot_classes(headline, outpath, gene_label):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    x = np.arange(len(headline))
    width = 0.35

    # clip to >=0 -- see plot_per_allele() for why float epsilon can make a
    # Wilson bound look infinitesimally negative
    t_err = [np.clip(headline["Tumor_pct"] - headline["Tumor_CI_low"], 0, None),
             np.clip(headline["Tumor_CI_high"] - headline["Tumor_pct"], 0, None)]
    n_err = [np.clip(headline["Normal_pct"] - headline["Normal_CI_low"], 0, None),
             np.clip(headline["Normal_CI_high"] - headline["Normal_pct"], 0, None)]

    ax.bar(x - width / 2, headline["Tumor_pct"], width, yerr=t_err, capsize=4,
           label="Tumor", color=TUMOR_COLOR)
    ax.bar(x + width / 2, headline["Normal_pct"], width, yerr=n_err, capsize=4,
           label="Normal", color=NORMAL_COLOR)

    ax.set_xticks(x)
    ax.set_xticklabels(headline["Event_class"])
    ax.set_ylabel("Carrier samples (%, 95% Wilson CI)")
    ax.set_title(f"{gene_label} alternative splicing: tumor vs normal")

    ymax = max(headline["Tumor_CI_high"].max(), headline["Normal_CI_high"].max(), 1)
    ax.set_ylim(0, ymax * 1.35)
    for i, row in headline.iterrows():
        top = max(row["Tumor_CI_high"], row["Normal_CI_high"])
        if gene_label == "HLA-C":
            # Direction ("enriched"/"depleted" in tumor) is RETIRED for HLA-C
            # (see MINI_PAPER.md Sec.4): intron5_retention's OR rests on a
            # single normal-side carrier and is non-significant, so a
            # directional claim isn't supported. Neutral caveat only --
            # Fisher p and OR values are unchanged.
            if row["Event_class"] == "intron5_retention":
                caveat = f"not significant (rests on {int(row['Normal_carriers'])} normal carrier)"
            else:
                caveat = "not significant"
            label = f"Fisher p={row['Fisher_p']:.2g}\nOR={row['Odds_ratio']:.2f}\n{caveat}"
        else:
            label = f"Fisher p={row['Fisher_p']:.2g}\nOR={row['Odds_ratio']:.2f}\n{row['Direction']}"
        ax.text(i, top + 0.04 * ymax, label, ha="center", va="bottom", fontsize=8)

    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  wrote {outpath}")


# ==============================================================================
# 12b. Figure: samples x alleles heatmap, colored by AS_pattern
# ==============================================================================
def plot_heatmap(events, sample_type, base_cmap, outpath, title, max_rows_per_panel=60):
    """Only samples carrying >=1 event of this Type are shown. Columns are
    Contig_rd (2-field allele); a cell is colored by the AS_pattern(s) seen
    in that sample x allele (joined with '; ' if more than one). Palette is
    a single hue family (Reds for tumor, Blues for normal) sampled across
    the distinct patterns actually observed, per the domain convention of
    tumor=red / normal=blue.

    Samples are split into separate files (each capped at max_rows_per_panel
    rows, suffixed _1/_2/_3/... when there's more than one) rather than one
    tall single-column heatmap -- a cohort the size of tumor (~100+
    carriers) renders as a figure many times taller than wide, which
    doesn't fit a slide even sideways. Splitting into standalone files
    (rather than side-by-side panels in one wide figure) means each one
    drops onto its own slide at a normal, readable size."""
    sub = events[events["Type"] == sample_type]
    if sub.empty:
        print(f"  no {sample_type} events -- skipping {outpath}")
        return

    cell = (
        sub.groupby(["Barcode", "Contig_rd"])["AS_pattern"]
        .agg(lambda s: "; ".join(sorted(set(s))))
    )
    wide = cell.unstack("Contig_rd").sort_index()

    patterns = sorted(p for p in pd.unique(wide.values.ravel()) if pd.notna(p))
    n = len(patterns)
    cmap_obj = plt.get_cmap(base_cmap)
    # sample from the middle-to-dark part of the ramp -- the pale end is
    # barely visible against a white background
    shades = [cmap_obj(0.35 + 0.55 * i / max(n - 1, 1)) for i in range(n)]
    code_of = {p: i for i, p in enumerate(patterns)}
    code_matrix = wide.replace(code_of).astype(float).to_numpy()

    n_samples, n_alleles = code_matrix.shape
    n_panels = int(np.ceil(n_samples / max_rows_per_panel))
    panel_bounds = np.array_split(np.arange(n_samples), n_panels)

    # If there are many patterns the legend gets unwieldy -- that's expected
    # for tumor (more distinct events); consider capping to top-N by
    # frequency if this becomes unreadable on a bigger cohort.
    handles = [Patch(facecolor=shades[i], edgecolor="black", label=p) for i, p in enumerate(patterns)]

    fig_w = max(6, 0.5 * n_alleles + 3)   # + room for legend
    written = []
    for panel_i, idx in enumerate(panel_bounds):
        part_outpath = outpath if n_panels == 1 else outpath.with_stem(f"{outpath.stem}_{panel_i + 1}")
        part_title = title if n_panels == 1 else f"{title} ({panel_i + 1}/{n_panels})"

        fig_h = max(4, 0.18 * len(idx))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        masked = np.ma.masked_invalid(code_matrix[idx, :])
        ax.pcolormesh(masked, cmap=ListedColormap(shades), edgecolors="black",
                      linewidth=0.4, vmin=-0.5, vmax=max(n - 0.5, 0.5))
        ax.set_xticks(np.arange(n_alleles) + 0.5)
        ax.set_xticklabels(wide.columns, rotation=90, fontsize=7)
        ax.set_yticks(np.arange(len(idx)) + 0.5)
        ax.set_yticklabels(wide.index[idx], fontsize=6)
        ax.invert_yaxis()  # first sample at top
        ax.set_xlabel("Contig_rd (2-field allele)")
        ax.set_ylabel("Barcode")
        ax.set_title(part_title)
        ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left",
                  fontsize=6, title="AS_pattern", frameon=False)

        fig.tight_layout()
        fig.savefig(part_outpath, bbox_inches="tight")
        plt.close(fig)
        written.append(part_outpath)

    names = ", ".join(str(p) for p in written)
    print(f"  wrote {names} ({len(wide.index)} samples x {len(wide.columns)} alleles, "
          f"{n} patterns, {n_panels} file{'s' if n_panels != 1 else ''})")


# ==============================================================================
# main
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tumor", default="data/export_BRCA_tumor_HLA-A.tsv", type=Path,
                         help="tumor SJ export tsv (13 cols, no header)")
    parser.add_argument("--normal", default="data/export_BRCA_normal_HLA-A.tsv", type=Path,
                         help="normal SJ export tsv (13 cols, no header)")
    parser.add_argument("--outdir", default="outputs", type=Path)
    parser.add_argument("--local-rate-threshold", type=float, default=0.1,
                         help="min share of local-bucket reads for a novel junction to count as an AS event")
    parser.add_argument("--exon-skip-min-jump", type=int, default=100,
                         help="min acceptor shift (bp) to call exon5_skip")
    parser.add_argument("--intron-retention-max-shift", type=int, default=50,
                         help="max donor shift (bp) to call intron5_retention")
    parser.add_argument("--min-allele-samples", type=int, default=20,
                         help="min total (tumor+normal) samples carrying an allele to include it in the per-allele view")
    parser.add_argument("--min-carriers", type=int, default=3,
                         help="min total carriers for a per-pattern comparison to be tested (Fisher p / FDR / volcano)")
    parser.add_argument("--volcano-label-n", type=int, default=6,
                         help="number of most-significant points to label on the volcano plot")
    parser.add_argument("--exclude-alleles", default="hla_a_11_01",
                         help="comma-separated Contig_rd values to drop from the tumor-vs-normal "
                              "headline and per-allele comparisons only (e.g. an allele diagnosed as a "
                              "constitutive annotation artifact by diagnostic_annotation_check.py -- "
                              "hla_a_11_01's intron5_retention is ~universal in both tumor and normal "
                              "carriers and isn't tumor-associated). They are NOT removed from the "
                              "descriptive/overview outputs (cohort overview, heatmaps, volcano, detailed "
                              "table) -- pass '' to disable exclusion entirely.")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    exclude_alleles = [a.strip() for a in args.exclude_alleles.split(",") if a.strip()]

    # ---- 1. load ----------------------------------------------------------
    print("=== Loading data ===")
    df = load_data([args.tumor, args.normal])
    print(f"  combined: {len(df):,} rows total")
    for t, n in df.groupby("Type")["Barcode"].nunique().items():
        print(f"  {t}: {n:,} distinct samples")

    # ---- 2. per-sample derived columns ------------------------------------
    print("\n=== Computing per-sample splice metrics ===")
    df = compute_derived_columns(df)
    df["Contig_rd"] = df["Contig"].str.slice(0, 11)

    # ---- 3-4. AS events + annotate -----------------------------------------
    events = find_as_events(df, local_rate_threshold=args.local_rate_threshold)
    print(f"  {len(events):,} candidate AS events (Annotated==0, Local_rate>{args.local_rate_threshold}, "
          f"Annotated_cumsum!=0) out of {len(df):,} rows")
    events = annotate_events(df, events)
    print(f"  {len(events):,} events retained after attaching annotated junction")

    # ---- 5. classify --------------------------------------------------------
    events = classify_events(
        events,
        exon_skip_min_jump=args.exon_skip_min_jump,
        intron_retention_max_shift=args.intron_retention_max_shift,
    )
    print("\n=== Event class counts (rows, not distinct carriers) ===")
    print(events.groupby(["Type", "Event_class"]).size().unstack(fill_value=0))

    # ---- 6. cohort overview ---------------------------------------------------
    overview, counts = cohort_overview(df, events, args.outdir)

    # ---- 7. headline comparison (with CI + effect sizes) -----------------------
    # Print BOTH versions so the effect of --exclude-alleles is visible: the
    # unfiltered number (every allele, including any diagnosed as constitutive
    # artifacts) and the presented number (those alleles' events excluded).
    # The denominator (total samples per Type) is unaffected either way --
    # only which events count as "this sample is a carrier" changes.
    headline_all = headline_comparison(df, events)
    print("\n--- Headline, ALL alleles (no exclusions) ---")
    print_headline_table(headline_all)
    headline_all.to_csv(args.outdir / "table_headline_classes_all_alleles.csv", index=False)

    if exclude_alleles:
        events_for_headline = events[~events["Contig_rd"].isin(exclude_alleles)]
        print(f"\n--- Headline, EXCLUDING {exclude_alleles} (this is the presented comparison) ---")
    else:
        events_for_headline = events
        print("\n--- Headline (--exclude-alleles disabled, same as ALL alleles above) ---")
    headline = headline_comparison(df, events_for_headline)
    print_headline_table(headline)
    headline.to_csv(args.outdir / "table_headline_classes.csv", index=False)
    plot_classes(headline, args.outdir / "fig_classes.pdf", gene_label="HLA-A")

    # ---- 8. per-allele view ---------------------------------------------------
    per_allele, common_alleles = per_allele_view(df, events, min_allele_samples=args.min_allele_samples,
                                                  exclude_alleles=exclude_alleles)
    per_allele.to_csv(args.outdir / "table_per_allele.csv", index=False)
    plot_per_allele(per_allele, args.outdir, gene_label="HLA-A")

    # "talk" subset: largest-n alleles union any allele with a nonzero
    # carrier bar -- the full grid above is unreadable projected on a slide,
    # this is the same numbers restricted to the panels that actually carry
    # a story. Panels/alleles are never dropped from the full grid, only
    # from this separate slide-sized figure.
    TALK_ALLELES_HLA_A = ["hla_a_02_01", "hla_a_03_01", "hla_a_01_01", "hla_a_24_02",
                           "hla_a_23_01", "hla_a_30_01"]
    talk_subset = per_allele[per_allele["Contig_rd"].isin(TALK_ALLELES_HLA_A)]
    plot_per_allele(talk_subset, args.outdir, gene_label="HLA-A", filename="fig_per_allele_talk.pdf",
                     title_suffix="\n(talk view: largest-n alleles + any with a nonzero carrier; "
                                   "full grid in fig_per_allele.pdf)")

    # ---- 9. detailed table (many patterns) + FDR --------------------------------
    detail = detailed_table(df, events, min_carriers=args.min_carriers)
    out_csv = args.outdir / "results_BRCA_tumor_vs_normal_HLA-A.csv"
    detail.to_csv(out_csv, index=False)
    print(f"  -> {out_csv}")

    # ---- 10. volcano plot -----------------------------------------------------
    volcano = volcano_data(detail, min_carriers=args.min_carriers)
    volcano.to_csv(args.outdir / "table_volcano.csv", index=False)
    plot_volcano(volcano, args.outdir, gene_label="HLA-A", label_n=args.volcano_label_n)

    # ---- 11. per-sample event burden --------------------------------------------
    plot_event_burden(counts, args.outdir)

    # ---- 12. heatmaps -------------------------------------------------------------
    print("\n=== Heatmaps ===")
    plot_heatmap(events, "tumor", "Reds", args.outdir / "fig_heatmap_tumor.pdf",
                 "HLA-A alternative splicing -- tumor")
    plot_heatmap(events, "normal", "Blues", args.outdir / "fig_heatmap_normal.pdf",
                 "HLA-A alternative splicing -- normal")

    print(f"\nDone. All outputs written to {args.outdir}/")


if __name__ == "__main__":
    main()
