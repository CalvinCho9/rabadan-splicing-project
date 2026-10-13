#!/usr/bin/env python3
"""
additional_analyses.py
================================================================================
Three analyses that round out the main pipeline (analyze_brca.py) for a
"have we done everything reasonable with this data" pass. Each is its own
function, its own CSV, and (where useful) its own figure in outputs/.

  1. Power analysis -- for the headline exon5_skip / intron5_retention test:
     what power did we actually have to detect the observed effect, and what
     effect size (minimum detectable effect, MDE) would 80% power have
     required at these sample sizes? Turns "not significant" into a
     quantified, honest "underpowered by how much."

  2. Threshold sensitivity -- the AS-event filter (Local_rate > 10%) and the
     two geometry cutoffs (exon_skip_min_jump, intron_retention_max_shift)
     are judgment calls. This sweeps each one-at-a-time around its default
     and re-runs the headline comparison, so we can show the conclusion
     doesn't hinge on having picked one specific threshold.

  3. Allele frequency spectrum -- descriptive background: which HLA-A
     alleles (Contig_rd) actually appear in this cohort, tumor vs normal,
     and how concentrated vs diverse the allele pool is. Not a tumor-vs-
     normal test; just cohort context for the audience.

Reuses analyze_brca's exact loading/bucketing/classification code so nothing
here can drift from the main analysis by reimplementing the logic
differently.

Usage: python3 additional_analyses.py
================================================================================
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.stats.power import NormalIndPower
from statsmodels.stats.proportion import proportion_effectsize

import analyze_brca as ab

OUTDIR = Path("outputs")
OUTDIR.mkdir(exist_ok=True)
EXCLUDE_ALLELES = ["hla_a_11_01"]  # keep in lockstep with analyze_brca.py's default
POWER = NormalIndPower()


def build_frame(local_rate_threshold=0.1, exon_skip_min_jump=100, intron_retention_max_shift=50):
    df = ab.load_data([Path("data/export_BRCA_tumor_HLA-A.tsv"), Path("data/export_BRCA_normal_HLA-A.tsv")])
    df = ab.compute_derived_columns(df)
    df["Contig_rd"] = df["Contig"].str.slice(0, 11)
    events = ab.find_as_events(df, local_rate_threshold=local_rate_threshold)
    events = ab.annotate_events(df, events)
    events = ab.classify_events(events, exon_skip_min_jump=exon_skip_min_jump,
                                 intron_retention_max_shift=intron_retention_max_shift)
    return df, events


# ==============================================================================
# 1. Power analysis
# ==============================================================================
def power_analysis(headline, outdir=OUTDIR):
    """For each event class in the (already exclude-alleles-filtered)
    headline table: achieved power to detect the OBSERVED effect at
    alpha=0.05 two-sided, and the minimum detectable effect (MDE, both as
    Cohen's h and as a tumor-rate range) for 80% power at the actual
    tumor:normal sample-size ratio. Two-proportion z-test power (Cohen's h
    effect size), via statsmodels -- standard for this kind of carrier-rate
    comparison."""
    rows = []
    for _, r in headline.iterrows():
        n_t, n_n = int(r["Tumor_total"]), int(r["Normal_total"])
        p_t, p_n = r["Tumor_pct"] / 100, r["Normal_pct"] / 100
        ratio = n_n / n_t

        h_obs = proportion_effectsize(p_t, p_n)
        achieved_power = POWER.power(effect_size=abs(h_obs), nobs1=n_t, ratio=ratio,
                                      alpha=0.05, alternative="two-sided")

        h_mde = POWER.solve_power(effect_size=None, nobs1=n_t, ratio=ratio,
                                   alpha=0.05, power=0.8, alternative="two-sided")
        # convert Cohen's h back to an actual tumor-rate range, holding the
        # observed normal rate fixed as baseline (arcsine transform inverse)
        phi_n = 2 * np.arcsin(np.sqrt(p_n))
        p_t_hi = np.sin(min((phi_n + h_mde) / 2, np.pi / 2)) ** 2
        p_t_lo = np.sin(max((phi_n - h_mde) / 2, 0)) ** 2

        rows.append({
            "Event_class": r["Event_class"],
            "Normal_pct_baseline": r["Normal_pct"], "Tumor_pct_observed": r["Tumor_pct"],
            "N_tumor": n_t, "N_normal": n_n,
            "Cohens_h_observed": h_obs,
            "Achieved_power_pct": 100 * achieved_power,
            "MDE_cohens_h_at_80pct_power": h_mde,
            "MDE_tumor_pct_low": 100 * p_t_lo, "MDE_tumor_pct_high": 100 * p_t_hi,
        })
    result = pd.DataFrame(rows)
    result.to_csv(outdir / "table_power_analysis.csv", index=False)
    return result


def print_power_analysis(power_table):
    print("\n=== 1. Power analysis (headline comparison, --exclude-alleles applied) ===")
    for _, r in power_table.iterrows():
        print(f"\n  {r['Event_class']}:")
        print(f"    observed: tumor {r['Tumor_pct_observed']:.1f}% (n={int(r['N_tumor'])}) vs "
              f"normal {r['Normal_pct_baseline']:.1f}% (n={int(r['N_normal'])})")
        print(f"    achieved power to detect the OBSERVED effect at alpha=0.05: {r['Achieved_power_pct']:.1f}%")
        print(f"    at 80% power, given these sample sizes, we could only detect a tumor rate "
              f"outside [{r['MDE_tumor_pct_low']:.1f}%, {r['MDE_tumor_pct_high']:.1f}%] "
              f"(holding the {r['Normal_pct_baseline']:.1f}% normal rate fixed)")
        if r["Achieved_power_pct"] < 80:
            print(f"    -> UNDERPOWERED: achieved power ({r['Achieved_power_pct']:.0f}%) is well below the "
                  f"conventional 80% target. A null result here is not informative either way.")


# ==============================================================================
# 2. Threshold sensitivity
# ==============================================================================
def threshold_sensitivity(outdir=OUTDIR):
    """One-at-a-time sweep around each of the three judgment-call
    thresholds, holding the other two at their analyze_brca.py default.
    Reports the intron5_retention headline (the class with enough carriers
    to move at all) at each setting, with --exclude-alleles applied."""
    sweeps = {
        "local_rate_threshold": [0.05, 0.10, 0.15, 0.20],
        "exon_skip_min_jump": [50, 75, 100, 150],
        "intron_retention_max_shift": [20, 35, 50, 70],
    }
    defaults = {"local_rate_threshold": 0.1, "exon_skip_min_jump": 100, "intron_retention_max_shift": 50}

    rows = []
    for param, values in sweeps.items():
        for v in values:
            kwargs = dict(defaults)
            kwargs[param] = v
            df, events = build_frame(**kwargs)
            events_f = events[~events["Contig_rd"].isin(EXCLUDE_ALLELES)]
            headline = ab.headline_comparison(df, events_f)
            r = headline.set_index("Event_class").loc["intron5_retention"]
            rows.append({
                "Swept_param": param, "Value": v, "Is_default": v == defaults[param],
                "Tumor_pct": r["Tumor_pct"], "Normal_pct": r["Normal_pct"],
                "Fisher_p": r["Fisher_p"], "Odds_ratio": r["Odds_ratio"],
                "OR_CI_low": r["OR_CI_low"], "OR_CI_high": r["OR_CI_high"],
                "Direction": r["Direction"],
            })
    result = pd.DataFrame(rows)
    result.to_csv(outdir / "table_threshold_sensitivity.csv", index=False)
    return result


def plot_threshold_sensitivity(sens, outdir=OUTDIR):
    params = sens["Swept_param"].unique()
    fig, axes = plt.subplots(1, len(params), figsize=(4.2 * len(params), 4), sharey=True)
    for ax, param in zip(axes, params):
        sub = sens[sens["Swept_param"] == param].sort_values("Value")
        yerr = [sub["Odds_ratio"] - sub["OR_CI_low"], sub["OR_CI_high"] - sub["Odds_ratio"]]
        yerr = [np.clip(e, 0, None) for e in yerr]
        colors = ["#c0392b" if d else "#555555" for d in sub["Is_default"]]
        ax.errorbar(sub["Value"], sub["Odds_ratio"], yerr=yerr, fmt="o", capsize=4,
                     ecolor="#888888", markerfacecolor="none")
        for x, y, c in zip(sub["Value"], sub["Odds_ratio"], colors):
            ax.scatter([x], [y], color=c, zorder=3, s=50)
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
        ax.set_xlabel(param, fontsize=9)
        ax.set_title(f"vary {param}\n(red = default)", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("intron5_retention odds ratio\n(tumor vs normal, 95% CI)")
    fig.suptitle("Threshold sensitivity: conclusion is stable across reasonable cutoffs", y=1.03)
    fig.tight_layout()
    outpath = outdir / "fig_threshold_sensitivity.pdf"
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outpath}")


def print_threshold_sensitivity(sens):
    print("\n=== 2. Threshold sensitivity (intron5_retention, --exclude-alleles applied) ===")
    with pd.option_context("display.width", 160, "display.max_columns", None):
        print(sens.drop(columns=["OR_CI_low", "OR_CI_high"]).to_string(index=False))
    or_range = (sens["Odds_ratio"].min(), sens["Odds_ratio"].max())
    p_range = (sens["Fisher_p"].min(), sens["Fisher_p"].max())
    all_same_direction = sens["Direction"].nunique() == 1
    print(f"\n  across all {len(sens)} threshold settings tested: OR ranges "
          f"{or_range[0]:.2f}-{or_range[1]:.2f}, Fisher p ranges {p_range[0]:.2g}-{p_range[1]:.2g}, "
          f"direction is {'CONSISTENT' if all_same_direction else 'NOT consistent'} "
          f"({sens['Direction'].unique().tolist()})")


# ==============================================================================
# 3. Allele frequency spectrum
# ==============================================================================
def allele_frequency_spectrum(df, outdir=OUTDIR, top_n=20):
    """Every HLA-A allele (Contig_rd) observed in this cohort -- not just
    ones above the --min-allele-samples cutoff used for the per-allele
    comparison -- with tumor/normal carrier counts. Descriptive context,
    not a tumor-vs-normal test."""
    spectrum = (
        df.groupby(["Contig_rd", "Type"])["Barcode"].nunique().unstack(fill_value=0)
    )
    for t in ("tumor", "normal"):
        if t not in spectrum.columns:
            spectrum[t] = 0
    spectrum["Total"] = spectrum["tumor"] + spectrum["normal"]
    spectrum = spectrum.sort_values("Total", ascending=False).reset_index()
    spectrum = spectrum.rename(columns={"tumor": "Tumor_carriers", "normal": "Normal_carriers"})
    spectrum.to_csv(outdir / "table_allele_frequency_spectrum.csv", index=False)

    n_alleles = len(spectrum)
    top_n_share = 100 * spectrum["Total"].head(top_n).sum() / spectrum["Total"].sum()
    print(f"\n=== 3. Allele frequency spectrum ===")
    print(f"  {n_alleles} distinct HLA-A alleles (Contig_rd) observed across both cohorts")
    print(f"  top {top_n} alleles account for {top_n_share:.0f}% of all allele-carrying samples")
    print(spectrum.head(top_n).to_string(index=False))

    top = spectrum.head(top_n)
    fig, ax = plt.subplots(figsize=(max(8, 0.4 * len(top)), 4.5))
    x = np.arange(len(top))
    width = 0.38
    ax.bar(x - width / 2, top["Tumor_carriers"], width, label="Tumor", color=ab.TUMOR_COLOR)
    ax.bar(x + width / 2, top["Normal_carriers"], width, label="Normal", color=ab.NORMAL_COLOR)
    ax.set_xticks(x)
    ax.set_xticklabels(top["Contig_rd"], rotation=90, fontsize=7)
    ax.set_ylabel("Samples carrying this allele")
    ax.set_title(f"HLA-A allele frequency spectrum (top {top_n} of {n_alleles})")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    outpath = outdir / "fig_allele_frequency_spectrum.pdf"
    fig.savefig(outpath)
    plt.close(fig)
    print(f"  wrote {outpath}")
    return spectrum


# ==============================================================================
# main
# ==============================================================================
def main():
    print("=== Loading data (analyze_brca defaults) ===")
    df, events = build_frame()
    events_for_headline = events[~events["Contig_rd"].isin(EXCLUDE_ALLELES)]
    headline = ab.headline_comparison(df, events_for_headline)

    power_table = power_analysis(headline)
    print_power_analysis(power_table)

    sens = threshold_sensitivity()
    print_threshold_sensitivity(sens)
    plot_threshold_sensitivity(sens)

    allele_frequency_spectrum(df)

    print("\nDone. New outputs: table_power_analysis.csv, table_threshold_sensitivity.csv, "
          "fig_threshold_sensitivity.pdf, table_allele_frequency_spectrum.csv, "
          "fig_allele_frequency_spectrum.pdf (all in outputs/).")


if __name__ == "__main__":
    main()
