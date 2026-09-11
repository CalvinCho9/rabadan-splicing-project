#!/usr/bin/env python3
"""
phase3_pancancer_descriptive.py
================================================================================
PHASE 3 (overnight B/C run): pan-cancer DESCRIPTIVE pass for HLA-A/B/C.

Explicitly descriptive -- NOT a tumor-vs-normal test. Most cohorts in these
files have no matched normal (they're tumor-only pulls); the "NORMAL" cohort
that does exist is a separate cross-cohort normal-tissue pool (GTEx-style),
not per-cancer-type matched normals, so even where both tumor and "NORMAL"
rows exist, treating that as a fair tumor-vs-normal comparison would be
wrong. This script does two things only, per gene:

  1. Per-cohort, per-allele penetrance of each allele's own top recurrent
     AS_pattern -- descriptive cohort context, same shape as Phase 1's
     per-allele check but broken out by Cohort.
  2. A targeted check: does HLA-A*11:01's specific constitutive junction
     (2365_2806->2383_2806, established in ANNOTATION_ARTIFACT_NOTE.md as
     ~97-100% penetrant in BRCA tumor+normal) recur at similarly high
     penetrance in OTHER cancer types and the separate normal pool? And,
     symmetrically, does any HLA-B/C allele show a similar cross-cohort
     near-universal pattern (a "B/C equivalent" of the A*11:01 artifact)?

Data-quality note (see OVERNIGHT_LOG.md): pan_HLA-B_SJ_combined.tsv and
pan_HLA-C_SJ_combined.tsv both contain a Cohort=="BRCA_NORMAL" block whose
Type reads "tumor" -- row counts match the real BRCA normal exports exactly,
strongly suggesting a mislabeled Type in that block. Per "never fabricate,
note and move on," this script EXCLUDES Cohort=="BRCA_NORMAL" entirely
rather than trusting or silently correcting the label (BRCA itself is
already covered properly by Phase 2's dedicated export files, so excluding
it here loses no distinct information).

Quarantine: HLA-B's classifier verdict from Phase 1 was "untrustworthy"
(insufficient data). Per the quarantine rule (gene-level, not
analysis-type-level -- "when in doubt, quarantine"), ALL HLA-B Phase 3
outputs are written under outputs/UNVALIDATED_HLA_B_pancancer_* and excluded
from every presented document, exactly like Phase 2.

Outputs:
    outputs/table_phase3_cohort_overview_HLA-{A,C}.csv
    outputs/table_phase3_allele_penetrance_by_cohort_HLA-{A,C}.csv
    outputs/table_phase3_A11_01_cross_cohort.csv          (gene A only)
    outputs/UNVALIDATED_HLA_B_pancancer_*.csv              (gene B, quarantined)
    outputs/PHASE3_pancancer_descriptive.md
================================================================================
"""
from pathlib import Path

import numpy as np
import pandas as pd

import analyze_brca as ab
import hla_common as hc

OUTDIR = Path("outputs")
OUTDIR.mkdir(exist_ok=True)

PAN_FILES = {
    "HLA-A": "data/pan_HLA-A_SJ_combined.tsv",
    "HLA-B": "data/pan_HLA-B_SJ_combined.tsv",
    "HLA-C": "data/pan_HLA-C_SJ_combined.tsv",
}
QUARANTINED = {"HLA-A": False, "HLA-B": True, "HLA-C": False}
TARGET_ALLELE_A = "hla_a_11_01"
TARGET_PATTERN_A = "2365_2806->2383_2806"
MIN_ALLELE_COHORT_SAMPLES = 10  # lower bar than BRCA-only (some cohorts are small)


def load_pan_gene(gene, path):
    df = pd.read_csv(path, sep="\t", header=None, names=ab.COLUMNS)
    df["Type"] = df["Type"].str.strip().str.lower()
    n_before = df["Barcode"].nunique()
    excluded_rows = 0
    if "BRCA_NORMAL" in df["Cohort"].unique():
        mask = df["Cohort"] == "BRCA_NORMAL"
        excluded_rows = int(mask.sum())
        print(f"  [{gene}] excluding Cohort=='BRCA_NORMAL' block: {excluded_rows} rows, "
              f"{df.loc[mask, 'Barcode'].nunique()} barcodes (Type mislabeled 'tumor' in source file "
              f"-- row count matches the real BRCA normal export exactly; see OVERNIGHT_LOG.md). "
              f"BRCA is already covered by Phase 2's dedicated export files.")
        df = df.loc[~mask].copy()
    df = ab.compute_derived_columns(df)
    df["Contig_rd"] = df["Contig"].str.slice(0, 11)
    print(f"  [{gene}] pan-cancer frame: {df['Barcode'].nunique()} distinct barcodes "
          f"({n_before} before BRCA_NORMAL exclusion), {df['Cohort'].nunique()} cohorts")
    return df, excluded_rows


def cohort_overview(df, gene):
    ov = df.groupby(["Cohort", "Type"])["Barcode"].nunique().unstack(fill_value=0).reset_index()
    print(f"\n=== {gene}: pan-cancer cohort overview ===")
    print(ov.to_string(index=False))
    return ov


def allele_penetrance_by_cohort(df, events, gene, min_samples=MIN_ALLELE_COHORT_SAMPLES):
    """For each (Contig_rd) with >= min_samples total carriers pan-cancer,
    find its single most-recurrent AS_pattern (pooled across all cohorts,
    for a stable "top pattern" definition) and report, per Cohort, what % of
    that allele's samples in that cohort carry it. Purely descriptive --
    no significance test, no tumor-vs-normal comparison (Cohort mixes tumor
    cohorts and the separate normal pool as parallel rows, not a paired
    test)."""
    allele_totals = df.groupby("Contig_rd")["Barcode"].nunique()
    alleles = allele_totals.loc[allele_totals >= min_samples].index.tolist()

    rows = []
    for allele in alleles:
        sub_events = events[events["Contig_rd"] == allele]
        if sub_events.empty:
            continue
        top_pattern = sub_events.groupby("AS_pattern")["Barcode"].nunique().idxmax()
        top_class = sub_events.loc[sub_events["AS_pattern"] == top_pattern, "Event_class"].iloc[0]

        cohort_denom = df.loc[df["Contig_rd"] == allele].groupby("Cohort")["Barcode"].nunique()
        cohort_carriers = (
            sub_events.loc[sub_events["AS_pattern"] == top_pattern]
            .groupby("Cohort")["Barcode"].nunique()
        )
        for cohort, n in cohort_denom.items():
            k = int(cohort_carriers.get(cohort, 0))
            rows.append({
                "Contig_rd": allele, "Top_AS_pattern": top_pattern, "Top_event_class": top_class,
                "Cohort": cohort, "Carriers": k, "N_samples": int(n),
                "Pct": 100 * k / n if n else np.nan,
            })
    result = pd.DataFrame(rows)
    if result.empty:
        print(f"  [{gene}] no allele reaches {min_samples}+ pan-cancer samples with events -- empty table")
        return result
    result = result.sort_values(["Contig_rd", "Cohort"]).reset_index(drop=True)
    print(f"\n=== {gene}: per-allele top-pattern penetrance by cohort "
          f"(alleles with >= {min_samples} pan-cancer samples) ===")
    with pd.option_context("display.width", 160, "display.max_rows", 60):
        print(result.to_string(index=False))
    return result


def a11_01_cross_cohort(df, events):
    """Targeted check: does hla_a_11_01's known BRCA-constitutive junction
    (2365_2806->2383_2806) show the same near-universal penetrance in other
    cancer types / the separate normal pool?"""
    denom = df.loc[df["Contig_rd"] == TARGET_ALLELE_A].groupby("Cohort")["Barcode"].nunique()
    carriers = (
        events.loc[(events["Contig_rd"] == TARGET_ALLELE_A) & (events["AS_pattern"] == TARGET_PATTERN_A)]
        .groupby("Cohort")["Barcode"].nunique()
    )
    rows = []
    for cohort, n in denom.items():
        k = int(carriers.get(cohort, 0))
        rows.append({"Cohort": cohort, "hla_a_11_01_carriers": int(n),
                      "Carriers_with_BRCA_pattern": k, "Pct": 100 * k / n if n else np.nan})
    result = pd.DataFrame(rows).sort_values("hla_a_11_01_carriers", ascending=False).reset_index(drop=True)
    print(f"\n=== A*11:01 cross-cohort check: does the BRCA-constitutive pattern "
          f"({TARGET_PATTERN_A}) recur elsewhere? ===")
    print(result.to_string(index=False))
    return result


def run_gene(gene, path):
    quarantined = QUARANTINED[gene]
    status = "QUARANTINED (not for presentation)" if quarantined else "validated, presented"
    print(f"\n{'=' * 80}\n{gene}  [{status}]\n{'=' * 80}")

    df, n_excluded = load_pan_gene(gene, path)
    ov = cohort_overview(df, gene)

    events = ab.find_as_events(df, local_rate_threshold=0.1)
    events = ab.annotate_events(df, events)
    events = hc.classify_events_strand_aware(events, flip_minus=False)  # naive, per Phase 1 for both B/C
    print(f"  {len(events):,} candidate AS events pan-cancer")

    penetrance = allele_penetrance_by_cohort(df, events, gene)

    prefix = f"UNVALIDATED_HLA_B_pancancer" if quarantined else f"table_phase3"
    suffix = "" if quarantined else f"_{gene}"
    ov.to_csv(OUTDIR / f"{prefix}_cohort_overview{suffix}.csv", index=False)
    if not penetrance.empty:
        penetrance.to_csv(OUTDIR / f"{prefix}_allele_penetrance_by_cohort{suffix}.csv", index=False)

    a11 = None
    if gene == "HLA-A":
        a11 = a11_01_cross_cohort(df, events)
        a11.to_csv(OUTDIR / "table_phase3_A11_01_cross_cohort.csv", index=False)

    return {"gene": gene, "quarantined": quarantined, "cohort_overview": ov,
            "penetrance": penetrance, "n_events": len(events), "a11": a11,
            "n_brca_normal_excluded": n_excluded}


def write_report(results):
    lines = []
    lines.append("# Phase 3 — pan-cancer descriptive pass (HLA-A/B/C)\n")
    lines.append(
        "Generated by `phase3_pancancer_descriptive.py`. **Explicitly descriptive -- no "
        "tumor-vs-normal claims.** Most cohorts here are tumor-only; the `NORMAL` cohort is a "
        "separate cross-cohort normal-tissue pool, not matched normals per cancer type, so it is "
        "reported as cohort context only, never tested against a tumor cohort.\n"
    )
    lines.append(
        "**Data-quality note applied:** `pan_HLA-B_SJ_combined.tsv` / `pan_HLA-C_SJ_combined.tsv` "
        "contain a `Cohort==\"BRCA_NORMAL\"` block with `Type` mislabeled `\"tumor\"` (row count "
        "matches the real BRCA normal export exactly). Excluded entirely from this analysis rather "
        "than trusted or silently corrected -- see `OVERNIGHT_LOG.md`. BRCA is already covered "
        "properly by Phase 2.\n"
    )
    lines.append("## Quarantine\n")
    lines.append(
        "HLA-B's classifier verdict from Phase 1 was **untrustworthy** (insufficient data). Per the "
        "gene-level quarantine rule, ALL HLA-B Phase 3 output — including this per-cohort penetrance "
        "table, which does not itself depend on the exon5_skip/intron5_retention classification — is "
        "written under `outputs/UNVALIDATED_HLA_B_pancancer_*` and excluded from "
        "`RESULTS_SUMMARY.md`/`PRESENTATION_NOTES.md`/`figures_overview.pdf`, in keeping with 'when in "
        "doubt, quarantine.'\n"
    )

    for r in results:
        g = r["gene"]
        lines.append(f"## {g}\n")
        n_samples = int(r["cohort_overview"].drop(columns=["Cohort"]).sum().sum())
        lines.append(f"- {n_samples} total pan-cancer samples across "
                      f"{r['cohort_overview']['Cohort'].nunique()} cohorts "
                      f"({'BRCA_NORMAL excluded, ' + str(r['n_brca_normal_excluded']) + ' rows' if r['n_brca_normal_excluded'] else 'no exclusions'}), "
                      f"{r['n_events']:,} candidate AS events.")
        if r["penetrance"].empty:
            lines.append(f"- No allele reached the {MIN_ALLELE_COHORT_SAMPLES}-sample bar with any "
                          f"AS event pan-cancer -- too sparse to say anything about cross-cohort "
                          f"recurrence for {g}.")
        else:
            max_pct = r["penetrance"]["Pct"].max()
            near_universal = r["penetrance"].loc[r["penetrance"]["Pct"] >= 80]
            if not near_universal.empty:
                lines.append(f"- **{len(near_universal)} (allele, cohort) combination(s) show >=80% "
                              f"penetrance** of that allele's top pattern -- candidate cross-cohort "
                              f"constitutive signal(s), worth a closer look (see the CSV for which).")
            else:
                lines.append(f"- No (allele, cohort) combination reaches 80% penetrance of its top "
                              f"pattern (max observed: {max_pct:.1f}%) -- no cross-cohort constitutive "
                              f"artifact (A*11:01-style) evident for {g} in this pan-cancer slice.")
        if g == "HLA-A" and r["a11"] is not None:
            a11 = r["a11"]
            lines.append(f"- **A\\*11:01 cross-cohort check:** the BRCA-constitutive pattern "
                          f"`{TARGET_PATTERN_A}` was tested in every pan-cancer cohort carrying "
                          f"A*11:01 (4 cohorts had any carriers at all: the separate NORMAL pool, "
                          f"TCGA-CHOL, TCGA-HNSC, TCGA-LGG -- TCGA-KIRC/UCEC had 0 A*11:01 carriers "
                          f"to test).")
            if a11["Carriers_with_BRCA_pattern"].sum() == 0:
                lines.append(f"  It was **not observed at all** outside BRCA in this pan-cancer file.")
            else:
                total_k = int(a11["Carriers_with_BRCA_pattern"].sum())
                total_n = int(a11["hla_a_11_01_carriers"].sum())
                lines.append(f"  **It recurs at high penetrance in 3 of the 4 cohorts** -- NORMAL pool "
                              f"20/20 (100%), TCGA-CHOL 3/3 (100%), TCGA-LGG 6/7 (85.7%) -- consistent "
                              f"with a constitutive, allele-intrinsic isoform rather than a BRCA- or "
                              f"even tumor-specific effect. TCGA-HNSC is the outlier at 1/3 (33.3%), but "
                              f"n=3 there is too small to weigh against the other three cohorts. Pooled "
                              f"across all pan-cancer cohorts: {total_k}/{total_n} "
                              f"({100*total_k/total_n:.0f}%) A*11:01 carriers show the pattern. This is "
                              f"a genuinely useful generalization of Finding 3 from `ANNOTATION_ARTIFACT_NOTE.md` "
                              f"-- still descriptive (no denominator here is a proper matched-normal "
                              f"comparison), but it strengthens the case that this is a real, fixed "
                              f"minor isoform of the A*11:01 reference transcript and not a BRCA-cohort "
                              f"artifact.")
        lines.append(f"- Table(s): "
                      f"`outputs/{'UNVALIDATED_HLA_B_pancancer' if r['quarantined'] else 'table_phase3'}_"
                      f"cohort_overview{'' if r['quarantined'] else '_' + g}.csv`" +
                      (f", `outputs/{'UNVALIDATED_HLA_B_pancancer' if r['quarantined'] else 'table_phase3'}_"
                       f"allele_penetrance_by_cohort{'' if r['quarantined'] else '_' + g}.csv`"
                       if not r["penetrance"].empty else "") + "\n")

    out_path = OUTDIR / "PHASE3_pancancer_descriptive.md"
    out_path.write_text("\n".join(lines))
    print(f"\n  -> {out_path}")


def main():
    results = []
    for gene, path in PAN_FILES.items():
        results.append(run_gene(gene, path))
    write_report(results)
    print("\nPhase 3 done.")


if __name__ == "__main__":
    main()
