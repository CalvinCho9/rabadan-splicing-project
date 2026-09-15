#!/usr/bin/env python3
"""
hla_common.py
================================================================================
Shared helpers for extending the HLA-A pipeline (analyze_brca.py) to HLA-B/C
and to the pan-cancer files, added for the overnight B/C run (see
OVERNIGHT_LOG.md, outputs/PHASE1_BC_strand_check.md).

Does NOT modify analyze_brca.py -- imports it as a library (`import
analyze_brca as ab`) and reuses its loading / bucketing / stats functions
unchanged, so HLA-A's existing outputs are untouched and every gene's numbers
come from the exact same math. Adds only what's new for B/C:

  - classify_events_strand_aware(): the exon5_skip/intron5_retention geometry
    classifier, generalized to optionally mirror coordinates for minus-strand
    junctions before classifying (see module docstring on that function for
    the geometric argument). HLA-A is unaffected either way (it's ~98%
    strand=1 already) -- this is only consequential for B/C.
  - build_gene_frame(): generic tumor+normal loader for any gene's two files.
  - top_event_penetrance() / add_verdict() / read_share_label(): generalized
    (allele-parametrized, not hardcoded to hla_a_11_01) versions of the
    constitutive-vs-variable-splicing check from diagnostic_annotation_check.py,
    reused as-is here (same thresholds, same logic) so the A*11:01-style
    artifact check is identical across genes.
================================================================================
"""
from pathlib import Path

import numpy as np
import pandas as pd

import analyze_brca as ab

# Same thresholds diagnostic_annotation_check.py uses for HLA-A -- kept
# identical here so artifact calls are comparable across genes.
HIGH_PENETRANCE = 80
LOW_PENETRANCE = 5
DOMINANT_SHARE = 50
MIN_ALLELE_SAMPLES = 20


def build_gene_frame(tumor_path, normal_path):
    """Load one gene's tumor+normal files and compute the shared derived
    columns (Annotated_cumsum, Splice_coord, Local_rate, Contig_rd). Strand-
    and gene-agnostic -- identical to what analyze_brca.main() does before
    finding events."""
    df = ab.load_data([Path(tumor_path), Path(normal_path)])
    df = ab.compute_derived_columns(df)
    df["Contig_rd"] = df["Contig"].str.slice(0, 11)
    return df


def classify_events_strand_aware(events, exon_skip_min_jump=100,
                                  intron_retention_max_shift=50, flip_minus=False):
    """Same exon5_skip / intron5_retention geometry as analyze_brca.classify_events,
    with an optional strand-aware coordinate mirror applied first.

    Why a mirror is the right transform: First_base <= Last_base always holds
    in this data (verified: 0 violations across A/B/C), i.e. junction
    coordinates are stored in genomic left-to-right order regardless of gene
    strand -- NOT already re-oriented to "transcript 5'->3'" order. The A
    classifier's geometry (donor = First_base, acceptor = Last_base;
    exon5_skip = donor fixed + acceptor jumps downstream; intron5_retention =
    acceptor fixed + donor shifts a few bp) is written in transcript-relative
    terms and is only correct as-is for plus-strand junctions, where genomic
    left-to-right IS transcript 5'->3'. For a minus-strand junction,
    transcript 5'->3' runs genomic right-to-left, so donor/acceptor swap:
    donor = Last_base, acceptor = First_base, and "downstream" means
    DEcreasing genomic coordinate.

    To reuse the exact same classify_events() geometry unchanged, mirror each
    minus-strand row's coordinates before classifying:
        W1' = -W2, W2' = -W1   (and same for V1/V2)
    This negate-and-swap preserves W1' <= W2' (since W1<=W2 => -W2<=-W1) and
    makes "W1' increases" correspond to "genomic Last_base decreases" i.e.
    "moves further downstream along a minus-strand transcript" -- exactly the
    plus-strand donor-at-W1 assumption, now correct for minus-strand rows too.
    Plus-strand rows (Strand==1) are left untouched.

    flip_minus=False reproduces the naive (strand-blind) classification, used
    in Phase 1 purely as the comparison baseline -- it is not expected to be
    correct for B/C given ~95-97% of their junctions are Strand==2.
    """
    events = events.copy()
    w = events["Annotated_coord"].str.split("_", expand=True).astype(int)
    v = events["Splice_coord"].str.split("_", expand=True).astype(int)
    w.columns = ["W1", "W2"]
    v.columns = ["V1", "V2"]
    events = pd.concat([events, w, v], axis=1)

    if flip_minus:
        is_minus = (events["Strand"] == 2).to_numpy()
        W1 = np.where(is_minus, -events["W2"].to_numpy(), events["W1"].to_numpy())
        W2 = np.where(is_minus, -events["W1"].to_numpy(), events["W2"].to_numpy())
        V1 = np.where(is_minus, -events["V2"].to_numpy(), events["V1"].to_numpy())
        V2 = np.where(is_minus, -events["V1"].to_numpy(), events["V2"].to_numpy())
    else:
        W1 = events["W1"].to_numpy()
        W2 = events["W2"].to_numpy()
        V1 = events["V1"].to_numpy()
        V2 = events["V2"].to_numpy()

    is_exon5_skip = (V1 == W1) & ((V2 - W2) >= exon_skip_min_jump)
    shift = V1 - W1
    is_intron5_ret = (V2 == W2) & (shift > 0) & (shift <= intron_retention_max_shift)

    events["Event_class"] = np.select(
        [is_exon5_skip, is_intron5_ret],
        ["exon5_skip", "intron5_retention"],
        default="other",
    )
    return events


# ==============================================================================
# Generalized constitutive-vs-variable-splicing check (from
# diagnostic_annotation_check.py, TARGET_ALLELE hardcoding removed -- same
# thresholds/logic, applies to every allele of any gene).
# ==============================================================================
def top_event_penetrance(df, events, min_allele_samples=MIN_ALLELE_SAMPLES):
    contig_denom = df.groupby(["Contig_rd", "Type"])["Barcode"].nunique().unstack(fill_value=0)
    for t in ("tumor", "normal"):
        if t not in contig_denom.columns:
            contig_denom[t] = 0
    contig_denom["Total"] = contig_denom["tumor"] + contig_denom["normal"]
    alleles = contig_denom.loc[contig_denom["Total"] >= min_allele_samples].index

    rows = []
    for allele in alleles:
        sub = events[events["Contig_rd"] == allele]
        if sub.empty:
            continue
        carrier_counts = sub.groupby("AS_pattern")["Barcode"].nunique().sort_values(ascending=False)
        top_pattern = carrier_counts.index[0]
        top_sub = sub[sub["AS_pattern"] == top_pattern]
        top_class = top_sub["Event_class"].iloc[0] if "Event_class" in top_sub.columns else "n/a"

        nt = int(contig_denom.loc[allele, "tumor"])
        nn = int(contig_denom.loc[allele, "normal"])
        carriers = top_sub.groupby("Type")["Barcode"].nunique()
        a = int(carriers.get("tumor", 0))
        c = int(carriers.get("normal", 0))

        rows.append({
            "Contig_rd": allele,
            "Top_AS_pattern": top_pattern,
            "Top_event_class": top_class,
            "Tumor_carriers": a, "Tumor_n": nt, "Tumor_pct": 100 * a / nt if nt else np.nan,
            "Normal_carriers": c, "Normal_n": nn, "Normal_pct": 100 * c / nn if nn else np.nan,
            "Min_group_n": min(nt, nn),
            "Avg_variant_local_rate_pct": 100 * top_sub["Local_rate"].mean(),
        })
    if not rows:
        # No allele reaches min_allele_samples with any qualifying event --
        # happens on small single-cohort slices (added for Phase 5's
        # per-cohort pan-cancer pass; never triggered on the full BRCA
        # cohort Phase 1/2 use this for). Return an empty-but-typed frame
        # so callers can check .empty rather than hit a KeyError from
        # sort_values() on a columnless DataFrame.
        return pd.DataFrame(columns=["Contig_rd", "Top_AS_pattern", "Top_event_class",
                                      "Tumor_carriers", "Tumor_n", "Tumor_pct",
                                      "Normal_carriers", "Normal_n", "Normal_pct",
                                      "Min_group_n", "Avg_variant_local_rate_pct"])
    result = pd.DataFrame(rows).sort_values(["Tumor_pct", "Normal_pct"], ascending=False).reset_index(drop=True)
    return result


def read_share_label(pct, dominant=DOMINANT_SHARE):
    if pd.isna(pct):
        return "n/a"
    return ("novel junction dominates local reads (>=%d%%) -- looks like a coordinate/GTF mismatch" % dominant
            if pct >= dominant else
            "annotated junction stays dominant -- novel is a minor, fairly fixed fraction")


def add_verdict(table, high=HIGH_PENETRANCE, low=LOW_PENETRANCE, min_n=10):
    def verdict(row):
        pairs = [(row["Tumor_pct"], row["Tumor_n"]), (row["Normal_pct"], row["Normal_n"])]
        pairs = [(p, n) for p, n in pairs if pd.notna(p) and n > 0]
        if not pairs:
            return "insufficient_data"
        pcts = [p for p, _ in pairs]
        if max(pcts) < low:
            return "rare_private_event"
        if min(pcts) >= high:
            return "constitutive_annotation_artifact"
        well_powered_pcts = [p for p, n in pairs if n >= min_n]
        if well_powered_pcts and max(well_powered_pcts) >= high and any(n < min_n for _, n in pairs):
            return "constitutive_annotation_artifact (one group underpowered)"
        return "real_variable_splicing_candidate"

    table = table.copy()
    table["Verdict"] = table.apply(verdict, axis=1)
    table["Read_share_pattern"] = table["Avg_variant_local_rate_pct"].apply(read_share_label)
    return table
