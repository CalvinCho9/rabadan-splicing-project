#!/usr/bin/env python3
"""
phase2_bc_headline.py
================================================================================
PHASE 2 (overnight B/C run): HLA-B/C tumor-vs-normal headline comparison,
per-allele view, and volcano -- the same statistics analyze_brca.py runs for
HLA-A (Fisher exact + Wilson CIs + odds ratio, BH-FDR across patterns),
applied to B and C, with the Phase 1 verdict and quarantine rule applied.

Verdicts from Phase 1 (outputs/PHASE1_BC_strand_check.md), reused here rather
than re-decided:
    HLA-C: valid_as_is             -> standard (unflipped) classifier,
                                       PRESENTED outputs.
    HLA-B: untrustworthy           -> standard (unflipped) classifier is still
                                       the best available (Phase 1 found naive
                                       weakly beats flipped, 1 vs 0 carriers --
                                       but that's on n=2 total events, nowhere
                                       near enough to validate), QUARANTINED:
                                       written only to outputs/UNVALIDATED_*,
                                       excluded from every presented document.

Constitutive-allele exclusion (the A*11:01-style check) is re-run per gene
here (not hardcoded) using the same penetrance/verdict logic as Phase 1 --
Phase 1 found none for either gene, but this script doesn't assume that; it
detects and excludes whatever's actually flagged, same as analyze_brca.py's
--exclude-alleles does for HLA-A.

Outputs (per gene; HLA-B's land under the UNVALIDATED_ prefix, all others
under the gene's own suffix):
    table_headline_classes_HLA-{B,C}.csv / fig_classes_HLA-{B,C}.pdf
    table_per_allele_HLA-{B,C}.csv       / fig_per_allele_HLA-{B,C}.pdf
    results_BRCA_tumor_vs_normal_HLA-{B,C}.csv  (detailed per-pattern + FDR)
    table_volcano_HLA-{B,C}.csv          / fig_volcano_HLA-{B,C}.pdf
================================================================================
"""
from pathlib import Path

import pandas as pd

import analyze_brca as ab
import hla_common as hc

OUTDIR = Path("outputs")
OUTDIR.mkdir(exist_ok=True)

GENES = {
    "HLA-B": ("data/export_BRCA_tumor_HLA-B.tsv", "data/export_BRCA_normal_HLA-B.tsv"),
    "HLA-C": ("data/export_BRCA_tumor_HLA-C.tsv", "data/export_BRCA_normal_HLA-C.tsv"),
}
# From Phase 1 (outputs/PHASE1_BC_strand_check.md) -- neither gene's flip
# beat naive, so both use the standard (unflipped) classifier; only the
# PRESENTATION status differs.
QUARANTINED = {"HLA-B": True, "HLA-C": False}

# "Talk" subset of the per-allele grid: story-carrying panels only (full grid
# is still written for every gene as the appendix figure). HLA-C: the
# Finding-4 real-variable-splicing candidate pair. HLA-B: largest-n alleles
# union the only two alleles with a nonzero carrier bar at all -- still
# quarantine-prefixed, audit convenience only, never referenced presented.
TALK_ALLELES = {
    "HLA-C": ["hla_c_03_03", "hla_c_03_04"],
    "HLA-B": ["hla_b_07_02", "hla_b_08_01", "hla_b_35_01", "hla_b_44_02", "hla_b_15_01", "hla_b_27_05"],
}


def gene_tag(gene):
    return gene.replace("-", "_")


def find_constitutive_alleles(df, events, gene):
    """Same A*11:01-style penetrance check as Phase 1 / diagnostic_annotation_check.py,
    re-run here so exclusion is derived, not assumed."""
    penetrance = hc.top_event_penetrance(df, events)
    penetrance = hc.add_verdict(penetrance)
    flagged = penetrance.loc[
        penetrance["Verdict"].str.startswith("constitutive_annotation_artifact"), "Contig_rd"
    ].tolist()
    print(f"  [{gene}] constitutive-artifact re-check: {len(flagged)} allele(s) flagged: {flagged or 'none'}")
    return flagged, penetrance


def run_gene(gene, tumor_path, normal_path):
    quarantined = QUARANTINED[gene]
    prefix = f"UNVALIDATED_{gene_tag(gene)}" if quarantined else None
    tag = gene_tag(gene)
    status = "QUARANTINED (not for presentation)" if quarantined else "validated, presented"
    print(f"\n{'=' * 80}\n{gene}  [{status}]\n{'=' * 80}")

    df = hc.build_gene_frame(tumor_path, normal_path)
    events = ab.find_as_events(df, local_rate_threshold=0.1)
    events = ab.annotate_events(df, events)
    events = hc.classify_events_strand_aware(events, flip_minus=False)  # both genes: naive wins Phase 1
    print(f"  {len(events):,} candidate AS events, "
          f"{df.loc[df['Type']=='tumor','Barcode'].nunique()} tumor / "
          f"{df.loc[df['Type']=='normal','Barcode'].nunique()} normal samples")
    print(events.groupby(["Type", "Event_class"]).size().unstack(fill_value=0))

    exclude_alleles, penetrance = find_constitutive_alleles(df, events, gene)
    pen_name = f"{prefix}_penetrance_recheck.csv" if quarantined else f"table_penetrance_recheck_HLA-{tag[4:]}.csv"
    penetrance.to_csv(OUTDIR / pen_name, index=False)

    def fname(base, ext):
        stem = f"{prefix}_{base}" if quarantined else f"{base}_HLA-{tag[4:]}"
        return OUTDIR / f"{stem}.{ext}"

    # ---- headline comparison ----------------------------------------------
    events_for_headline = events[~events["Contig_rd"].isin(exclude_alleles)] if exclude_alleles else events
    headline = ab.headline_comparison(df, events_for_headline)
    ab.print_headline_table(headline)
    headline.to_csv(fname("table_headline_classes", "csv"), index=False)
    ab.plot_classes(headline, fname("fig_classes", "pdf"), gene_label=gene)

    # ---- per-allele view ----------------------------------------------------
    # plot_per_allele()/plot_volcano() always write to the FIXED filenames
    # <outdir>/fig_per_allele.pdf / <outdir>/fig_volcano.pdf -- calling them
    # with OUTDIR directly would clobber HLA-A's own canonical
    # outputs/fig_per_allele.pdf / outputs/fig_volcano.pdf (whether or not
    # this gene actually produces a figure). Route through an isolated
    # scratch subdirectory instead and move the result into the gene-tagged
    # name in OUTDIR -- outputs/fig_per_allele.pdf and outputs/fig_volcano.pdf
    # (HLA-A's) are never touched by this script.
    scratch = OUTDIR / "_phase2_scratch"
    scratch.mkdir(exist_ok=True)

    per_allele, common_alleles = ab.per_allele_view(df, events, min_allele_samples=20,
                                                     exclude_alleles=exclude_alleles)
    per_allele.to_csv(fname("table_per_allele", "csv"), index=False)
    ab.plot_per_allele(per_allele, scratch, gene_label=gene)
    scratch_fig = scratch / "fig_per_allele.pdf"
    if scratch_fig.exists():
        scratch_fig.rename(fname("fig_per_allele", "pdf"))
    else:
        print(f"  [{gene}] no alleles meet --min-allele-samples -- no fig_per_allele for this gene")

    # "talk" subset -- same full-grid-is-unreadable-on-a-slide fix as HLA-A's
    # fig_per_allele_talk.pdf, restricted to this gene's story-carrying panels.
    talk_alleles = TALK_ALLELES.get(gene, [])
    talk_subset = per_allele[per_allele["Contig_rd"].isin(talk_alleles)]
    ab.plot_per_allele(talk_subset, scratch, gene_label=gene, filename="fig_per_allele_talk.pdf",
                        title_suffix=f"\n(talk view: story-carrying alleles only; full grid in "
                                      f"{fname('fig_per_allele', 'pdf').name})")
    scratch_talk = scratch / "fig_per_allele_talk.pdf"
    if scratch_talk.exists():
        scratch_talk.rename(fname("fig_per_allele_talk", "pdf"))

    # ---- detailed table (FDR) + volcano --------------------------------------
    detail = ab.detailed_table(df, events, min_carriers=3)
    detail.to_csv(fname("results_BRCA_tumor_vs_normal", "csv"), index=False)
    volcano = ab.volcano_data(detail, min_carriers=3)
    volcano.to_csv(fname("table_volcano", "csv"), index=False)
    ab.plot_volcano(volcano, scratch, gene_label=gene, label_n=6)
    scratch_volcano = scratch / "fig_volcano.pdf"
    if scratch_volcano.exists():
        scratch_volcano.rename(fname("fig_volcano", "pdf"))
    else:
        print(f"  [{gene}] no patterns meet --min-carriers -- no fig_volcano for this gene")

    scratch.rmdir()  # empty by now -- both files moved out (or never created)

    if quarantined:
        note_path = OUTDIR / f"{prefix}_README.md"
        note_path.write_text(
            f"# {gene} tumor-vs-normal — NOT FOR PRESENTATION\n\n"
            f"Quarantined per Phase 1 (`outputs/PHASE1_BC_strand_check.md`): the "
            f"exon5_skip/intron5_retention classifier could not be validated for {gene} -- "
            f"only 2 candidate AS events were found in the entire 935-tumor/106-normal "
            f"cohort at the standard 10%-local-rate threshold (vs 60 for HLA-C, ~200+ for "
            f"HLA-A), far too few to tell whether the naive or strand-flipped geometry is "
            f"correct, or whether the exon5_skip/intron5_retention labels mean anything at "
            f"all for this gene.\n\n"
            f"These files are the same computations Phase 2 runs for validated genes, kept "
            f"here for audit/reference ONLY. **Do not fold any number from these files into "
            f"RESULTS_SUMMARY.md, PRESENTATION_NOTES.md, or figures_overview.pdf** unless a "
            f"manual strand/geometry review specifically validates {gene}.\n\n"
            f"Pending manual review: (1) confirm with the alignment/reference-construction "
            f"pipeline (outside this repo) whether personalized HLA-B contigs are built "
            f"the same way as HLA-C's (Phase 1 found HLA-C's coordinates are already "
            f"transcript-oriented, contradicting the initial reverse-strand-flip hypothesis) "
            f"and (2) investigate why HLA-B's maximum observed Local_rate in this cohort is "
            f"only 12.6% (barely above the 10% detection floor) vs 52-60% for A/C -- real "
            f"biology or a technical/reference-quality issue specific to HLA-B.\n"
        )
        print(f"  -> {note_path} (quarantine note)")

    return {
        "gene": gene, "quarantined": quarantined, "n_events": len(events),
        "exclude_alleles": exclude_alleles, "headline": headline,
    }


def main():
    results = []
    for gene, (tp, np_) in GENES.items():
        results.append(run_gene(gene, tp, np_))

    print("\n\n=== Phase 2 summary ===")
    for r in results:
        status = "QUARANTINED" if r["quarantined"] else "PRESENTED"
        print(f"  {r['gene']}: {status}, {r['n_events']} candidate AS events, "
              f"excluded alleles: {r['exclude_alleles'] or 'none'}")
    print("\nPhase 2 done.")


if __name__ == "__main__":
    main()
