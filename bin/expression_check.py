#!/usr/bin/env python3
"""Classify cells as B-like, T-like, both or neither from gene expression, flag
doublets, and cross-tabulate against each tool's VDJ calls.

B and T scores use only non-receptor marker genes (assets/markers.tsv, sets B
and T), so the expression evidence is independent of the VDJ transcripts. Ig and
TR constant-region expression is reported alongside but never used for the call.

Thresholds: a score is called positive above the crossing point of a two-component
Gaussian mixture when the score is clearly bimodal (BIC); otherwise above
--score-min. A B-enriched sample has no T population to fit, so the fallback is
the normal case for the T score there.

"both" cells are split by Scrublet's doublet call into `both_doublet` and
`both_singlet`, the latter being candidate B/T hybrid (dual-expresser) cells.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture

CLASS_ORDER = ["B", "T", "both_singlet", "both_doublet", "neither"]


def mixture_threshold(scores, score_min):
    """Return (threshold, method) for one marker score."""
    x = np.asarray(scores).reshape(-1, 1)
    one = GaussianMixture(1, random_state=0).fit(x)
    two = GaussianMixture(2, random_state=0).fit(x)
    lo, hi = np.argsort(two.means_.ravel())
    sd = np.sqrt(two.covariances_.ravel())
    separated = two.means_.ravel()[hi] - two.means_.ravel()[lo] > 2 * max(sd)
    if two.bic(x) < one.bic(x) - 10 and separated and min(two.weights_) > 0.01:
        grid = np.linspace(x.min(), x.max(), 2000).reshape(-1, 1)
        post_hi = two.predict_proba(grid)[:, hi]
        between = (grid.ravel() > two.means_.ravel()[lo]) & (grid.ravel() < two.means_.ravel()[hi])
        crossing = grid.ravel()[between][np.argmin(np.abs(post_hi[between] - 0.5))]
        return float(crossing), "mixture"
    return float(score_min), "fallback"


def gene_set(markers, name, var_names):
    genes = markers.loc[markers["set"] == name, "gene"].tolist()
    present = [g for g in genes if g in var_names]
    return present, sorted(set(genes) - set(present))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5", required=True, help="Cell Ranger sample_filtered_feature_bc_matrix.h5")
    parser.add_argument("--vdj", required=True, help="vdj_cell_classes.tsv from vdj_compare.py")
    parser.add_argument("--markers", required=True)
    parser.add_argument("--score-min", type=float, default=0.1)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-pct-mt", type=float, default=20.0)
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    out = Path(args.outdir)

    adata = sc.read_10x_h5(args.h5, gex_only=True)
    adata.var_names_make_unique()
    adata.obs_names = adata.obs_names.str.replace(r"-\d+$", "", regex=True)

    # QC is a flag, not a filter: low-quality cells still get a class.
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    adata.obs["qc_pass"] = (adata.obs["n_genes_by_counts"] >= args.min_genes) & (
        adata.obs["pct_counts_mt"] <= args.max_pct_mt
    )

    # Scrublet on raw counts. Its automatic threshold needs scikit-image; fall back
    # to a fixed one if that is missing, and record which was used.
    try:
        sc.pp.scrublet(adata)
        doublet_threshold = "auto"
    except ImportError:
        sc.pp.scrublet(adata, threshold=0.25)
        doublet_threshold = "0.25 (fixed: scikit-image missing)"

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    markers = pd.read_csv(args.markers, sep="\t", comment=None)
    summary = []
    thresholds = {}
    for name in ["B", "T"]:
        genes, missing = gene_set(markers, name, adata.var_names)
        sc.tl.score_genes(adata, genes, score_name=f"{name}_score", random_state=0)
        ref = adata.obs.loc[adata.obs["qc_pass"], f"{name}_score"]
        thresholds[name], method = mixture_threshold(ref, args.score_min)
        adata.obs[f"{name}_pos"] = adata.obs[f"{name}_score"] > thresholds[name]
        summary += [
            {"metric": f"{name}_genes_used", "value": ",".join(genes)},
            {"metric": f"{name}_genes_missing", "value": ",".join(missing)},
            {"metric": f"{name}_threshold", "value": f"{thresholds[name]:.3f} ({method})"},
        ]
    for name in ["Ig_const", "TR_const"]:
        genes, _ = gene_set(markers, name, adata.var_names)
        adata.obs[f"{name}_log_expr"] = np.asarray(adata[:, genes].X.sum(axis=1)).ravel() if genes else 0.0

    b, t, dbl = adata.obs["B_pos"], adata.obs["T_pos"], adata.obs["predicted_doublet"].astype(bool)
    adata.obs["expr_class"] = np.select(
        [b & ~t, t & ~b, b & t & ~dbl, b & t & dbl],
        ["B", "T", "both_singlet", "both_doublet"],
        default="neither",
    )

    cells = adata.obs[[
        "qc_pass", "n_genes_by_counts", "pct_counts_mt", "doublet_score", "predicted_doublet",
        "B_score", "T_score", "Ig_const_log_expr", "TR_const_log_expr", "expr_class",
    ]].rename_axis("cell_id").reset_index()
    cells.to_csv(out / "expression_cells.tsv", sep="\t", index=False, float_format="%.4f")

    # Cross-tab per tool. Cells a tool called VDJ for but that are not in the
    # expression matrix are counted as `not_in_gex` rather than dropped.
    vdj = pd.read_csv(args.vdj, sep="\t")
    crosstabs = []
    for tool, calls in vdj.groupby("tool"):
        merged = cells[["cell_id", "expr_class"]].merge(
            calls[["cell_id", "vdj_class"]], on="cell_id", how="outer"
        )
        merged["expr_class"] = merged["expr_class"].fillna("not_in_gex")
        merged["vdj_class"] = merged["vdj_class"].fillna("none")
        table = pd.crosstab(merged["expr_class"], merged["vdj_class"]).reset_index()
        table.insert(0, "tool", tool)
        crosstabs.append(table)
    pd.concat(crosstabs, ignore_index=True).fillna(0).to_csv(out / "expression_vs_vdj.tsv", sep="\t", index=False)

    counts = cells["expr_class"].value_counts()
    summary += [{"metric": "cells_total", "value": len(cells)},
                {"metric": "cells_qc_pass", "value": int(cells["qc_pass"].sum())},
                {"metric": "doublet_threshold", "value": doublet_threshold},
                {"metric": "predicted_doublets", "value": int(cells["predicted_doublet"].sum())}]
    summary += [{"metric": f"expr_{c}", "value": int(counts.get(c, 0))} for c in CLASS_ORDER]
    for tool, calls in vdj.groupby("tool"):
        with_vdj = cells["cell_id"].isin(calls["cell_id"])
        for c in CLASS_ORDER:
            in_class = cells["expr_class"] == c
            n = int(in_class.sum())
            frac = with_vdj[in_class].mean() if n else float("nan")
            summary.append({"metric": f"{tool}_vdj_fraction_{c}", "value": f"{frac:.3f}"})
    pd.DataFrame(summary).to_csv(out / "expression_summary.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, name in zip(axes[:2], ["B", "T"]):
        ax.hist(cells[f"{name}_score"], bins=80, color="0.4")
        ax.axvline(thresholds[name], color="crimson")
        ax.set(title=f"{name} score", xlabel="score", ylabel="cells")
    for c in CLASS_ORDER:
        sel = cells["expr_class"] == c
        axes[2].scatter(cells.loc[sel, "B_score"], cells.loc[sel, "T_score"], s=4, label=f"{c} ({sel.sum()})")
    axes[2].axvline(thresholds["B"], color="0.6", lw=0.8)
    axes[2].axhline(thresholds["T"], color="0.6", lw=0.8)
    axes[2].set(xlabel="B score", ylabel="T score", title="Expression class")
    axes[2].legend(markerscale=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "marker_scores.png", dpi=120)


if __name__ == "__main__":
    main()
