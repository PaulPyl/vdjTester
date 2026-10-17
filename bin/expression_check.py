#!/usr/bin/env python3
"""Classify cells as B-like, T-like, both or neither from gene expression, flag
doublets, and cross-tabulate against each tool's VDJ calls.

B and T scores use only non-receptor marker genes (assets/markers.tsv, sets B
and T), so the expression evidence is independent of the VDJ transcripts. Ig and
TR constant-region expression is reported alongside but never used for the call.

Scores are UCell (Andreatta & Carmona 2021) by default: within each cell, genes are
ranked by expression (ties averaged) and genes ranked at or beyond maxRank count as
not expressed; score = 1 - U / (n * maxRank - n(n+1)/2), where U is the
Mann-Whitney U statistic of the n signature genes' ranks (as in UCell 2.x). It depends only on the cell itself, not on the
rest of the dataset, so thresholds transfer between samples. maxRank is set once
per dataset from the data: the median number of detected genes per cell (QC-pass
cells), so that genes a typical cell does detect are not discarded. The 10th
percentile is computed too, and the Spearman correlation between the two settings
is reported as a sensitivity check. `--score-method scanpy` keeps the earlier
sc.tl.score_genes scores for comparison.

A cell is positive for a set only if at least --min-markers different genes of
that set are detected AND its score clears the threshold. On the 10x 1k B-cell
sample, ~120 B cells crossed the T threshold on LCK alone (no CD3D/E/G or CD2),
which a score threshold by itself would have called B/T hybrids.

Thresholds: a score is called positive above the crossing point of a two-component
Gaussian mixture when the score is clearly bimodal (BIC); otherwise above
--score-min, which is also a floor under the mixture threshold. A B-enriched sample
has no T population to fit, so the fallback is the normal case for the T score there.
Default floors for UCell: 0.2 for T (1k B cells: LCK-only cells score T ~0.06, max
0.14) and 0.15 for B (--score-min-b); 0.1 for both with scanpy scores.

The B set has 19 genes. The first 7 (CD79A/B, MS4A1, CD19, CD22, PAX5, BANK1) missed
CLL cells, which are dim for CD20, CD22 and CD79B: only 56-88 % of BCR-only cells
in cll_1 and cll_5 were called B. With the 12 added B-lineage genes and the B floor
lowered from 0.2 to 0.15 (a longer signature gives lower UCell scores for the same
evidence), 93-99 % are, while at most 1.3 % of TCR-only cells, 1 of 520 myeloid-like
and 1 of 88 NK-like cells are called B; healthy 1k B cells stay at 100 % (median
score 0.57). Scratch analysis 2026-10-17, see docs/PLAN.md.

"both" cells are split by Scrublet's doublet call into `both_doublet` and
`both_singlet`, the latter being candidate B/T hybrid (dual-expresser) cells.
"""

import argparse
import importlib.util
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
        crossing = float(grid.ravel()[between][np.argmin(np.abs(post_hi[between] - 0.5))])
        if crossing < score_min:
            # A mixture can split exact zeros from a sliver of near-zero scores
            # (e.g. single-UMI LCK in B cells); never go below the floor.
            return float(score_min), f"mixture at {crossing:.3f}, floored"
        return crossing, "mixture"
    return float(score_min), "fallback"


def ucell(counts, sig_cols, max_rank):
    """UCell scores from a cells x genes CSR matrix of raw counts.

    Ranks are within-cell and descending, ties averaged; genes absent from a cell
    tie at the bottom. Matches R UCell 2.16 (checked to 1e-12 on the 1k B-cell sample).
    Normalisation and log transform do not change within-cell ranks, so raw
    counts give the same scores as normalised data.
    """
    counts = counts.tocsr()
    n_genes = counts.shape[1]
    scores = {name: np.zeros(counts.shape[0]) for name in sig_cols}
    for i in range(counts.shape[0]):
        row = counts.getrow(i)
        nonzero = np.sort(row.data)
        nnz = len(nonzero)
        for name, cols in sig_cols.items():
            n = len(cols)
            values = row[:, cols].toarray().ravel()
            ranks = np.empty(n)
            for j, v in enumerate(values):
                if v > 0:
                    greater = nnz - np.searchsorted(nonzero, v, side="right")
                    equal = np.searchsorted(nonzero, v, side="right") - np.searchsorted(nonzero, v, side="left")
                else:
                    greater, equal = nnz, n_genes - nnz
                ranks[j] = greater + (equal + 1) / 2
            # As in UCell 2.x (u_stat): ranks >= maxRank count as maxRank, and the
            # score is scaled so a cell expressing none of the genes scores 0.
            ranks = np.minimum(ranks, max_rank)
            rank_sum_min = n * (n + 1) / 2
            scores[name][i] = 1 - (ranks.sum() - rank_sum_min) / (n * max_rank - rank_sum_min)
    return scores


def gene_set(markers, name, var_names):
    genes = markers.loc[markers["set"] == name, "gene"].tolist()
    present = [g for g in genes if g in var_names]
    return present, sorted(set(genes) - set(present))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5", required=True, help="Cell Ranger sample_filtered_feature_bc_matrix.h5")
    parser.add_argument("--vdj", required=True, help="vdj_cell_classes.tsv from vdj_compare.py")
    parser.add_argument("--markers", required=True)
    parser.add_argument("--score-method", choices=["ucell", "scanpy"], default="ucell")
    parser.add_argument("--max-rank", type=int, default=None,
                        help="UCell maxRank; default: median detected genes per QC-pass cell")
    parser.add_argument("--score-min", type=float, default=None,
                        help="T score floor / fallback threshold; default 0.2 (ucell) or 0.1 (scanpy)")
    parser.add_argument("--score-min-b", type=float, default=None,
                        help="B score floor / fallback threshold; default 0.15 (ucell) or 0.1 (scanpy)")
    parser.add_argument("--min-markers", type=int, default=3,
                        help="distinct marker genes of a set that must be detected")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-pct-mt", type=float, default=20.0)
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    out = Path(args.outdir)
    if args.score_min is None:
        args.score_min = 0.2 if args.score_method == "ucell" else 0.1
    if args.score_min_b is None:
        args.score_min_b = 0.15 if args.score_method == "ucell" else 0.1
    floors = {"B": args.score_min_b, "T": args.score_min}

    adata = sc.read_10x_h5(args.h5, gex_only=True)
    adata.var_names_make_unique()
    adata.obs_names = adata.obs_names.str.replace(r"-\d+$", "", regex=True)

    # QC is a flag, not a filter: low-quality cells still get a class.
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    adata.obs["qc_pass"] = (adata.obs["n_genes_by_counts"] >= args.min_genes) & (
        adata.obs["pct_counts_mt"] <= args.max_pct_mt
    )

    # Scrublet on raw counts. Its automatic threshold needs scikit-image (scanpy
    # raises ValueError without it); fall back to a fixed one and record which.
    if importlib.util.find_spec("skimage"):
        sc.pp.scrublet(adata)
        doublet_threshold = "auto"
    else:
        sc.pp.scrublet(adata, threshold=0.25)
        doublet_threshold = "0.25 (fixed: scikit-image missing)"

    markers = pd.read_csv(args.markers, sep="\t", comment=None)
    for name in ["B", "T"]:
        genes, _ = gene_set(markers, name, adata.var_names)
        adata.obs[f"{name}_markers_detected"] = np.asarray((adata[:, genes].X > 0).sum(axis=1)).ravel()

    summary = []
    if args.score_method == "ucell":
        detected = adata.obs.loc[adata.obs["qc_pass"], "n_genes_by_counts"]
        max_rank = args.max_rank or int(round(detected.median()))
        max_rank_low = int(round(detected.quantile(0.10)))
        sig_cols = {name: [adata.var_names.get_loc(g) for g in gene_set(markers, name, adata.var_names)[0]]
                    for name in ["B", "T"]}
        if max(len(c) for c in sig_cols.values()) >= max_rank_low:
            raise SystemExit("a signature is longer than the 10th percentile of detected genes")
        main_scores = ucell(adata.X, sig_cols, max_rank)
        low_scores = ucell(adata.X, sig_cols, max_rank_low)
        for name in ["B", "T"]:
            adata.obs[f"{name}_score"] = main_scores[name]
            rho = pd.Series(main_scores[name]).corr(pd.Series(low_scores[name]), method="spearman")
            summary.append({"metric": f"{name}_spearman_maxrank_{max_rank}_vs_{max_rank_low}", "value": f"{rho:.4f}"})
        summary = [{"metric": "score_method", "value": "UCell"},
                   {"metric": "max_rank", "value": f"{max_rank} (median detected genes; 10th pct {max_rank_low})"}] + summary

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    thresholds = {}
    for name in ["B", "T"]:
        genes, missing = gene_set(markers, name, adata.var_names)
        if args.score_method == "scanpy":
            sc.tl.score_genes(adata, genes, score_name=f"{name}_score", random_state=0)
        ref = adata.obs.loc[adata.obs["qc_pass"], f"{name}_score"]
        thresholds[name], method = mixture_threshold(ref, floors[name])
        adata.obs[f"{name}_pos"] = (adata.obs[f"{name}_score"] > thresholds[name]) & (
            adata.obs[f"{name}_markers_detected"] >= args.min_markers
        )
        above_only = (adata.obs[f"{name}_score"] > thresholds[name]) & ~adata.obs[f"{name}_pos"]
        summary += [
            {"metric": f"{name}_genes_used", "value": ",".join(genes)},
            {"metric": f"{name}_genes_missing", "value": ",".join(missing)},
            {"metric": f"{name}_threshold", "value": f"{thresholds[name]:.3f} ({method})"},
            {"metric": f"{name}_above_threshold_but_too_few_markers", "value": int(above_only.sum())},
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

    # UMAP for the report. Immunoglobulin and TCR genes are kept out of the variable
    # genes, so a dominant clone's own transcripts do not shape the embedding.
    adata.obs["umap1"] = np.nan
    adata.obs["umap2"] = np.nan
    if adata.n_obs >= 50:
        emb = adata.copy()
        receptor = emb.var_names.str.match(r"^(IG[HKL][VDJC]|TR[ABDG][VDJC])")
        sc.pp.highly_variable_genes(emb, n_top_genes=2000)
        emb.var.loc[receptor, "highly_variable"] = False
        sc.pp.pca(emb, n_comps=30, mask_var="highly_variable", random_state=0)
        sc.pp.neighbors(emb, n_neighbors=15, random_state=0)
        sc.tl.umap(emb, random_state=0)
        adata.obs[["umap1", "umap2"]] = emb.obsm["X_umap"]
        summary.append({"metric": "umap", "value": "2000 HVGs without IG/TR genes, 30 PCs, 15 neighbours"})

    cells = adata.obs[[
        "qc_pass", "n_genes_by_counts", "pct_counts_mt", "doublet_score", "predicted_doublet",
        "B_score", "T_score", "B_markers_detected", "T_markers_detected", "Ig_const_log_expr", "TR_const_log_expr", "expr_class",
        "umap1", "umap2",
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
