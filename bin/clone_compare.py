#!/usr/bin/env python3
"""Clone-level comparison of the VDJ methods for one sample.

The cell-level agreement in vdj_pairwise_agreement.tsv weights every cell equally,
so in a CLL sample one clone of thousands of cells decides it. Here every clonotype
counts once, and the question is whether two methods find the same clones, give
them the same identity, and group the same cells into them.

Input: vdj_calls_long.tsv from vdj_compare.py. A cell's clonotype comes from its top
chain (most UMIs) of the anchor locus, IGH for BCR and TRB for TCR, as in
clone_tracking.py. Four clonotype definitions:
  strict   V gene | J gene | CDR3 aa (without the conserved C and W/F)
  paralog  as strict, with duplicated V genes collapsed (IGHV1-69D = IGHV1-69)
  cdr3     CDR3 only
  lineage  clonal family: same paralog-collapsed V, same J, same CDR3 length, and
           CDR3s joined by single linkage at >= --lineage-identity amino-acid
           identity (0.85), clustered over all methods of the sample together so
           the labels are shared. In CLL samples the exact definitions split the
           CLL clone into dozens of one-residue CDR3 variants (sequencing errors,
           or real subclones when several methods report them); this one does not.

Outputs (prefixed with the sample by the pipeline):
  clone_summary.tsv   per method: cells, clonotypes, expanded clonotypes (>= 2 cells),
                      share of cells in expanded clonotypes
  clone_pairwise.tsv  per method pair:
    clones_a/_b, shared_clones, jaccard   clonotype sets, by identity
    expanded_a_found_in_b (and b in a)    share of a method's expanded clonotypes
                                          the other also reports, at any size
    cells_both                            cells with a clonotype from both methods
    same_clone_cells                      cell-weighted: share of those cells with the
                                          same clonotype (as in the cell-level table)
    same_clone_clones                     clone-weighted: per clonotype of either
                                          method, the share of its cells the other
                                          method gives the same clonotype, averaged
                                          over clonotypes
    purity_a_in_b, purity_b_in_a          grouping, independent of naming: per
                                          clonotype with >= 2 shared cells, the share
                                          of its cells in the other method's largest
                                          matching clonotype, averaged over clonotypes.
                                          Low purity_a_in_b: b splits a's clones.
    ari                                   adjusted Rand index of the two groupings of
                                          the shared cells (pair-weighted, so large
                                          clones dominate it)
  clone_table.tsv     the 30 largest clonotypes (by any method), cells per method
"""

import argparse
import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd

ANCHOR = {"BCR": "IGH", "TCR": "TRB"}
DEFINITIONS = ["strict", "paralog", "cdr3", "lineage"]
TOP_N = 30


def collapse_duplicate(gene):
    """As in vdj_compare.py: IGKV1D-39 -> IGKV1-39, IGHV1-69D -> IGHV1-69,
    TRAV14/DV4 -> TRAV14."""
    g = re.sub(r"^(IG[HKL]V\d+)D-", r"\1-", gene)
    g = re.sub(r"(-\d+)D$", r"\1", g)
    return re.sub(r"/?DV\d+$", "", g)


def lineages(keys, min_identity):
    """Map each (v, j, cdr3) key to a lineage label: single linkage within the same
    V, J and CDR3 length, joining CDR3s with >= min_identity amino-acid identity.
    keys: one row per key with its cell count over all methods ("cells")."""
    label = {}
    for (v, j, length), group in keys.groupby(["v", "j", "length"]):
        cells = group.groupby("cdr3")["cells"].sum()
        cdr3s = sorted(cells.index)
        parent = list(range(len(cdr3s)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        max_mismatch = int(np.floor((1 - min_identity) * length + 1e-9))
        for x, y in itertools.combinations(range(len(cdr3s)), 2):
            if sum(a != b for a, b in zip(cdr3s[x], cdr3s[y])) <= max_mismatch:
                parent[find(x)] = find(y)
        members = {}
        for i in range(len(cdr3s)):
            members.setdefault(find(i), []).append(cdr3s[i])
        for family in members.values():
            # Label: the family's most common CDR3 (ties: first in sort order), "~"
            # marking a lineage rather than an exact clonotype.
            top = max(family, key=lambda c: (cells[c], [-ord(ch) for ch in c]))
            for c in family:
                label[(v, j, c)] = f"{v}|{j}|{top}~"
    return label


def clonotypes(calls, min_identity):
    """One row per method, receptor and cell, with the four clonotype keys."""
    rows = []
    for receptor, locus in ANCHOR.items():
        top = calls[(calls["rank"] == 1) & (calls["locus"] == locus)].copy()
        top["receptor"] = receptor
        rows.append(top)
    cells = pd.concat(rows, ignore_index=True).fillna({"v_gene": "", "j_gene": "", "cdr3_core": ""})
    cells = cells[cells["cdr3_core"] != ""].copy()
    cells["strict"] = cells["v_gene"] + "|" + cells["j_gene"] + "|" + cells["cdr3_core"]
    cells["paralog"] = cells["v_gene"].map(collapse_duplicate) + "|" + cells["j_gene"] + "|" + cells["cdr3_core"]
    cells["cdr3"] = cells["cdr3_core"]
    cells["v"] = cells["v_gene"].map(collapse_duplicate)
    keys = (cells.groupby(["v", "j_gene", "cdr3_core"]).size().rename("cells").reset_index()
            .rename(columns={"j_gene": "j", "cdr3_core": "cdr3"}))
    keys["length"] = keys["cdr3"].str.len()
    label = lineages(keys, min_identity)
    cells["lineage"] = [label[k] for k in zip(cells["v"], cells["j_gene"], cells["cdr3_core"])]
    return cells[["tool", "receptor", "cell_id"] + DEFINITIONS].drop_duplicates(["tool", "receptor", "cell_id"])


def comb2(n):
    n = np.asarray(n, dtype=float)
    return (n * (n - 1) / 2).sum()


def ari(table):
    """Adjusted Rand index from a contingency table (rows: a's clones, cols: b's)."""
    n = table.to_numpy().sum()
    if n < 2:
        return float("nan")
    index = comb2(table.to_numpy().ravel())
    a, b = comb2(table.sum(axis=1)), comb2(table.sum(axis=0))
    expected = a * b / comb2([n])
    top = (a + b) / 2
    return float("nan") if top == expected else (index - expected) / (top - expected)


def purity(table):
    """Clone-weighted share of each row clone's cells in its largest column clone,
    over row clones with >= 2 cells."""
    sizes = table.sum(axis=1)
    keep = sizes >= 2
    if not keep.any():
        return float("nan")
    return float((table[keep].max(axis=1) / sizes[keep]).mean())


def pair_metrics(a, b):
    """a, b: Series cell_id -> clonotype for one receptor and definition."""
    clones_a, clones_b = set(a), set(b)
    shared = clones_a & clones_b
    size_a, size_b = a.value_counts(), b.value_counts()
    exp_a, exp_b = set(size_a[size_a >= 2].index), set(size_b[size_b >= 2].index)
    row = {
        "clones_a": len(clones_a), "clones_b": len(clones_b), "shared_clones": len(shared),
        "jaccard": len(shared) / len(clones_a | clones_b) if clones_a | clones_b else float("nan"),
        "expanded_a": len(exp_a), "expanded_b": len(exp_b),
        "expanded_a_found_in_b": len(exp_a & clones_b) / len(exp_a) if exp_a else float("nan"),
        "expanded_b_found_in_a": len(exp_b & clones_a) / len(exp_b) if exp_b else float("nan"),
    }
    both = a.index.intersection(b.index)
    row["cells_both"] = len(both)
    if len(both) == 0:
        return row
    pa, pb = a[both], b[both]
    same = pa.values == pb.values
    row["same_clone_cells"] = float(same.mean())
    # Clone-weighted identity: each clonotype of either method counts once.
    per_clone = pd.concat([pd.Series(same, index=pa.values).groupby(level=0).mean(),
                           pd.Series(same, index=pb.values).groupby(level=0).mean()])
    row["same_clone_clones"] = float(per_clone.mean())
    table = pd.crosstab(pa.values, pb.values)
    row["purity_a_in_b"] = purity(table)
    row["purity_b_in_a"] = purity(table.T)
    row["ari"] = ari(table)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, help="vdj_calls_long.tsv")
    ap.add_argument("--lineage-identity", type=float, default=0.85,
                    help="CDR3 amino-acid identity for joining clonotypes into a lineage")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    out = Path(args.outdir)

    calls = pd.read_csv(args.calls, sep="\t", dtype={"cell_id": str})
    cells = clonotypes(calls, args.lineage_identity)
    tools = sorted(cells["tool"].unique())

    summary, pairwise, tables = [], [], []
    for receptor, by_receptor in cells.groupby("receptor"):
        keyed = {t: g.set_index("cell_id") for t, g in by_receptor.groupby("tool")}
        for definition in DEFINITIONS:
            for tool, g in keyed.items():
                sizes = g[definition].value_counts()
                summary.append({
                    "receptor": receptor, "definition": definition, "tool": tool,
                    "cells": int(sizes.sum()), "clonotypes": len(sizes),
                    "expanded_clonotypes": int((sizes >= 2).sum()),
                    "cells_in_expanded": float(sizes[sizes >= 2].sum() / sizes.sum()),
                })
            for ta, tb in itertools.combinations([t for t in tools if t in keyed], 2):
                row = {"receptor": receptor, "definition": definition, "tool_a": ta, "tool_b": tb}
                row.update(pair_metrics(keyed[ta][definition], keyed[tb][definition]))
                pairwise.append(row)
            # Largest clonotypes by any method, cells per method.
            counts = pd.DataFrame({t: g[definition].value_counts() for t, g in keyed.items()}).fillna(0).astype(int)
            counts = counts.loc[counts.max(axis=1).sort_values(ascending=False).index[:TOP_N]]
            counts.index.name = "clonotype"
            counts = counts.reset_index()
            counts.insert(0, "definition", definition)
            counts.insert(0, "receptor", receptor)
            tables.append(counts)

    pd.DataFrame(summary).to_csv(out / "clone_summary.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(pairwise).to_csv(out / "clone_pairwise.tsv", sep="\t", index=False, float_format="%.4f")
    pd.concat(tables, ignore_index=True).to_csv(out / "clone_table.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
