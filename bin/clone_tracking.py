#!/usr/bin/env python3
"""Clone sizes per timepoint and their change over time, per tool, for one patient.

Input: one vdj_calls_long.tsv (from vdj_compare.py) per sample, with its timepoint:
  --sample <sample_id> <timepoint> <vdj_calls_long.tsv>   (repeat per timepoint)

A cell's clonotype comes from its top chain (most UMIs) of the receptor's
"anchor" locus: IGH for BCR, TRB for TCR. Two clonotype definitions are reported:
  strict  V gene | J gene | CDR3 (aa, without conserved C and W/F)
  cdr3    CDR3 only, which is robust to the tools' different V-gene naming
          (e.g. IGHV1-69 vs IGHV1-69D) and so is the one to compare tools on

Outputs (all long-format TSV):
  clonality.tsv                 cells, clonotypes, top-1 / top-10 share, Simpson
                                index, Shannon clonality (1 - H / log N)
  top_clones.tsv                the 10 largest clonotypes per tool, sample, receptor
  clone_tracking.tsv            share of each clonotype that is top-10 at any
                                timepoint, at every timepoint (0 when absent)
  timepoint_overlap.tsv         Morisita-Horn overlap between timepoint pairs
  dominant_clone_agreement.tsv  per sample and receptor, each tool's largest clone,
                                how many clonotypes tie for largest, and whether
                                all tools picked the same one
"""

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

ANCHOR = {"BCR": "IGH", "TCR": "TRB"}
TOP_N = 10


def load(samples):
    frames = []
    for sample_id, timepoint, path in samples:
        df = pd.read_csv(path, sep="\t", dtype={"cell_id": str})
        df.insert(0, "timepoint", int(timepoint))
        df.insert(0, "sample_id", sample_id)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def clonotypes(calls):
    """One row per tool, sample, receptor and cell, with both clonotype keys."""
    rows = []
    for receptor, locus in ANCHOR.items():
        top = calls[(calls["rank"] == 1) & (calls["locus"] == locus)].copy()
        top["receptor"] = receptor
        rows.append(top)
    cells = pd.concat(rows, ignore_index=True).fillna({"v_gene": "", "j_gene": "", "cdr3_core": ""})
    cells = cells[cells["cdr3_core"] != ""]
    cells["strict"] = cells["v_gene"] + "|" + cells["j_gene"] + "|" + cells["cdr3_core"]
    cells["cdr3"] = cells["cdr3_core"]
    return cells[["tool", "sample_id", "timepoint", "receptor", "cell_id", "strict", "cdr3"]]


def shares(cells, key):
    counts = cells.groupby(["tool", "sample_id", "timepoint", "receptor", key]).size().rename("cells").reset_index()
    totals = counts.groupby(["tool", "sample_id", "receptor"])["cells"].transform("sum")
    counts["share"] = counts["cells"] / totals
    return counts.rename(columns={key: "clonotype"})


def clonality(counts):
    rows = []
    for (tool, sid, tp, receptor), df in counts.groupby(["tool", "sample_id", "timepoint", "receptor"]):
        p = df["share"].to_numpy()
        n = len(p)
        h = -(p * np.log(p)).sum()
        top = np.sort(p)[::-1]
        rows.append({
            "tool": tool, "sample_id": sid, "timepoint": tp, "receptor": receptor,
            "cells": int(df["cells"].sum()), "clonotypes": n,
            "top1_share": top[0], "top10_share": top[:TOP_N].sum(),
            "simpson": (p ** 2).sum(),
            "shannon_clonality": 1 - h / np.log(n) if n > 1 else 1.0,
        })
    return pd.DataFrame(rows)


def top_clones(counts):
    counts = counts.sort_values("cells", ascending=False)
    top = counts.groupby(["tool", "sample_id", "receptor"]).head(TOP_N).copy()
    top["rank"] = top.groupby(["tool", "sample_id", "receptor"]).cumcount() + 1
    return top.sort_values(["tool", "receptor", "timepoint", "rank"])


def tracking(counts, top):
    tracked = top[["tool", "receptor", "clonotype"]].drop_duplicates()
    timepoints = counts[["sample_id", "timepoint"]].drop_duplicates()
    grid = tracked.merge(timepoints, how="cross")
    out = grid.merge(counts[["tool", "receptor", "clonotype", "sample_id", "timepoint", "cells", "share"]],
                     on=["tool", "receptor", "clonotype", "sample_id", "timepoint"], how="left")
    return out.fillna({"cells": 0, "share": 0.0}).sort_values(["tool", "receptor", "clonotype", "timepoint"])


def morisita_horn(a, b):
    keys = a.index.union(b.index)
    x = a.reindex(keys, fill_value=0).to_numpy(float)
    y = b.reindex(keys, fill_value=0).to_numpy(float)
    X, Y = x.sum(), y.sum()
    if X == 0 or Y == 0:
        return float("nan")
    dx, dy = (x ** 2).sum() / X ** 2, (y ** 2).sum() / Y ** 2
    return 2 * (x * y).sum() / ((dx + dy) * X * Y)


def overlap(counts):
    rows = []
    for (tool, receptor), df in counts.groupby(["tool", "receptor"]):
        by_tp = {tp: d.set_index("clonotype")["cells"] for tp, d in df.groupby("timepoint")}
        for t1, t2 in itertools.combinations(sorted(by_tp), 2):
            shared = len(by_tp[t1].index.intersection(by_tp[t2].index))
            rows.append({"tool": tool, "receptor": receptor, "timepoint_a": t1, "timepoint_b": t2,
                         "shared_clonotypes": shared, "morisita_horn": morisita_horn(by_tp[t1], by_tp[t2])})
    return pd.DataFrame(rows)


def dominant(counts):
    rows = []
    for (sid, tp, receptor), df in counts.groupby(["sample_id", "timepoint", "receptor"]):
        best = df.sort_values(["cells", "clonotype"], ascending=[False, True]).groupby("tool").head(1).set_index("tool")
        same = best["clonotype"].nunique() == 1
        for tool, r in best.iterrows():
            # Clonotypes tied for largest: with ties, "dominant" and agreement mean little.
            tied = int((df.loc[df["tool"] == tool, "cells"] == r["cells"]).sum())
            rows.append({"sample_id": sid, "timepoint": tp, "receptor": receptor, "tool": tool,
                         "dominant_clonotype": r["clonotype"], "cells": int(r["cells"]), "share": r["share"],
                         "tied_for_largest": tied, "all_tools_agree": same})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", nargs=3, action="append", required=True,
                        metavar=("SAMPLE_ID", "TIMEPOINT", "CALLS"))
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    out = Path(args.outdir)

    cells = clonotypes(load(args.sample))
    tables = {name: [] for name in ["clonality", "top_clones", "clone_tracking", "timepoint_overlap",
                                    "dominant_clone_agreement"]}
    for key in ["strict", "cdr3"]:
        counts = shares(cells, key)
        top = top_clones(counts)
        for name, df in [("clonality", clonality(counts)), ("top_clones", top),
                         ("clone_tracking", tracking(counts, top)), ("timepoint_overlap", overlap(counts)),
                         ("dominant_clone_agreement", dominant(counts))]:
            if len(df):
                df.insert(0, "clonotype_def", key)
                tables[name].append(df)
    for name, frames in tables.items():
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        df.to_csv(out / f"{name}.tsv", sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
