#!/usr/bin/env python3
"""Harmonise VDJ calls from several tools and compare them cell by cell.

Input files are recognised by their name, `<tool>.<sample>.<library>.<ext>`:
  cellranger.*.csv     Cell Ranger filtered_contig_annotations.csv
  immcantation.*.tsv   Change-O MakeDb output (AIRR rearrangement)
  mixcr.*.tsv          MiXCR exportClones --split-by-tags Cell

Every tool is reduced to the same columns, keeping productive chains only:
  tool, cell_id, locus, v_gene, j_gene, c_gene, cdr3_core, umis
For each tool, cell and locus the chain with most UMIs is the call.

Normalisation, so that differences are biological rather than notational:
  cell_id    barcode without the GEM-well suffix (`-1`)
  genes      first call only, allele (`*01`) and scores dropped
  cdr3_core  CDR3 amino acids without the conserved C and W/F, since tools differ
             on whether those are part of the CDR3
"""

import argparse
import itertools
import re
from pathlib import Path

import pandas as pd

LOCI = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]
IG_LOCI = {"IGH", "IGK", "IGL"}
TR_LOCI = {"TRA", "TRB", "TRD", "TRG"}

CODONS = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        itertools.product("TCAG", repeat=3),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}


def translate(nt):
    if not isinstance(nt, str):
        return ""
    nt = nt.upper()
    return "".join(CODONS.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def clean_cell(barcode):
    barcode = re.sub(r"_contig_\d+$", "", str(barcode))
    return re.sub(r"-\d+$", "", barcode)


def clean_gene(call):
    if not isinstance(call, str) or call.strip() in ("", "None", "nan"):
        return ""
    first = re.split(r"[,;]", call)[0].strip()
    first = re.sub(r"\(.*?\)", "", first)  # MiXCR hit scores
    return first.split("*")[0].strip()


def gene_set(call):
    """All tied calls (IgBLAST lists equally good genes), allele-stripped, as 'a,b'."""
    if not isinstance(call, str) or call.strip() in ("", "None", "nan"):
        return ""
    genes = {re.sub(r"\(.*?\)", "", g).split("*")[0].strip() for g in re.split(r"[,;]", call)}
    return ",".join(sorted(g for g in genes if g))


def cdr3_core(aa):
    if not isinstance(aa, str):
        return ""
    aa = aa.strip().upper()
    if len(aa) > 2 and aa[0] == "C" and aa[-1] in "WF":
        aa = aa[1:-1]
    return aa


def locus_from_gene(*genes):
    for gene in genes:
        if isinstance(gene, str) and gene[:3] in LOCI:
            return gene[:3]
    return ""


def as_bool(values):
    return values.astype(str).str.strip().str.upper().isin(["TRUE", "T"])


def find_col(df, *patterns, exclude=()):
    """First column whose name matches any pattern (case-insensitive)."""
    for pattern in patterns:
        for col in df.columns:
            name = col.lower()
            if re.search(pattern, name) and not any(x in name for x in exclude):
                return col
    return None


def read_cellranger(path):
    df = pd.read_csv(path)
    df = df[as_bool(df["productive"])]
    if "is_cell" in df:
        df = df[as_bool(df["is_cell"])]
    return pd.DataFrame({
        "cell_id": df["barcode"].map(clean_cell),
        "locus": df["chain"],
        "v_gene": df["v_gene"].map(clean_gene),
        "j_gene": df["j_gene"].map(clean_gene),
        "c_gene": df["c_gene"].map(clean_gene),
        "cdr3_core": df["cdr3"].map(cdr3_core),
        "umis": df["umis"],
    })


def read_immcantation(path):
    df = pd.read_csv(path, sep="\t")
    df = df[as_bool(df["productive"])]
    if "junction_aa" in df:
        cdr3 = df["junction_aa"]
    else:
        cdr3 = df["junction"].map(translate)
    umis = df["umi_count"] if "umi_count" in df else df.get("consensus_count", 0)
    c_call = df["c_call"] if "c_call" in df else ""
    return pd.DataFrame({
        "cell_id": df["cell_id"].map(clean_cell),
        "locus": df["locus"],
        "v_gene": df["v_call"].map(clean_gene),
        "j_gene": df["j_call"].map(clean_gene),
        "c_gene": pd.Series(c_call, index=df.index).map(clean_gene),
        "v_set": df["v_call"].map(gene_set),
        "j_set": df["j_call"].map(gene_set),
        "cdr3_core": cdr3.map(cdr3_core),
        "umis": umis,
    })


def read_mixcr(path):
    df = pd.read_csv(path, sep="\t")
    # Column names differ between MiXCR versions, so match them loosely.
    cell = find_col(df, r"cell")
    v = find_col(df, r"^(best)?vgene", r"vhit", r"^v", exclude=("family",))
    j = find_col(df, r"^(best)?jgene", r"jhit", r"^j", exclude=("family",))
    c = find_col(df, r"^(best)?cgene", r"chit", r"^c(gene|hit)", exclude=("family",))
    aa = find_col(df, r"aaseq.*cdr3", r"cdr3.*aa", r"aa.*cdr3")
    umis = find_col(df, r"molecule", r"umi")
    reads = find_col(df, r"readcount")
    missing = [n for n, col in [("cell", cell), ("V", v), ("J", j), ("CDR3 aa", aa)] if col is None]
    if missing:
        raise SystemExit(f"{path}: cannot find MiXCR column(s) for {missing}; columns are {list(df.columns)}")
    out = pd.DataFrame({
        "cell_id": df[cell].map(clean_cell),
        "v_gene": df[v].map(clean_gene),
        "j_gene": df[j].map(clean_gene),
        "c_gene": df[c].map(clean_gene) if c else "",
        "cdr3_core": df[aa].map(cdr3_core),
        "umis": df[umis] if umis else (df[reads] if reads else 0),
    })
    out["locus"] = [locus_from_gene(vg, jg) for vg, jg in zip(out["v_gene"], out["j_gene"])]
    # MiXCR reports no productive flag: drop stop codons (*) and frameshifts (_).
    productive = ~out["cdr3_core"].str.contains(r"[*_]", regex=True) & (out["cdr3_core"] != "")
    return out[productive & (out["cell_id"] != "")]


READERS = {"cellranger": read_cellranger, "immcantation": read_immcantation, "mixcr": read_mixcr}


def load(paths):
    frames = []
    for path in paths:
        tool = Path(path).name.split(".")[0]
        if tool not in READERS:
            raise SystemExit(f"Don't know which tool produced {path}")
        if Path(path).stat().st_size == 0:
            continue
        df = READERS[tool](path)
        df.insert(0, "tool", tool)
        frames.append(df)
    calls = pd.concat(frames, ignore_index=True)
    # Tools without tied calls: the set is just the single call. MiXCR's hit list
    # is ranked by score rather than tied, so it also uses its top hit only.
    for f in ["v", "j"]:
        if f"{f}_set" not in calls:
            calls[f"{f}_set"] = ""
        calls[f"{f}_set"] = calls[f"{f}_set"].fillna("").where(calls[f"{f}_set"].fillna("") != "", calls[f"{f}_gene"])
    calls["umis"] = pd.to_numeric(calls["umis"], errors="coerce").fillna(0)
    calls = calls[calls["locus"].isin(LOCI)]
    # Rank chains within tool, cell and locus; rank 1 is the call.
    calls = calls.sort_values(["tool", "cell_id", "locus", "umis"], ascending=[True, True, True, False])
    calls["rank"] = calls.groupby(["tool", "cell_id", "locus"]).cumcount() + 1
    return calls.reset_index(drop=True)


def cell_classes(calls):
    per_cell = calls.groupby(["tool", "cell_id"])["locus"].agg(set).reset_index()
    has_ig = per_cell["locus"].map(lambda loci: bool(loci & IG_LOCI))
    has_tr = per_cell["locus"].map(lambda loci: bool(loci & TR_LOCI))
    per_cell["vdj_class"] = [
        "both" if ig and tr else ("BCR" if ig else "TCR") for ig, tr in zip(has_ig, has_tr)
    ]
    for locus in ["IGH", "IGK", "IGL", "TRA", "TRB"]:
        per_cell[f"has_{locus}"] = per_cell["locus"].map(lambda loci, l=locus: l in loci)
    multi = calls.groupby(["tool", "cell_id", "locus"]).size().unstack(fill_value=0)
    per_cell = per_cell.merge(
        (multi.get("IGH", 0) > 1).rename("multiple_IGH").reset_index(), on=["tool", "cell_id"], how="left"
    )
    return per_cell.drop(columns="locus")


def tool_summary(classes):
    rows = []
    for tool, df in classes.groupby("tool"):
        rows.append({
            "tool": tool,
            "cells_with_vdj": len(df),
            "BCR_only": (df["vdj_class"] == "BCR").sum(),
            "TCR_only": (df["vdj_class"] == "TCR").sum(),
            "BCR_and_TCR": (df["vdj_class"] == "both").sum(),
            "BCR_paired": (df["has_IGH"] & (df["has_IGK"] | df["has_IGL"])).sum(),
            "TCR_paired": (df["has_TRA"] & df["has_TRB"]).sum(),
            "multiple_IGH": df["multiple_IGH"].fillna(False).sum(),
        })
    return pd.DataFrame(rows)


def pairwise(calls):
    top = calls[calls["rank"] == 1].set_index(["cell_id", "locus"])
    tools = sorted(calls["tool"].unique())
    rows = []
    for a, b in itertools.combinations(tools, 2):
        ta = top[top["tool"] == a]
        tb = top[top["tool"] == b]
        for locus in sorted(set(ta.index.get_level_values("locus")) | set(tb.index.get_level_values("locus"))):
            la = ta.xs(locus, level="locus", drop_level=True) if locus in ta.index.get_level_values("locus") else ta.iloc[0:0]
            lb = tb.xs(locus, level="locus", drop_level=True) if locus in tb.index.get_level_values("locus") else tb.iloc[0:0]
            both = la.join(lb, how="inner", lsuffix="_a", rsuffix="_b")
            n = len(both)
            agree = {f: (both[f"{f}_a"] == both[f"{f}_b"]) for f in ["v_gene", "j_gene", "cdr3_core"]}
            all_agree = agree["v_gene"] & agree["j_gene"] & agree["cdr3_core"]
            # Compatible: the two tools' sets of tied calls overlap.
            compat = {
                f: pd.Series(
                    [bool(set(x.split(",")) & set(y.split(","))) for x, y in zip(both[f"{f}_set_a"], both[f"{f}_set_b"])],
                    index=both.index, dtype=bool,
                )
                for f in ["v", "j"]
            }
            rows.append({
                "tool_a": a, "tool_b": b, "locus": locus,
                "cells_a": len(la), "cells_b": len(lb), "cells_both": n,
                "only_a": len(la) - n, "only_b": len(lb) - n,
                "v_agree": agree["v_gene"].mean() if n else float("nan"),
                "v_compatible": compat["v"].mean() if n else float("nan"),
                "j_agree": agree["j_gene"].mean() if n else float("nan"),
                "j_compatible": compat["j"].mean() if n else float("nan"),
                "cdr3_agree": agree["cdr3_core"].mean() if n else float("nan"),
                "all_agree": all_agree.mean() if n else float("nan"),
                "all_compatible": (compat["v"] & compat["j"] & agree["cdr3_core"]).mean() if n else float("nan"),
            })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("calls", nargs="+", help="per-tool call files")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    out = Path(args.outdir)

    calls = load(args.calls)
    classes = cell_classes(calls)
    calls.to_csv(out / "vdj_calls_long.tsv", sep="\t", index=False)
    classes.to_csv(out / "vdj_cell_classes.tsv", sep="\t", index=False)
    tool_summary(classes).to_csv(out / "vdj_tool_summary.tsv", sep="\t", index=False)
    pairwise(calls).to_csv(out / "vdj_pairwise_agreement.tsv", sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
