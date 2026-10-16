#!/usr/bin/env python3
"""Harmonise VDJ calls from several tools and compare them cell by cell.

Input files are recognised by their name, `<tool>.<sample>.<library>.<ext>`:
  cellranger.*.csv     Cell Ranger filtered_contig_annotations.csv
  immcantation.*.tsv   Change-O MakeDb output (AIRR rearrangement)
  mixcr.*.tsv          MiXCR exportClones --split-by-tags Cell
  trust4.*.tsv         TRUST4 <prefix>_barcode_airr.tsv (AIRR, one row per chain)

TRUST4 has no cell calling, so with --gex-h5 its barcodes are restricted to the
cells in Cell Ranger's filtered GEX matrix. It also reports ambient chains: on
cll_1 (timepoint 3) TCR chains in Cell Ranger's B cells had a median of 1 UMI
(99 % <= 2) against 25 for real BCR chains. A TRUST4 chain is therefore kept only
with >= --trust4-min-umis UMIs and >= --trust4-min-frac of the strongest chain in
its barcode (across BCR and TCR libraries). That took agreement with Cell
Ranger's cell classes from 38 % to 99 %. TRUST4 run on a GEX library (no VDJ
enrichment, --gex-libraries) is reported as tool `trust4_gex`, with its own UMI
minimum (--trust4-gex-min-umis).

Every tool is reduced to the same columns, keeping productive chains only:
  tool, cell_id, locus, v_gene, j_gene, c_gene, cdr3_core, umis
For each tool, cell and locus the chain with most UMIs is the call.

Normalisation, so that differences are biological rather than notational:
  cell_id    barcode without the GEM-well suffix (`-1`)
  genes      first call only, allele (`*01`) and scores dropped
  cdr3_core  CDR3 amino acids without the conserved C and W/F, since tools differ
             on whether those are part of the CDR3

V-gene agreement is reported at four levels (vdj_pairwise_agreement.tsv):
  v_agree       same first call
  v_compatible  overlapping tied calls (IgBLAST lists equally good genes)
  v_paralog     the same after collapsing duplicated genes: IGKV1D-39 = IGKV1-39,
                IGHV1-69D = IGHV1-69, TRAV14/DV4 = TRAV14
  v_germline    paralog, or germline V sequences >= --germline-identity identical
                (IMGT V-REGIONs from --germline-fasta, else the Cell Ranger VDJ
                reference, --vdj-ref)
vdj_v_disagreements.tsv lists every remaining V disagreement by gene pair, with
the germline identity and which level, if any, resolves it.
"""

import argparse
import difflib
import itertools
import re
from functools import lru_cache
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


def read_trust4(path):
    """TRUST4 barcode AIRR. Same layout as Change-O's AIRR output, but `locus` and
    count columns are not guaranteed, so fall back where needed."""
    df = pd.read_csv(path, sep="\t")
    if "productive" in df:
        df = df[as_bool(df["productive"])]
    cdr3 = df["junction_aa"] if "junction_aa" in df else df["junction"].map(translate)
    umis = next((df[c] for c in ["umi_count", "consensus_count", "duplicate_count"] if c in df), 0)
    out = pd.DataFrame({
        "cell_id": df["cell_id"].map(clean_cell),
        "v_gene": df["v_call"].map(clean_gene),
        "j_gene": df["j_call"].map(clean_gene),
        "c_gene": df["c_call"].map(clean_gene) if "c_call" in df else "",
        "v_set": df["v_call"].map(gene_set),
        "j_set": df["j_call"].map(gene_set),
        "cdr3_core": cdr3.map(cdr3_core),
        "umis": umis,
    })
    if "locus" in df:
        out["locus"] = df["locus"].fillna("")
    else:
        out["locus"] = [locus_from_gene(v, j) for v, j in zip(out["v_gene"], out["j_gene"])]
    out["locus"] = [l if l in LOCI else locus_from_gene(v, j) for l, v, j in zip(out["locus"], out["v_gene"], out["j_gene"])]
    return out[(out["cdr3_core"] != "") & ~out["cdr3_core"].str.contains(r"[*_]", regex=True)]


READERS = {"cellranger": read_cellranger, "immcantation": read_immcantation, "mixcr": read_mixcr,
           "trust4": read_trust4}


def gex_barcodes(h5_path):
    import h5py
    with h5py.File(h5_path, "r") as f:
        return {clean_cell(b.decode()) for b in f["matrix"]["barcodes"][:]}


# --- Paralog-aware V-gene comparison ---------------------------------------------

def collapse_duplicate(gene):
    """Map duplicated genes onto their partner: IGKV1D-39 -> IGKV1-39,
    IGHV1-69D -> IGHV1-69, TRAV14/DV4 or TRAV14DV4 -> TRAV14."""
    g = re.sub(r"^(IG[HKL]V\d+)D-", r"\1-", gene)
    g = re.sub(r"(-\d+)D$", r"\1", g)
    g = re.sub(r"/?DV\d+$", "", g)
    return g


def load_imgt_v(fasta):
    """Gene -> V-REGION sequence (lowest allele) from an IMGT-style FASTA with
    headers like >IGHV3-30-3*01 (TRUST4's human_IMGT+C.fa)."""
    seqs, name, buf = {}, None, []
    with open(fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name and name not in seqs:
                    seqs[name] = "".join(buf)
                gene = line[1:].split("*")[0]
                name = gene if re.match(r"^(IG[HKL]|TR[ABDG])V", gene) else None
                buf = []
            elif name:
                buf.append(line.upper())
    if name and name not in seqs:
        seqs[name] = "".join(buf)
    return seqs


def load_v_germlines(vdj_ref, imgt_fasta=None):
    """Germline V sequences per gene: IMGT first (complete gene set), then the Cell
    Ranger VDJ reference (fasta/regions.fa, L-REGION+V-REGION) for anything else.
    Ensembl-based 10x references lack some IMGT genes, e.g. IGHV3-30-3."""
    seqs = load_imgt_v(imgt_fasta) if imgt_fasta else {}
    if vdj_ref:
        for gene, seq in load_10x_v(vdj_ref).items():
            seqs.setdefault(gene, seq)
    return seqs


def load_10x_v(vdj_ref):
    seqs, name, buf = {}, None, []
    with open(Path(vdj_ref) / "fasta" / "regions.fa") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name and name not in seqs:
                    seqs[name] = "".join(buf)
                fields = line[1:].split("|")
                name = fields[2] if len(fields) > 3 and fields[3] == "L-REGION+V-REGION" else None
                buf = []
            elif name:
                buf.append(line)
    if name and name not in seqs:
        seqs[name] = "".join(buf)
    return seqs


class VComparer:
    def __init__(self, germlines, min_identity):
        self.germlines = germlines
        self.min_identity = min_identity

    @lru_cache(maxsize=None)
    def identity(self, a, b):
        sa = self.germlines.get(a) or self.germlines.get(collapse_duplicate(a))
        sb = self.germlines.get(b) or self.germlines.get(collapse_duplicate(b))
        if not sa or not sb:
            return float("nan")
        return difflib.SequenceMatcher(None, sa, sb, autojunk=False).ratio()

    def level(self, a, b):
        if a == b:
            return "same"
        if collapse_duplicate(a) == collapse_duplicate(b):
            return "duplicate"
        ident = self.identity(a, b)
        if ident == ident and ident >= self.min_identity:
            return "similar_germline"
        return "different"


def filter_trust4(df, min_umis, min_frac):
    """Drop ambient chains: too few UMIs, or weak relative to the barcode's best."""
    umis = pd.to_numeric(df["umis"], errors="coerce").fillna(0)
    best = umis.groupby(df["cell_id"]).transform("max")
    return df[(umis >= min_umis) & (umis >= min_frac * best)]


def load(paths, cells=None, gex_libraries=(), trust4_min_umis=3, trust4_min_frac=0.1, trust4_gex_min_umis=1):
    frames = []
    for path in paths:
        parts = Path(path).name.split(".")
        tool = parts[0]
        if tool not in READERS:
            raise SystemExit(f"Don't know which tool produced {path}")
        if Path(path).stat().st_size == 0:
            continue
        df = READERS[tool](path)
        if tool == "trust4" and len(parts) > 2 and parts[2] in gex_libraries:
            tool = "trust4_gex"
        if tool.startswith("trust4") and cells is not None:
            df = df[df["cell_id"].isin(cells)]
        df.insert(0, "tool", tool)
        frames.append(df)
    calls = pd.concat(frames, ignore_index=True)
    # TRUST4 filters work per barcode across its BCR and TCR files, so apply after concat.
    for tool, min_umis in [("trust4", trust4_min_umis), ("trust4_gex", trust4_gex_min_umis)]:
        sel = calls["tool"] == tool
        if sel.any():
            kept = filter_trust4(calls[sel], min_umis, trust4_min_frac)
            calls = pd.concat([calls[~sel], kept], ignore_index=True)
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


def pairwise(calls, vcmp=None):
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
            if vcmp is not None and n:
                levels = pd.Series([vcmp.level(x, y) for x, y in zip(both["v_gene_a"], both["v_gene_b"])], index=both.index)
                paralog = levels.isin(["same", "duplicate"]) | compat["v"]
                germline = paralog | (levels == "similar_germline")
            else:
                paralog = germline = pd.Series(dtype=bool)
            rows.append({
                "tool_a": a, "tool_b": b, "locus": locus,
                "cells_a": len(la), "cells_b": len(lb), "cells_both": n,
                "only_a": len(la) - n, "only_b": len(lb) - n,
                "v_agree": agree["v_gene"].mean() if n else float("nan"),
                "v_compatible": compat["v"].mean() if n else float("nan"),
                "v_paralog": paralog.mean() if n and vcmp is not None else float("nan"),
                "v_germline": germline.mean() if n and vcmp is not None else float("nan"),
                "j_agree": agree["j_gene"].mean() if n else float("nan"),
                "j_compatible": compat["j"].mean() if n else float("nan"),
                "cdr3_agree": agree["cdr3_core"].mean() if n else float("nan"),
                "all_agree": all_agree.mean() if n else float("nan"),
                "all_compatible": (compat["v"] & compat["j"] & agree["cdr3_core"]).mean() if n else float("nan"),
                "all_germline": (germline & compat["j"] & agree["cdr3_core"]).mean() if n and vcmp is not None else float("nan"),
            })
    return pd.DataFrame(rows)


def v_disagreements(calls, vcmp):
    """Every V-gene disagreement between two tools on the same cell and locus,
    counted by gene pair, with germline identity and the level resolving it."""
    top = calls[calls["rank"] == 1]
    wide = top.pivot_table(index=["cell_id", "locus"], columns="tool", values="v_gene", aggfunc="first")
    rows = []
    for a, b in itertools.combinations(sorted(wide.columns), 2):
        pair = wide[[a, b]].dropna()
        pair = pair[(pair[a] != pair[b]) & (pair[a] != "") & (pair[b] != "")]
        for (locus, ga, gb), n in pair.reset_index().groupby(["locus", a, b]).size().items():
            rows.append({"tool_a": a, "tool_b": b, "locus": locus, "v_gene_a": ga, "v_gene_b": gb, "cells": n,
                         "germline_identity": vcmp.identity(ga, gb), "level": vcmp.level(ga, gb)})
    cols = ["tool_a", "tool_b", "locus", "v_gene_a", "v_gene_b", "cells", "germline_identity", "level"]
    return pd.DataFrame(rows, columns=cols).sort_values("cells", ascending=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("calls", nargs="+", help="per-tool call files")
    parser.add_argument("--gex-h5", help="Cell Ranger filtered GEX matrix; restricts TRUST4 to these cells")
    parser.add_argument("--vdj-ref", help="Cell Ranger VDJ reference dir, for germline V identity")
    parser.add_argument("--germline-fasta", help="IMGT-style V germlines, e.g. TRUST4's human_IMGT+C.fa")
    parser.add_argument("--germline-identity", type=float, default=0.97)
    parser.add_argument("--gex-libraries", default="", help="comma-separated fastq_ids of GEX libraries")
    parser.add_argument("--trust4-min-umis", type=int, default=3)
    parser.add_argument("--trust4-min-frac", type=float, default=0.1)
    parser.add_argument("--trust4-gex-min-umis", type=int, default=1)
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    out = Path(args.outdir)

    cells = gex_barcodes(args.gex_h5) if args.gex_h5 else None
    have_germlines = args.vdj_ref or args.germline_fasta
    vcmp = VComparer(load_v_germlines(args.vdj_ref, args.germline_fasta), args.germline_identity) if have_germlines else None
    calls = load(args.calls, cells, set(filter(None, args.gex_libraries.split(","))),
                 args.trust4_min_umis, args.trust4_min_frac, args.trust4_gex_min_umis)
    classes = cell_classes(calls)
    calls.to_csv(out / "vdj_calls_long.tsv", sep="\t", index=False)
    classes.to_csv(out / "vdj_cell_classes.tsv", sep="\t", index=False)
    tool_summary(classes).to_csv(out / "vdj_tool_summary.tsv", sep="\t", index=False)
    pairwise(calls, vcmp).to_csv(out / "vdj_pairwise_agreement.tsv", sep="\t", index=False, float_format="%.4f")
    if vcmp is not None:
        v_disagreements(calls, vcmp).to_csv(out / "vdj_v_disagreements.tsv", sep="\t", index=False, float_format="%.4f")
    else:
        (out / "vdj_v_disagreements.tsv").write_text("")


if __name__ == "__main__":
    main()
