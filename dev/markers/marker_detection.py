"""Per-gene detection in VDJ-defined B (BCR only) vs T (TCR only) vs no-VDJ cells."""
import sys
import numpy as np, pandas as pd, scanpy as sc

R = "/Volumes/vdj/vdjTester"
samples = {f"cll_1_t{i}": f"{R}/results" for i in (1, 2, 3)} | {f"cll_5_t{i}": f"{R}/results_cll_5" for i in (1, 2, 3)}
CAND = """CD19 CD79A CD79B MS4A1 CD22 PAX5 BANK1 CD37 CD74 BLK FCRL1 FCRL2 FCRL3 FCRL5 FCRLA FCER2
TCL1A VPREB3 CD24 CD72 CD180 TNFRSF13B TNFRSF13C CR2 SPIB POU2AF1 IRF8 BLNK CD200 ROR1 LEF1 CD5
EBF1 BACH2 AFF3 MEF2C CIITA CD40 SWAP70 OSBPL10 COBLL1 RALGPS2 ADAM28 KLHL14 PLEKHG1 LINC00926
HLA-DQA1 HLA-DRA FCMR CD52 LTB SELL CXCR4 CD83 MZB1 IGLL5 TNFRSF17 CD27 CD38 BTK SYK LYN PTPRC
CD3D CD3E CD3G CD2 CD247 LCK IL7R CD7 LYZ CST3 NKG7 GNLY S100A9""".split()
rows, uni = [], []
for s, root in samples.items():
    a = sc.read_10x_h5(f"{root}/cellranger/{s}/per_sample_outs/sample_filtered_feature_bc_matrix.h5")
    a.var_names_make_unique()
    a.obs_names = [c.split("-")[0] for c in a.obs_names]
    e = pd.read_csv(f"{root}/expression/{s}/{s}.expression_cells.tsv", sep="\t", index_col=0)
    v = pd.read_csv(f"{root}/compare/{s}/{s}.vdj_cell_classes.tsv", sep="\t")
    v = v[v.tool == "cellranger"].set_index("cell_id")["vdj_class"]
    a = a[a.obs_names.isin(e.index[e.qc_pass])].copy()
    grp = pd.Series("noVDJ", index=a.obs_names)
    grp[v.reindex(a.obs_names).eq("BCR").values] = "vdjB"
    grp[v.reindex(a.obs_names).eq("TCR").values] = "vdjT"
    det = pd.DataFrame((a.X > 0).toarray() if hasattr(a.X, "toarray") else (a.X > 0), columns=a.var_names, index=a.obs_names) if False else None
    X = (a.X > 0).tocsc()
    for g in ["vdjB", "vdjT", "noVDJ"]:
        m = (grp == g).values
        frac = np.asarray(X[m].mean(axis=0)).ravel()
        uni.append(pd.DataFrame({"sample": s, "group": g, "n": m.sum(), "gene": a.var_names, "frac": frac}))
    print(s, grp.value_counts().to_dict(), "median genes", int(e.loc[e.qc_pass, "n_genes_by_counts"].median()), file=sys.stderr)
d = pd.concat(uni)
w = d.pivot_table(index="gene", columns=["group", "sample"], values="frac")
out = pd.DataFrame({
    "B_min": w["vdjB"].min(axis=1), "B_mean": w["vdjB"].mean(axis=1),
    "T_max": w["vdjT"].max(axis=1), "noVDJ_mean": w["noVDJ"].mean(axis=1)})
out["B_minus_T"] = out.B_min - out.T_max
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200)
print("== candidates (fraction of cells detecting the gene; B_min = worst CLL timepoint)")
print(out.reindex([g for g in CAND if g in out.index]).round(3).to_string())
print("\n== unbiased: top 40 by B_min - T_max, excluding IG/TR genes")
keep = ~out.index.str.match(r"^(IG[HKL][VDJC]|TR[ABDG][VDJC]|IGLL|MT-|RP[SL])")
print(out[keep].sort_values("B_minus_T", ascending=False).head(40).round(3).to_string())
w.to_csv("/Volumes/vdj/scratch_markers/detection_by_sample.tsv", sep="\t")
