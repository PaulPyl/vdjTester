"""Compare B marker sets: B call rate in VDJ-B cells vs false calls in T, myeloid, NK."""
import sys
sys.path.insert(0, "/Volumes/vdj/vdjTester/bin")
import numpy as np, pandas as pd, scanpy as sc
from expression_check import ucell, mixture_threshold

R = "/Volumes/vdj/vdjTester"
samples = {f"cll_1_t{i}": (f"{R}/results/cellranger/cll_1_t{i}/per_sample_outs", f"{R}/results/compare/cll_1_t{i}/cll_1_t{i}.vdj_cell_classes.tsv") for i in (1, 2, 3)}
samples |= {f"cll_5_t{i}": (f"{R}/results_cll_5/cellranger/cll_5_t{i}/per_sample_outs", f"{R}/results_cll_5/compare/cll_5_t{i}/cll_5_t{i}.vdj_cell_classes.tsv") for i in (1, 2, 3)}
samples["healthy_1k_B"] = (f"{R}/results/cellranger/per_sample_outs", f"{R}/results/compare/vdj_cell_classes.tsv")
OLD = "CD79A CD79B MS4A1 CD19 CD22 PAX5 BANK1".split()
EXTRA = "FCER2 BLK FCRLA FCRL1 VPREB3 LINC00926 TNFRSF13C POU2AF1 AFF3 OSBPL10 RALGPS2 BLNK".split()
SETS = {"new19": OLD + EXTRA}
res, det = [], []
for s, (d, vf) in samples.items():
    a = sc.read_10x_h5(f"{d}/sample_filtered_feature_bc_matrix.h5", gex_only=True)
    a.var_names_make_unique()
    a.obs_names = a.obs_names.str.replace(r"-\d+$", "", regex=True)
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    a = a[(a.obs.n_genes_by_counts >= 200) & (a.obs.pct_counts_mt <= 20)].copy()
    v = pd.read_csv(vf, sep="\t")
    v = v[v.tool == "cellranger"].set_index("cell_id")["vdj_class"].reindex(a.obs_names)
    X = a.X.tocsc()
    pos = lambda g: np.asarray((X[:, a.var_names.get_loc(g)] > 0).todense()).ravel()
    t_free = ~pos("CD3E") & ~pos("CD3D") & ~pos("CD79A")
    myeloid = v.isna().values & t_free & ((pos("LYZ").astype(int) + pos("CST3") + pos("S100A9")) >= 2)
    nk = v.isna().values & t_free & pos("NKG7") & pos("GNLY") & ~myeloid
    groups = {"vdjB": v.eq("BCR").values, "vdjT": v.eq("TCR").values, "myeloid": myeloid, "NK": nk}
    max_rank = int(round(a.obs.n_genes_by_counts.median()))
    Xr = a.X.tocsr()
    sig = {name: [a.var_names.get_loc(g) for g in genes if g in a.var_names] for name, genes in SETS.items()}
    allscores = ucell(Xr, sig, max_rank)
    for name, cols in sig.items():
        score = allscores[name]
        n_det = np.asarray((Xr[:, cols] > 0).sum(axis=1)).ravel()
        for floor in (0.2, 0.15, 0.1):
            thr, how = mixture_threshold(score, floor)
            called = (score >= thr) & (n_det >= 3)
            row = {"sample": s, "floor": floor, "thr": f"{thr:.2f} {how.split()[0]}"}
            for g, m in groups.items():
                row[f"{g}_n"] = int(m.sum()); row[g] = called[m].mean() if m.sum() else np.nan
                row[f"{g}_count"] = int(called[m].sum())
            res.append(row)
    for g, m in groups.items():
        if m.sum():
            det.append(pd.Series({x: pos(x)[m].mean() for x in OLD + EXTRA if x in a.var_names}, name=(s, g)))
    print(s, {g: int(m.sum()) for g, m in groups.items()}, "maxRank", max_rank, file=sys.stderr)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200); pd.set_option("display.max_columns", 30)
r = pd.DataFrame(res)
print(r.drop(columns=[c for c in r if c.endswith("_n")]).round(4).to_string(index=False)); sys.exit(0)
print("\ngroup sizes"); print(r.drop_duplicates("sample")[["sample"] + [c for c in r if c.endswith("_n")]].to_string(index=False))
D = pd.DataFrame(det); D.index = pd.MultiIndex.from_tuples(D.index)
print("\n== per-gene detection, mean over CLL samples (and healthy 1k)")
cll = D.loc[[i for i in D.index if i[0].startswith("cll")]].groupby(level=1).mean().T
cll["healthy_vdjB"] = D.loc[("healthy_1k_B", "vdjB")]
print(cll[["vdjB", "vdjT", "myeloid", "NK", "healthy_vdjB"]].round(3).to_string())
