"""PBMC 3k single-cell RNA-seq: QC, clustering, annotation and a stability check.

Public data (10x Genomics PBMC 3k, via scanpy). Run:
    micromamba run -n binf-7000-2026-sc python src/analysis.py
Writes figures/*.png and results.json (every number quoted in the report comes from there).
"""
import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation
from sklearn.metrics import adjusted_rand_score

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
SEED = 0
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=150, frameon=False)
INK, RUST, GREEN, BROWN = "#2A2521", "#8C4B34", "#53624A", "#5C4736"
res = {}

# ---- load ----
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
res["n_cells_raw"], res["n_genes_raw"] = int(adata.n_obs), int(adata.n_vars)

# ---- QC: adaptive (MAD) thresholds, written before looking at the filtered result ----
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
def is_outlier(x, nmads, both=True):
    m, mad = np.median(x), median_abs_deviation(x)
    lo = x < m - nmads * mad
    hi = x > m + nmads * mad
    return lo | hi if both else hi
log_counts = np.log1p(adata.obs["total_counts"])
log_genes = np.log1p(adata.obs["n_genes_by_counts"])
out = is_outlier(log_counts, 5) | is_outlier(log_genes, 5) | is_outlier(adata.obs["pct_counts_mt"], 3, both=False) | (adata.obs["pct_counts_mt"] > 20)
res["n_removed_qc"] = int(out.sum())
fixed = (adata.obs["n_genes_by_counts"] < 200) | (adata.obs["n_genes_by_counts"] > 2500) | (adata.obs["pct_counts_mt"] > 5)
res["n_removed_fixed_thresholds"] = int(fixed.sum())
res["n_removed_both"] = int((out & fixed).sum())
fig, ax = plt.subplots(1, 3, figsize=(9, 2.8))
for a, col, ttl in zip(ax, ["total_counts", "n_genes_by_counts", "pct_counts_mt"], ["UMI counts", "genes detected", "% mitochondrial"]):
    a.scatter(np.random.default_rng(1).uniform(-.3, .3, adata.n_obs), adata.obs[col], s=2, c=np.where(out, RUST, "#8a8174"), alpha=.6)
    a.set_xticks([]); a.set_title(ttl, fontsize=9); [s.set_visible(False) for s in (a.spines["top"], a.spines["right"])]
plt.tight_layout(); plt.savefig(FIG / "qc.png", dpi=170); plt.close()
adata = adata[~out].copy()
sc.pp.filter_genes(adata, min_cells=3)
res["n_cells"], res["n_genes"] = int(adata.n_obs), int(adata.n_vars)

# ---- normalise, HVGs, PCA ----
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
adata.raw = adata
ad = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(ad, max_value=10)
sc.tl.pca(ad, n_comps=50, random_state=SEED)
vr = ad.uns["pca"]["variance_ratio"]
res["pc_var_first10_pct"] = round(float(vr[:10].sum() * 100), 1)
fig, a = plt.subplots(figsize=(3.6, 2.6)); a.plot(range(1, 31), vr[:30] * 100, "o-", ms=3, c=BROWN)
a.set_xlabel("principal component"); a.set_ylabel("% variance"); [s.set_visible(False) for s in (a.spines["top"], a.spines["right"])]
plt.tight_layout(); plt.savefig(FIG / "elbow.png", dpi=170); plt.close()
NPC = 15; res["n_pcs_used"] = NPC

# ---- graph, clusters at several resolutions ----
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=NPC, random_state=SEED)
sc.tl.umap(ad, random_state=SEED)
for r in (0.3, 0.5, 0.8, 1.2):
    sc.tl.leiden(ad, resolution=r, key_added=f"leiden_{r}", random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
res["clusters_by_resolution"] = {str(r): int(ad.obs[f"leiden_{r}"].nunique()) for r in (0.3, 0.5, 0.8, 1.2)}
CHOSEN = "leiden_0.5"; res["chosen_resolution"] = 0.5

# ---- stability: re-cluster 80% subsamples, compare on the shared cells ----
rng = np.random.default_rng(SEED); aris = []
base = ad.obs[CHOSEN].astype(int).values
for i in range(20):
    idx = np.sort(rng.choice(ad.n_obs, int(.8 * ad.n_obs), replace=False))
    sub = ad[idx].copy()
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=NPC, random_state=i)
    sc.tl.leiden(sub, resolution=0.5, key_added="l", random_state=i, flavor="igraph", n_iterations=2, directed=False)
    aris.append(adjusted_rand_score(base[idx], sub.obs["l"].astype(int).values))
res["subsample_ari_mean"], res["subsample_ari_sd"] = round(float(np.mean(aris)), 3), round(float(np.std(aris)), 3)
res["ari_between_resolutions_0.5_vs_0.8"] = round(float(adjusted_rand_score(ad.obs["leiden_0.5"], ad.obs["leiden_0.8"])), 3)

# ---- markers and annotation from canonical genes (annotation written down BEFORE comparing with any reference) ----
sc.tl.rank_genes_groups(ad, CHOSEN, method="wilcoxon", use_raw=False, pts=True)
markers = {"CD4 T": ["IL7R", "CCR7"], "CD8 T": ["CD8A", "CD8B"], "NK": ["GNLY", "NKG7"], "B": ["MS4A1", "CD79A"],
           "CD14 monocyte": ["CD14", "LYZ"], "FCGR3A monocyte": ["FCGR3A", "MS4A7"], "Dendritic": ["FCER1A", "CST3"], "Megakaryocyte": ["PPBP", "PF4"]}
adr = adata.raw.to_adata(); adr = adr[ad.obs_names]
# Annotation is an explicit, ordered rule table on the fraction of a cluster's cells expressing canonical markers.
# Written down before comparing with any reference. Two earlier attempts failed and are kept in the report:
#  (1) z-scoring cluster means labelled the CD16+ monocyte cluster "Megakaryocyte" (9% PPBP-positive);
#  (2) an argmax over marker-panel scores labelled the CD8 T cluster "NK" (CD8A/CD8B are poorly detected in 3k PBMC data).
Xg = lambda g: np.asarray(adr[:, g].X.todense()).ravel() if hasattr(adr.X, "todense") else np.asarray(adr[:, g].X).ravel()
cl = ad.obs[CHOSEN].values
need = ["MS4A1", "FCER1A", "FCGR3A", "MS4A7", "CD14", "LYZ", "GNLY", "CD3E", "CD8A", "IL7R", "PPBP"]
frac = pd.DataFrame({g: pd.Series(Xg(g) > 0).groupby(cl).mean() for g in need})
rules = [("B", lambda f: f.MS4A1 > 0.5), ("Dendritic", lambda f: f.FCER1A > 0.5), ("Megakaryocyte", lambda f: f.PPBP > 0.5),
         ("FCGR3A monocyte", lambda f: f.FCGR3A > 0.5 and f.MS4A7 > 0.5), ("CD14 monocyte", lambda f: f.CD14 > 0.4 and f.LYZ > 0.9),
         ("NK", lambda f: f.GNLY > 0.8), ("CD8 T", lambda f: f.CD3E > 0.5 and f.CD8A > 0.3), ("CD4 T", lambda f: f.CD3E > 0.5 and f.IL7R > 0.4 and f.CD8A <= 0.3)]
assign = {}
for k, f in frac.iterrows():
    hits = [name for name, rule in rules if rule(f)]
    assign[k] = hits[0] if len(hits) == 1 else ("AMBIGUOUS:" + "/".join(hits) if hits else "UNASSIGNED")
assert not any(v.startswith(("AMBIGUOUS", "UNASSIGNED")) for v in assign.values()), assign
res["cluster_assignment"] = {str(k): v for k, v in assign.items()}
res["marker_fractions"] = {str(k): {g: round(float(v), 2) for g, v in row.items()} for k, row in frac.iterrows()}
res["top_genes_by_cluster"] = {str(g): [str(x) for x in ad.uns["rank_genes_groups"]["names"][g][:6]] for g in ad.obs[CHOSEN].cat.categories}
ad.obs["cell_type"] = ad.obs[CHOSEN].map(assign).astype("category")
res["cells_per_type"] = {k: int(v) for k, v in ad.obs["cell_type"].value_counts().items()}

# ---- compare with the scanpy tutorial's annotated version of the same dataset (a reference, not ground truth) ----
ref = sc.datasets.pbmc3k_processed()
common = ad.obs_names.intersection(ref.obs_names)
res["n_cells_in_reference"] = int(len(common))
res["ari_vs_tutorial_louvain"] = round(float(adjusted_rand_score(ref.obs.loc[common, "louvain"].astype(str), ad.obs.loc[common, CHOSEN].astype(str))), 3)
ct = pd.crosstab(ad.obs.loc[common, "cell_type"].astype(str), ref.obs.loc[common, "louvain"].astype(str))
res["crosstab"] = {r: {c: int(v) for c, v in row.items() if v > 0} for r, row in ct.iterrows()}

# ---- figures ----
pal = {"CD4 T": BROWN, "CD8 T": "#B8923A", "NK": RUST, "B": GREEN, "CD14 monocyte": "#7A8C9E", "FCGR3A monocyte": "#3f5566", "Dendritic": "#9a6b8a", "Megakaryocyte": "#2A2521"}
fig, a = plt.subplots(figsize=(5.2, 4.2))
um = ad.obsm["X_umap"]
for t in ad.obs["cell_type"].cat.categories:
    m = (ad.obs["cell_type"] == t).values; a.scatter(um[m, 0], um[m, 1], s=4, c=pal.get(t, "#999"), label=f"{t} ({m.sum()})", alpha=.85)
a.set_xticks([]); a.set_yticks([]); a.set_xlabel("UMAP 1"); a.set_ylabel("UMAP 2"); [s.set_visible(False) for s in (a.spines["top"], a.spines["right"])]
a.legend(fontsize=6.5, markerscale=2.5, frameon=False, loc="center left", bbox_to_anchor=(1, .5))
plt.tight_layout(); plt.savefig(FIG / "umap.png", dpi=170); plt.close()

genes = [g for gs in markers.values() for g in gs if g in adr.var_names]
order = list(ad.obs["cell_type"].cat.categories)
X = pd.DataFrame(adr[:, genes].X.toarray() if hasattr(adr.X, "toarray") else adr[:, genes].X, index=ad.obs_names, columns=genes)
X["t"] = ad.obs["cell_type"].values
mean = X.groupby("t").mean().loc[order]; frac = (X[genes].gt(0).groupby(X["t"]).mean()).loc[order]
fig, a = plt.subplots(figsize=(6.4, 3.2))
for i, t in enumerate(order):
    for j, g in enumerate(genes):
        a.scatter(j, i, s=frac.loc[t, g] * 90, c=mean.loc[t, g], cmap="Greens", vmin=0, vmax=mean.values.max(), edgecolors=INK, linewidths=.3)
a.set_xticks(range(len(genes))); a.set_xticklabels(genes, rotation=90, fontsize=7); a.set_yticks(range(len(order))); a.set_yticklabels(order, fontsize=7)
a.set_xlim(-.6, len(genes) - .4); a.invert_yaxis(); [s.set_visible(False) for s in a.spines.values()]
plt.tight_layout(); plt.savefig(FIG / "dotplot.png", dpi=170); plt.close()

fig, a = plt.subplots(figsize=(3.6, 2.6)); a.hist(aris, bins=10, color=GREEN, edgecolor=INK); a.axvline(np.mean(aris), c=RUST)
a.set_xlabel("ARI vs full-data clusters"); a.set_ylabel("subsamples"); [s.set_visible(False) for s in (a.spines["top"], a.spines["right"])]
plt.tight_layout(); plt.savefig(FIG / "stability.png", dpi=170); plt.close()
(ROOT / "results.json").write_text(json.dumps(res, indent=2))
print(json.dumps({k: v for k, v in res.items() if k != "crosstab"}, indent=1))
