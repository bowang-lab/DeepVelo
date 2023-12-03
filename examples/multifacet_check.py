# %%
import os

import numpy as np
from anndata import AnnData
import scvelo as scv

datasets = [
    "Dentate Gyrus",
    "Pancreas",
    "Hindbrain",
    "Hippocampus",
    "Chondrocyte Organogenesis",
    "Gastrulation Erythroid",
]

# %%
results = {}
for data in datasets:
    if data == "Dentate Gyrus":
        adata = scv.datasets.dentategyrus()
        groupby = "clusters"
    if data == "Pancreas":
        adata = scv.datasets.pancreas()
        groupby = "clusters"
    if data == "Hindbrain":
        adata = scv.read(
            "deepvelo_data/h5ad_files/Hindbrain_GABA_Glio.h5ad", cache=True
        )
        groupby = "Celltype"
    if data == "Hippocampus":
        adata = scv.datasets.dentategyrus_lamanno()
        groupby = "clusters"
    if data == "Chondrocyte Organogenesis":
        adata = scv.read("../data/cao_organogenesis_chondrocyte.h5ad", cache=True)
        adata = adata[np.random.choice(adata.obs_names, 30000, replace=False)]
        adata.X = (
            adata.layers["spliced"]
            + adata.layers["unspliced"]
            + adata.layers["ambiguous"]
        )
        groupby = "Main_cell_type"
    if data == "Mouse Gastrulation":
        adata = scv.datasets.gastrulation_erythroid()
        groupby = "celltype"

    scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata)

    var_names = adata[:, adata.var["velocity_genes"]].var_names
    # var_names = adata[:, adata.var["velocity_genes"]].var_names[:10]
    scv.tl.recover_dynamics(adata, n_jobs=8)
    scv.tl.differential_kinetic_test(adata, var_names=var_names, groupby=groupby)

    res_df = scv.get_df(
        adata[:, var_names], ["fit_diff_kinetics", "fit_pval_kinetics"], precision=2
    )

    num_total_genes = len(var_names)
    num_multifacet_genes = len(res_df)
    ratio = num_multifacet_genes / num_total_genes

    results[data] = {
        "num_total_genes": num_total_genes,
        "num_multifacet_genes": num_multifacet_genes,
        "ratio": ratio,
    }

# %% save results and visualize in sns barplot
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

df = pd.DataFrame(results).T
df = df.reset_index()
df.to_csv("saved/multifacet_check.csv", index=False)

# set default font size
sns.set(font_scale=1.3)

fig = plt.figure(figsize=(12, 9))
# plot with larger font size and make sure x-axis labels are well separated
sns.barplot(x="index", y="ratio", data=df, palette="Blues_d")
plt.xticks(rotation=45, ha="right")
plt.xlabel("")
plt.ylabel("")
plt.title(f"Ratio of multifaceted genes per dataset, avg {df['ratio'].mean():.2f}")
plt.tight_layout()
plt.savefig("saved/multifacet_check.png", dpi=300)

# %%
