# %%
from pathlib import Path
import shutil
import sys
import time
import warnings

import numpy as np
import pandas as pd
import scvelo as scv
import torch
from umap import UMAP
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import wandb

sys.path.append("..")

import deepvelo as dv
from deepvelo.utils import continuity_confidence, update_dict, velocity
from deepvelo.utils.scatter import scatter
from deepvelo.utils.preprocess import autoset_coeff_s
from deepvelo import train, Constants

warnings.resetwarnings()

hyperparameter_defaults = dict(
    do_sweep=False,
    seed=123,
    pearson_scale=18.0,
    pp_hvg=3000,
    pp_neighbors=30,
    velocity_genes=False,
    use_scaled_u=False,
    grad_clip=True,
    stop_pearson_after=10,
    new_gene_selector=False,
    remove_extreme_expr_genes=True,
)
run = wandb.init(
    entity="deepvelo-team",
    config=hyperparameter_defaults,
    project="deepvelo-sweep-gastrulation",
    reinit=True,
)
wargs = wandb.config

# %%
# fix random seeds for reproducibility
SEED = wargs.seed
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)

save_dir = (
    Path(f"saved/gastrulation/sweep")
    if wargs.do_sweep
    else Path(f"saved/gastrulation/{time.strftime('%b-%d-%H%M')}")
)
save_dir.mkdir(parents=True, exist_ok=True)
shutil.copy(__file__, save_dir)

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params(
    "scvelo", transparent=False
)  # for beautified visualization
scv.settings.figdir = str(save_dir)
scv.settings.plot_prefix = ""
use_methods = ["DeepVelo", "Dynamical", "Steady-state"]

remove_murk_genes = False
additional_run = False

# %% [markdown]
# # Load gastrulation data and preprocess

# %%
adata = scv.datasets.gastrulation_erythroid()
adata.obs["clusters"] = adata.obs["celltype"]
print(f"expression value stats, mean {adata.X.mean()}, max {adata.X.max()}")
print(
    f"spliced value stats, mean {adata.layers['spliced'].mean()}, max {adata.layers['spliced'].max()}"
)
# murk_gene_file = "data/Gastrulation/MURK_genes.csv"
# murk_genes = pd.read_csv(murk_gene_file, index_col=0).index.tolist()
# if remove_murk_genes:
#     ori_n_genes = adata.n_vars
#     adata = adata[:, ~adata.var_names.isin(murk_genes)]
#     print(f"remove {ori_n_genes - adata.n_vars} murk genes")


# %%
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=wargs.pp_hvg)
if wargs.remove_extreme_expr_genes:
    # remove genes have extreme high expression
    # 98% quantile of expression in adata.X.mean
    mean_expr = adata.X.mean(axis=0)
    top95 = np.quantile(mean_expr, 0.95)
    adata = adata[:, mean_expr < top95]

# %%
scv.pp.moments(adata, n_neighbors=wargs.pp_neighbors, n_pcs=30)

print(f"processed expression, mean {adata.X.mean()}, max {adata.X.max()}")
print(
    f"processed spliced stats, mean {adata.layers['spliced'].mean()}, max {adata.layers['spliced'].max()}"
)
print(f"Ms stats, mean {adata.layers['Ms'].mean()}, max {adata.layers['Ms'].max()}")
print(f"Ms stats, mean {adata.layers['Mu'].mean()}, max {adata.layers['Mu'].max()}")

adata_raw = adata.copy()
result_adatas = {}

# %% [markdown]
# # DeepVelo
method = "DeepVelo"
adata = adata_raw.copy()

# %%
# specific configs to overide the default configs
configs = {
    "name": "DeepVelo",  # name of the experiment
    "data_loader": {
        "args": {
            "velocity_genes": wargs.velocity_genes,
            "use_scaled_u": wargs.use_scaled_u,
        },
    },
    "loss": {
        "args": {
            "coeff_s": autoset_coeff_s(adata),
            "pearson_scale": wargs.pearson_scale,
            "stop_pearson_after": wargs.stop_pearson_after,
        },
    },
    "trainer": {
        "verbosity": 1,
        "grad_clip": wargs.grad_clip,
    },  # increase verbosity to show training progress
}
configs = update_dict(Constants.default_configs, configs)


# %%
# initial velocity
if wargs.new_gene_selector:
    dv.tl.get_velo_genes(adata)
else:
    velocity(adata)
trainer = train(adata, configs)

# %%
# adata.var["velocity_genes"] = True  # use all genes for plotting
scv.tl.velocity_graph(adata, n_jobs=8)

# %%
# velocity plot
if wargs.do_sweep:
    ax = scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="clusters",
        legend_fontsize=9,
        dpi=150,  # increase dpi for higher resolution
        show=False,
    )
    wandb.log({"velocity_embedding_stream": wandb.Image(ax.figure)})
else:
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="clusters",
        legend_fontsize=9,
        dpi=150,  # increase dpi for higher resolution
        save=f"{method}_velocity.png",
    )

cluster_edges = [
    ("Blood progenitors 1", "Blood progenitors 2"),
    ("Blood progenitors 2", "Erythroid1"),
    ("Erythroid1", "Erythroid2"),
    ("Erythroid2", "Erythroid3"),
]

from deepvelo.utils import cross_boundary_correctness

cbcs, avg_cbc = cross_boundary_correctness(
    adata,
    "clusters",
    "velocity",
    cluster_edges,
    x_emb_key="umap",  # or Ms
)
print(f"Average cross-boundary correctness of DeepVelo: {avg_cbc:.2f}")

wandb.log({"direction_score": avg_cbc})

if wargs.do_sweep:
    sys.exit(0)
# %%
continuity_confidence(adata, trainer)

# %%
# show histogram of cell continuity
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="cell_continuity",
    # color_map="heat",
    legend_fontsize=6,
    # perc=[2, 98],
    dpi=150,  # increase dpi for higher resolution
    save=f"{method}_continuity.png",
)

# %%
from deepvelo.utils import latent_time

latent_time(adata)
scv.pl.scatter(
    adata,
    color="velocity_pseudotime",
    color_map="gnuplot",
    save=f"{method}_pseudotime.png",
)

# %%
if "velocity_unspliced" not in adata.layers:
    adata.layers["velocity_unspliced"] = np.zeros_like(adata.layers["velocity"])
scatter(
    adata,
    basis=[
        "Blvrb",
        "Alad",
        "Clta",
        "Cpox",
        "Gmpr",
        "Gnai2",
        "Tpm4",
        "Ccnd2",
        "Rbms1",
        "Fech",
    ],
    add_quiver=True,
    # dpi=150,
    legend_loc_lines="none",
    ncols=2,
)

# %% [markdown]
# ### plot MURK genes

# %%
if not remove_murk_genes:
    scatter(
        adata,
        basis=["Gclm", "Abcg2", "Hemgn", "Hebp1", "Smim1", "Hba-x"],
        add_quiver=True,
        # dpi=150,
        legend_loc_lines="none",
        ncols=2,
    )

# %% show genes with direction evaluation
cluster_edges = [
    ("Blood progenitors 1", "Blood progenitors 2"),
    ("Blood progenitors 2", "Erythroid1"),
    ("Erythroid1", "Erythroid2"),
    ("Erythroid2", "Erythroid3"),
]

from deepvelo.utils import genewise_cross_boundary_correctness

genewise_cross_boundary_correctness(
    adata,
    cluster_key="clusters",
    velocity_key="velocity",
    cluster_edges=cluster_edges,
)

# visualize adta.var genewise direction scores
gene_ds_sorted = adata.var[adata.var["velocity_genes"]][
    "gene_direction_scores"
].sort_values(ascending=False)
assert not gene_ds_sorted.isna().any()

# histogram of gene direction scores
scv.pl.hist(
    gene_ds_sorted,
    labels=f"gene direction score, avg {gene_ds_sorted.mean():.2f}",
    kde=True,
    normed=True,
    bins=100,
    fontsize=18,
    legend_fontsize=12,
)

# Top 10 velocity genes with highest direction scores
scatter(
    adata,
    basis=gene_ds_sorted.index[:10].tolist(),
    add_quiver=True,
    legend_loc_lines="none",
    ncols=2,
)

# Last 10 velocity genes with lowest direction scores
scatter(
    adata,
    basis=gene_ds_sorted.index[-10:].tolist(),
    add_quiver=True,
    legend_loc_lines="none",
    ncols=2,
)

# %% plot MURK genes in the order of their direction scores
if not remove_murk_genes:
    scatter(
        adata,
        basis="Smim1",
        add_quiver=True,
        quiver_size=0.3,
        legend_loc_lines="none",
        ncols=2,
        save=f"{method}_Smim1.png",
    )
    scatter(
        adata,
        basis="Gypa",
        add_quiver=True,
        quiver_size=0.1,
        legend_loc_lines="none",
        ncols=2,
        save=f"{method}_Gypa.png",
    )

# %%
scv.pl.velocity_embedding(
    adata,
    basis="umap",
    arrow_length=6,
    arrow_size=1.2,
    dpi=150,
)

# %%
scv.pl.velocity_embedding_grid(
    adata,
    basis="umap",
    arrow_length=4,
    # alpha=0.1,
    arrow_size=2,
    arrow_color="tab:blue",
    dpi=150,
)

result_adatas["DeepVelo"] = adata.copy()

# %% [markdown]
# # scVelo (dynamical) or steady-state
for method in [m for m in use_methods if m != "DeepVelo"]:
    adata = adata_raw.copy()
    if method == "Dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=8)
        scv.tl.velocity(adata, mode="dynamical")
    elif method == "Steady-state":
        scv.tl.velocity(adata, mode="stochastic")  # deterministic
    scv.tl.velocity_graph(adata, n_jobs=8)

    # velocity plot
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="clusters",
        legend_fontsize=9,
        legend_loc="none",
        dpi=150,  # increase dpi for higher resolution
        save=f"{method}_velocity_embedding_stream.png",
    )
    scv.pl.velocity_embedding(
        adata,
        basis="umap",
        arrow_length=6,
        arrow_size=1.2,
        dpi=150,
        save=f"{method}_velocity_embedding.png",
    )
    latent_time(adata)
    scv.pl.scatter(
        adata,
        color="velocity_pseudotime",
        color_map="gnuplot",
        save=f"{method}_pseudotime.png",
    )
    if not remove_murk_genes:
        if "velocity_unspliced" not in adata.layers:
            adata.layers["velocity_unspliced"] = np.zeros_like(adata.layers["velocity"])
        scatter(
            adata,
            basis="Smim1",
            add_quiver=True,
            # quiver_size=0.3,
            legend_loc_lines="none",
            ncols=2,
            save=f"{method}_Smim1.png",
        )
        scatter(
            adata,
            basis="Gypa",
            add_quiver=True,
            # quiver_size=0.1,
            legend_loc_lines="none",
            ncols=2,
            save=f"{method}_Gypa.png",
        )
    # save adata for next steps
    result_adatas[method] = adata.copy()

# %% [markdown]
# # Evaluation
dv.pipe.evaluate(
    result_adatas,
    cluster_edges=cluster_edges,
    cluster_key="clusters",
    vkey="velocity",
    save_dir=save_dir,
)

# %% skip the following cells by default
if not additional_run:
    sys.exit(0)

# %%
adata = result_adatas["DeepVelo"]
continuity_confidence(adata, trainer)

# %%
# show histogram of cell continuity
cell_continuity_mse = adata.layers["continuity_relative_error"].mean(axis=1) / (
    adata.layers["Ms"].mean(axis=1) + 1e-6
)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="cell_continuity",
    # color_map="heat",
    legend_fontsize=6,
    # perc=[2, 98],
    dpi=150,  # increase dpi for higher resolution
)

# %%
# get kinetic_rates
if "cell_specific_alpha" in adata.layers:
    all_rates = np.concatenate(
        [
            adata.layers["cell_specific_beta"],
            adata.layers["cell_specific_gamma"],
            adata.layers["cell_specific_alpha"],
        ],
        axis=1,
    )
    gene_wise_rates = np.concatenate(
        [
            adata.layers["cell_specific_beta"],
            adata.layers["cell_specific_gamma"],
            adata.layers["cell_specific_alpha"],
        ],
        axis=0,
    ).T
else:
    all_rates = np.concatenate(
        [
            adata.layers["cell_specific_beta"],
            adata.layers["cell_specific_gamma"],
        ],
        axis=1,
    )
    gene_wise_rates = np.concatenate(
        [
            adata.layers["cell_specific_beta"],
            adata.layers["cell_specific_gamma"],
        ],
        axis=0,
    ).T

# pca and umap of all rates
rates_pca = PCA(n_components=30, random_state=SEED).fit_transform(all_rates)
adata.obsm["X_rates_pca"] = rates_pca

rates_umap = UMAP(
    n_neighbors=60,
    min_dist=0.6,
    spread=0.9,
    random_state=SEED,
).fit_transform(rates_pca)
adata.obsm["X_rates_umap"] = rates_umap

# pca and umap of gene-wise rates
rates_pca_gene_wise = PCA(n_components=30, random_state=SEED).fit_transform(
    adata.layers["Ms"].T
)
adata.varm["rates_pca"] = rates_pca_gene_wise

rates_umap_gene_wise = UMAP(
    n_neighbors=60,
    min_dist=0.6,
    spread=0.9,
    random_state=SEED,
).fit_transform(rates_pca_gene_wise)
adata.varm["rates_umap"] = rates_umap_gene_wise


# %%
# plot kinetic rates umap
scv.pl.scatter(
    adata,
    basis="rates_umap",
    # omit_velocity_fit=True,
    add_outline="Granule mature, Granule immature, Neuroblast",
    outline_width=(0.15, 0.3),
    title="umap of cell-specific kinetic rates",
    legend_loc="none",
    dpi=150,
)

# %%
num_show_genes = 10
top_ranked_genes_corr = adata.var.sort_values("gene_corr", ascending=False).index[
    :num_show_genes
]
least_ranked_genes_corr = adata.var.sort_values("gene_corr", ascending=False).index[
    -num_show_genes:
]

ax = scv.pl.velocity(
    adata,
    var_names=top_ranked_genes_corr,
    basis="umap",
    show=False,
)
ax.get_figure().suptitle("top ranked genes by gene_corr")

ax = scv.pl.velocity(
    adata,
    var_names=least_ranked_genes_corr,
    basis="umap",
    show=False,
)
ax.get_figure().suptitle("least ranked genes by gene_corr")


# %%
top_ranked_genes_continuity = adata.var.sort_values(
    "gene_continuity", ascending=False
).index[:num_show_genes]
least_ranked_genes_continuity = adata.var.sort_values(
    "gene_continuity", ascending=False
).index[-num_show_genes:]

ax = scv.pl.velocity(
    adata,
    var_names=top_ranked_genes_continuity,
    basis="umap",
    show=False,
)
ax.get_figure().suptitle("top ranked genes by gene_continuity")
# NOTE: we have Ppp3ca listed here nicely
# And nicely it is explainable in the way, in most of the top ranked ones, the
# clear difference of spliced gene expression horisontally can be observed
# among the large portion of neuralblasts to granule immature and mature cells.

ax = scv.pl.velocity(
    adata,
    var_names=least_ranked_genes_continuity,
    basis="umap",
    show=False,
)
ax.get_figure().suptitle("least ranked genes by gene_continuity")

# %%
adata.var["gene_confidence"] = adata.var["gene_continuity"]  # + adata.var["gene_corr"]
# select the top 300 genes by gene_confidence
top_ranked_genes_confidence = adata.var.sort_values(
    "gene_confidence", ascending=False
).index[:300]

adata.var["velocity_genes"] = False
adata.var.loc[top_ranked_genes_confidence, "velocity_genes"] = True

# %%
# recompute velocity with selected genes only

# %% [markdown]
# # scVelo stochastic

adata = adata_raw.copy()

scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata, n_jobs=8)

# velocity plot
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="clusters",
    legend_fontsize=9,
    dpi=150,  # increase dpi for higher resolution
)

if "velocity_unspliced" not in adata.layers:
    adata.layers["velocity_unspliced"] = np.zeros_like(adata.layers["velocity"])
scatter(
    adata,
    basis=["Blvrb", "Alad", "Clta", "Cpox"],
    add_quiver=True,
    # dpi=150,
    legend_loc_lines="none",
    ncols=2,
)

# %% [markdown]
# # scVelo dynamical

adata = adata_raw.copy()

scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=8)

# velocity plot
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="clusters",
    legend_fontsize=9,
    dpi=150,  # increase dpi for higher resolution
)

if "velocity_unspliced" not in adata.layers:
    adata.layers["velocity_unspliced"] = np.zeros_like(adata.layers["velocity"])
scatter(
    adata,
    basis=["Blvrb", "Alad", "Clta", "Cpox"],
    add_quiver=True,
    # dpi=150,
    legend_loc_lines="none",
    ncols=2,
)

# %%
