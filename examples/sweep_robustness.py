# %%
import numpy as np
import scvelo as scv
import torch
from umap import UMAP
from sklearn.decomposition import PCA
from scipy.stats import mannwhitneyu
import wandb

from deepvelo.utils import velocity, velocity_confidence, update_dict
from deepvelo.utils.preprocess import autoset_coeff_s
from deepvelo.utils.plot import statplot, compare_plot
from deepvelo import train, Constants

hyperparameter_defaults = dict(
    seed=123,
    layers=[64, 64],
    topC=30,
    topG=20,
    lr=0.001,
    pearson_scale=18.0,
    pp_hvg=2000,
    pp_neighbors=30,
    pp_pcs=30
    # NOTE: add any hyperparameters you want to sweep here
)
run = wandb.init(config=hyperparameter_defaults, project="scFormer", reinit=True)
wargs = wandb.config


# fix random seeds for reproducibility
SEED = wargs.seed
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params(
    "scvelo", transparent=False
)  # for beautified visualization

# %% [markdown]
# # Load DG data and preprocess

# %%
adata = scv.datasets.dentategyrus()
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=wargs.pp_hvg)
scv.pp.moments(adata, n_neighbors=wargs.pp_neighbors, n_pcs=wargs.pp_pcs)

# %% [markdown]
# # DeepVelo

# %%
# specific configs to overide the default configs, #NOTE: see train.py for complete args
configs = {
    "name": "DeepVelo",  # name of the experiment
    "arch": {
        "args": {
            "layers": wargs.layers,
        },
    },
    "data_loader": {
        "args": {
            "topC": wargs.topC,
            "topG": wargs.topG,
        },
    },
    "optimizer": {
        "args": {
            "lr": wargs.lr,
        },
    },
    "loss": {
        "args": {
            "pearson_scale": wargs.pearson_scale,
            "coeff_s": autoset_coeff_s(adata),
        },
    },
    "trainer": {"verbosity": 0},  # increase verbosity to show training progress
}
configs = update_dict(Constants.default_configs, configs)

# %%
# initial velocity
velocity(adata, mask_zero=False)
trainer = train(adata, configs)


# %%
scv.tl.velocity_graph(adata, n_jobs=8)

# %%
# velocity plot
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="clusters",
    legend_fontsize=9,
    dpi=150,  # increase dpi for higher resolution
    show=False,
)
# NOTE: may log the plot to wandb using wandb.log({"velocity_embedding_stream": wandb.Image(plt)})


# %%
scv.pl.velocity_embedding(
    adata,
    basis="umap",
    arrow_length=6,
    arrow_size=1.2,
    dpi=150,
    show=False,
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
    show=False,
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
else:
    all_rates = np.concatenate(
        [
            adata.layers["cell_specific_beta"],
            adata.layers["cell_specific_gamma"],
        ],
        axis=1,
    )
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
    show=False,
)

# %%
# plot genes
scv.pl.velocity(
    adata,
    var_names=["Tmsb10", "Ppp3ca"],
    basis="umap",
    show=False,
)

# %%
# save adata for next steps
deepvelo_adata = adata.copy()


# %% [markdown]
# # Compare consistency score

# %%
vkey = "velocity"
method = "cosine"
velocity_confidence(deepvelo_adata, vkey=vkey, method=method)
deepvelo_adata.obs["overall_consistency"] = deepvelo_adata.obs[
    f"{vkey}_confidence_{method}"
].copy()

# %%
vkey = "velocity"
method = "cosine"
scope_key = "clusters"
# 3. cosine similarity, compute within Celltype
velocity_confidence(deepvelo_adata, vkey=vkey, method=method, scope_key=scope_key)
deepvelo_adata.obs["celltype_consistency"] = deepvelo_adata.obs[
    f"{vkey}_confidence_{method}"
].copy()


# NOTE: example of logging metrics
wandb.log(
    {
        "celltype_consistency": deepvelo_adata.obs["celltype_consistency"].mean(),
        "overall_consistency": deepvelo_adata.obs["overall_consistency"].mean(),
    }
)
