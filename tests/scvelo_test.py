# Panels refer to https://docs.google.com/document/d/1byHTSTjPRaXJ9KKFndWORH3fzccE-TcTWjBesUZjsZA/edit

# Update Feb 2022
# 1. add unspliced prediction
# %%
import os
from time import time

import pickle
import torch
import numpy as np
import scvelo as scv
from scvelo.plotting.utils import get_obs_vector
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from umap import UMAP

import deepvelo.data_loader.data_loaders as module_data
import deepvelo.model.loss as module_loss
import deepvelo.model.metric as module_metric
import deepvelo.model.model as module_arch
from deepvelo.parse_config import ConfigParser
from deepvelo.trainer import Trainer
from deepvelo.data_loader.data_loaders import VeloDataset
from deepvelo.utils.temporal import latent_time
from deepvelo.utils.confidence import velocity_confidence
from deepvelo.utils.scatter import scatter
from deepvelo.utils.util import update_dict
from deepvelo.utils.velocity import velocity
from deepvelo import Constants

# fix random seeds for reproducibility
SEED = 123
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(SEED)

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params(
    "scvelo", transparent=False
)  # for beautified visualization
meta = {
    "algo": "DeepVelo",  # choice of {DeepVelo, dynamical}
    "mode": "run",  # choice of {run, show-target}
    "targets": "s",  # choice of {u, s, both}, works for show-target mode only
    "mask_zero": False,
    "monitor": {
        "verbose": True,  # the overall cotrol of the verbose
        "kinetic_rates": True,  # show kinetic rates during training
        "confidence": True,  # show confidence score distribution during training
    },
    "fig_dir": "online",
}
MASK_ZERO = meta["mask_zero"]
DEEPVELO = True if meta["algo"] == "DeepVelo" else False
DYNAMICAL = (
    True if meta["algo"] == "dynamical" else False
)  # whether use the dynamical mode of scvelo and compute latent time
DEEPVELO_FILE = "scvelo_mat.npz"
data = "GABAInterneuraons"
fig_dir = meta["fig_dir"]
# make fig_dir if not exist
os.makedirs(fig_dir, exist_ok=True)
SURFIX = "[dynamical]" if DYNAMICAL else ""
SURFIX += "[deep_velo]" if DEEPVELO else ""

# specific configs to overide the default configs
if data == "DG":
    configs = {
        "name": "DeepVelo",
        # "arch": {"args": {"pred_unspliced": True}},
        "data_loader": {"args": {"type": "pca, t"}},
        "trainer": {"epochs": 210},
    }
elif data == "GABAInterneuraons":
    configs = {
        "name": "DeepVelo",
        "arch": {"args": {"layers": [128, 128, 128], "pred_unspliced": False}},
        "data_loader": {"args": {"type": "pca, t"}},
        "trainer": {"epochs": 400},
    }
configs = update_dict(Constants.default_configs, configs)


# %%loading and cleaningup data
if data == "DG":
    adata = scv.datasets.dentategyrus()
elif data == "EP":
    adata = scv.datasets.pancreas()
elif data == "GABAInterneuraons":
    adata = scv.read("../data/MDT_GABAInterneuraons_v2.h5ad", cache=True)
else:
    raise ValueError(
        "choose data from \{'EP', 'DG', 'velocyto_dg', 'velocyto_hg',"
        "'E9M2_Glial', 'E9-11F1_Glial', 'E9-11M2_Glial', 'E9-11F1_Gluta'\}"
    )


# %% Preprocessing Data
# here we have the size normalization
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# here comes the NN graph and dynamic estimations
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)
# %%


def train(adata):
    batch_size, n_genes = adata.layers["velocity"].shape
    configs["arch"]["args"]["n_genes"] = n_genes
    configs["data_loader"]["args"]["batch_size"] = batch_size
    print(configs)
    config = ConfigParser(configs)
    logger = config.get_logger("train")

    # setup data_loader instances
    data_loader = config.init_obj("data_loader", module_data)
    valid_data_loader = data_loader.split_validation()

    # build model architecture, then print to console
    if config["arch"]["type"] in ["VeloGCN", "VeloGIN"]:
        model = config.init_obj("arch", module_arch, g=data_loader.dataset.g)
    else:
        model = config.init_obj("arch", module_arch)
    logger.info(model)

    # get function handles of loss and metrics
    criterion = getattr(module_loss, configs["loss"])
    metrics = [getattr(module_metric, met) for met in configs["metrics"]]

    # build optimizer, learning rate scheduler. delete every lines containing lr_scheduler for disabling scheduler
    trainable_params = filter(lambda p: p.requires_grad, model.parameters())
    optimizer = config.init_obj("optimizer", torch.optim, trainable_params)
    # lr_scheduler = config.init_obj('lr_scheduler', torch.optim.lr_scheduler, optimizer)

    trainer = Trainer(
        model,
        criterion,
        metrics,
        optimizer,
        config=config,
        data_loader=data_loader,
        valid_data_loader=valid_data_loader,
    )

    def callback(epoch):
        # evaluate all and return the velocity matrix (cells, features)
        config_copy = configs["data_loader"]["args"].copy()
        config_copy.update(shuffle=False, training=False)
        eval_loader = getattr(module_data, configs["data_loader"]["type"])(
            **config_copy
        )
        velo_mat, velo_mat_u, kinetic_rates = trainer.eval(
            eval_loader, return_kinetic_rates=meta["monitor"]["kinetic_rates"]
        )
        # _plot_kinetic_rates(
        #     kinetic_rates,
        #     select_gene="Tmsb10",
        #     fig_surfix=f"E{epoch}",
        # )
        run_basic_velocity(velo_mat, velo_mat_u, fig_surfix=f"E{epoch}")

    trainer.train_with_epoch_callback(callback=callback, freq=30)

    # evaluate all and return the velocity matrix (cells, features)
    config_copy = configs["data_loader"]["args"].copy()
    config_copy.update(shuffle=False, training=False)
    eval_loader = getattr(module_data, configs["data_loader"]["type"])(**config_copy)
    velo_mat, velo_mat_u, kinetic_rates = trainer.eval(
        eval_loader, return_kinetic_rates=meta["monitor"]["kinetic_rates"]
    )
    _plot_kinetic_rates(
        kinetic_rates,
        select_gene="Tmsb10",
    )

    print("velo_mat shape:", velo_mat.shape)
    # np.savez(f"./data/{configs['online_test']}", velo_mat=velo_mat)

    return velo_mat, velo_mat_u


# %% Compute velocity and velocity graph
def run_basic_velocity(velo_mat=None, velo_mat_u=None, fig_surfix=None):
    """
    :param velo_mat: spliced velocity matrix (cells, features)
    :param velo_mat_u: unspliced velocity matrix (cells, features)
    :param fig_surfix: str, suffix of figure name
    """
    if velo_mat is not None:
        # TODO: careful with this line, delete or whether set mode to deterministic
        velocity(adata, mask_zero=MASK_ZERO)
        adata.layers["velocity"] = velo_mat  # (cells, genes)
        if velo_mat_u is not None and len(velo_mat_u) > 0:
            adata.layers["velocity_unspliced"] = velo_mat_u
    else:
        if DYNAMICAL:
            scv.tl.recover_dynamics(adata)
            velocity(adata, mode="dynamical", mask_zero=MASK_ZERO)
        else:
            velocity(adata, mask_zero=MASK_ZERO)

        # %% output and change the velocity
        if DEEPVELO:
            now = time()
            velo_mat = train(adata)
            print(f"finished in {time()-now:.2f}s")

            if isinstance(velo_mat, np.ndarray):
                assert adata.layers["velocity"].shape == velo_mat.shape
                adata.layers["velocity"] = velo_mat  # (cells, genes)
            elif isinstance(velo_mat, tuple) and len(velo_mat) == 2:
                adata.layers["velocity"] = velo_mat[0]
                adata.layers["velocity_unspliced"] = velo_mat[1]
            else:
                raise ValueError("velocity matrix should be numpy array or tuple")

    scv.tl.velocity_graph(adata, n_jobs=8)

    # %% generate umap if need
    if not ("X_umap" in adata.obsm or "X_tsne" in adata.obsm):
        scv.tl.umap(adata)  # this will add adata.obsm: 'X_umap'

    # %% plot panel a
    _plot(adata, fig_surfix)

    if meta["monitor"]["confidence"]:
        # _plot_confidence(adata, fig_surfix)
        _plot_confidence(
            adata,
            fig_surfix,
            scope_key="clusters" if data != "GABAInterneuraons" else "Celltype",
        )


def _plot(adata, fig_surfix=None):
    fig_surfix = fig_surfix if fig_surfix else SURFIX
    if data == "GABAInterneuraons":
        scv.pl.velocity_embedding_stream(
            adata,
            basis="tsne",
            color="Celltype",
            legend_fontsize=9,
            dpi=300,
            save=f"{fig_dir}/velo_emb_stream{fig_surfix}[Celltype].png",
        )
        scv.tl.velocity_pseudotime(adata)
        scv.tl.paga(adata, groups="Celltype")
        scv.pl.paga(
            adata,
            basis="tsne",
            size=30,
            alpha=0.1,
            min_edge_width=2,
            node_size_scale=1,
            dpi=150,
            save=f"{fig_dir}/trajectory{fig_surfix}.png",
        )
    else:
        scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            color="clusters",
            legend_fontsize=6,
            dpi=300,
            save=f"{fig_dir}/velo_emb_stream{fig_surfix}.png",
        )
        scv.tl.velocity_pseudotime(adata)
        scv.tl.paga(adata, groups="clusters")
        scv.pl.paga(
            adata,
            basis="umap",
            size=30,
            alpha=0.1,
            min_edge_width=2,
            node_size_scale=1,
            dpi=150,
            save=f"{fig_dir}/trajectory{fig_surfix}.png",
        )
        # scatter(
        #     adata,
        #     basis=["Tmsb10"],
        #     add_quiver=True,
        #     dpi=300,
        #     save=f"{fig_dir}/phase{SURFIX}{fig_surfix}.png",
        #     legend_loc_lines="none",
        # )
        # scv.pl.velocity_embedding(
        #     adata,
        #     basis="umap",
        #     arrow_length=6,
        #     arrow_size=1.2,
        #     dpi=150,
        #     save=f"{fig_dir}/velo_emb{fig_surfix}.png",
        # )
    # scv.pl.velocity_embedding_grid(adata, basis='umap', arrow_length=3, arrow_size=2, arrow_color='tab:blue',
    #                                dpi=300, save=f'{fig_dir}/velo_emb_grid{fig_surfix}.png')


def _plot_confidence(adata, fig_surfix=None, scope_key=None):
    vkey = "velocity"
    method = "cosine"
    velocity_confidence(adata, vkey=vkey, method=method, scope_key=scope_key)
    scope = scope_key if scope_key else "Overall"
    print(
        f"{scope} Consistency mean score: {adata.obs[vkey + '_confidence_' + method].mean():.4f}"
        f" std: {adata.obs[vkey + '_confidence_' + method].std():.4f}"
    )
    fig, ax = plt.subplots(figsize=None, dpi=None)
    scv.pl.hist(
        adata.obs[vkey + "_confidence_" + method].values,
        ax=ax,
        normed=True,
        bins=200,
        xlim=[0, 1],
        fontsize=18,
        legend_fontsize=16,
    )
    ax.set_title(
        f"{vkey}_confidence, mean: {adata.obs[vkey + '_confidence_' + method].mean():.4f},"
        f" std: {adata.obs[vkey + '_confidence_' + method].std():.4f}"
    )
    fig.savefig(f"{fig_dir}/{scope}confidence{SURFIX}{fig_surfix}.png")


def _plot_kinetic_rates(kinetic_rates, select_gene=None, fig_surfix=None):
    for k, v in kinetic_rates.items():
        if v is not None:
            adata.layers["cell_specific_" + k] = v

    all_rates = np.concatenate([v for k, v in kinetic_rates.items()], axis=1)

    # pca and umap of all rates
    rates_pca = PCA(n_components=30, random_state=SEED).fit_transform(all_rates)
    adata.obsm["X_rates_pca"] = rates_pca

    rates_umap = UMAP(
        n_neighbors=60,
        min_dist=0.9,
        spread=1.0,
        random_state=SEED,
    ).fit_transform(rates_pca)
    adata.obsm["X_rates_umap"] = rates_umap

    # pca plot or umap plot
    scv.pl.scatter(
        adata,
        basis="rates_pca",
        # omit_velocity_fit=True,
        dpi=300,
        add_outline="Endothelial",
        save=f"{fig_dir}/kinetics_pca{SURFIX}{fig_surfix}.png",
        title="pca of cell-specific kinetic rates",
        legend_loc="right margin",  # whether and how to show the cluster legend
    )
    scv.pl.scatter(
        adata,
        basis="rates_umap",
        # omit_velocity_fit=True,
        dpi=300,
        add_outline="Granule mature, Granule immature, Neuroblast",
        outline_width=(0.15, 0.3),
        save=f"{fig_dir}/kinetics_umap{SURFIX}{fig_surfix}.png",
        title="umap of cell-specific kinetic rates",
        legend_loc="right margin",  # whether and how to show the cluster legend
    )

    if select_gene:  # gene to look at, e.g. "Tmsb10"
        scatter(
            adata,
            x="cell_specific_beta",  # also support direct input ndarray to scatter
            y="cell_specific_gamma",
            basis=select_gene,
            omit_velocity_fit=True,
            dpi=300,
            save=f"{fig_dir}/kinetics_{select_gene}{SURFIX}{fig_surfix}.png",
            legend_loc="right margin",
        )
        # # another way to plot the gene kinetic rate
        # betas = get_obs_vector(adata, select_gene, layer="cell_specific_beta")
        # gammas = get_obs_vector(adata, select_gene, layer="cell_specific_gamma")
        # scatter(
        #     adata,
        #     x=betas,  # also support direct input ndarray to scatter
        #     y=gammas,
        #     omit_velocity_fit=True,
        #     dpi=300,
        #     save=f"{fig_dir}/kinetics_{select_gene}{SURFIX}{fig_surfix}.png",
        #     legend_loc="right margin",  # whether and how to show the cluster legend
        # )


def _plot_data(
    basis="umap",
    color="clusters",
    color_map=None,
    legend_loc="right margin",
    title="",
):
    # umap plot
    # assign color_map to categorical labels is not working properly,
    # see https://github.com/theislab/scvelo/issues/720
    scatter(
        adata,
        basis=basis,
        color=color,
        color_map=color_map,
        dpi=300,
        title=title,
        save=f"{fig_dir}/{basis}_{color}.png",
        legend_loc=legend_loc,  # whether and how to show the cluster legend
    )


def show_target():
    if DYNAMICAL:
        scv.tl.recover_dynamics(adata)
        velocity(adata, mode="dynamical", mask_zero=MASK_ZERO)
    else:
        velocity(adata, mode="deterministic", mask_zero=MASK_ZERO)

    print("computing target velocities")

    # VeloDataset computes the targets
    ds = VeloDataset(adata, basis=configs["data_loader"]["args"]["basis"])
    if meta["targets"] == "s":
        velo_mat = ds.velo.numpy()
        xkey = "Ms"
    elif meta["targets"] == "u":
        velo_mat = ds.velo_u.numpy()
        xkey = "Mu"
    elif meta["targets"] == "both":
        raise NotImplementedError("to be implemented for target both")
    else:
        raise ValueError("meta['targets'] should be 's' or 'u' or 'both'")
    assert adata.layers["velocity"].shape == velo_mat.shape
    adata.layers["velocity"] = velo_mat  # (2930 cells, 1999 genes)

    # compute velocity graph
    scv.tl.velocity_graph(adata, vkey="velocity", xkey=xkey)

    # plot
    _plot(adata, fig_surfix="[target]")


# %% actual run
if meta["mode"] == "run":
    run_basic_velocity()

elif meta["mode"] == "show-target":
    show_target()

# %%
scv.pl.velocity(
    adata,
    var_names=["Dlg2", "Tmsb4x", "Mybl1", "Eif2s3y"],
    basis="tsne",
    dpi=300,
    save=f"{fig_dir}/phase_velo_exp{SURFIX}.png",
)


# %%
latent_time(adata)
scv.pl.scatter(
    adata,
    color="latent_time",
    color_map="gnuplot",
    size=80,
    dpi=300,
    save=f"{fig_dir}/latent_time{SURFIX}.png",
)

top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
# scv.pl.heatmap(adata, var_names=top_genes, tkey='latent_time', n_convolve=100, col_color='clusters')
scv.pl.heatmap(
    adata, var_names=top_genes, tkey="latent_time", n_convolve=100, col_color="clusters"
)
