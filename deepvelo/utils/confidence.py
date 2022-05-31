import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

from scvelo import logging as logg
from scvelo.core import l2_norm, prod_sum
from scvelo.preprocessing.neighbors import get_neighs
from .util import get_indices

# Code modified from function velocity_confidence from
# https://github.com/theislab/scvelo/blob/master/scvelo/tools/velocity_confidence.py
def velocity_confidence(
    data,
    vkey="velocity",
    method="corr",
    scope_key=None,
    copy=False,
    only_velocity_genes=False,
    only_high_spearman=False,
):
    """Computes confidences of velocities.

    .. code:: python

        scv.tl.velocity_confidence(adata)
        scv.pl.scatter(adata, color='velocity_confidence', perc=[2,98])

    .. image:: https://user-images.githubusercontent.com/31883718/69626334-b6df5200-1048-11ea-9171-495845c5bc7a.png
       :width: 600px


    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    vkey: `str` (default: `'velocity'`)
        Name of velocity estimates to be used.
    method: `str` (default: `'corr'`), choice of `'corr'` or `'cosine'`
        Method to use for computing confidence, whether to use correlation or
        cosine similarity.
    scope_key: `str` (default: `None`)
        For each cell, cells with in scope_key are used for computing confidence.
        If `None`, use cell neighbors. Else, pick the cells in scope_key. Valid
        scope_key has to be in adata.obs.
    copy: `bool` (default: `False`)
        Return a copy instead of writing to adata.
    only_velocity_genes: `bool` (default: `False`)
        Only use velocity genes.
    only_high_spearman: `bool` (default: `False`)
        Only use high spearman.

    Returns
    -------
    velocity_length: `.obs`
        Length of the velocity vectors for each individual cell
    velocity_confidence: `.obs`
        Confidence for each cell
    """  # noqa E501

    adata = data.copy() if copy else data
    if vkey not in adata.layers.keys():
        raise ValueError("You need to run `tl.velocity` first.")
    if method not in ["corr", "cosine"]:
        raise ValueError("Method must be either 'corr' or 'cosine'.")
    if scope_key is not None and method == "corr":
        raise ValueError("Cannot use scope_key with method 'corr'.")

    V = np.array(adata.layers[vkey])

    # filter genes if needed
    tmp_filter = np.invert(np.isnan(np.sum(V, axis=0)))
    if only_velocity_genes and (f"{vkey}_genes" in adata.var.keys()):
        tmp_filter &= np.array(adata.var[f"{vkey}_genes"], dtype=bool)
    if only_high_spearman and ("spearmans_score" in adata.var.keys()):
        tmp_filter &= adata.var["spearmans_score"].values > 0.1
    V = V[:, tmp_filter]

    # zero mean, only need for correlation
    if method == "corr":
        V -= V.mean(1)[:, None]
    V_norm = l2_norm(V, axis=1)
    # normalize, only need for cosine similarity
    if method == "cosine":
        V /= V_norm[:, None]
    R = np.zeros(adata.n_obs)

    indices = (
        get_indices(dist=get_neighs(adata, "distances"))[0]
        if not scope_key
        else adata.obs[scope_key]  # the scope_key (e.g. cluster) of each cell
    )
    Vi_neighs_avg_cache = {}
    for i in range(adata.n_obs):
        if not scope_key:
            # use the neighbors of the cell
            Vi_neighs = V[indices[i]]
        else:
            # use the cells in scope_key
            if indices[i] not in Vi_neighs_avg_cache:
                Vi_neighs = V[indices == indices[i]]
                Vi_neighs_avg_cache[indices[i]] = Vi_neighs.mean(0)
        if method == "corr":
            Vi_neighs -= Vi_neighs.mean(1)[:, None]
            R[i] = np.mean(
                np.einsum("ij, j", Vi_neighs, V[i])
                / (l2_norm(Vi_neighs, axis=1) * V_norm[i])[None, :]
            )
        elif method == "cosine":
            # could compute mean first, because V has been normed
            Vi_neighs_avg = (
                Vi_neighs_avg_cache[indices[i]] if scope_key else Vi_neighs.mean(0)
            )
            R[i] = np.inner(V[i], Vi_neighs_avg)

    adata.obs[f"{vkey}_length"] = V_norm.round(2)
    adata.obs[f"{vkey}_confidence_{method}"] = R

    logg.hint(f"added '{vkey}_length' (adata.obs)")
    logg.hint(f"added '{vkey}_confidence_{method}' (adata.obs)")

    # if f"{vkey}_confidence_transition" not in adata.obs.keys():
    #     velocity_confidence_transition(adata, vkey)

    return adata if copy else None


def split_by_cluster():
    """split scores by cluster."""
    pass


# Code modified from function inner_cluster_coh from
# https://github.com/qiaochen/VeloAE/blob/main/veloproj/eval_util.py
def velocity_cosine(
    data,
    vkey="velocity",
    copy=False,
):
    pass


def inner_cluster_coh(adata, k_cluster, k_velocity, return_raw=False):
    """In-cluster Coherence Score.

    Args:
        adata (Anndata): Anndata object.
        k_cluster (str): key to the cluster column in adata.obs DataFrame.
        k_velocity (str): key to the velocity matrix in adata.obsm.
        return_raw (bool): return aggregated or raw scores.

    Returns:
        dict: all_scores indexed by cluster_edges
        or
        dict: mean scores indexed by cluster_edges
        float: averaged score over all cells.

    """
    clusters = np.unique(adata.obs[k_cluster])
    scores = {}
    all_scores = {}
    for cat in clusters:
        sel = adata.obs[k_cluster] == cat
        nbs = adata.uns["neighbors"]["indices"][sel]
        same_cat_nodes = map(lambda nodes: keep_type(adata, nodes, cat, k_cluster), nbs)
        velocities = adata.layers[k_velocity]
        cat_vels = velocities[sel]
        cat_score = [
            cosine_similarity(cat_vels[[ith]], velocities[nodes]).mean()
            for ith, nodes in enumerate(same_cat_nodes)
            if len(nodes) > 0
        ]
        all_scores[cat] = cat_score
        scores[cat] = np.mean(cat_score)

    if return_raw:
        return all_scores

    return scores, np.mean([sc for sc in scores.values()])


def keep_type(adata, nodes, target, k_cluster):
    """Select cells of targeted type

    Args:
        adata (Anndata): Anndata object.
        nodes (list): Indexes for cells
        target (str): Cluster name.
        k_cluster (str): Cluster key in adata.obs dataframe
    Returns:
        list: Selected cells.
    """
    return nodes[adata.obs[k_cluster][nodes].values == target]
