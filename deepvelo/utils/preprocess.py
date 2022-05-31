from typing import Optional, Tuple

import numpy as np
import scvelo as scv
from scvelo import logging as logg

from deepvelo.utils.plot import dist_plot


def clip_and_norm_Ms_Mu(
    adata,
    do_clip: bool = True,
    do_norm: bool = True,
    target_mean: float = 0.4,
    replace: bool = False,
    save_fig: Optional[str] = None,
    plot: bool = True,
    print_summary: bool = True,
) -> Tuple[float, float]:
    """
    Normalize using the mean and standard deviation of the gene expression matrix.

    Args:
        adata (Anndata): Anndata object.
        target_mean (float): target mean.
        replace (bool): replace the original data.
        save_fig (str): directory to save figures.
        plot (bool): plot the distribution of the normalized data.
        print_summary (bool): print the summary of the normalized data.

    Returns:
        Tupel[float, float]: scale factor for Ms and Mu.
    """
    non_zero_Ms = adata.layers["Ms"][adata.layers["Ms"] > 0]
    non_zero_Mu = adata.layers["Mu"][adata.layers["Mu"] > 0]
    if print_summary:
        print(
            f"Raw Ms: mean {adata.layers['Ms'].mean():.2f},"
            f" max {adata.layers['Ms'].max():.2f},"
            f" std {adata.layers['Ms'].std():.2f},"
            f" 99.5% quantile {np.percentile(adata.layers['Ms'], 99.5):.2f}"
            f" 99.5% of non-zero: {np.percentile(non_zero_Ms, 99.5):.2f}"
        )
        print(
            f"Raw Mu: mean {adata.layers['Mu'].mean():.2f},"
            f" max {adata.layers['Mu'].max():.2f},"
            f" std {adata.layers['Mu'].std():.2f},"
            f" 99.5% quantile {np.percentile(adata.layers['Mu'], 99.5):.2f}"
            f" 99.5% of non-zero: {np.percentile(non_zero_Mu, 99.5):.2f}"
        )

    if do_clip:
        # clip the max value to 99.5% quantile
        adata.layers["NMs"] = np.clip(
            adata.layers["Ms"], None, np.percentile(non_zero_Ms, 99.5)
        )
        adata.layers["NMu"] = np.clip(
            adata.layers["Mu"], None, np.percentile(non_zero_Mu, 99.5)
        )
    else:
        adata.layers["NMs"] = adata.layers["Ms"]
        adata.layers["NMu"] = adata.layers["Mu"]
    logg.hint(f"added 'NMs' (adata.layers)")
    logg.hint(f"added 'NMu' (adata.layers)")

    if plot:
        dist_plot(
            adata.layers["NMs"].flatten(),
            adata.layers["NMu"].flatten(),
            bins=20,
            labels=["NMs", "NMu"],
            title="Distribution of Ms and Mu",
            save=f"{save_fig}/hist-Ms-Mu.png" if save_fig is not None else None,
        )

    scale_Ms, scale_Mu = 1.0, 1.0
    if do_norm:
        scale_Ms = adata.layers["NMs"].mean() / target_mean
        scale_Mu = adata.layers["NMu"].mean() / target_mean
        adata.layers["NMs"] = adata.layers["NMs"] / scale_Ms
        adata.layers["NMu"] = adata.layers["NMu"] / scale_Mu
        print(f"Normalized Ms and Mu to mean of {target_mean}")
        if plot:
            ax = scv.pl.hist(
                [adata.layers["NMs"].flatten(), adata.layers["NMu"].flatten()],
                labels=["NMs", "NMu"],
                kde=False,
                normed=False,
                bins=20,
                # xlim=[0, 1],
                fontsize=18,
                legend_fontsize=16,
                show=False,
            )
            if save_fig is not None:
                ax.get_figure().savefig(f"{save_fig}/hist-normed-Ms-Mu.png")

    if print_summary:
        print(
            f"New Ms: mean {adata.layers['NMs'].mean():.2f},"
            f" max {adata.layers['NMs'].max():.2f},"
            f" std {adata.layers['NMs'].std():.2f},"
            f" 99.5% quantile {np.percentile(adata.layers['NMs'], 99.5):.2f}"
        )
        print(
            f"New Mu: mean {adata.layers['NMu'].mean():.2f},"
            f" max {adata.layers['NMu'].max():.2f},"
            f" std {adata.layers['NMu'].std():.2f},"
            f" 99.5% quantile {np.percentile(adata.layers['NMu'], 99.5):.2f}"
        )

    if replace:
        adata.layers["Ms"] = adata.layers["NMs"]
        adata.layers["Mu"] = adata.layers["NMu"]
        logg.hint(f"replaced 'Ms' (adata.layers) with 'NMs'")
        logg.hint(f"replaced 'Mu' (adata.layers) with 'NMu'")

    return scale_Ms, scale_Mu
