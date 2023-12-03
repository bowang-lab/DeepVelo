# This file process the rds and csv files to generate the anndata object for the
# cao et al organogenesis dataset
# %%
from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np
import pyreadr
import rpy2.robjects as robjects
import anndata2ri
from scipy import sparse

data_folder = Path("/cluster/projects/bwanggroup/cao_et_al/alevin/runs")
cell_annotation_file = data_folder / "cell_annotations_updated_aligned.rds"
gene_annotation_file = data_folder / "gene_col_idx.rds"
matrix_file = data_folder / "quants.rds"
quants_col_names_file = data_folder / "quants_col_names.rds"
quants_row_names_file = data_folder / "quants_row_names.rds"

output_folder = Path("/cluster/home/haotianc/process_cao_organogenesis")

#%% Read in the rds file
cell_annotation = pyreadr.read_r(cell_annotation_file)[None]
gene_annotation = pyreadr.read_r(gene_annotation_file)[None]

gene_names = pyreadr.read_r(quants_col_names_file)[None]
gene_names = gene_names.iloc[:, 0].to_list()
assert len(gene_names) == len(gene_annotation) * 3
gene_names = gene_names[: len(gene_annotation)]
cell_barcodes = pyreadr.read_r(quants_row_names_file)[None]
cell_barcodes = cell_barcodes.iloc[:, 0].to_list()
assert len(cell_barcodes) == len(cell_annotation)
assert all(
    [a == b for a, b in zip(cell_barcodes, cell_annotation["sample_name"])]
), "cell barcodes need to be aligned"

# # align the cell annotation and cell barcodes
# barcode_dict = cell_annotation.set_index("cb")["id"].to_dict()
# assert all([x in barcode_dict for x in cell_barcodes])
# # reorder the cell annotation
# cell_annotation.set_index("cb", inplace=True)

readRDS = robjects.r["readRDS"]
df = readRDS(str(matrix_file))  # dgCMatrix mapped to RS4
# summary = robjects.r["summary"]
# as_dataframe = robjects.r["as.data.frame"]
# as_matrix = robjects.r["as.matrix"]
# matrix_ = as_matrix(df)
anndata2ri.activate()
matrix = anndata2ri.scipy2ri.rpy2py(df)  # anndata2ri.ri2py
matrix = matrix.tocsr()
assert isinstance(matrix, sparse.csr_matrix)

#%% Wrap into an anndata object
# meta info
num_genes = len(gene_annotation)
num_cells = len(cell_annotation)

assert (gene_annotation["spliced"].to_numpy() - np.arange(num_genes) - 1).any() == 0
spliced = matrix[:, :num_genes]
unspliced = matrix[:, num_genes : 2 * num_genes]
ambiguous = matrix[:, 2 * num_genes : 3 * num_genes]

#%% Create an anndata object
adata = sc.AnnData(
    X=spliced,
    obs=cell_annotation.set_index("sample_name"),
    var=pd.DataFrame(index=gene_names),
    layers={"spliced": spliced, "unspliced": unspliced, "ambiguous": ambiguous},
)

# add tsne and other coordinates
adata.obsm["X_tsne"] = cell_annotation[["tsne_1", "tsne_2"]].to_numpy()
adata.obs.drop(columns=["tsne_1", "tsne_2"], inplace=True)
adata.obsm["sub_tsne"] = cell_annotation[["sub_tsne_1", "sub_tsne_2"]].to_numpy()
adata.obs.drop(columns=["sub_tsne_1", "sub_tsne_2"], inplace=True)
adata.obsm["Main_cluster_tsne"] = cell_annotation[
    ["Main_cluster_tsne_1", "Main_cluster_tsne_2"]
].to_numpy()
adata.obs.drop(columns=["Main_cluster_tsne_1", "Main_cluster_tsne_2"], inplace=True)
adata.obsm["Sub_cluster_tsne"] = cell_annotation[
    ["Sub_cluster_tsne_1", "Sub_cluster_tsne_2"]
].to_numpy()
adata.obs.drop(columns=["Sub_cluster_tsne_1", "Sub_cluster_tsne_2"], inplace=True)
adata.obsm["Main_trajectory_umap"] = cell_annotation[
    ["Main_trajectory_umap_1", "Main_trajectory_umap_2", "Main_trajectory_umap_3"]
].to_numpy()
adata.obs.drop(
    columns=[
        "Main_trajectory_umap_1",
        "Main_trajectory_umap_2",
        "Main_trajectory_umap_3",
    ],
    inplace=True,
)
adata.obsm["Main_trajectory_refined_umap"] = cell_annotation[
    [
        "Main_trajectory_refined_umap_1",
        "Main_trajectory_refined_umap_2",
        "Main_trajectory_refined_umap_3",
    ]
].to_numpy()
adata.obs.drop(
    columns=[
        "Main_trajectory_refined_umap_1",
        "Main_trajectory_refined_umap_2",
        "Main_trajectory_refined_umap_3",
    ],
    inplace=True,
)
adata.obsm["Sub_trajectory_umap"] = cell_annotation[
    ["Sub_trajectory_umap_1", "Sub_trajectory_umap_2"]
].to_numpy()
adata.obs.drop(
    columns=[
        "Sub_trajectory_umap_1",
        "Sub_trajectory_umap_2",
    ],
    inplace=True,
)
# remove other unused columns
adata.obs.drop(
    columns=[
        "id.x",
        "sample",
        "id.y",
    ],
    inplace=True,
)

# Save
# convert obs columns to string if contains nan
for col in adata.obs.columns:
    if adata.obs[col].isnull().any():
        adata.obs[col] = adata.obs[col].astype(str)
# save the anndata object
adata.write(output_folder / "cao_organogenesis.h5ad", compression="gzip")

# %% subset the data to only include Chondrocyte trajectory
chondrocyte_idx = adata.obs["Sub_trajectory_name"] == "Chondrocyte trajectory"
adata_chondrocyte = adata[chondrocyte_idx, :]

print(adata_chondrocyte.obs["Sub_trajectory_name"].unique())
print(adata_chondrocyte.obs["Main_cell_type"].value_counts())

# Save
# convert obs columns to string if contains nan
for col in adata_chondrocyte.obs.columns:
    if adata_chondrocyte.obs[col].isnull().any():
        adata_chondrocyte.obs[col] = adata_chondrocyte.obs[col].astype(str)
# save the anndata object
adata_chondrocyte.write(
    output_folder / "cao_organogenesis_chondrocyte.h5ad", compression="gzip"
)

# %% subset the data to only include Mesenchymal trajectory
mesenchymal_idx = (
    adata.obs["Main_trajectory_refined_by_cluster"] == "Mesenchymal trajectory"
)
adata_mesenchymal = adata[mesenchymal_idx, :]

print(adata_mesenchymal.obs["Main_trajectory_refined_by_cluster"].unique())
print(adata_mesenchymal.obs["Main_cell_type"].value_counts())

# Save
# convert obs columns to string if contains nan
for col in adata_mesenchymal.obs.columns:
    if adata_mesenchymal.obs[col].isnull().any():
        adata_mesenchymal.obs[col] = adata_mesenchymal.obs[col].astype(str)
# save the anndata object
adata_mesenchymal.write(
    output_folder / "cao_organogenesis_mesenchymal.h5ad", compression="gzip"
)

# filter in the celltypes as in the figure 4 of cao et al. paper
celltypes = [
    "Chondrocytes & osteoblasts",
    "Connective tissue progenitors",
    "Intermediate Mesoderm",
    "Early mesenchyme",
    "Myocytes",
    "Chondroctye progenitors",
    "Limb mesenchyme",
]
adata_mesenchymal_filtered = adata_mesenchymal[
    adata_mesenchymal.obs["Main_cell_type"].isin(celltypes), :
]
print(adata_mesenchymal_filtered.obs["Main_cell_type"].value_counts())

# Save
adata_mesenchymal_filtered.write(
    output_folder / "cao_organogenesis_mesenchymal_filtered.h5ad",
    compression="gzip",
)

# %% draw tsnes and umaps using the saved coordinates
def plots(adata, dir_name=None):
    if dir_name is not None:
        default_figdir = sc.settings.figdir
        sc.settings.figdir = Path(dir_name)

    sc.pl.scatter(
        adata,
        basis="tsne",
        color=["Main_cell_type", "development_stage"],
        # color="Main_cell_type",
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=1,
        show=False,
        save="_cao_organogenesis_tsne.png",
    )

    # sub_tsne
    adata.obsm["X_sub_tsne"] = adata.obsm["sub_tsne"]
    sc.pl.scatter(
        adata,
        basis="sub_tsne",
        color=["Main_cell_type", "development_stage"],
        legend_fontsize=8,
        legend_fontoutline=1,
        show=False,
        save="_cao_organogenesis.png",
    )

    # Main_cluster_tsne
    adata.obsm["X_Main_cluster_tsne"] = adata.obsm["Main_cluster_tsne"]
    sc.pl.scatter(
        adata,
        basis="Main_cluster_tsne",
        color=["Main_cell_type", "development_stage"],
        legend_fontsize=8,
        legend_fontoutline=1,
        show=False,
        save="_cao_organogenesis.png",
    )

    # Sub_cluster_tsne
    adata.obsm["X_Sub_cluster_tsne"] = adata.obsm["Sub_cluster_tsne"]
    sc.pl.scatter(
        adata,
        basis="Sub_cluster_tsne",
        color=["Main_cell_type", "development_stage"],
        legend_fontsize=8,
        legend_fontoutline=1,
        show=False,
        save="_cao_organogenesis.png",
    )

    # Main_trajectory_umap
    adata.obsm["X_Main_trajectory_umap"] = adata.obsm["Main_trajectory_umap"]
    sc.pl.scatter(
        adata,
        basis="Main_trajectory_umap",
        color=["Main_cell_type", "development_stage"],
        legend_fontsize=8,
        legend_fontoutline=1,
        show=False,
        save="_cao_organogenesis.png",
    )

    # NOTE: this one is good, Main_trajectory_refined_umap
    adata.obsm["X_Main_trajectory_refined_umap"] = adata.obsm[
        "Main_trajectory_refined_umap"
    ]
    sc.pl.scatter(
        adata,
        basis="Main_trajectory_refined_umap",
        color=["Main_cell_type", "development_stage"],
        legend_fontsize=8,
        legend_fontoutline=1,
        show=False,
        save="_cao_organogenesis.png",
    )

    # Sub_trajectory_umap
    adata.obsm["X_Sub_trajectory_umap"] = adata.obsm["Sub_trajectory_umap"]
    sc.pl.scatter(
        adata,
        basis="Sub_trajectory_umap",
        color=["Main_cell_type", "development_stage"],
        legend_fontsize=8,
        legend_fontoutline=1,
        show=False,
        save="_cao_organogenesis.png",
    )

    if dir_name is not None:
        sc.settings.figdir = default_figdir


plots(adata_chondrocyte, dir_name=output_folder / "chondrocyte")
plots(adata_mesenchymal, dir_name=output_folder / "mesenchymal")
plots(adata_mesenchymal_filtered, dir_name=output_folder / "mesenchymal_filtered")
