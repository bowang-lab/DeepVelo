# The multi-lineages figure
import anndata as ad
import pandas as pd
import pickle
import scvelo as scv
import numpy as np
from time import time
from deepvelo.utils.temporal import latent_time
from deepvelo.utils.scatter import scatter
from deepvelo.utils.velocity import velocity
import os
# %%
# ===================================
# run this block to process the data for MDT gaba interneuron data
# exactly following the original Extended Figure 3 in the nature paper.
# ===================================
filePath = './data/concat.loom'
adata = scv.read(filePath, cache=False)
obs = adata.obs

metaLabels = pd.read_csv('data/MDT_ExtFig3_tsne.csv')
cellID2metaLabels = {x.split('-')[0]: x for x in metaLabels['CellID'].tolist()}

gabaInterneurons = []
cellIDs_str = [x.split('_')[1] for x in obs.index.tolist()]
for cellID in cellIDs_str:  # loop over sequenced cells
    if cellID2metaLabels.get(cellID):  # if found in the matadata of lineages
        if metaLabels.loc[metaLabels['CellID'] == cellID2metaLabels.get(cellID)].ident.values[0] in [4, 17, 12, 15, 10, 9]:
            gabaInterneurons.append(True)
        else:
            gabaInterneurons.append(False)
    else:
        gabaInterneurons.append(False)

# filter cells
cell_mask = np.array(gabaInterneurons, dtype=bool)
assert len(adata.uns) == 0
assert len(adata.obsm) == 0
assert len(adata.varm) == 0
gabaInterneurons_adata = ad.AnnData(
    X=adata.X[cell_mask],
    obs=adata.obs[cell_mask],
    var=adata.var,
    layers={k: v[cell_mask] for k, v in adata.layers.items()}
)

ident2info = {
    4: {'Celltype': 'Neural stem cells', 'lineages': 'Early Progenitors'},
    17: {'Celltype': 'Proliferating VZ progenitors', 'lineages': 'GABAergic'},
    12: {'Celltype': 'VZ progenitors', 'lineages': 'GABAergic'},
    15: {'Celltype': 'Differentiating GABA interneurons', 'lineages': 'GABAergic'},
    10: {'Celltype': 'GABA interneurons', 'lineages': 'GABAergic'},
    9: {'Celltype': 'Gliogenic progenitors', 'lineages': 'Glial'}
}

# add obs columns
obs = gabaInterneurons_adata.obs
cell_types = []
lineages = []
timepoints = [x.split('_')[0] for x in obs.index.tolist()]
cellIDs_str = [x.split('_')[1] for x in obs.index.tolist()]
for cellID in cellIDs_str:
    tmp = metaLabels.loc[metaLabels['CellID'] == cellID2metaLabels.get(cellID)].ident.values[0]
    cell_types.append(ident2info[tmp]['Celltype'])
    lineages.append(ident2info[tmp]['lineages'])
gabaInterneurons_adata.obs['Celltype'] = cell_types
gabaInterneurons_adata.obs['Lineage'] = lineages
gabaInterneurons_adata.obs['Timepoint'] = timepoints

# add obsm columns
X_tsne = np.zeros([len(cellIDs_str), 2], dtype='float64')
for i, cellID in enumerate(cellIDs_str):
    tmp = metaLabels.loc[metaLabels['CellID'] == cellID2metaLabels.get(cellID)]
    X_tsne[i, 0] = tmp.tSNE_1.values[0]
    X_tsne[i, 1] = tmp.tSNE_2.values[0]
gabaInterneurons_adata.obsm['X_tsne'] = X_tsne

gabaInterneurons_adata.write('data/MDT_GABAInterneuraons_v2.h5ad', compression='gzip')
# ===================================
# %% If run for the first time. Prepare the glutamatergic and gabaergic lineage data into separate adata file.
filePath = "/cluster/projects/bwanggroup/for_haotian/loom/concat.loom"
adata = scv.read(filePath, cache=True)

"""
AnnData object with n_obs × n_vars = 91182 × 55421 
    obs: 'TotalUMIs'
    var: 'Accession', 'AccessionVersion', 'Aliases', 'CcdsID', 'Chromosome', 'ChromosomeEnd', 'ChromosomeStart', 'CosmicID', 'DnaBindingDomain', 'FullName', 'GeneType', 'HgncID', 'IsTFi (TcoF-DB)', 'Location', 'LocationSortable', 'LocusGroup', 'LocusType', 'MgdID', 'MirBaseID', 'OmimID', 'PubmedID', 'RefseqID', 'Regulates (TRRUST)', 'RgdID', 'Strand', 'UcscID', 'UniprotID', 'VegaID'
    layers: 'matrix', 'spliced', 'unspliced'
"""
obs = adata.obs

# lineage file path
lineageFilePath = "data/concat_loom_metadata.tsv"
metaLabels = pd.read_csv(lineageFilePath, sep='\t')

to_filter = 'gaba interneurons'
# %% filter for gaba interneurons
if to_filter == 'gaba interneurons':
    gabaInterneurons = []
    cellIDs_str = obs.index if isinstance(obs.index[0], str) else [x_byte.decode('utf-8') for x_byte in obs.index]
    for cellID in cellIDs_str:  # loop over sequenced cells
        if cellID in metaLabels['CellID'].values:  # if found in the matadata of lineages
            if metaLabels.loc[metaLabels['CellID'] == cellID].Celltype.values[0].strip() in ['Neural stem cells', 'Proliferating VZ progenitors', 'VZ progenitors', 'Differentiating GABA interneurons', 'GABA interneurons', 'Gliogenic progenitors']:
                # this cell is either gaba or gluta
                gabaInterneurons.append(True)
            else:
                gabaInterneurons.append(False)
        else:
            gabaInterneurons.append(False)

    # filter cells
    cell_mask = np.array(gabaInterneurons, dtype=bool)
    assert len(adata.uns) == 0
    assert len(adata.obsm) == 0
    assert len(adata.varm) == 0
    gabaInterneurons_adata = ad.AnnData(
        X=adata.X[cell_mask],
        obs=adata.obs[cell_mask],
        var=adata.var,
        layers={k: v[cell_mask] for k, v in adata.layers.items()}
    )

    # add obs columns
    cell_types = []
    lineages = []
    timepoints = []
    obs = gabaInterneurons_adata.obs
    cellIDs_str = obs.index if isinstance(obs.index[0], str) else [x_byte.decode('utf-8') for x_byte in obs.index]
    for cellID in cellIDs_str:
        tmp = metaLabels.loc[metaLabels['CellID'] == cellID]
        cell_types.append(tmp.Celltype.values[0].strip())
        lineages.append(tmp.Lineage.values[0].strip())
        timepoints.append(tmp.Timepoint.values[0].strip())
    gabaInterneurons_adata.obs['Celltype'] = cell_types
    gabaInterneurons_adata.obs['Lineage'] = lineages
    gabaInterneurons_adata.obs['Timepoint'] = timepoints

    # save data
    gabaInterneurons_adata.write('data/MDT_GABAInterneuraons.h5ad', compression='gzip')


# %% filter for gluta and gaba lineages
# make the cell filter list
glutaAndGABACells = []
cellIDs_str = obs.index if isinstance(obs.index[0], str) else [x_byte.decode('utf-8') for x_byte in obs.index]
for cellID in cellIDs_str:  # loop over sequenced cells
    if cellID in metaLabels['CellID'].values:  # if found in the matadata of lineages
        if metaLabels.loc[metaLabels['CellID'] == cellID].Lineage.values[0] in ['GABAergic', 'Glutamatergic']:
            # this cell is either gaba or gluta
            glutaAndGABACells.append(True)
        else:
            glutaAndGABACells.append(False)
    else:
        glutaAndGABACells.append(False)


# filter cells
cell_mask = np.array(glutaAndGABACells, dtype=bool)
assert len(adata.uns) == 0
assert len(adata.obsm) == 0
assert len(adata.varm) == 0
gluta_gaba_cells_adata = ad.AnnData(
    X=adata.X[cell_mask],
    obs=adata.obs[cell_mask],
    var=adata.var,
    layers={k: v[cell_mask] for k, v in adata.layers.items()}
)

# add obs columns
cell_types = []
lineages = []
timepoints = []
obs = gluta_gaba_cells_adata.obs
cellIDs_str = obs.index if isinstance(obs.index[0], str) else [x_byte.decode('utf-8') for x_byte in obs.index]
for cellID in cellIDs_str:
    tmp = metaLabels.loc[metaLabels['CellID'] == cellID]
    cell_types.append(tmp.Celltype.values[0].strip())
    lineages.append(tmp.Lineage.values[0].strip())
    timepoints.append(tmp.Timepoint.values[0].strip())
gluta_gaba_cells_adata.obs['Celltype'] = cell_types
gluta_gaba_cells_adata.obs['Lineage'] = lineages
gluta_gaba_cells_adata.obs['Timepoint'] = timepoints

# sub sampling to 6000 cells
shuffled_index = np.random.permutation(len(gluta_gaba_cells_adata))
gluta_gaba_cells_adata[shuffled_index[:6000]].write('data/MDT_gluta_gaba_cells_adata_6000.h5ad', compression='gzip')

# save data
gluta_gaba_cells_adata.write('data/MDT_gluta_gaba_cells_adata.h5ad', compression='gzip')
# %% it has 89893 cells
# use the first 6000 cells to dry
adataAll = adata
adata = adata[:6000]

# %% settings
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo', transparent=False)  # for beautified visualization
MASK_ZERO = False
DEEPVELO = False  # choice of {True, False, 'ShowTarget'}
DYNAMICAL = False  # whether use the dynamical mode of scvelo and compute latent time
DEEPVELO_FILE = 'scvelo_mat.npz'
data = 'ML'
SURFIX = '[dynamical]' if DYNAMICAL else ''
SURFIX += '[deep_velo]' if DEEPVELO else ''

# # %% found all unspliced reads are zero
# adataAll.layers['unspliced'].max()
# # try read using vlocyto again
# filePath = "data/concat.loom"
# vlm = vcy.VelocytoLoom(filePath)
# %% Preprocessing Data
# here we have the size normalization
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# here comes the NN graph and dynamic estimations
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)
# %% Compute velocity and velocity graph
# import pudb; pudb.set_trace()
if DYNAMICAL:
    scv.tl.recover_dynamics(adata)
    velocity(adata, mode='dynamical', mask_zero=MASK_ZERO)
else:
    velocity(adata, mask_zero=MASK_ZERO)
# %% output and change the velocity
to_save = {
    'Ux_sz': adata.layers['Mu'].T,
    'Sx_sz': adata.layers['Ms'].T,
    'velo': adata.layers['velocity'].T,
    'conn': adata.obsp['connectivities'].T  # (features, cells)
}  # have to input in dimmention order (1999 genes, 2930 cells)
with open('./data/scveloDG.npz', 'wb') as f:
    pickle.dump(to_save, f)
if DEEPVELO == 'ShowTarget':
    print('computing target velocities')
    n_genes, batch_size = adata.layers['velocity'].T.shape
    from data_loader.data_loaders import VeloDataset
    ds = VeloDataset(data_dir='./data/scveloDG.npz')
    velo_mat = ds.velo.numpy()
    assert adata.layers['velocity'].shape == velo_mat.shape
    adata.layers['velocity'] = velo_mat  # (2930 cells, 1999 genes)
elif DEEPVELO:
    n_genes, batch_size = adata.layers['velocity'].T.shape
    now = time()
    os.system(
        f'python train.py -c config_figure3.json --ng {n_genes} --bs {batch_size} --ot {DEEPVELO_FILE} --dd ./data/scveloDG.npz')
    print(f'finished in {time()-now:.2f}s')

    # load
    velo_mat = np.load(f'./data/{DEEPVELO_FILE}')
    assert adata.layers['velocity'].shape == velo_mat['velo_mat'].shape
    adata.layers['velocity'] = velo_mat['velo_mat']  # (2930 cells, 1999 genes)
    # adata.layers['velocity'] = - adata.layers['velocity']
scv.tl.velocity_graph(adata)
# %% generate umap if need
if not ('X_umap' in adata.obsm or 'tsne' in adata.obsm):
    scv.tl.umap(adata)  # this will add adata.obsm: 'X_umap'

# %% plot panel a
scv.pl.velocity_embedding_stream(adata, basis='umap', dpi=300, save=f'figure5/[{data}]velo_emb_stream{SURFIX}.png',
                                 legend_fontsize=9)
