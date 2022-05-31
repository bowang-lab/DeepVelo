# DeepVelo - A Deep Learning-based velocity estimation tool with cell-specific kinetic rates

This is the official implementation of the [DeepVelo](https://www.biorxiv.org/content/10.1101/2022.04.03.486877) method.
DeepVelo employs cell-specific kinetic rates and provides more accurate RNA velocity estimates for complex differentiation and lineage decision events in heterogeneous scRNA-seq data. Please check out the paper for more details.

![alt text](https://user-images.githubusercontent.com/11674033/171066682-a899377f-fae1-452a-8b67-8bc8c244b641.png)

## Installation

```bash
pip install deepvelo
```

The `dgl` package is required, the cpu version is installed by default. Feel free to install the [dgl cuda](https://www.dgl.ai/pages/start.html) version for GPU acceleration.

```bash
pip install dgl-cu101>=0.4.3 # an example for CUDA 10.1

```

### Install the development version

We use poetry to manage dependencies.

```bash
poetry install
```

This will install the exact versions in the provided [poetry.lock](poetry.lock) file. If you want to install the latest version for all dependencies, use the following command.

```bash
poetry update
```

## Usage

We provide a number of notebooks in the [exmaples](examples) folder to help you get started. DeepVelo fullly integrates with [scanpy](https://scanpy.readthedocs.io/en/latest/) and [scVelo](https://scvelo.readthedocs.io/). The basic usage is as follows:

```python
import deepvelo as dv
import scvelo as scv

adata = ... # load your data in AnnData format

# preprocess the data
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

# run DeepVelo using the default configs
trainer = dv.train(adata, dv.Constants.default_configs)
# this will train the model and predict the velocity vectore. The result is stored in adata.layers['velocity']. You can use trainer.model to access the model.
```
