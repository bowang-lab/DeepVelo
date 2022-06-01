# Preprocessing of the mouse hindbrain development data for analysis

### This directory contains scripts for preprocessing mouse hindbrain developmental data from Vladoiu et al. for the subsequent Velocity analysis with scVelo and DeepVelo

## Instructions:

### Downloading data and kallisto files 
To be added after larger data repository prepared. 

### Running Snakemake pipeline

1) Assuming conda is installed, resolve and install environment 
```
conda env create -f environment.yaml
```

2) Ensure files paths are correct in `config.json` by setting script and working directory

3) Test Snakemake (dry run) to ensure validity of DAG and files 
```
conda activate hindbrain_velocity
snakemake -np
```

4) Run snakemake pipeline either locally (a) or using a custom configuration on HPC

(Note for this step, creating a cluster.json or profile will be necessary. See Snakemake documentation for details - https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html)
```
a) snakemake --cores=8
```
