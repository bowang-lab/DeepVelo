# Preprocessing of the Cao et al. organogenesis data


### This directory contains scripts for preprocessing mouse hindbrain developmental data from Vladoiu et al. for the subsequent Velocity analysis with scVelo and DeepVelo

## Instructions:

### Install the R processing conda environment
```
conda env create -f environment.yaml
```

### Downloading the Cao et al. organogenesis data
Download the SRR files corresponding to the fastqs from the SRA archive - https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA490754&o=acc_s%3Aa.
Create a `data` directory and save the SRR files in that directory.

Convert the SRA files to fastq using `scripts/sra_to_fastq.sh`

### Downloading the alevin fry quantification data 
The data and scripts used to reprocess the Cao et al. data are taken from the Alevin-Fry quantification pipeline (https://combine-lab.github.io/alevin-fry-tutorials/2021/sci-rna-seq3/), with two exceptions:

1) All of the samples are processed 
2) The spliced and unspliced reads are kept separate for RNA velocity analysis

To download the necessary data:
```
sh scripts/alevin_fry_data_dl.sh
```

### Running the Alevin-fry quantification pipeline 
To quantify the fastqs and obtain `.rds` files with the spliced and unspliced counts, with matched barcodes, run:
```
Rscript cao_et_al_get_quant.R
Rscript cao_et_al_update_anno.R
```

### Converting the counts to h5ad and running the RNA velocity analysis

Set the `data_folder` path to the output data in the previous step in `process2anndata.py` and run
```python
python process2anndata.py
```
The h5ad files will be stored in the `output_folder`.
