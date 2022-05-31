# Driver-gene analysis of mouse hindbrain developmental data

This folder contains the necessary scripts and environment to the run the R-based analysis for the mouse hindbrain developmental data, and reproduce the 
figures related to this data in the main manuscript. 

## Instructions 

1. Install the conda environment that contains the necessary R-build and libraries for running the analysis:
```
# This will create an env named 'deepvelo_r_analysis'
conda env create -f env.yaml 
```

2. Create outs and data folders, for output of scripts and necessary data, respectively
```
mkdir data
mkdir -p outs/figures
```

3. Copy the necessary files to run the scripts from the DeepVelo resource folder into the data dir
```
cp ../deepvelo_data/deepvelo_outputs/MDT_driver_genes[DYNAMICAL].csv data # scVelo driver genes
cp ../deepvelo_data/deepvelo_outputs/MDT_driver_genes.csv data # DeepVelo driver genes 
cp ../deepvelo_data/metadata_files/41586_2019_1158_MOESM4_ESM.csv data # Vladiou et al. marker genes
cp ../deepvelo_data/metadata_files/HGNC_AllianceHomology.rpt data # Full list of mouse genes  
cp ../deepvelo_data/metadata_files/mm_go_mf_bp_reac_feb_2022.gmt data # Pathway GMT files for pathway analysis
cp ../deepvelo_data/metadata_files/02_human_mouse_tfs_matched.tsv data # List of matched human-mouse TFs 
```
4. Run the R-rscripts in order to reproduce analysis 
```
conda activate deepvelo_r_analysis
cd scripts
Rscript --verbose 00_extra_package_installs.R
Rscript --verbose 01_deepvelo_marker_analysis.R
Rscript --verbose 02_deepvelo_scvelo_pathway_analysis.R
Rscript --verbose 03_deepvelo_scvelo_activepathways_results_analysis.R
Rscript --verbose 04_deepvelo_tf_driver_analysis.R
Rscript --verbose 05_supplementary_table_formatting.R
```
