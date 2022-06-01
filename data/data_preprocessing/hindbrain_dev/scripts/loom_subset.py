import argparse

import numpy as np
import pandas as pd
import scanpy as sc 
import anndata as ann
import loompy as lp

# Define main function to subset loom file (after importing)
# as anndata for the relevant celltypes/lineages used in analysis
def loom_subset(full_loomfile, loom_metadata, out_file):
    # Read in loom file using scanpy 
    full_data = sc.read_loom(
        full_loomfile,
        sparse = True
    )

    # Rename index    
    full_data.obs.index.rename("id", inplace = True) 

    # Add new column to obs using CellIDs from index
    full_data.obs["CellID"] = full_data.obs.index   
 
    # Read in loom metadata file
    loom_metadata = pd.read_csv(
        loom_metadata,
        sep = "\t"
    )
    
    # Pull out obs columns from full data 
    full_data_obs = full_data.obs.copy()
    
    # Merge anndata obs and loom metadata 
    full_data_obs_meta = full_data_obs.merge(
        loom_metadata,
        how = "inner",
        on = "CellID"
    )
    
    # Subset original anndata object for only matched cellids and combine information 
    full_data = full_data[full_data.obs["CellID"].isin(full_data_obs_meta["CellID"])]
    
    # Subset metadata for only matched cellids
    loom_metadata_sub = loom_metadata[loom_metadata["CellID"].isin(full_data_obs_meta["CellID"])]
    
    # Append subset loom metadata to original anndata object and reset index
    full_data.obs = full_data.obs.merge(
        loom_metadata_sub,
        how = "inner",
        on = "CellID"
    )
    full_data.obs.set_index("CellID", inplace = True)
    
    # Subset data to only include celltypes/lineages used in analysis
    full_data_sub = full_data[
        full_data.obs["Celltype"].isin([
            "Differentiating GABA interneurons", 
            "GABA interneurons", 
            "Gliogenic progenitors",
            "Neural stem cells", 
            "Proliferating VZ progenitors", 
            "VZ progenitors"
        ])
    ]
    
    # Return subsetted data
    full_data_sub.write(
        out_file
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Output file and input files " +
        "for subsetting of loom data"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Path of output for subset data as anndata object"
    )
    parser.add_argument(
        "--loomfile",
        type = str,
        help = "Path of concatenated loom file for developmental data"
    )
    parser.add_argument(
        "--metadata",
        type = str,
        help = "Path of inner joined metadata for loom file"
    )
    args = parser.parse_args()
    loom_subset(
        full_loomfile = args.loomfile,
        loom_metadata = args.metadata,
        out_file = args.outfile
    )
    
