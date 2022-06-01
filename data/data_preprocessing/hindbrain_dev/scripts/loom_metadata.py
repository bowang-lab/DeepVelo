import argparse

import loompy
import numpy as np
import pandas as pd

# Function to get barcodes from loom file and subset
# barcode metadata from vladoiu et al and return
def loom_metadata(loom_file, barcode_file, cluster_file, out_file):
    # Connect to loom file
    ds = loompy.connect(loom_file)
    
    # Load barcode and cluster files
    barcode_meta = pd.read_csv(barcode_file)
    cluster_meta = pd.read_csv(cluster_file, header = None)
    
    # Format data for both
    barcode_meta.columns = ["CellID", "Cluster"]
    cluster_meta.columns = ["Cluster", "Celltype", "Lineage"]
    
    # Create dataframe of CellID and UMI data from loom file
    loom_cellid = ds.ca["CellID"]
    loom_umi = ds.ca["TotalUMIs"]
    loom_meta = pd.DataFrame({
        "CellID": loom_cellid,
        "TotalUMIs": loom_umi
    })
    
    # Consecutively merge three datasets 
    # (loom meta -> barcode meta -> cluster meta)
    loom_barcode_merge = loom_meta.merge(
        barcode_meta,
        how = "inner",
        on = "CellID"
    )
    loom_barcode_cluster_merge = loom_barcode_merge.merge(
        cluster_meta,
        how = "inner",
        on = "Cluster"
    )
    
    # Add timepoint information
    barcode_list = list(loom_barcode_cluster_merge["CellID"].values)
    timepoints = [i.split("_", 1)[0] for i in barcode_list]
    loom_barcode_cluster_merge["Timepoint"] = timepoints
    
    # Save dataframe to tsv 
    loom_barcode_cluster_merge.to_csv(
        out_file,
        sep = "\t",
        header = True,
        index = False
    )
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Output file and input files " +
        "for concatenation of loom data"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Path of output for merged loom metadata"
    )
    parser.add_argument(
        "--loomfile",
        type = str,
        help = "Path of concatenated loom file for developmental data"
    )
    parser.add_argument(
        "--barcodefile",
        type = str,
        help = "Path of barcode-cluster metadata file"
    )
    parser.add_argument(
        "--clusterfile",
        type = str,
        help = "Path of cluster-celltype-lineage metadata file"
    )
    args = parser.parse_args()
    loom_metadata(
        loom_file = args.loomfile,
        barcode_file = args.barcodefile,
        cluster_file = args.clusterfile,
        out_file = args.outfile
    )