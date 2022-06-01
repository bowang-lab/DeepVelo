import sys
import os
import numpy as np
import pandas as pd

def write_metadata(save_loc):
    # Define subsets, tech, and target cells
    # Target cell number defined from
    # https://www.nature.com/articles/s41586-019-1158-7 
    # rounded up to the nearest thousand
    samples = [
        "E10", 
        "E12", 
        "E14", 
        "E16", 
        "E18", 
        "P0", 
        "P5", 
        "P7", 
        "P14"
    ]
    technology = "10xv2"
    target_cells = [
        "8000",
        "8000",
        "7000",
        "8000",
        "6000",
        "5000",
        "12000",
        "8000",
        "5000"
    ]
    
    # Create and save df as tsv
    meta_df = pd.DataFrame({
        "name" : samples,
        "technology" : technology,
        "targetnumcells": target_cells
    })
    meta_df.to_csv(
        save_loc,
        sep = "\t",
        index = False
    )

if __name__ == '__main__':
    write_metadata(sys.argv[1])