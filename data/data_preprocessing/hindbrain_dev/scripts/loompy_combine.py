import loompy
import argparse

# Function to take in arbitrary number of loompy
# files, return and save concatenated 
def loompy_combine(out_file, loom_file_list):
    # Define params
    files = loom_file_list
    output_file = out_file
    key = "Accession"
    
    # Concat and save loom file
    loompy.combine(
        files = files,
        output_file = output_file,
        key = key
    )
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Output file and input files " +
        "for concatenation of loom data"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Path of output for concatenated loom file"
    )
    parser.add_argument(
        "--loomfiles",
        type = str,
        nargs = "*",
        help = "Any number of loom file paths"
    )
    args = parser.parse_args()
    loompy_combine(
        out_file = args.outfile,
        loom_file_list = args.loomfiles
    )