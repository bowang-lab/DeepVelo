library(data.table)

# Change to working directory of alevin runs
setwd("/cluster/projects/bwanggroup/cao_et_al/alevin/runs")

# Load the cell annotations, with and without the cell-types
cell_annos <- readRDS("cell_annotations.rds")
cell_annos_with_type <- fread(
    "cell_annotations_with_labels.csv",
    sep = ","
)

# Merge the two together (left merge on previously saved cell annotations)
cell_annos_merged_with_type <- merge(
    cell_annos,
    cell_annos_with_type,
    by = "sample",
    all.x = TRUE
)

# Reload the rows corresponding to the data cell barcodes
quants_rows <- readRDS("quants_row_names.rds")

# Create a dataframe of the quants rows and merge with the cell annotations
quants_rows_df <- data.frame(
    sample_name = quants_rows,
    stringsAsFactors = FALSE
)
quants_rows_df$id <- 1:nrow(quants_rows_df)

# Merge the quants rows with the cell annotations
quants_rows_cell_annos_merged <- merge(
    quants_rows_df,
    cell_annos_merged_with_type,
    by = "sample_name",
    all.x = TRUE
)

# Get index for reordering based on sample name
reorder_idx <- match(
    quants_rows_df$sample_name, quants_rows_cell_annos_merged$sample_name
)

# Order the cell annotations by the order of the quants rows
quants_rows_cell_annos_merged <- quants_rows_cell_annos_merged[reorder_idx,]

# Ensure that ids from quants rows cell annos merged are in same order as quants rows
if (all(quants_rows_df$sample_name == quants_rows_cell_annos_merged$sample_name)){
    print("Cell annotations are in same order as quants rows")
} else {
    stop("Cell annotations are not in same order as quants rows")
}

# Re-save the merged cell annotations
saveRDS(
    quants_rows_cell_annos_merged,
    "cell_annotations_updated_aligned.rds"
)