library(ggplot2)
library(ggthemes)
library(data.table)
library(rjson)
library(plyr)
library(stringr)
library(tibble)
library(RColorBrewer)
library(dplyr)
library(scales)
library(ggimage)
library(cowplot)

# Change to data directory
setwd("../data/")

# Load all of the json files
json_files <- list.files(pattern = ".json")

# Map the json file names to dataset names 
json_names <- stringr::str_split_fixed(
    json_files, 
    pattern = "_eval_results", 
    n = 3
)[,1]
dataset_names <- plyr::mapvalues(
    json_names,
    from = c(
        "figure2_dentategyrus",
        "figure3_pancreas",
        "figure4_hindbrain",
        "organogenesis_chondrocyte",
        "la_manno_hippocampus"
    ),
    to = c(
        "Dentate gyrus",
        "Pancreas",
        "Hindbrain",
        "Chondrocyte",
        "Hippocampus"
    )
)

# Load all of the files into a dataframe object
jsons_loaded <- lapply(json_files, function(x) {
  json <- fromJSON(file = x)
  return(json)
})

dfs_per_dataset <- lapply(jsons_loaded, function(x) {
  dir_scores <- unlist(x$direction_score[1:4])
  oc_scores <- unlist(x$overall_consistency[1:4])
  csc_scores <- unlist(x$celltype_consistency[1:4])
  methods <- c("DeepVelo", "scVelo (dynamical)", "Velocyto (steady state)", "CellDancer")
  results_df <- data.frame(
    "Method" = methods, 
    "Direction score" = dir_scores, 
    "Overall consistency" = oc_scores, 
    "Cell-type wise consistency" = csc_scores
  )
  return(results_df)
})

# Add dataset names to the dataframes
dfs_per_dataset <- mapply(
    FUN = function(x, y) {
        x$Dataset <- y
        return(x)
    },
    x = dfs_per_dataset,
    y = dataset_names,
    SIMPLIFY = FALSE
)

# Combine all of the dataframes into one
dfs_combined <- Reduce(rbind, dfs_per_dataset)

# Save the processed dataframe
fwrite(
    dfs_combined, 
    file = "full_scores_processed.csv",
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

# Source the scib knit table function
source("../helpers/scIB_knit_table.R")

# Scale the scores between 0 and 1 
relevant_columns = c("Direction score", "Overall consistency", "Cell-type wise consistency")
colnames(dfs_combined) <- c(
    "Method", 
    "Direction score", 
    "Overall consistency", 
    "Cell-type wise consistency", 
    "Dataset"
)
dfs_combined[, relevant_columns] <- lapply(
    dfs_combined[, relevant_columns], 
    function(x) {
        return((x - min(x)) / (max(x) - min(x)))
    }
)

# Create an overall score column that uses the formula - 0.5*direction_score + 0.25*overall_consistency + 0.25*celltype_consistency
dfs_combined$`Overall score` <- 
    0.5*dfs_combined$`Direction score` + 
    0.25*dfs_combined$`Overall consistency` + 
    0.25*dfs_combined$`Cell-type wise consistency`

# Prepare the data, column info, row info, and palettes 
data <- dfs_combined
data$Embedding <- rep("graph", nrow(data))
data <- data[, c("Method", "Dataset", "Overall score", "Overall consistency", "Cell-type wise consistency", "Direction score")]
row_info <- data.frame("id" = data$Method, "group" = NA)
column_info <- data.frame(
    "id" = colnames(data), 
    "group" = c("Text", "Text", "S0", "S1", "S2", "S3"),
    "geom" = c("text", "text", "bar", "bar", "bar", "bar"),
    "width" = c(6, 2.5, 2, 2, 2, 2),
    overlay = FALSE
)
palettes <- list(
    "S0" = "YlOrRd",
    "S1" = "YlGnBu",
    "S2" = "BuPu",
    "S3" = "RdPu"
)
g <- scIB_knit_table(data = data, column_info = column_info, row_info = row_info, palettes = palettes, usability = F)
now <- Sys.time()
outdir <- "../outs/figures/"
ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), "velocity_summary_metrics.pdf"), g, device = cairo_pdf, width = 297, height = 420, units = "mm")
ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), "velocity_summary_metrics.tiff"), g, device = "tiff", dpi = "retina", width = 297, height = 420, units = "mm")
ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), "velocity_summary_metrics.png"), g, device = "png", dpi = "retina", width = 297, height = 420, units = "mm")