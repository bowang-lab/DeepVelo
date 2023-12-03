library(Matrix)
library(data.table)

# Change to working directory of alevin runs
setwd("/cluster/projects/bwanggroup/cao_et_al/alevin/runs")

# Get all of the directories that start with "res"
dirs <- list.files(pattern = "res")

# Extract the run name from the directory name
srr_names <- lapply(dirs, function(fp) {
  strsplit(fp, "_")[[1]][2]
})

# Load quant files from each directory
quantFiles <- lapply(dirs, function(fp) {
  readMM(file.path(fp, "alevin/quants_mat.mtx"))
})
names(quantFiles) <- srr_names

# Load the cell files from each directory
cellFiles <- lapply(dirs, function(fp) {
  fread(
    file.path(fp, "alevin/quants_mat_rows.txt"),
    header = F,
    col.names = "cellbarcode"
  )
})
names(cellFiles) <- srr_names

# Load one file for the gene names
genes <- fread(
    file.path(dirs[1], "alevin/quants_mat_cols.txt"),
    header = F,
    col.names = "geneId"
)

# Change to top level dir (alevin)
setwd("..")

# Load cell annotations
cell_anno <- fread("cell_annotations.txt")
cell_anno <- cell_anno[run %in% names(quantFiles)]
cellb <- sapply(cell_anno$cb,function(cb){
  if(nchar(cb)==20){
    return(paste0(cb,"A"))
  } else {
    return(paste0(cb,"AC"))
  }
})
cell_anno[,cb:=cellb]

# Check if the dimensions of quantFiles and cellFiles match up
for (i in 1:length(quantFiles)){
  if (dim(quantFiles[[i]])[1] != dim(cellFiles[[i]])[1]){
    print(paste0("Mismatch in dimensions for ", names(quantFiles)[i]))
  }
}

# subset the quantFiles to include cells used in the study
quantFiles <- lapply(1:length(quantFiles), function(i){
  quantfile_sub <- quantFiles[[i]]
  cellfile_sub <- cellFiles[[i]]
  rownames(quantfile_sub) <- cellfile_sub$cellbarcode
  colnames(quantfile_sub) <- genes$geneId
  cbi <- cell_anno[run == names(quantFiles)[i], cb]
  rws <- cellfile_sub$cellbarcode %in% cbi
  quantfile_sub <- quantfile_sub[rws,]
  return(quantfile_sub)
})
names(quantFiles) <- srr_names

# get matching annotations for each quantFile
cell_annotations <- lapply(1:length(quantFiles), function(i){
  c_anno <- cell_anno[run == names(quantFiles)[i],]
  cb_matches <- match(rownames(quantFiles[[i]]), c_anno$cb)
  cb_matches <- cb_matches[!is.na(cb_matches)]
  return(c_anno[cb_matches])
})

gene_annotations <- fread("gene_annotation_table.txt")
g_matches <- match(colnames(quantFiles[[1]]),gene_annotations$Geneid)
gene_annotations <- gene_annotations[g_matches]

# Don't add the spliced, unspliced and ambiguous counts, 
# but keep an index of which columns they correspond to
gene_col_idx = data.frame(
  "spliced" = 1:55401, 
  "unspliced" = 55402:110802, 
  "ambiguous" = 110803:166203
)

# For each quantfile add the full correct rowname based on the
# cell barcode and the gene id
quantFiles <- lapply(1:length(quantFiles), function(i){
  quantfile_sub <- quantFiles[[i]]
  quantfile_name <- names(quantFiles)[i]
  quantfile_rownames <- rownames(quantfile_sub)
  rownames(quantfile_sub) <- paste0(
    "sci3-me-", quantfile_name, ".", quantfile_rownames
  )
  return(quantfile_sub)
})
names(quantFiles) <- srr_names

# Similarly for the cell annotations add a column based on this procedure
cell_annotations <- lapply(1:length(quantFiles), function(i) {
  quantfile_name <- names(quantFiles)[i]
  cell_annotations[[i]]$sample_name <- paste0(
    "sci3-me-", quantfile_name, ".", cell_annotations[[i]]$cb
  )
  return(cell_annotations[[i]])
})

# Concatenate the quant files and cell annotations
quants <- do.call(rbind, quantFiles)
cell_annotations <- do.call(rbind, cell_annotations)

# Save quants row and column names
saveRDS(rownames(quants), "runs/quants_row_names.rds")
saveRDS(colnames(quants), "runs/quants_col_names.rds")

# Save the quants, cell annotations, and gene indices
saveRDS(quants, "runs/quants.rds")
saveRDS(cell_annotations, "runs/cell_annotations.rds")
saveRDS(gene_col_idx, "runs/gene_col_idx.rds")
