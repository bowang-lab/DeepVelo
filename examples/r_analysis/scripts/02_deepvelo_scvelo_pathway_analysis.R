library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)
library(ActivePathways)

# Define not in 
`%ni%` <- Negate(`%in%`)

######################## DeepVelo pathway analysis ########################

# Deepvelo genes and associated p-values
dvelo <- fread("../data/MDT_driver_genes.csv")

# Get drivers ordered by gliogenic and gabaergic lineages 
gaba_dvelo <- dvelo[
  order(dvelo$`GABA interneurons_pval`, decreasing = FALSE), 
]
glio_dvelo <- dvelo[
  order(dvelo$`Gliogenic progenitors_pval`, decreasing = FALSE),
]
colnames(gaba_dvelo)[1] <- "gene"
colnames(glio_dvelo)[1] <- "gene"

# Perform more strict bonferroni correction on p-values
gaba_dvelo$gaba_bonfer_p <- p.adjust(
  gaba_dvelo$`GABA interneurons_pval`,
  method = "bonferroni"
)
glio_dvelo$glio_bonfer_p <- p.adjust(
  glio_dvelo$`Gliogenic progenitors_pval`,
  method = "bonferroni"
)

# Create dataframe for gaba and glio significant values and convert to numeric
# matrix
gaba_glio_signif <- merge(gaba_dvelo, glio_dvelo)
gaba_glio_mat_signif <- as.matrix(
  gaba_glio_signif[, c("gaba_bonfer_p", "glio_bonfer_p")]
)
rownames(gaba_glio_mat_signif) <- gaba_glio_signif$gene
diff_exp_mat <- gaba_glio_mat_signif
colnames(diff_exp_mat) <- c("GABAergic", "Gliogenic")

# Run activepathways 
dir.create("../outs/02_activepathways_dvelo_gaba_glio_out", recursive = TRUE)
gaba_glio_ap <- ActivePathways(
  diff_exp_mat,
  "../data/mm_go_mf_bp_reac_feb_2022.gmt",
  geneset.filter = c(5, 2000),
  significant = 0.05, 
  cytoscape.file.tag = "../outs/02_activepathways_dvelo_gaba_glio_out"
)

# Run activepathways per subset independantly 
gaba_mat_signif <- as.matrix(gaba_glio_signif[, c("gaba_bonfer_p")])
glio_mat_signif <- as.matrix(gaba_glio_signif[, c("glio_bonfer_p")])
rownames(gaba_mat_signif) <- gaba_glio_signif$gene
rownames(glio_mat_signif) <- gaba_glio_signif$gene

dir.create("../outs/02_activepathways_dvelo_gaba_out", recursive = TRUE)
gaba_ap <- ActivePathways(
  gaba_mat_signif,
  "../data/mm_go_mf_bp_reac_feb_2022.gmt",
  geneset.filter = c(5, 2000),
  significant = 0.05, 
  cytoscape.file.tag = "../outs/02_activepathways_dvelo_gaba_out"
)

dir.create("../outs/02_activepathways_dvelo_glio_out", recursive = TRUE)
glio_ap <- ActivePathways(
  glio_mat_signif,
  "../data/mm_go_mf_bp_reac_feb_2022.gmt",
  geneset.filter = c(5, 2000),
  significant = 0.05, 
  cytoscape.file.tag = "../outs/02_activepathways_dvelo_glio_out"
)

# Run ActivePathway on top 100, 250, 500, and full list subsets for both gaba and glio
# independantly 
subset_list <- list(100, 250, 500)
gaba_subset_mats <- lapply(subset_list, function(x) {
  subset <- gaba_glio_signif[order(gaba_glio_signif$gaba_bonfer_p, decreasing = FALSE), c("gene", "gaba_bonfer_p")]
  subset[(x+1):nrow(subset), c("gaba_bonfer_p")] <- 1
  colnames(subset)[2] <- paste0("GABAergic ", x)
  return(subset)
})
glio_subset_mats <- lapply(subset_list, function(x) {
  subset <- gaba_glio_signif[order(gaba_glio_signif$glio_bonfer_p, decreasing = FALSE), c("gene", "glio_bonfer_p")]
  subset[(x+1):nrow(subset), c("glio_bonfer_p")] <- 1
  colnames(subset)[2] <- paste0("Gliogenic ", x)
  return(subset)
})
glio_mat_subset <- as.data.frame(glio_mat_signif)
gaba_mat_subset <- as.data.frame(gaba_mat_signif)
glio_mat_subset$gene <- rownames(glio_mat_subset)
gaba_mat_subset$gene <- rownames(gaba_mat_subset)
colnames(glio_mat_subset)[1] <- "Gliogenic all"
colnames(gaba_mat_subset)[1] <- "GABAergic all"

# Merge all together
glio_subsets_merged <- Reduce(merge, glio_subset_mats)
gaba_subsets_merged <- Reduce(merge, gaba_subset_mats)
subsets_list <- list(
  glio_subsets_merged,
  gaba_subsets_merged, 
  glio_mat_subset,
  gaba_mat_subset
)
subsets_merged <- Reduce(merge, subsets_list)
subsets_merged_copy <- subsets_merged

# Remove gene column and perform ActivePathways enrichment
subsets_merged <- subsets_merged[, -c("gene")]
subsets_merged_mat <- as.matrix(subsets_merged)
rownames(subsets_merged_mat) <- subsets_merged_copy$gene

# Perform ActivePathways enrichment analysis 
dir.create("../outs/02_activepathways_dvelo_glio_gaba_top_100_250_500_all_out", recursive = TRUE)
glio_ap <- ActivePathways(
  subsets_merged_mat,
  "../data/mm_go_mf_bp_reac_feb_2022.gmt",
  geneset.filter = c(5, 2000),
  significant = 0.05, 
  cytoscape.file.tag = "../outs/02_activepathways_dvelo_glio_gaba_top_100_250_500_all_out"
)

# Run activepathways on only the top 100 - based on the correlation values instead of 
# the p-values
glio_top_100 <- gaba_glio_signif[
  order(gaba_glio_signif$`Gliogenic progenitors_corr`, decreasing = TRUE)
][0:100]
glio_top_100 <- glio_top_100[, c("gene", "glio_bonfer_p")]
glio_top_100_mat <- as.matrix(glio_top_100[, -1])
rownames(glio_top_100_mat) <- glio_top_100$gene
colnames(glio_top_100_mat) <- "Gliogenic"

dir.create("../outs/02_activepathways_dvelo_glio_top_100_corr_out", recursive = TRUE)
glio_ap <- ActivePathways(
  glio_top_100_mat,
  "../data/mm_go_mf_bp_reac_feb_2022.gmt",
  geneset.filter = c(5, 2000),
  significant = 0.05, 
  cytoscape.file.tag = "../outs/02_activepathways_dvelo_glio_top_100_corr_out"
)

gaba_top_100 <- gaba_glio_signif[
  order(gaba_glio_signif$`GABA interneurons_corr`, decreasing = TRUE)
][0:100]
gaba_top_100 <- gaba_top_100[, c("gene", "gaba_bonfer_p")]
gaba_top_100_mat <- as.matrix(gaba_top_100[, -1])
rownames(gaba_top_100_mat) <- gaba_top_100$gene
colnames(gaba_top_100_mat) <- "GABAergic"

dir.create("../outs/02_activepathways_dvelo_gaba_top_100_corr_out", recursive = TRUE)
gaba_ap <- ActivePathways(
  gaba_top_100_mat,
  "../data/mm_go_mf_bp_reac_feb_2022.gmt",
  geneset.filter = c(5, 2000),
  significant = 0.05, 
  cytoscape.file.tag = "../outs/02_activepathways_dvelo_gaba_top_100_corr_out"
)

######################## scVelo pathway analysis ########################

# scVelo genes and associated p-values
scvelo <- fread("../data/MDT_driver_genes[DYNAMICAL].csv")

# Get drivers ordered by gliogenic and gabaergic lineages 
gaba_scvelo <- scvelo[
  order(scvelo$`GABA interneurons_pval`, decreasing = FALSE), 
]
glio_scvelo <- scvelo[
  order(scvelo$`Gliogenic progenitors_pval`, decreasing = FALSE),
]
colnames(gaba_scvelo)[1] <- "gene"
colnames(glio_scvelo)[1] <- "gene"

# Perform more strict bonferroni correction on p-values
gaba_scvelo$gaba_bonfer_p <- p.adjust(
  gaba_scvelo$`GABA interneurons_pval`,
  method = "bonferroni"
)
glio_scvelo$glio_bonfer_p <- p.adjust(
  glio_scvelo$`Gliogenic progenitors_pval`,
  method = "bonferroni"
)

# Create dataframe for gaba and glio significant values and convert to numeric
# matrix
gaba_glio_signif <- merge(gaba_scvelo, glio_scvelo)
gaba_glio_mat_signif <- as.matrix(
  gaba_glio_signif[, c("gaba_bonfer_p", "glio_bonfer_p")]
)
rownames(gaba_glio_mat_signif) <- gaba_glio_signif$gene
diff_exp_mat <- gaba_glio_mat_signif
colnames(diff_exp_mat) <- c("GABAergic", "Gliogenic")

# Run activepathways on only the top 100 - based on the correlation values instead of 
# the p-values
glio_top_100 <- gaba_glio_signif[
  order(gaba_glio_signif$`Gliogenic progenitors_corr`, decreasing = TRUE)
][0:100]
glio_top_100 <- glio_top_100[, c("gene", "glio_bonfer_p")]
glio_top_100_mat <- as.matrix(glio_top_100[, -1])
rownames(glio_top_100_mat) <- glio_top_100$gene
colnames(glio_top_100_mat) <- "Gliogenic"

dir.create("../outs/02_activepathways_scvelo_glio_top_100_corr_out", recursive = TRUE)
glio_ap <- ActivePathways(
  glio_top_100_mat,
  "../data/mm_go_mf_bp_reac_feb_2022.gmt",
  geneset.filter = c(5, 2000),
  significant = 0.05, 
  cytoscape.file.tag = "../outs/02_activepathways_scvelo_glio_top_100_corr_out"
)

gaba_top_100 <- gaba_glio_signif[
  order(gaba_glio_signif$`GABA interneurons_corr`, decreasing = TRUE)
][0:100]
gaba_top_100 <- gaba_top_100[, c("gene", "gaba_bonfer_p")]
gaba_top_100_mat <- as.matrix(gaba_top_100[, -1])
rownames(gaba_top_100_mat) <- gaba_top_100$gene
colnames(gaba_top_100_mat) <- "GABAergic"

dir.create("../outs/02_activepathways_scvelo_gaba_top_100_corr_out", recursive = TRUE)
gaba_ap <- ActivePathways(
  gaba_top_100_mat,
  "../data/mm_go_mf_bp_reac_feb_2022.gmt",
  geneset.filter = c(5, 2000),
  significant = 0.05, 
  cytoscape.file.tag = "../outs/02_activepathways_scvelo_gaba_top_100_corr_out"
)
