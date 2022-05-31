library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)
library(ggvenn)

# Define not in 
`%ni%` <- Negate(`%in%`)

# Dynamical genes 
scvelo <- fread("../data/MDT_driver_genes[DYNAMICAL].csv")

# Deepvelo genes
dvelo <- fread("../data/MDT_driver_genes.csv")

colnames(scvelo)[1] <- "gene"
colnames(dvelo)[1] <- "gene"

# Load transcription factor dataframe 
tf_df <- fread("../data/02_human_mouse_tfs_matched.tsv")

# Subset for top 100 driver genes for GABAergic and Gliogenic for 
# both deepvelo and scvelo 
scvelo_top_100_gaba <- scvelo[
  order(scvelo$`GABA interneurons_corr`, decreasing = TRUE)
][0:100]
scvelo_top_100_glio <- scvelo[
  order(scvelo$`Gliogenic progenitors_corr`, decreasing = TRUE)
][0:100]

dvelo_top_100_gaba <- dvelo[
  order(dvelo$`GABA interneurons_corr`, decreasing = TRUE)
][0:100]
dvelo_top_100_glio <- dvelo[
  order(dvelo$`Gliogenic progenitors_corr`, decreasing = TRUE)
][0:100]

scvelo_top_100_gaba <- scvelo_top_100_gaba[,c("gene", "GABA interneurons_corr")]
scvelo_top_100_glio <- scvelo_top_100_glio[,c("gene", "Gliogenic progenitors_corr")]
dvelo_top_100_gaba <- dvelo_top_100_gaba[,c("gene", "GABA interneurons_corr")]
dvelo_top_100_glio <- dvelo_top_100_glio[,c("gene", "Gliogenic progenitors_corr")]

# Append TF info 
scvelo_top_100_gaba$tf <- ifelse(
  scvelo_top_100_gaba$gene %in% tf_df$mm_gene_name, 
  "Yes",
  "No"
)
scvelo_top_100_glio$tf <- ifelse(
  scvelo_top_100_glio$gene %in% tf_df$mm_gene_name, 
  "Yes",
  "No"
)
dvelo_top_100_gaba$tf <- ifelse(
  dvelo_top_100_gaba$gene %in% tf_df$mm_gene_name, 
  "Yes",
  "No"
)
dvelo_top_100_glio$tf <- ifelse(
  dvelo_top_100_glio$gene %in% tf_df$mm_gene_name, 
  "Yes",
  "No"
)

# Compare TF lengths across subsets
tf_overlaps <- list(
  "scVelo \nGABAergic" = scvelo_top_100_gaba$gene[which(scvelo_top_100_gaba$tf %in% "Yes")],
  "DeepVelo \nGABAergic" = dvelo_top_100_gaba$gene[which(dvelo_top_100_gaba$tf %in% "Yes")],
  "scVelo \nGliogenic" = scvelo_top_100_glio$gene[which(scvelo_top_100_glio$tf %in% "Yes")],
  "DeepVelo \nGliogenic" = dvelo_top_100_glio$gene[which(dvelo_top_100_glio$tf %in% "Yes")]
)

# Get VennDiagram of comparison
ggvenn(
  tf_overlaps,
  set_name_size = 5
)
ggsave(
  "../outs/figures/06_dvelo_scvelo_top_100_gaba_glio_tf_venn.pdf",
  height = 6,
  width = 12
)

# Pick out transcription factors found unique to deepvelo for both lineages
dvelo_top_100_gaba$gene_unique <- ifelse(
  dvelo_top_100_gaba$gene %in% scvelo_top_100_gaba$gene,
  "No",
  "Yes"
)
dvelo_top_100_glio$gene_unique <- ifelse(
  dvelo_top_100_glio$gene %in% scvelo_top_100_glio$gene,
  "No",
  "Yes"
)

# Create barplot of TF overlap between two methods for driver genes 
tf_overlap_deepvelo <- data.frame(
  "GABAergic" = length(which(dvelo_top_100_gaba$tf == "Yes")),
  "Gliogenic" = length(which(dvelo_top_100_glio$tf == "Yes")),
  "Method" = "DeepVelo"
)
tf_overlap_scvelo <- data.table(
  "GABAergic" = length(which(scvelo_top_100_gaba$tf == "Yes")),
  "Gliogenic" = length(which(scvelo_top_100_glio$tf == "Yes")),
  "Method" = "scVelo"    
)
tf_overlap <- rbind(tf_overlap_deepvelo, tf_overlap_scvelo)
tf_overlap_melted <- reshape2::melt(tf_overlap)
colnames(tf_overlap_melted) <- c("Method", "Lineage", "TF_Overlap")

ggplot(data = tf_overlap_melted, aes(
  x = Lineage, 
  y = TF_Overlap,
  fill = factor(Method, levels = c("scVelo", "DeepVelo"))
)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +
  theme_few() +
  labs(
    fill = "Method",
    x = "Lineage",
    y = "TF overlap in top 100 driver genes"
  ) +
  geom_text(
    aes(Lineage, label = TF_Overlap), 
    position = position_dodge(width = 1), 
    vjust = -0.2,
    size = 5.5
  ) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 20)) +
  theme(aspect.ratio = 1)
ggsave(
  "../outs/figures/04_deepvelo_glio_gaba_tf_overlap_comparison.pdf",
  height = 7, 
  width = 7
)