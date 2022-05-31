library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)

# Define not in 
`%ni%` <- Negate(`%in%`)

# Dynamical genes 
dyn <- fread("../data/MDT_driver_genes[DYNAMICAL].csv")
colnames(dyn)[1] <- "gene"

# Deepvelo genes
dvelo <- fread("../data/MDT_driver_genes.csv")
colnames(dvelo)[1] <- "gene"

# Load marker list from Vladiou et al.
vla_markers <- fread("../data/41586_2019_1158_MOESM4_ESM.csv")

# Fill down values
vla_markers$Annotations <- na_if(vla_markers$Annotations, "")
vla_markers <- vla_markers %>% tidyr::fill(Annotations, .direction = "down")

# Subset for GABAergic and gliogenic markers
gaba_markers <- vla_markers[vla_markers$Annotations %in% c(
  "GABA interneurons",
  "Differentiating GABA interneurons"
)]
glio_markers <- vla_markers[vla_markers$Annotations %in% c(
  "Gliogenic progenitors"
)]

# Load mouse gene list
mouse_genes <- fread("../data/HGNC_AllianceHomology.rpt")
mouse_gene_list <- unique(mouse_genes$`MGI Accession ID`)

### Gliogenic 

# Perform hypergeom test between markers and dynamical
dvelo_list <- dvelo[order(dvelo$`Gliogenic progenitors_corr`, decreasing = TRUE)]$gene[0:100]
dyn_list <- dyn[order(dyn$`Gliogenic progenitors_corr`, decreasing = TRUE)]$gene[0:100]
vla_list <- unique(glio_markers$gene)
vla_mouse_diff <- vla_list[which(vla_list %ni% mouse_gene_list)]

# Test using dynamical model
A1 <- length(which(dyn_list %in% vla_list))
A2 <- length(which(dyn_list %ni% vla_list))
A3 <- length(which(dyn_list %in% vla_mouse_diff))
A4 <- length(which(dyn_list %ni% vla_mouse_diff))

dyn_chiseq_table <- matrix(c(A1, A2, A3, A4), nrow = 2, ncol = 2)
dyn_fishers <- fisher.test(dyn_chiseq_table, y = NULL, workspace = 200000, hybrid = FALSE,
                           hybridPars = c(expect = 5, percent = 80, Emin = 1),
                           control = list(), or = 1, alternative = "two.sided",
                           conf.int = TRUE, conf.level = 0.95,
                           simulate.p.value = FALSE, B = 2000)
dyn_fishers

# Test using deepvelo model
B1 <- length(which(dvelo_list %in% vla_list))
B2 <- length(which(dvelo_list %ni% vla_list))
B3 <- length(which(dvelo_list %in% vla_mouse_diff))
B4 <- length(which(dvelo_list %ni% vla_mouse_diff))
dvelo_chiseq_table <- matrix(c(B1, B2, B3, B4), nrow = 2, ncol = 2)
dvelo_fishers <- fisher.test(dvelo_chiseq_table, y = NULL, workspace = 200000, hybrid = FALSE,
                           hybridPars = c(expect = 5, percent = 80, Emin = 1),
                           control = list(), or = 1, alternative = "two.sided",
                           conf.int = TRUE, conf.level = 0.95,
                           simulate.p.value = FALSE, B = 2000)
dvelo_fishers

# Save contingency table results
glio_table <- data.frame(
  "scVelo" = c(A1, A2),
  "DeepVelo" = c(B1, B2)
)
rownames(glio_table) <- c("Markers", "Non-markers")
fwrite(
  glio_table,
  "../outs/01_deepvelo_top_100_glio_drivers_marker_overlap.tsv",
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# Save deepvelo and scvelo top 100 list
fwrite(
  dvelo[order(dvelo$`Gliogenic progenitors_corr`, decreasing = TRUE)][0:100],
  "../outs/01_deepvelo_top_100_glio_drivers_corr.tsv",
  sep = "\t"
)
fwrite(
  dyn[order(dyn$`Gliogenic progenitors_corr`, decreasing = TRUE)][0:100],
  "../outs/01_scvelo_top_100_glio_drivers_corr.tsv",
  sep = "\t"
)

### Test for gabergic progenitors - test with top 100

dvelo_list <- dvelo[order(dvelo$`GABA interneurons_corr`, decreasing = TRUE)]$gene[0:100]
dyn_list <- dyn[order(dyn$`GABA interneurons_corr`, decreasing = TRUE)]$gene[0:100]
vla_list <- unique(gaba_markers$gene)
vla_mouse_diff <- vla_list[which(vla_list %ni% mouse_gene_list)]

# Test using dynamical model
A1 <- length(which(dyn_list %in% vla_list))
A2 <- length(which(dyn_list %ni% vla_list))
A3 <- length(which(dyn_list %in% vla_mouse_diff))
A4 <- length(which(dyn_list %ni% vla_mouse_diff))

dyn_chiseq_table <- matrix(c(A1, A2, A3, A4), nrow = 2, ncol = 2)
dyn_fishers <- fisher.test(dyn_chiseq_table, y = NULL, workspace = 200000, hybrid = FALSE,
                           hybridPars = c(expect = 5, percent = 80, Emin = 1),
                           control = list(), or = 1, alternative = "two.sided",
                           conf.int = TRUE, conf.level = 0.95,
                           simulate.p.value = FALSE, B = 2000)
dyn_fishers

# Test using deepvelo model
B1 <- length(which(dvelo_list %in% vla_list))
B2 <- length(which(dvelo_list %ni% vla_list))
B3 <- length(which(dvelo_list %in% vla_mouse_diff))
B4 <- length(which(dvelo_list %ni% vla_mouse_diff))
dvelo_chiseq_table <- matrix(c(B1, B2, B3, B4), nrow = 2, ncol = 2)
dvelo_fishers <- fisher.test(dvelo_chiseq_table, y = NULL, workspace = 200000, hybrid = FALSE,
                             hybridPars = c(expect = 5, percent = 80, Emin = 1),
                             control = list(), or = 1, alternative = "two.sided",
                             conf.int = TRUE, conf.level = 0.95,
                             simulate.p.value = FALSE, B = 2000)
dvelo_fishers

# Save contingency table results
gaba_table <- data.frame(
  "scVelo" = c(A1, A2),
  "DeepVelo" = c(B1, B2)
)
rownames(gaba_table) <- c("Markers", "Non-markers")
fwrite(
  gaba_table,
  "../outs/01_deepvelo_top_100_gaba_drivers_marker_overlap.tsv",
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# Save deepvelo and scvelo top 100 list
fwrite(
  dvelo[order(dvelo$`GABA interneurons_corr`, decreasing = TRUE)][0:100],
  "../outs/01_deepvelo_top_100_gaba_drivers_corr.tsv",
  sep = "\t"
)
fwrite(
  dyn[order(dyn$`GABA interneurons_corr`, decreasing = TRUE)][0:100],
  "../outs/01_scvelo_top_100_gaba_drivers_corr.tsv",
  sep = "\t"
)


# Subset data for GABA markers from each and plot correlations 
dyn_gaba_subset <- dyn[
  (dyn$gene %in% gaba_markers$gene) &
  (dyn$`GABA interneurons_corr` > 0)
]
dvelo_gaba_subset <- dvelo[
  (dvelo$gene %in% gaba_markers$gene) &
  (dvelo$`GABA interneurons_corr` > 0) 
]

# std normalize the correlation for each gene
dyn_gaba_subset$`GABA interneurons_corr` <- (
  dyn_gaba_subset$`GABA interneurons_corr` / sd(dyn$`GABA interneurons_corr`)
)
dvelo_gaba_subset$`GABA interneurons_corr` <- (
  dvelo_gaba_subset$`GABA interneurons_corr` / sd(dvelo$`GABA interneurons_corr`)
)

dyn_corrs <- data.frame(
  "Subset" = "GABAergic",
  "Method" = "scVelo",
  "Corr" = dyn_gaba_subset$`GABA interneurons_corr`  
)
dvelo_corrs <- data.frame(
  "Subset" = "GABAergic",
  "Method" = "DeepVelo",
  "Corr" = dvelo_gaba_subset$`GABA interneurons_corr`
)

dvelo_dyn_corrs_melted <- rbind(
  dyn_corrs,
  dvelo_corrs
)

kst_gaba <- ks.test(
  dvelo_corrs$Corr, 
  dyn_corrs$Corr
)
kst_gaba_df <- data.frame(
  "group1" = "DeepVelo",
  "group2" = "scVelo",
  "p-value" = kst_gaba$p.value,
  "subset" = "GABAergic"
)

ggplot(dvelo_dyn_corrs_melted, aes(x = Corr)) +
  geom_density(aes(
    fill = factor(Method, levels = c("scVelo", "DeepVelo"))), 
    stat = "density", 
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +
  theme_few() +
  labs(
    x = "GABAergic driver gene correlation", 
    y = "Density", 
    fill = "Method"
  ) + 
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 20)) +
  theme(aspect.ratio = 1)
ggsave(
  "../outs/figures/01_deepvelo_gaba_subset_all_hvgs_corr_comparison.pdf",
  height = 7, 
  width = 7
)
fwrite(
  kst_gaba_df,
  "../outs/01_deepvelo_gaba_subset_all_hvgs_corr_ks_test.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# Subset genes from glio markers for each and plot correlations 
dyn_glio_subset <- dyn[
  (dyn$gene %in% glio_markers$gene) &
  (dyn$`Gliogenic progenitors_corr` > 0)
]
dvelo_glio_subset <- dvelo[
  (dvelo$gene %in% glio_markers$gene) &
  (dvelo$`Gliogenic progenitors_corr` > 0)
]

# std normalize the correlation for each gene
dyn_glio_subset$`Gliogenic progenitors_corr` <- (
  dyn_glio_subset$`Gliogenic progenitors_corr` / sd(dyn$`Gliogenic progenitors_corr`)
)
dvelo_glio_subset$`Gliogenic progenitors_corr` <- (
  dvelo_glio_subset$`Gliogenic progenitors_corr` / sd(dvelo$`Gliogenic progenitors_corr`)
)

dyn_corrs <- data.frame(
  "Subset" = "Gliogenic",
  "Method" = "scVelo",
  "Corr" = dyn_glio_subset$`Gliogenic progenitors_corr`
)
dvelo_corrs <- data.frame(
  "Subset" = "Gliogenic",
  "Method" = "DeepVelo",
  "Corr" = dvelo_glio_subset$`Gliogenic progenitors_corr`
)

dvelo_dyn_corrs_melted <- rbind(
  dyn_corrs,
  dvelo_corrs
)

kst_glio <- ks.test(
  dvelo_corrs$Corr, 
  dyn_corrs$Corr
)
kst_glio_df <- data.frame(
  "group1" = "DeepVelo",
  "group2" = "scVelo",
  "p-value" = kst_glio$p.value,
  "subset" = "Gliogenic"
)

ggplot(dvelo_dyn_corrs_melted, aes(x = Corr)) +
  geom_density(
    aes(fill = factor(Method, levels = c("scVelo", "DeepVelo"))), 
    stat = "density", 
    alpha = 0.5
  ) +
  theme_few() +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +
  labs(
    x = "Gliogenic driver gene correlation", 
    y = "Density",
    fill = "Method"
  ) + 
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 20)) +
  theme(aspect.ratio = 1)
ggsave(
  "../outs/figures/01_deepvelo_glio_subset_all_hvgs_corr_comparison.pdf",
  height = 7, 
  width = 7
)
fwrite(
  kst_glio_df,
  "../outs/01_deepvelo_glio_subset_all_hvgs_corr_ks_test.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# Create ggplot of driver/marker gene combinations for lineages subset
# by method 
glio_table_melted <- glio_table
glio_table_melted$Subset <- rownames(glio_table_melted)
glio_table_melted <- reshape2::melt(glio_table_melted, id.vars = "Subset")
colnames(glio_table_melted) <- c("Subset", "Method", "Overlap")
glio_table_melted$Lineage <- "Gliogenic"

gaba_table_melted <- gaba_table
gaba_table_melted$Subset <- rownames(gaba_table_melted)
gaba_table_melted <- reshape2::melt(gaba_table_melted, id.vars = "Subset")
colnames(gaba_table_melted) <- c("Subset", "Method", "Overlap")
gaba_table_melted$Lineage <- "GABAergic"

glio_gaba_table <- rbind(glio_table_melted, gaba_table_melted)
glio_gaba_table_markers <- glio_gaba_table[
  glio_gaba_table$Subset %in% c("Markers"),
]

ggplot(data = glio_gaba_table_markers, aes(
    x = Lineage, 
    y = Overlap,
    fill = factor(Method, levels = c("scVelo", "DeepVelo"))
  )) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +
  theme_few() +
  labs(
    fill = "Method",
    x = "Lineage",
    y = "Marker overlap in top 100 driver genes"
  ) +
  geom_text(
    aes(Lineage, label = Overlap), 
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
  "../outs/figures/01_deepvelo_glio_gaba_marker_overlap_comparison.pdf",
  height = 7, 
  width = 7
)

### Comparison of DeepVelo and scVelo marker overlap based on rankings ###

# Sort dyn (scVelo) and DeepVelo dataframes by the specific lineage 
dyn_glio_sorted <- dyn[order(dyn$`Gliogenic progenitors_corr`, decreasing = TRUE)]
dyn_gaba_sorted <- dyn[order(dyn$`GABA interneurons_corr`, decreasing = TRUE)]
dvelo_glio_sorted <- dvelo[order(dvelo$`Gliogenic progenitors_corr`, decreasing = TRUE)]
dvelo_gaba_sorted <- dvelo[order(dvelo$`GABA interneurons_corr`, decreasing = TRUE)]

# Get the order of marker gene overlap in the driver gene estimations, for
# both methods and both lineages
dyn_glio_ranking <- which(dyn_glio_sorted$gene %in% glio_markers$gene)
dyn_gaba_ranking <- which(dyn_gaba_sorted$gene %in% gaba_markers$gene)

dvelo_glio_ranking <- which(dvelo_glio_sorted$gene %in% glio_markers$gene)
dvelo_gaba_ranking <- which(dvelo_gaba_sorted$gene %in% gaba_markers$gene)

# Plot and visualize and the differences between the GABAergic rankings 
gaba_ranking_dyn <- data.frame(
  "Method" = "scVelo",
  "Lineage" = "GABAergic",
  "Ranking" = dyn_gaba_ranking
)
gaba_ranking_dvelo <- data.frame(
  "Method" = "DeepVelo",
  "Lineage" = "GABAergic",
  "Ranking" = dvelo_gaba_ranking
)
gaba_ranking_combined <- rbind(
  gaba_ranking_dyn,
  gaba_ranking_dvelo
)

ggplot(data = gaba_ranking_combined, aes(x = Ranking)) + 
  geom_density(aes(
    fill = factor(Method, levels = c("scVelo", "DeepVelo"))), 
    stat = "density", 
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +
  theme_few() +
  labs(
    x = "GABAergic marker gene ranking", 
    y = "Density", 
    fill = "Method"
  ) + 
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 20)) +
  theme(aspect.ratio = 1)
ggsave(
  "../outs/figures/01_deepvelo_gaba_subset_all_hvgs_marker_ranking_comparison.pdf",
  height = 7, 
  width = 7
)

# Plot and visualize and the differences between the Gliogenic rankings 
glio_ranking_dyn <- data.frame(
  "Method" = "scVelo",
  "Lineage" = "Gliogenic",
  "Ranking" = dyn_glio_ranking
)
glio_ranking_dvelo <- data.frame(
  "Method" = "DeepVelo",
  "Lineage" = "Gliogenic",
  "Ranking" = dvelo_glio_ranking
)
glio_ranking_combined <- rbind(
  glio_ranking_dyn,
  glio_ranking_dvelo
)

ggplot(data = glio_ranking_combined, aes(x = Ranking)) + 
  geom_density(aes(
    fill = factor(Method, levels = c("scVelo", "DeepVelo"))), 
    stat = "density", 
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +
  theme_few() +
  labs(
    x = "Gliogenic marker gene ranking", 
    y = "Density", 
    fill = "Method"
  ) + 
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 20)) +
  theme(aspect.ratio = 1)
ggsave(
  "../outs/figures/01_deepvelo_glio_subset_all_hvgs_marker_ranking_comparison.pdf",
  height = 7, 
  width = 7
)

# Perform Wilcoxon rank-sum tests for testing differences in ranking
# between the scVelo and DeepVelo methods, for both the GABAergic 
# and gliogenic lineages - save the results 
wrst_gaba <- wilcox.test(dvelo_gaba_ranking, dyn_gaba_ranking)
wrst_glio <- wilcox.test(dvelo_glio_ranking, dyn_glio_ranking)

wrst_glio_df <- data.frame(
  "group1" = "DeepVelo",
  "group2" = "scVelo",
  "p-value" = wrst_glio$p.value,
  "subset" = "Gliogenic"
)
wrst_gaba_df <- data.frame(
  "group1" = "DeepVelo",
  "group2" = "scVelo",
  "p-value" = wrst_gaba$p.value,
  "subset" = "GABAergic"
)

fwrite(
  wrst_glio_df,
  "../outs/01_deepvelo_glio_subset_all_hvgs_marker_ranking_wrs_test.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)
fwrite(
  wrst_gaba_df,
  "../outs/01_deepvelo_gaba_subset_all_hvgs_marker_ranking_wrs_test.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

