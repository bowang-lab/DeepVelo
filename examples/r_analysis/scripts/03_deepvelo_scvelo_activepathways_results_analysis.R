library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(stringr)
library(ggpubr)
library(ggvenn)
library(VennDiagram)
library(GSA)
library(ActivePathways)

######################## DeepVelo pathway analysis ########################

## Load result of ActivePathways list for all 100, 250, 500 GABA/Glio genes
ap_result_all <- fread(
  "../outs/02_activepathways_dvelo_glio_gaba_top_100_250_500_all_outsubgroups.txt"
)
ap_result_all_df <- as.data.frame(ap_result_all)

# Get list of set intersections 
col_list = list(
  "Gliogenic 100",
  "Gliogenic 250",
  "GABAergic 100",
  "GABAergic 250"
)
set_lists <- lapply(col_list, function(x) {
  ap_result_sub <- ap_result_all_df[which(ap_result_all_df[x] == 1), ]
  ap_result_sub_values <- ap_result_sub$term.id
  return(ap_result_sub_values)
})

# Get Venn Diagram of intersection between sets
names(set_lists) <- col_list
ggvenn(
  set_lists,
  set_name_size = 4,
)

# This doesn't really seem to be working - non-significant unique pathways
# for the gabaergic/glutamatergic lineages

## Load result of Deepvelo ActivePathways for GABA/Glio top 100 corr subsets
dvelo_gaba_top_100 <- fread(
  "../outs/01_deepvelo_top_100_gaba_drivers_corr.tsv"
)
dvelo_glio_top_100 <- fread(
  "../outs/01_deepvelo_top_100_glio_drivers_corr.tsv"
)
colnames(dvelo_gaba_top_100)[1] <- "gene"
colnames(dvelo_glio_top_100)[1] <- "gene"

dvelo_gaba_ap <- fread(
  "../outs/02_activepathways_dvelo_gaba_top_100_corr_outpathways.txt"
)
dvelo_glio_ap <- fread(
  "../outs/02_activepathways_dvelo_glio_top_100_corr_outpathways.txt"
)
dvelo_gaba_gmt <- GSA.read.gmt(
  "../outs/02_activepathways_dvelo_gaba_top_100_corr_outpathways.gmt"
)
dvelo_glio_gmt <- GSA.read.gmt(
  "../outs/02_activepathways_dvelo_glio_top_100_corr_outpathways.gmt"
)

# Get set and intersection sizes using venn diagrams 
gaba_enriched_dvelo <- dvelo_gaba_ap$term.name
glio_enriched_dvelo <- dvelo_glio_ap$term.name
enriched_list <- list(
  "DeepVelo \nGABAergic \nenriched \npathways" = gaba_enriched_dvelo,
  "DeepVelo \nGliogenic \nenriched \npathways" = glio_enriched_dvelo
)
ggvenn(
  enriched_list,
  fill_color = c("#E66100", "#5D3A9B"),
  set_name_size = 7,
  text_size = 7
)
ggsave(
  "../outs/figures/03_dvelo_top_100_gaba_glio_pathway_venn.pdf",
  height = 6.2,
  width = 6.2
)

# Do regex matching for relevant catagories in both subsets
dvelo_gaba_ap$Lineage <- "GABAergic"
dvelo_gaba_ap$rel <- ifelse(
  str_detect(dvelo_gaba_ap$term.name, regex("neur", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("brain", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("axon", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("cerebell", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("synap", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("nerv", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("gaba", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("myelin", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("metenceph", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("dendrite", ignore_case = TRUE)) |
  str_detect(dvelo_gaba_ap$term.name, regex("projection", ignore_case = TRUE)),
  "Neurogenesis", ifelse(
    str_detect(dvelo_gaba_ap$term.name, regex("devel", ignore_case = TRUE)) |
    str_detect(dvelo_gaba_ap$term.name, regex("growth", ignore_case = TRUE)) |
    str_detect(dvelo_gaba_ap$term.name, regex("prolifer", ignore_case = TRUE)) |
    str_detect(dvelo_gaba_ap$term.name, regex("growth", ignore_case = TRUE)) |
    str_detect(dvelo_gaba_ap$term.name, regex("morpho", ignore_case = TRUE)) |
    str_detect(dvelo_gaba_ap$term.name, regex("differen", ignore_case = TRUE)),
    "Developmental non-neuronal",
    "Non-specific"
  )
)

dvelo_glio_ap$Lineage <- "Gliogenic"
dvelo_glio_ap$rel <- ifelse(
  str_detect(dvelo_glio_ap$term.name, regex("neur", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("brain", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("axon", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("cerebell", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("synap", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("nerv", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("gli", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("metenceph", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("astro", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("oligo", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("myelin", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("dendrite", ignore_case = TRUE)) |
  str_detect(dvelo_glio_ap$term.name, regex("projection", ignore_case = TRUE)),
  "Neurogenesis", ifelse(
    str_detect(dvelo_glio_ap$term.name, regex("devel", ignore_case = TRUE)) |
      str_detect(dvelo_glio_ap$term.name, regex("growth", ignore_case = TRUE)) |
      str_detect(dvelo_glio_ap$term.name, regex("prolifer", ignore_case = TRUE)) |
      str_detect(dvelo_glio_ap$term.name, regex("growth", ignore_case = TRUE)) |
      str_detect(dvelo_glio_ap$term.name, regex("morpho", ignore_case = TRUE)) |
      str_detect(dvelo_glio_ap$term.name, regex("differen", ignore_case = TRUE)),
    "Developmental non-neuronal",
    "Non-specific"
  )
)

# Plot results by percentage 
glio_gaba_ap <- rbind(dvelo_glio_ap, dvelo_gaba_ap)
ggplot(data = glio_gaba_ap, aes(x = Lineage, fill = rel)) +
  theme_few() +
  geom_bar(stat = "count", position = "fill") + 
  scale_fill_brewer(palette = "Dark2") + 
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18, face = "bold")) +
  theme(axis.title.y = element_text(size = 18, face = "bold")) +
  theme(legend.title = element_text(size = 18, face = "bold")) +
  theme(legend.text = element_text(size = 18)) +
  theme(aspect.ratio = 1) +
  labs(
    fill = "Pathway type",
    x = "Lineage",
    y = "Percentage of DeepVelo driver \nenriched pathways"
) 
ggsave(
  "../outs/figures/03_dvelo_top_100_pathway_neurogenesis_enrichment.pdf",
  height = 8,
  width = 8
)  
  
# Subset for top 20 pathways for each lineage (based on p-value)
dvelo_glio_ap_top_20 <- dvelo_glio_ap[order(dvelo_glio_ap$adjusted.p.val)][1:20]
dvelo_gaba_ap_top_20 <- dvelo_gaba_ap[order(dvelo_gaba_ap$adjusted.p.val)][1:20]

# Merge with geneset sizes + add overlap information
glio_genesets <- data.frame(
  "term.id" = dvelo_glio_gmt$geneset.names,
  "term.name" = dvelo_glio_gmt$geneset.descriptions,
  "set_size" = sapply(dvelo_glio_gmt$genesets, length)
)
glio_genesets$enrich_pct <- sapply(
  dvelo_glio_gmt$genesets, function(x) {
    length(which(dvelo_glio_top_100$gene %in% x))/length(dvelo_glio_top_100$gene)
  }
)

gaba_genesets <- data.frame(
  "term.id" = dvelo_gaba_gmt$geneset.names,
  "term.name" = dvelo_gaba_gmt$geneset.descriptions,
  "set_size" = sapply(dvelo_gaba_gmt$genesets, length)
)
gaba_genesets$enrich_pct <- sapply(
  dvelo_gaba_gmt$genesets, function(x) {
    length(which(dvelo_gaba_top_100$gene %in% x))/length(dvelo_gaba_top_100$gene)
  }
)

dvelo_glio_ap_top_20 <- merge(
  dvelo_glio_ap_top_20,
  glio_genesets,
  by = "term.name",
  join = "inner"
)
dvelo_gaba_ap_top_20 <- merge(
  dvelo_gaba_ap_top_20,
  gaba_genesets,
  by = "term.name",
  join = "inner"
)

# Plot pathway enrichment results for top 20 in each lineage
ggplot(data = dvelo_gaba_ap_top_20, aes(x = enrich_pct, y = term.name)) +
  theme_few() +
  geom_point(aes(color = rel, size = set_size), stat = "identity") +
  scale_color_brewer(palette = "Dark2") + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 18, face = "bold")) +
  theme(axis.title.y = element_text(size = 18, face = "bold")) +
  theme(legend.title = element_text(size = 18, face = "bold")) +
  theme(legend.text = element_text(size = 18)) +
  theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_size(
    breaks = c(500, 750, 1000, 1500, 2000)
  ) +
  labs(
    size = "Gene set size",
    color = "Pathway type",
    x = "Enrichment Percentage",
    y = "GABAergic pathways"
  ) +
  theme(aspect.ratio = 1)
ggsave(
  "../outs/figures/03_dvelo_top_20_gabaergic_pathways_full_plot.pdf",
  height = 13,
  width = 13
)

ggplot(data = dvelo_glio_ap_top_20, aes(x = enrich_pct, y = term.name)) +
  theme_few() +
  geom_point(aes(color = rel, size = set_size), stat = "identity") +
  scale_color_brewer(palette = "Dark2") + 
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(size = 20, face = "bold")) +
  theme(legend.text = element_text(size = 20)) +
  theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_size(
    breaks = c(500, 750, 1000, 1500, 2000)
  ) +
  labs(
    size = "Gene set size",
    color = "Pathway type",
    x = "Enrichment Percentage",
    y = "Gliogenic pathways"
  ) +
  theme(aspect.ratio = 1)
ggsave(
  "../outs/figures/03_dvelo_top_20_gliogenic_pathways_full_plot.pdf",
  height = 12.5,
  width = 12.5
) 

# Get contingency table for deepvelo gaba and glio functional enrichment
dvelo_glio_table <- table(dvelo_glio_ap$rel)
dvelo_gaba_table <- table(dvelo_gaba_ap$rel)

######################## scVelo pathway analysis ########################

## Load result of scVelo ActivePathways for GABA/Glio top 100 corr subsets
scvelo_gaba_top_100 <- fread(
  "../outs/01_scvelo_top_100_gaba_drivers_corr.tsv"
)
scvelo_glio_top_100 <- fread(
  "../outs/01_scvelo_top_100_glio_drivers_corr.tsv"
)
colnames(scvelo_gaba_top_100)[1] <- "gene"
colnames(scvelo_glio_top_100)[1] <- "gene"

scvelo_gaba_ap <- fread(
  "../outs/02_activepathways_scvelo_gaba_top_100_corr_outpathways.txt"
)
scvelo_glio_ap <- fread(
  "../outs/02_activepathways_scvelo_glio_top_100_corr_outpathways.txt"
)
scvelo_gaba_gmt <- GSA.read.gmt(
  "../outs/02_activepathways_scvelo_gaba_top_100_corr_outpathways.gmt"
)
scvelo_glio_gmt <- GSA.read.gmt(
  "../outs/02_activepathways_scvelo_glio_top_100_corr_outpathways.gmt"
)

# Get set and intersection sizes using venn diagrams 
gaba_enriched_scvelo <- scvelo_gaba_ap$term.name
glio_enriched_scvelo <- scvelo_glio_ap$term.name
enriched_list <- list(
  "scVelo \nGABAergic \nenriched \npathways" = gaba_enriched_scvelo,
  "scVelo \nGliogenic \nenriched \npathways" = glio_enriched_scvelo
)
ggvenn(
  enriched_list,
  fill_color = c("#E66100", "#5D3A9B"),
  set_name_size = 7,
  text_size = 7
)
ggsave(
  "../outs/figures/03_scvelo_top_100_gaba_glio_pathway_venn.pdf",
  height = 6.2,
  width = 6.2
)

# Do regex matching for relevant catagories in both subsets
scvelo_gaba_ap$Lineage <- "GABAergic"
scvelo_gaba_ap$rel <- ifelse(
  str_detect(scvelo_gaba_ap$term.name, regex("neur", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("brain", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("axon", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("cerebell", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("synap", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("nerv", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("gaba", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("myelin", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("metenceph", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("dendrite", ignore_case = TRUE)) |
    str_detect(scvelo_gaba_ap$term.name, regex("projection", ignore_case = TRUE)),
  "Neurogenesis", ifelse(
    str_detect(scvelo_gaba_ap$term.name, regex("devel", ignore_case = TRUE)) |
      str_detect(scvelo_gaba_ap$term.name, regex("growth", ignore_case = TRUE)) |
      str_detect(scvelo_gaba_ap$term.name, regex("prolifer", ignore_case = TRUE)) |
      str_detect(scvelo_gaba_ap$term.name, regex("growth", ignore_case = TRUE)) |
      str_detect(scvelo_gaba_ap$term.name, regex("morpho", ignore_case = TRUE)) |
      str_detect(scvelo_gaba_ap$term.name, regex("differen", ignore_case = TRUE)),
    "Developmental non-neuronal",
    "Non-specific"
  )
)

scvelo_glio_ap$Lineage <- "Gliogenic"
scvelo_glio_ap$rel <- ifelse(
  str_detect(scvelo_glio_ap$term.name, regex("neur", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("brain", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("axon", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("cerebell", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("synap", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("nerv", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("gli", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("metenceph", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("astro", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("oligo", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("myelin", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("dendrite", ignore_case = TRUE)) |
    str_detect(scvelo_glio_ap$term.name, regex("projection", ignore_case = TRUE)),
  "Neurogenesis", ifelse(
    str_detect(scvelo_glio_ap$term.name, regex("devel", ignore_case = TRUE)) |
      str_detect(scvelo_glio_ap$term.name, regex("growth", ignore_case = TRUE)) |
      str_detect(scvelo_glio_ap$term.name, regex("prolifer", ignore_case = TRUE)) |
      str_detect(scvelo_glio_ap$term.name, regex("growth", ignore_case = TRUE)) |
      str_detect(scvelo_glio_ap$term.name, regex("morpho", ignore_case = TRUE)) |
      str_detect(scvelo_glio_ap$term.name, regex("differen", ignore_case = TRUE)),
    "Developmental non-neuronal",
    "Non-specific"
  )
)

# Plot results by percentage 
glio_gaba_ap <- rbind(scvelo_glio_ap, scvelo_gaba_ap)
ggplot(data = glio_gaba_ap, aes(x = Lineage, fill = rel)) +
  theme_few() +
  geom_bar(stat = "count", position = "fill") + 
  scale_fill_brewer(palette = "Dark2") + 
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18, face = "bold")) +
  theme(axis.title.y = element_text(size = 18, face = "bold")) +
  theme(legend.title = element_text(size = 18, face = "bold")) +
  theme(legend.text = element_text(size = 18)) +
  theme(aspect.ratio = 1) +
  labs(
    fill = "Pathway type",
    x = "Lineage",
    y = "Percentage of scVelo driver \nenriched pathways"
  ) 
ggsave(
  "../outs/figures/03_scvelo_top_100_pathway_neurogenesis_enrichment.pdf",
  height = 8,
  width = 8
)  

# Get contingency table for scvelo gaba and glio functional enrichment
scvelo_glio_table <- table(scvelo_glio_ap$rel)
scvelo_gaba_table <- table(scvelo_gaba_ap$rel)

# Create contingency table for scvelo and dvelo comparison and run
# fisher's exact test 
scvelo_dvelo_glio_table <- cbind(
  as.matrix(scvelo_glio_table), as.matrix(dvelo_glio_table)
)
colnames(scvelo_dvelo_glio_table) <- c("scVelo", "DeepVelo")
scvelo_dvelo_glio_fishers <- fisher.test(scvelo_dvelo_glio_table, y = NULL, 
                           workspace = 200000, hybrid = FALSE,
                           hybridPars = c(expect = 5, percent = 80, Emin = 1),
                           control = list(), or = 1, alternative = "two.sided",
                           conf.int = TRUE, conf.level = 0.95,
                           simulate.p.value = FALSE, B = 2000)
scvelo_dvelo_glio_fishers

scvelo_dvelo_gaba_table <- cbind(
  as.matrix(scvelo_gaba_table), as.matrix(dvelo_gaba_table)
)
colnames(scvelo_dvelo_gaba_table) <- c("scVelo", "DeepVelo")
scvelo_dvelo_gaba_fishers <- fisher.test(scvelo_dvelo_gaba_table, y = NULL, 
                                         workspace = 200000, hybrid = FALSE,
                                         hybridPars = c(expect = 5, percent = 80, Emin = 1),
                                         control = list(), or = 1, alternative = "two.sided",
                                         conf.int = TRUE, conf.level = 0.95,
                                         simulate.p.value = FALSE, B = 2000)
scvelo_dvelo_gaba_fishers

### GABA significant at 2.919e-06, Glio non-siginficant

### Supplementary figures ###

## Get VennDiagram of scVelo GABA/Glio and DeepVelo GABA/Glio predicted
## pathways - jointly
enriched_list_all <- list(
  "scVelo \nGABAergic \nenriched \npathways" = gaba_enriched_scvelo,
  "DeepVelo \nGABAergic \nenriched \npathways" = gaba_enriched_dvelo,
  "scVelo \nGliogenic \nenriched \npathways" = glio_enriched_scvelo,
  "DeepVelo \nGliogenic \nenriched \npathways" = glio_enriched_dvelo
)
venn.diagram(
  enriched_list_all,
  filename = "../outs/figures/03_dvelo_scvelo_top_100_gaba_glio_path_venn.tiff",
  cat.just = list(
    c(0.5, 0.2),
    c(0.5, 0.2),
    c(0.5, 0.2),
    c(0.5, 0.2)
  ),
  fill = c(
    "#0C7BDC",
    "#FFC20A", 
    "#0C7BDC",
    "#FFC20A"
  ),
  height = 3000, 
  width = 3000, 
  resolution = 500, 
  imagetype = "tiff"
)
