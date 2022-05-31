library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggpubr)

# Load marker list from Vladiou et al.
vla_markers <- fread("../data/41586_2019_1158_MOESM4_ESM.csv")
vla_markers <- vla_markers[, -1]

# Fill down values
vla_markers$Annotations <- na_if(vla_markers$Annotations, "")
vla_markers <- vla_markers %>% tidyr::fill(Annotations, .direction = "down")

# Subset and save for GABAergic and gliogenic markers
gaba_glio_markers <- vla_markers[vla_markers$Annotations %in% c(
  "GABA interneurons",
  "Differentiating GABA interneurons",
  "Gliogenic progenitors"
)]
fwrite(
  gaba_glio_markers,
  "../outs/05_marker_gaba_glio_subset.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

