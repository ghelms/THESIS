library(tidyverse)
library(Seurat)

# IMPORT DISCO
disco = readRDS("../disco_data/disco_combined_normalized.rds")
disco <- JoinLayers(disco)

disco@meta.data <- disco@meta.data %>% 
  mutate(tissue_cell = paste(tissue, cell_type, sep= "_")) %>% 
  group_by(tissue_cell) %>% 
  slice_head(n = 1000) %>% 
  ungroup()


# Raw count matrix (genes × cells)
X <- GetAssayData(disco[["RNA"]], slot = "counts")

# Cell metadata (cells × variables)
meta <- disco@meta.data

library(Matrix)

writeMM(X, "disco_counts.mtx")
write.csv(meta_cluster, "disco_meta.csv")
write.csv(data.frame(gene = rownames(X)), "disco_var.csv", row.names = FALSE)
