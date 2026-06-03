library(dplyr)
library(Seurat)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_dir <- args[2]
# creates seurat obj
pbmc.data <- Read10X(data.dir = data_dir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
saveRDS(pbmc, file = file.path(out_dir, "pbmc_seurat_obj.rds"))
