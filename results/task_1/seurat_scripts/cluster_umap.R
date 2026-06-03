library(Seurat)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]

pbmc <- readRDS(file.path(out_dir, "pbmc_pca.rds"))

# build KNN graph using first 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# cluster cells with Louvain algorithm, resolution 0.4-1.2 good for ~3k cells
pbmc <- FindClusters(pbmc, resolution = 0.7)

# print cluster IDs of first 5 cells for debug
# print(head(Idents(pbmc), 5))

# run UMAP using same 10 PCs as clustering
pbmc <- RunUMAP(pbmc, dims = 1:10)

# plot UMAP with cluster labels
png(file.path(out_dir, "umap.png"), width = 800, height = 600)
DimPlot(pbmc, reduction = "umap", label = TRUE)
dev.off()

# save
saveRDS(pbmc, file = file.path(out_dir, "pbmc_clustered.rds"))
