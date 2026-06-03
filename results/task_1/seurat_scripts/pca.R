library(Seurat)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_dir <- args[2]

pbmc <- readRDS(file.path(data_dir, "pbmc_filtered_qc.rds"))

# run PCA on variable features
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# print top 5 genes for first 5 PCs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# visualize gene loadings for PC1 and PC2
png(file.path(out_dir, "pca_loadings.png"), width = 1200, height = 600)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
dev.off()

# plot cells in PCA space
png(file.path(out_dir, "pca_dimplot.png"), width = 800, height = 600)
DimPlot(pbmc, reduction = "pca") + NoLegend()
dev.off()

# heatmap for PC1
png(file.path(out_dir, "pca_heatmap_pc1.png"), width = 800, height = 600)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# heatmap for PC1-15
png(file.path(out_dir, "pca_heatmap_pc1_15.png"), width = 1200, height = 1200)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# elbow plot to determine dimensionality
png(file.path(out_dir, "elbow_plot.png"), width = 800, height = 600)
ElbowPlot(pbmc)
dev.off()

# save
saveRDS(pbmc, file = file.path(out_dir, "pbmc_pca.rds"))