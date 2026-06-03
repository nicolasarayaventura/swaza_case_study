library(Seurat)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
out_dir <- args[2]

pbmc <- readRDS(file.path(data_dir, "pbmc_seurat_obj.rds"))

# add mitochondrial % to metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# save violin plot
png(file.path(out_dir, "vlnplot_qc.png"), width = 1200, height = 400)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# save scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(file.path(out_dir, "scatter_qc.png"), width = 1000, height = 500)
plot1 + plot2
dev.off()

# filter cells between 200 and 2500 genes
# more than 200 probable droplets or dead cells
# less than 2500 probable doublets that have 2 cells in 1 droplet
# less than 5% are for dying/ damaged cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) #normalization with default paramaters for precaution

# find top 2000 variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# plot top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(file.path(out_dir, "variable_features.png"), width = 1200, height = 600)
plot1 + plot2
dev.off()

# scale data across all genes
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# save
saveRDS(pbmc, file = file.path(out_dir, "pbmc_filtered_qc.rds"))