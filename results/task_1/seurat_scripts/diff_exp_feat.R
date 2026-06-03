library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]

pbmc <- readRDS(file.path(out_dir, "pbmc_clustered.rds"))

# find markers for every cluster compared to all remaining cells
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

# violin plot of known marker genes
png(file.path(out_dir, "markers_vlnplot.png"), width = 1200, height = 600)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
dev.off()

# raw counts violin plot
png(file.path(out_dir, "markers_vlnplot_raw.png"), width = 1200, height = 600)
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
dev.off()

# feature plot of canonical markers on UMAP
png(file.path(out_dir, "markers_featureplot.png"), width = 1200, height = 1000)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
dev.off()

# heatmap of top 10 markers per cluster
top10 <- pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup()
png(file.path(out_dir, "markers_heatmap.png"), width = 1200, height = 1000)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()

# assign cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# final UMAP with cell type labels
png(file.path(out_dir, "umap_celltypes.png"), width = 1200, height = 700)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
    guides(colour = guide_legend(override.aes = list(size = 10)))
dev.off()

# save final object
saveRDS(pbmc, file = file.path(out_dir, "pbmc_final.rds"))