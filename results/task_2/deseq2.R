library(DESeq2)
library(airway)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]

# load airway dataset
data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# filter low count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# run DESeq2 - handles normalization of library size and composition
dds <- DESeq(dds)

# get results treated vs untreated
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# ranked DE table by log2FoldChange
res_ordered <- as.data.frame(res[order(-res$log2FoldChange), ])
write.csv(res_ordered, file = file.path(out_dir, "deseq2_results.csv"))

# PCA plot
vsd <- vst(dds, blind = FALSE)
png(file.path(out_dir, "pca_samples.png"), width = 800, height = 600)
print(plotPCA(vsd, intgroup = c("dex", "cell")))
dev.off()

# volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

png(file.path(out_dir, "volcano_plot.png"), width = 800, height = 600)
print(ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot: Treated vs Untreated"))
dev.off()

# MA plot
png(file.path(out_dir, "ma_plot.png"), width = 800, height = 600)
plotMA(res, main = "MA Plot: Treated vs Untreated")
dev.off()