# Swaza_case_study

## Task 1
### Introduction & Background
  For this small case study I set up a simple mock repo as if I were working on a normal project under a laboratory. I set up a mock scratch directory under the ./scratch directory for any files that would utilize too much storage under a work directory to save storage and be able to save space on an HPC. I kept files organized and worked with shell to maintain an easy to read and friendly debugging method for each step during each task for reproducibility and tracking steps. I will elaborate further in my methods for packages used and why certain parameters were used.
### Methods
#### Pipeline Organization
The pipeline was written in bash (run.sh) and structured as modular functions, each calling a dedicated R script. This design keeps each step isolated, easier to debug, and reproducible. The script was run with set -ex to print each command before execution and exit immediately on any error. A shared output directory under ../../scratch/results/task_1/seurat_scripts/ was used across all steps to pass intermediate .rds objects between scripts for any debugging issues as well as easy submissions in slurm. I also used a dockerfile to have a safe enviroment of reproduicble experiment.
##### Why Seurat?
I used Seurat because of familiarity of use in previous research settings with well documentated pipeline and explanations on parameters as well as acceseible small training data set. 
#### Data
Raw 10x Genomics PBMC data was downloaded using wget via a dataset link file (dataset_link.txt), which contained the dataset label and download URL. The data was saved to a structured scratch directory to avoid storing large files in the main working directory.
#### Seurat Object Creation (makeobj.R)
The raw 10x data was loaded using Read10X() from the Seurat package and converted into a Seurat object using CreateSeuratObject() with a minimum cell threshold of 3 and a minimum feature threshold of 200. The object was saved as pbmc_seurat_obj.rds for use in downstream steps.
#### Quality Control & Normalization (qc_norm_plot.R)
Mitochondrial gene percentage was calculated using PercentageFeatureSet() with the pattern ^MT- and added to cell metadata. Cells were filtered using the following thresholds: fewer than 200 unique genes (likely empty droplets or dead cells), more than 2,500 unique genes (likely doublets), and greater than 5% mitochondrial reads (likely dying or damaged cells).
<img width="1000" height="500" alt="scatter_qc" src="https://github.com/user-attachments/assets/a2c49670-a977-4790-bc63-8436b5a8e66b" />

#### QC metrics before filtering (qc.R)
Following filtering, data was normalized using NormalizeData() with the LogNormalize method and a scale factor of 10,000, which scales each cell to a total count of 10,000 before log-transforming. The top 2,000 highly variable features were identified using FindVariableFeatures() with the VST selection method, and all genes were scaled using ScaleData() to ensure equal contribution in downstream dimensionality reduction.
PCA was performed on the scaled data using RunPCA() with the previously identified variable features as input. Gene loadings were visualized using VizDimLoadings() and DimHeatmap() across the first 15 PCs. An elbow plot was generated using ElbowPlot() to determine dataset dimensionality, with an elbow observed around PC 9-10, suggesting the first 10 PCs capture the majority of biological signal.
<img width="1200" height="600" alt="pca_loadings" src="https://github.com/user-attachments/assets/ce7c1ab8-6877-4042-8bcb-ad345ab6cd0a" />

#### Clustering & UMAP (cluster_umap.R)
A K-nearest neighbor graph was constructed in PCA space using FindNeighbors() with the first 10 PCs. Cells were clustered using FindClusters() with the Louvain algorithm at a resolution of 0.7, which is within the recommended 0.4-1.2 range for datasets of approximately 3,000 cells, producing 10 clusters.
##### Marker Gene Identification & Cell Type Annotation (diff_exp_feat.R)
Differentially expressed marker genes were identified for each cluster using FindAllMarkers(), retaining only positive markers with a log2 fold change greater than 1. The top 10 markers per cluster were visualized in a heatmap using DoHeatmap(). A FeaturePlot() of canonical marker genes was generated to visualize expression across clusters on the UMAP. Cell type identities were assigned to 10 clusters using RenameIdents(), including Naive CD4+ T cells, CD14+ Monocytes, Memory CD4+ T cells, B cells, CD8+ T cells, FCGR3A+ Monocytes, NK cells, Dendritic Cells, Platelets, and one cluster labeled Unknown due to futher resolution being used at 0.7 as median between given recommendation in documentation.
<img width="1200" height="700" alt="umap_celltypes" src="https://github.com/user-attachments/assets/51c13365-0139-44b9-b8d0-6ebafd5ed611" />

#### Two-Cluster Comparison: Naive CD4+ T vs Memory CD4+ T (cluster_umap.R cont.)
To directly compare two biologically related clusters, FindMarkers() was used to perform differential expression between cluster 0 (Naive CD4+ T) and cluster 2 (Memory CD4+ T). This comparison is biologically meaningful as Naive and Memory CD4+ T cells share lineage but differ in activation and differentiation state, making their transcriptional differences informative for understanding T cell heterogeneity in the PBMC dataset.
<img width="1000" height="500" alt="cluster0_vs_2_dotplot" src="https://github.com/user-attachments/assets/1179f68a-b721-4c05-ba4e-a8970d2408bf" />
## Task 2
### Introduction & Background
For Task 2 I performed a bulk RNA-seq differential expression analysis using the Bioconductor airway dataset (GEO: GSE52778), which contains RNA-seq counts from dexamethasone-treated and untreated human airway smooth muscle cells across four cell lines. This is a well-characterized two-group design making it ideal for demonstrating a standard DE analysis workflow. Similar to Task 1, I used a dedicated Dockerfile to ensure a reproducible environment.
### Methods
#### Pipeline Organization
The pipeline was written in bash (run.sh) and structured as a single function calling deseq2.R. The script was run with set -ex for debugging and reproducibility. A shared output directory under results/task_2/ was used for all outputs. A separate Dockerfile was built using the bioconductor/bioconductor_docker base image, which is more appropriate for DESeq2 work than the Seurat image used in Task 1.
#### Why DESeq2?
DESeq2 was chosen for its well-documented and widely adopted negative binomial model for count data, robust normalization procedure, and accessible Bioconductor documentation. It is the standard tool for bulk RNA-seq differential expression analysis in the field.
#### Data
The airway dataset was loaded directly from Bioconductor using data(airway), requiring no external download. It contains 8 samples across 4 cell lines (SRR1039508, SRR1039509, SRR1039512, SRR1039513, SRR1039516, SRR1039517, SRR1039520, SRR1039521) with two conditions: dexamethasone treated (trt) and untreated (untrt).
#### Normalization & Design Model (deseq2.R)
Low count genes were removed by filtering out genes with fewer than 10 total counts across all samples. DESeq2's median-of-ratios normalization was applied automatically via DESeq(), which corrects for both library size and composition effects. The design model ~ cell + dex was used to control for cell line as a confounder while testing the effect of dexamethasone treatment. Including cell in the design accounts for the fact that the four cell lines may have baseline expression differences unrelated to treatment.
<img width="800" height="600" alt="ma_plot" src="https://github.com/user-attachments/assets/2464ad7c-4713-4df1-a6a2-a73723f7e723" />

#### Differential Expression
Results were extracted using results() with the contrast trt vs untrt, ranked by log2FoldChange. Genes with an adjusted p-value below 0.05 and absolute log2 fold change greater than 1 were considered significant.
#### Diagnostic Plots
A PCA plot of variance-stabilized counts was generated using plotPCA() to assess sample quality and batch structure. A volcano plot was produced to visualize the overall distribution of fold changes and significance. An MA plot was generated using plotMA() to assess the relationship between mean expression and fold change across all genes.
<img width="800" height="600" alt="volcano_plot" src="https://github.com/user-attachments/assets/6f269857-36bf-499b-82ba-0f2a34f4bed5" />
<img width="800" height="600" alt="pca_samples" src="https://github.com/user-attachments/assets/2947e891-7ec9-4223-be34-561b43fa01e6" />

## Results
### Task 1
The Seurat pipeline successfully processed the PBMC 3k dataset from raw 10x Genomics count matrices through to annotated cell type clusters. QC filtering removed low quality cells, empty droplets, and doublets, with the scatter plot confirming expected relationships between total counts and gene counts with no major outliers. PCA loadings showed biologically meaningful signal, with PC1 separating myeloid from lymphoid lineages and the elbow plot confirming the first 10 PCs as the dimensionality cutoff. Louvain clustering at resolution 0.7 produced 10 distinct clusters with clear separation on the UMAP. Canonical marker genes confirmed 9 known PBMC cell types with one cluster labeled Unknown, likely reflecting a rare or transitional population warranting further investigation. Direct comparison of Naive CD4+ T vs Memory CD4+ T cells using FindMarkers() revealed expected transcriptional differences consistent with known biology.
### Task 2
The DESeq2 pipeline successfully identified differentially expressed genes between dexamethasone-treated and untreated airway smooth muscle cells. The PCA plot showed clear separation between treated and untreated samples along PC1, confirming the treatment effect dominates the variance. The cell line effect was visible along PC2, validating the decision to include cell in the design model as a confounder. The volcano and MA plots showed a well-distributed set of significant genes, indicating reliable DE results with no major quality issues.

