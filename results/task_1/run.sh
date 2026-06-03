set -x -e
scratch="../../scratch"
data_10x="${scratch}/data/task_1/pbmc_10x"
output="${scratch}/task_1/results/seurat_scripts"

function seurat_object {
    mkdir -p "${scratch}/task_1/results/seurat_scripts"
    Rscript analysis.R "${data_10x}" "${output}"
}
function qc_norm {
    mkdir -p "${output}"
    Rscript qc_norm_plot.R "${data_10x}" "${output}"
}
function pca {
    mkdir -p "${output}"
    Rscript pca.R "${output}" "${output}"
}

function cluster_umap {
    mkdir -p "${output}"
    Rscript cluster_umap.R "${output}"
}
function diff_exp_feat {
    mkdir -p "${output}"
    Rscript diff_exp_feat.R "${output}"
}
seurat_object
qc_norm
pca
cluster_umap
diff_exp_feat