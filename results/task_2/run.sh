set -x -e
output="/swaza_case_study/results/task_2/graphs"

function deseq2 {
    mkdir -p "${output}"
    Rscript deseq2.R "${output}"
}

deseq2