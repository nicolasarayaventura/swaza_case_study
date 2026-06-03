set -x -e
scratch="../../scratch/"

function data_download {
    cat dataset_link.txt | while read -r dataset url ; do
        output="${scratch}/data/task_1/${dataset}"
        mkdir -p "${output}"
        wget -P "${output}" "${url}"
    done
}
