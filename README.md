# Swaza_case_study
Interview Case Study 

## Write Up
### Introduction & Background
  For this small case study I set up a simple mock repo as if I were working on a normal project under a laboratory. I set up a mock scratch directory under the ./scratch directory for any files that would utilize too much storage under a work directory to save storage and be able to save space on an HPC. I kept files organized and worked with shell to maintain an easy to read and friendly debugging method for each step during each task for reproducibility and tracking steps. I will elaborate further in my methods for packages used and why certain parameters were used.
### Methods
#### Pipeline Organization
The pipeline was written in bash (run.sh) and structured as modular functions, each calling a dedicated R script. This design keeps each step isolated, easier to debug, and reproducible. The script was run with set -ex to print each command before execution and exit immediately on any error. A shared output directory under ../../scratch/results/task_1/seurat_scripts/ was used across all steps to pass intermediate .rds objects between scripts for any debugging issues as well as easy submissions in slurm.
##### Why Seurat?
I used Seurat because of familiarity of use in previous research settings with well documentated pipeline and explanations on parameters as well as acceseible small training data set. 
