#!/usr/bin/env bash

# 1. Set variables to make the script works from any directory
BASH_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
SM_DIR=$BASH_DIR/../..
RMD_DIR=$BASH_DIR/../rmd/belhocine2021

# 3. Find Snakemake targets in the Rmd files
SMI=`grep '^smi <- "' $RMD_DIR/*.Rmd | sed 's/^.*smi <- "//' | tr '"' ' ' | tr '\r\n' ' '` 

# 4. Run Snakemake on all targets in the Rmd files
eval "$(conda shell.bash hook)"
conda activate snakemake
cd $SM_DIR

mkdir -p $SM_DIR/out/snakemake
snakemake -prk --rerun-incomplete --cores 2 --use-conda $SMI --dag > $SM_DIR/out/snakemake/all_dag.dot
snakemake -prk --rerun-incomplete --cores 2 --use-conda $SMI --rulegraph > $SM_DIR/out/snakemake/all_rulegraph.dot
snakemake -prk --rerun-incomplete --cores 2 --use-conda $SMI -n > $SM_DIR/out/snakemake/all_tasks.txt

cd $SM_DIR/out/snakemake
dot -Tpdf -o all_dag.pdf all_dag.dot
dot -Tpdf -o all_rulegraph.pdf all_rulegraph.dot


# 5. Compile the Bookdown report
# cd $RMD_DIR
# conda activate r_rmd_belhocine2021
# ./compile.R
