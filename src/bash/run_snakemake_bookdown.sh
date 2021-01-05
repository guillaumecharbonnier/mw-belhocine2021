#!/usr/bin/env bash

# 1. Set variables to make the script works from any directory
BASH_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
SM_DIR=$BASH_DIR/../..
RMD_DIR=$BASH_DIR/../rmd/belhocine2021

# 2. Create required conda environnements
conda env create -f $SM_DIR/../mw-lib/src/snakemake/envs/snakemake.yaml
conda env create -f $SM_DIR/src/snakemake/envs/r_rmd_belhocine2021.yaml

# 3. Find Snakemake targets in the Rmd files
SMI=`grep '^smi <- "' $RMD_DIR/*.Rmd | sed 's/^.*smi <- "//' | tr '"' ' ' | tr '\r\n' ' '` 
echo $SMI

# 4. Run Snakemake on all targets in the Rmd files
eval "$(conda shell.bash hook)"
conda activate snakemake
cd $SM_DIR
snakemake -prk --rerun-incomplete --cores 2 --use-conda $SMI

# 5. Compile the Bookdown report
# cd $RMD_DIR
# conda activate r_rmd_belhocine2021
# ./compile.R
