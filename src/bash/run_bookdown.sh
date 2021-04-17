#!/usr/bin/env bash

# 1. Set variables to make the script works from any directory
BASH_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
SM_DIR=$BASH_DIR/../..
RMD_DIR=$BASH_DIR/../rmd/belhocine2021

# 2. Create required conda environnements
conda env create -f $SM_DIR/src/snakemake/envs/r_rmd_belhocine2021.yaml

# 3. Compile the Bookdown report
cd $RMD_DIR
conda activate r_rmd_belhocine2021
./compile.R
