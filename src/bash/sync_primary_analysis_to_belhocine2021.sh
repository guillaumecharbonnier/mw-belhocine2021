#!/usr/bin/env bash

# 1. Set variables to make the script works from any directory
BASH_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
SM_DIR=$BASH_DIR/../..
RMD_DIR=$BASH_DIR/../rmd/belhocine2021

cd $SM_DIR

# 2. Get the latest version of the repo
git pull

# 3. Copy from private mw-tall repository to public mw-belhocine2021 one
grep '^smi <- "' $RMD_DIR/*.Rmd | sed -e 's/^.*smi <- "//' -e 's/"//' > to_sync.txt

rsync -zarmvP --files-from to_sync.txt $SM_DIR/../mw-tall  $SM_DIR

rm -f to_sync.txt
