# mw-belhocine2021

This repository contains code for the analyses from the article "Dynamic of broad H3K4me3 domains uncover an epigenetic switch between cell identity and cancer-related genes".
The workflow is written as a Snakemake wokflow which processes raw files into a Bookdown report. Since some of the raw files come from the DAC-restricted [Blueprint](https://www.blueprint-epigenome.eu/index.cfm?p=9CE408D2-B4A2-846D-C3608A21096EACF4) project, the whole analysis can not be reproduced without access to these files. 

The workflow has been compiled on a Linux environment with [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html) installed, following these instructions:

```
git clone git@github.com:guillaumecharbonnier/mw-lib.git
git clone git@github.com:guillaumecharbonnier/mw-belhocine2021.git
mw-belhocine2021/src/bash/run_snakemake_bookdown.sh
```

The [rulegaph](out/snakemake/all_rulegraph.pdf), [DAG](out/snakemake/all_dag.pdf) and listing of [all tasks](out/snakemake/all_tasks.txt) may help giving an overview of the Snakemake workflow without having to dive into the code.

Please open an issue on [github](https://github.com/guillaumecharbonnier/mw-belhocine2021/issues) if you have any questions.
