# INSTRUCTIONS
----

## Software
#### Missing software
if any of the following tools are missing, please install it/them
- install git
- install conda (anaconda or miniconda)
- install snakemake, via conda or pip (easier via pip; version 6.15.5 works fine with the current snakemake-file of the pipeline)
```
pip install snakemake==6.15.5
```

----


## Preparational steps

#### Screen session
create a linux-screen-session and change to base environment
```
screen -S session_name
conda deactivate
```

#### Folder
create an empty folder for your read-in-task and change into it
```
mkdir my_folder
cd my_folder
```

#### Working directory
change into the main folder "ngs_pipeline" to set it as working directory
```
cd ngs_pipeline
```

#### File/Sample overview for the NGS preprocessing pipeline
- create a file/smp overview that shows which files should be preprocessed; in the same time smpIDs for the files are given
- when the input-folder is somewhere else, then edit the variable "input_folder" and run then the script
```
time Rscript src/create_sample_sheet_for_core_preprocessing.R
```

#### edit folders in the config.yaml
edit following variables:
- sample_sheet_for_core_preprocessing: <PATH>/sample_sheet_for_core_preprocessing.txt
- dir_results: <PATH>/preprocessed/

----

## ACTUAL NGS preprocessing pipeline

#### NGS Preprocessing pipeline
- runtime around 1-3 hours
- editing the config-file of the actual pipeline (snakemake pipeline) is not needed: configs/config.yaml
  - you can edit the parameters for trim_and_collapse if needed (when different adapters are used, etc)
- run in terminal following command
```
time snakemake -s=Snakefile -rp -j 4 --use-conda --reason --configfile=configs/config.yaml
```
- -j 4 means using 4 cores (use maximum 4 cores, minimum 1)

