# MPRA_lib_python

# Sequencing Analysis Pipeline

This is an internal pipeline for the sequencing analysis of the data obtained from MPRA experiments (transMPRA/scMPRA/bulkMPRA) based on the following workflow:

![Workflow Image](path/to/workflow_image.png)

## Directory Structure

- `results/`: This directory will store the results locally (output CSV file with counts, Rmarkdown output)
- `data/`: You can add locally your data in this directory if you want or specify full paths otherwise
- `data/example_data`: This directory contains small example data extracted from the real experiments, as well as the corresponding reference files
- `scripts/`: This directory contains Python and R files for processing


# Installation

1. Clone this repository:
```shell
git clone git@github.com:veltenlab/MPRA_lib_python.git
```
2. Change directory to the repository: `cd MPRA_lib_python`
3. Install the conda environment. It is recommended to use mamba [mamba](https://mamba.readthedocs.io/en/latest/index.html) or [miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)

```shell
conda env create -n MPRA_env -f MPRA_env.yaml 
```
4. Activate the environment: `mamba activate MPRA_env`

## Configuration

The pipeline is applicable for 3 modes: TRANS, SC and BULK, which need to be specified in the config file prior to running the pipeline. Additionally, full paths to the FASTQ files must be provided in the config file.

## Quick start: Customize pipeline

Before running a snakefile, config.yaml file has to be customized:
- choose the type of the performed MPRA experiment (trans/sc/bulk)
- add data paths for sequencing data
- adjust the number of threads used for alignment







