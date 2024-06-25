# MPRA_lib_python

# Sequencing Analysis Pipeline

This is an internal pipeline for the sequencing analysis of the data obtained from MPRA experiments (transMPRA/ssMPRA/bulkMPRA) based on the following workflow:

![Workflow Image](path/to/workflow_image.png)

## Directory Structure

- `results/`: This directory will store the results locally (output CSV file with counts, Rmarkdown output).
- `data/`: This directory contains small example data extracted from the real experiments, as well as the corresponding reference files.

## Configuration

The pipeline is applicable for 3 modes: trans, ss and bull, which need to be specified in the config file prior to running the pipeline. Additionally, full paths to the FASTQ files must be provided in the config file.

