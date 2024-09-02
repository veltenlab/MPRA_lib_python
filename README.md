# MPRA_lib_python

# Sequencing Analysis Pipeline

This is an internal pipeline for the sequencing analysis of the data obtained from MPRA experiments (transMPRA/scMPRA/bulkMPRA) based on the following workflow:

![Workflow Image](path/to/workflow_image.png)

## Directory Structure

- `results/`: This directory will store the results locally (output CSV file with counts, Rmarkdown output).
- `data/example_data`: This directory contains small example data extracted from the real experiments, as well as the corresponding reference files.

## Configuration

The pipeline is applicable for 3 modes: TRANS, SC and BULK, which need to be specified in the config file prior to running the pipeline. Additionally, full paths to the FASTQ files must be provided in the config file.

## Prerequisites

Before you can run this pipeline, you will need to have the following installed on your system:

## Installation and Setup

Follow the steps below to set up and run the pipeline.

### 1. Clone the Repository

First, you need to clone the Git repository containing the pipeline to your local machine.

```bash
git clone our repo
cd your-repo
```

### 2. Create local environment from the given in this folder file.

conda env create -f MPRA_lib.yml
conda activate MPRA_lib






