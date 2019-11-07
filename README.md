# HMM-scRNA-seq
Author: Michael Wang (fw262@cornell.edu)
Link to preprint: 

## Outline
HMM scRNA-seq is a workflow that enables the discovery of transcripts beyond those listed in gene annotations in single-cell analysis. This workflow uses the Drop-seq package (https://github.com/broadinstitute/Drop-seq/releases/tag/v2.3.0) to align single-cell tagged sequencing reads to a genome without gene annotations. The alignment files are then used to define *de-novo* features using a modified version of groHMM (https://github.com/dankoc/groHMM) generating a peak feature annotation file. These *de-novo* features are compared against an existing set of gene annotations to label the de-novo features as in-gene when they overlap with existing gene annotations or out-gene when they do not. These de-novo features are then used to generate a feature expression matrix in parallel with a gene expression matrix.

## Required Software
1. [snakemake](https://snakemake.readthedocs.io/en/stable/)

Please install the Snakemake workflow managment system as instructed on their website.

2. [Drop-seq computational tools](https://github.com/broadinstitute/Drop-seq/releases)

Download the latest release of the Drop-seq tools (i.e. the Drop-seq_tools-2.3.0.zip package), unzip the folder and add a link to this folder in the config.yaml file.

3. 
