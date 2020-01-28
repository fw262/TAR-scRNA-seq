# HMM-scRNA-seq
Author: Michael Wang (fw262@cornell.edu)
Link to preprint: 

## Outline
HMM scRNA-seq is a workflow that enables the discovery of transcripts beyond those listed in gene annotations in single-cell analysis. This workflow uses the Drop-seq package (https://github.com/broadinstitute/Drop-seq/releases/tag/v2.3.0) to align single-cell tagged sequencing reads to a genome without gene annotations. The alignment files are then used to define *de-novo* features using a modified version of groHMM (https://github.com/dankoc/groHMM) generating a peak feature annotation file. These *de-novo* features are compared against an existing set of gene annotations to label the *de-novo* features as in-gene when they overlap with existing gene annotations or out-gene when they do not. These *de-novo* features are then used to generate a feature expression matrix in parallel with a gene expression matrix.

## Required Software
This workflow requires the following packages listed below. Most are standard common bioinformatics tools.

1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)

Please install the Snakemake workflow managment system as instructed on their website.

2. [Drop-seq Computational Tools](https://github.com/broadinstitute/Drop-seq/releases)

Download the Drop-seq tools v1.13 (i.e. conda install -c bioconda dropseq_tools=1.13).

3. [Picard Tools](https://broadinstitute.github.io/picard/)

4. [STAR Aligner](https://github.com/alexdobin/STAR/releases)

5. [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

6. [R](https://www.r-project.org/)

7. [Samtools](http://www.htslib.org/)

8. [groHMM](https://bioconductor.org/packages/release/bioc/html/groHMM.html)

## Procedure
### 1. Download required software.

Please ensure to include all required software before starting.

### 2. Download Snakefile, config.yaml, and the scripts folder.

Please download all of the files in this repository including the scripts folder.

### 2.1 Edit the SingleCellHMM_MW.bash file to include full path of SingleCellHMM.R script file

### 3. Move or create a symbolic link to fastq files in the same directory.

Please move all your fastq sequencing files to one data directory. Please ensure the sequence files end in "{sample}\_R1.fastq.gz" and "{sample}\_R2.fastq.gz" in your data directory. You can also create a symbolic link of your sequence files into a data directory. 

### 4. Edit the config.yaml file as needed for your experiment.

Please change the variable names in the config.yaml as required for your analysis. This includes the following changes:
- **Samples**: Samples prefix (before the \_R1.fastq.gz)
- **DATADIR**: Path to where the sequencing samples are stored.
- **GENOMEREF**: Path to your genome assembly file.
- **REFGTF**: Path to your gene annotations in gtf format.
- **GLOBAL**: Variables used to define cell barcode and unique molecule identifier ranges in your read 1 fastq sequencing file.
- **expectedCells**: Expected number of cells in each scRNA-seq experiment.
- **MINCOV**: Minimum number of aligned reads for the region to be considered as a feature for further scRNA-seq analysis.
- **CORES**: Number of cores used in the pipeline.
Please also include paths to the required software packages. Please note that the following tools are included in the scripts folder:
- **SingleCellHMM**: scripts/SingleCellHMM_2.bash
- **generate_refFlat_script**: scripts/generate_refFlat_script_both.R # included in scripts folder
- **GTFTOGENEPRED**: scripts/gtfToGenePred 

### 5. Run snakemake with the command "snakemake".

Please ensure the Snakefile and the edited config.yaml file is contained in the output directory.
