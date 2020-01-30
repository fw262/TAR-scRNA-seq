# TAR-scRNA-seq
Author: Michael Wang (fw262@cornell.edu)
Link to preprint: 

## Outline
TAR-scRNA-seq (Transcriptionally Active Region single-cell RNA-seq) is a workflow that enables the discovery of transcripts beyond those listed in gene annotations in scRNA-seq analysis. The workflow aligns single-cell tagged sequencing reads to a genome without gene annotations using [STAR](https://github.com/alexdobin/STAR). The alignment files are then used to define *de-novo* transcriptionally active regions (TARs) using a modified version of [groHMM](https://github.com/dankoc/groHMM) based on regions along the genome with signficant read coverage. TARs are compared against an existing set of gene annotations to label them as annotated (aTARs) when they overlap with existing gene annotations or unanotated (uTARs) when they do not. These labeled TARs are then used to generate a TAR feature expression matrix in parallel with a gene expression matrix.

## Required Software
This workflow requires the following packages listed below.

1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)

2. [Drop-seq Computational Tools v1.13](https://github.com/broadinstitute/Drop-seq/releases)

Download the Drop-seq tools v1.13 (i.e. conda install -c bioconda dropseq_tools=1.13).

3. [Picard Tools](https://broadinstitute.github.io/picard/)

4. [STAR Aligner](https://github.com/alexdobin/STAR/releases)

6. [R](https://www.r-project.org/)

7. [Samtools](http://www.htslib.org/)

8. [groHMM](https://bioconductor.org/packages/release/bioc/html/groHMM.html)

9. [gtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)

This tool is used to convert gtf annotation files to refFlat format.

## Procedure
### 1. Download required software.

Please ensure to include all required software before starting.

### 2. Download Snakefile, config.yaml, and the scripts folder.

Please download all of the files in this repository including the scripts folder and its content.

### 3. Store or link paired end sequencing files.

Please combine and move the fastq files for each experiment into one data directory. Please ensure the sequence files end in "{sample}\_R1.fastq.gz" and "{sample}\_R2.fastq.gz" in your data directory.

### 4. Edit the config.yaml file as needed for your experiment.

Please change the variable names in the config.yaml as required for your analysis. This includes the following changes:
- **Samples**: Samples prefix (before the \_R1.fastq.gz)
- **GENOMEREF**: Path to your genome assembly file.
- **REFGTF**: Path to your gene annotations in gtf format.
- **SAMPWDIR**: Directory where summary information is stored.
- **DATADIR**: Path to where the sequencing samples ({sample}\_R1.fastq.gz) are stored.
- **TMPDIR**: Directory to store temporary files.
- **PIPELINE_MAJOR**: Directory where the expression matrices are stored.
- **GLOBAL**: Define global variables for pipeline including number of mismatches allowed in STAR, cell barcode base pair range in read 1, and UMI base pair range in read 1.
- **PICARD**: Path to the picard tools .jar file.
- **DROPSEQ**: Path to the Dropseq tools folder.
- **GTFTOGENEPRED**: Path to gtfToGenePred tool. *NOTE* If your reference file is already in the refFlat format, please rename the annotation file to "refFlat.refFlat" and store in the same directory as the Snakefile and config.yaml files.
- **STAREXEC**: Path to the STAR aligner tool.
- **CORES**: Number of cores used in each step of the pipeline. To run multiple samples in parallel, please specify total number of cores in the snakemake command (i.e. "snakemake -j {total cores}").
- **expectedCells**: Expected number of cells in each scRNA-seq experiment.
- **MERGEBP**: Number of bases to merge in groHMM. Smaller numbers creates more TARs but takes longer to run. We recommend keeping the default value of 500.
- **THRESH**: Used to set TARs coverage threshold. This is sequence depth dependent. Default coverage threshold set at 1 in 10,000,000 uniquely aligned reads. For example, if there are 500,000,000 total aligned reads, TARs with at least 50 reads are kept when *THRESH* is set to 10,000,000. A higher *THRESH* value increases the number of TARs kept after filtering. We recommend keeping the default value of 10000000.

### 5. Run snakemake with the command "snakemake".

Please ensure the Snakefile and config.yaml files as well as the scripts folder are in the directory where you intend to run the pipeline.

## Output

- Digital expression matrix for gene features is stored in "results_out/{sample}/{sample}\_gene_expression_matrix.txt.gz".
- Digital expression matrix for TAR features, without consideration of TAR directionality relative to annotated gene features, is stored in "results_out/{sample}/{sample}\_TAR_expression_matrix_noDir.txt.gz".
- Digital expression matrix for TAR features, with consideration of TAR directionality relative to annotated gene features, is stored in "results_out/{sample}/{sample}\_TAR_expression_matrix_withDir.txt.gz".
- refFlat format of TAR features with and without consideration of directionality stored in "TAR_reads.bed.gz.withDir.refFlat.refFlat" and "TAR_reads.bed.gz.noDir.refFlat.refFlat".

### Format of TAR feature label

TAR features, listed in the refFlat and expression matrix files, are named based on their position, total coverage, and whether they overlap with an existing gene annotation. Examples listed below

- *chr3_40767549_40767699_+_187_0* means that this TAR feature is located at chr3:40767549-40767699 on the positive strand with a total read coverage of 187. The "\_0" means that this is a uTAR feature, no overlap with an existing gene annotation.
- *chr6_42888199_42888349_-_983_RPS6KA2_+_1* means that this TAR feature is located at chr6:42888199-42888349 on the negative strand with a total read coverage of 983. This feature overlaps in genomic position with the gene annotated as RPS6KA2, which is annotated on the positive strand. The "\_1" means that this is an aTAR feature overlapping an existing gene annotation without considering directionality.
