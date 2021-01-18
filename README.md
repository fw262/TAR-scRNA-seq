# TAR-scRNA-seq
Author: Michael Wang (fw262@cornell.edu)
Link to preprint: 

## Outline
TAR-scRNA-seq (Transcriptionally Active Region single-cell RNA-seq) is a workflow that enables the discovery of transcripts beyond those listed in gene annotations in scRNA-seq analysis. The workflow aligns single-cell tagged sequencing reads to a genome without gene annotations using [STAR](https://github.com/alexdobin/STAR). The alignment files are then used to define *de-novo* transcriptionally active regions (TARs) using a modified version of [groHMM](https://github.com/dankoc/groHMM) based on regions along the genome with signficant read coverage. TARs are compared against an existing set of gene annotations to label them as annotated (aTARs) when they overlap with existing gene annotations or unanotated (uTARs) when they do not. These labeled TARs are then used to generate a TAR feature expression matrix in parallel with a gene expression matrix. Seurat tools are used to find differentially expressed uTARs which are labeled based on BLASTn sequence homology.

## Required Software
This workflow requires the following packages listed below. Please ensure that tool can be called from the command line (i.e. the paths to each tool is in your path variable).

### 1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)

### 2. [Drop-seq Computational Tools v2.3.0](https://github.com/broadinstitute/Drop-seq/releases)

Run the following command in your command line.
```
conda install -c bioconda dropseq_tools
```
### 3. [Picard Tools](https://broadinstitute.github.io/picard/)

### 4. [STAR Aligner](https://github.com/alexdobin/STAR/releases)

### 5. [R, version 3.6 or greater](https://www.r-project.org/)

Please also ensure that you have downloaded the following R packages. They will be used throughout the pipeline. 
- [BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
- [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [groHMM](https://www.bioconductor.org/packages/release/bioc/html/groHMM.html)
- [Seurat, version 3.2.1](https://satijalab.org/seurat/install.html)
- [data.table](https://github.com/Rdatatable/data.table)
- [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
- [stringr](https://cran.r-project.org/web/packages/stringr/readme/README.html)

### 6. [Samtools](http://www.htslib.org/)

### 7. [GtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)

This tool is used to convert gtf annotation files to refFlat format.

### 8. [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)

Please make sure this tool is available in your working environment.

### 9. [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

Please also download the nt database.

## Procedure

### 1. Clone this repository.

Run the following command in your command line.
```
git clone https://github.com/fw262/TAR-scRNA-seq
```

### 2. Download required software listed above.

Please ensure to include all required software before starting.

### 3. Store or link paired end sequencing files.

Please combine and move the fastq files for each experiment into one data directory. Please ensure the sequence files end in "{sample}\_R1.fastq.gz" and "{sample}\_R2.fastq.gz" in your data directory.

### 4. Edit the config.yaml file for your experiment.

Please change the variable names in the config.yaml as required for your analysis. This includes the following changes:
- **Samples**: Samples prefix (before the \_R1.fastq.gz)
- **GENOMEREF**: Path to your genome assembly file.
- **REFGTF**: Path to your gene annotations in gtf format. **NOTE** If your gene annotation file is already in the refFlat format, please rename the annotation file to "refFlat.refFlat" and store in the same directory as the Snakefile and config.yaml files.
- **SAMPWDIR**: Directory where summary information is stored.
- **DATADIR**: Path to where the sequencing samples ({sample}\_R1.fastq.gz) are stored.
- **TMPDIR**: Directory to store temporary files.
- **PIPELINE_MAJOR**: Directory where the outputs (expression matrices, differentially expressed uTARs) are stored.
- **GLOBAL**: Define global variables for pipeline including number of mismatches allowed in STAR, cell barcode base pair range in read 1, and UMI base pair range in read 1.
- **PICARD**: Path to the picard tools .jar file.
- **DROPSEQ**: Path to the Dropseq tools folder.
- **GTFTOGENEPRED**: Path to gtfToGenePred tool. **NOTE** If your gene annotation file is already in the refFlat format, please rename the annotation file to "refFlat.refFlat" and store in the same directory as the Snakefile and config.yaml files.
- **STAREXEC**: Path to the STAR aligner tool.
- **BLASTDB**: Directory in which the nt BLAST database is stored.
- **CORES**: Number of cores used in each step of the pipeline. To run multiple samples in parallel, please specify total number of cores in the snakemake command (i.e. "snakemake -j {total cores}").
- **expectedCells**: Expected number of cells in each scRNA-seq experiment.
- **MERGEBP**: Number of bases to merge in groHMM. Smaller numbers creates more TARs but takes longer to run. We recommend keeping the default value of 500.
- **THRESH**: Used to set TARs coverage threshold. This is sequence depth dependent. Default coverage threshold set at 1 in 10,000,000 uniquely aligned reads. For example, if there are 500,000,000 total aligned reads, TARs with at least 50 reads are kept when **THRESH** is set to 10,000,000. A higher **THRESH** value increases the number of TARs kept after filtering. We recommend keeping the default value of 10000000.

### 5. Run snakemake with the command "snakemake".

Please ensure the Snakefile and config.yaml files as well as the scripts folder are in the directory where you intend to run the pipeline.

## Test datset

A subset of the chicken embryonic heart development sequencing data is attached in the **testData_small** folder. To download the corresponding references for this chicken dataset, please visit https://useast.ensembl.org/Gallus_gallus/Info/Index. The full chicken embryonic heart development dataset is available at [GSE149457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149457).


To generate expression matrices (gene, stranded and unstranded TAR), please run the snakemake rule "getMats" with the following command (filling in # cores):
```
snakemake -R --until getMats -j [# cores]
```
Please note that the test data will generate small expression matrices. 

To run the full pipeline including generating a list of labeled differentially expressed uTARs, run the default snakemake command:
```
snakemake -j [# cores]
```
After making sure that the pipeline generates expression matrices, you may test the labeling of uTARs 

You should be able to run through this small dataset without errors, generating gene and TAR expression matrices as well as a list of labeled differentially expressed uTARs. The entire pipeline should take less than 15 minutes to complete on the test dataset with 1 core.


## Output

- RefFlat format of TAR features with and without consideration of directionality stored in "**TAR_reads.bed.gz.withDir.refFlat.refFlat**" and "**TAR_reads.bed.gz.noDir.refFlat.refFlat**".
- Digital expression matrix for gene features is stored in "**results_out/{sample}/{sample}\_gene_expression_matrix.txt.gz**".
- Digital expression matrix for TAR features, without consideration of TAR directionality relative to annotated gene features, is stored in "**results_out/{sample}/{sample}\_TAR_expression_matrix_noDir.txt.gz**".
- Digital expression matrix for TAR features, with consideration of TAR directionality relative to annotated gene features, is stored in "**results_out/{sample}/{sample}\_TAR_expression_matrix_withDir.txt.gz**".
- A list of differentially expressed genes and uTARs in "**results_out/{sample}/{sample}\_diffMarkers.txt**".
- A list of differentially expressed uTARs and their labels based on BLASTn results in "**results_out/{sample}/{sample}\_diffuTARMarkersLabeled.txt**".
- Results of the BLASTn analysis for differentially expressed uTARs in "**results_out/{sample}/{sample}\_blastResults.txt**".

### Format of TAR feature label

TAR features, listed in the refFlat and expression matrix files, are named based on their position, total coverage, and whether they overlap with an existing gene annotation. Examples listed below

- **chr3_40767549_40767699_+\_187_0** means that this TAR feature is located at chr3:40767549-40767699 on the positive strand with a total read coverage of 187. The "\_0" means that this is a **uTAR** feature, no overlap with an existing gene annotation.
- **chr6_42888199_42888349_-\_983_RPS6KA2_+\_1** means that this TAR feature is located at chr6:42888199-42888349 on the negative strand with a total read coverage of 983. This feature overlaps in genomic position with the gene annotated as RPS6KA2, which is annotated on the positive strand. The "\_1" means that this is an **aTAR** feature overlapping an existing gene annotation without considering directionality.

## Data Availability

The following publicly available datasets are used in the manuscript. 

- [Human PBMC data](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k?)
- [Mouse Atlas](https://tabula-muris.ds.czbiohub.org/), [GSE109774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109774)
- Naked mole rat spleen data [GSE132642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132642)
- Sea urchin embryo data [GSE134350](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134350)
- Chicken embryonic heart [GSE149457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149457)
- Mouse lemur

## Frequently Asked Questions (FAQs)

### 1. Error in rule calcHMMrefFlat stating "Error in read.table(file = file, header = header, sep = sep, quote = quote,  : no lines available in input". 
  
Please refer to the "SingleCellHMM_Run_combined_bam_HMM_features.log" created in the same directory as your Snakefile. You will likely see the following error:
```
cannot read: chr*_HMM.bed: No such file or directory
```
Please manually install [groHMM](https://bioconductor.org/packages/release/bioc/html/groHMM.html) and [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) and make sure you can load these packages in R. This error indicates that the groHMM R script did not finish running to generate groHMM bed files.
