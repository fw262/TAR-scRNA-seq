# TAR-scRNA-seq cellranger workflow
Authors: Michael Wang (fw262@cornell.edu) & David McKellar (dwm269@cornell.edu)

[Nature Communications link](https://www.nature.com/articles/s41467-021-22496-3)

## Outline
This workflow for TAR-scRNA-seq was written for datasets which have already been aligned with [cellranger count] (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct), a widely used pipeline from 10x Genomics. All output files are saved into the cellranger count output directory (```sampleID_cellranger_count/TAR```).  Note that this pipeline will only be compatible with data generated by the 10x Genomics Chromium platform.

TAR-scRNA-seq (Transcriptionally Active Region single-cell RNA-seq) is a workflow that enables the discovery of transcripts beyond those listed in gene annotations in scRNA-seq analysis. The workflow aligns single-cell tagged sequencing reads to a genome without gene annotations using [STAR](https://github.com/alexdobin/STAR). The alignment files are then used to define *de-novo* transcriptionally active regions (TARs) using a modified version of [groHMM](https://github.com/dankoc/groHMM) based on regions along the genome with signficant read coverage. TARs are compared against an existing set of gene annotations to label them as annotated (aTARs) when they overlap with existing gene annotations or unanotated (uTARs) when they do not. These labeled TARs are then used to generate a TAR feature expression matrix in parallel with a gene expression matrix. Seurat tools are used to find differentially expressed uTARs which are labeled based on BLASTn sequence homology.

## Required Software
This workflow requires the following packages listed below. Please ensure that tool can be called from the command line (i.e. the paths to each tool is in your path variable). We recommend using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to organize your packages. We also included a .yml (```scTAR_cellranger.yml```) file which can be used to initialize a conda environment with all of the required dependencies.

### 1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)

### 2. [Drop-seq Computational Tools v2.3.0](https://github.com/broadinstitute/Drop-seq/releases)

If using conda, run the following command in your command line.
```
conda install -c bioconda dropseq_tools
```
### 3. [Picard Tools, version 2.18.29-0 or greater](https://broadinstitute.github.io/picard/)

### 4. [R, version 3.6 or greater](https://www.r-project.org/)

Please also ensure that you have downloaded the following R packages. They will be used throughout the pipeline.
- [BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
- [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [groHMM](https://www.bioconductor.org/packages/release/bioc/html/groHMM.html)
- [Seurat, version >= 4.0](https://satijalab.org/seurat/install.html)
- [data.table](https://github.com/Rdatatable/data.table)
- [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
- [stringr](https://cran.r-project.org/web/packages/stringr/readme/README.html)

### 5. [Samtools](http://www.htslib.org/)

### 6. [GtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)

This tool is used to convert gtf annotation files to refFlat format.

### 7. [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)

Please make sure this tool is available in your working environment.

### 8. [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

**Please also download the nt database.**

## Procedure #TODO- update

### 1. Clone this repository.

Run the following command in your command line.
```
git clone https://github.com/fw262/TAR-scRNA-seq
```

### 2. Download required software listed above.

Please ensure to include all required software before starting. If using conda, you can create a new environment with the included ```scTAR_cellranger.yml``` file. Be sure to change the path in the last line of the .yml file so that it points to your miniconda3 installation. Install everything using the following command:
```
conda create -f scTAR_cellranger.yml
```
By default, the new environment will be named "scTAR_cellranger", but you can rename it by changing the .yml file prior to environment initialization if you'd like (see the last line of the file).

### 3. Align .fastq files with cellranger count and name the output directory "sampleID_cellranger_count".

Our snakemake pipeline uses this naming convention to simplify tracking each sample. Also please place each cellranger count output into the same directory ("DATADIR" in the config file).

### 4. Edit the config.yaml file for your experiment. TODO

Please change the variable names in the config.yaml as required for your analysis. This includes the following changes:
- **Samples**: Samples prefix (before \_cellranger_count)
- **CR_REF**: Path to your genomic reference, which you used to align your samples (use ```cellranger mkref``` to generate, or download from [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references)).
- **DATADIR**: Path to where the cellranger count outputs ({sample}\_cellranger_count) are stored. Also note that outputs will be stored in each cellranger_count directory, individually.
- **TMPDIR**: Directory to store temporary files.
- **PICARD**: Path to the picard tools .jar file.
- **DROPSEQ**: Path to the Dropseq tools folder.
- **GTFTOGENEPRED**: Path to gtfToGenePred tool. **NOTE** In this workflow, we use the cellranger mkref format. This step will only need to run once for each reference genome you use, and the resulting .refFlat file ("genes.refFlat") will be stored in ""{CR_REF}/genes".
- **BLASTDB**: Directory in which the nt BLAST database is stored.
- **CORES**: Number of cores used in each step of the pipeline. To run multiple samples in parallel, please specify total number of cores in the snakemake command (i.e. "snakemake -j {total cores}").
- **MERGEBP**: Number of bases to merge in groHMM. Smaller numbers creates more TARs but takes longer to run. We recommend keeping the default value of 500.
- **THRESH**: Used to set TARs coverage threshold. This is sequence depth dependent. Default coverage threshold set at 1 in 10,000,000 uniquely aligned reads. For example, if there are 500,000,000 total aligned reads, TARs with at least 50 reads are kept when **THRESH** is set to 10,000,000. A higher **THRESH** value increases the number of TARs kept after filtering. We recommend keeping the default value of 10000000.

### 5. Run snakemake with the command ```snakemake -j [# total cores]```.

Please ensure the Snakefile and config.yaml files as well as the scripts folder are in the directory where you intend to run the pipeline.

## Test datset - TODO- make toy dataset

A subset of the chicken embryonic heart development sequencing data is attached in the **testData_small** folder. To download the corresponding references for this chicken dataset, please visit https://useast.ensembl.org/Gallus_gallus/Info/Index. The full chicken embryonic heart development dataset is available at [GSE149457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149457). Please note that full expression matrices (gene and unstranded TAR) are also included for the day 4 and day 7 datasets in the **testData_small** folder.


To generate expression matrices (gene, stranded and unstranded TAR), please run the snakemake rule "getMats" with the following command:
```
snakemake -R --until getMats -j [# total cores]
```

To test the labeling of differentially expressed uTARs through scRNA-seq and BLASTn analysis (after issuing the command above), please move the expression matrices in the **testData_small** folder to the corresponding results folder with the following commands:
```
cp testData_small/day7_0.25m_*expression_matrix* results_chicken/day7_0.25m/
cp testData_small/day4_0.25m_*expression_matrix* results_chicken/day4_0.25m/
```

To run the full pipeline including generating a list of labeled differentially expressed uTARs, run the default snakemake command:
```
snakemake -j [# cores]
```

Assuming the expression matrices are available, the test dataset should take less than 30 minutes to generate a list of labeled uTARs with 12 cores and 16GB of RAM.

## Output - TODO- update

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

## Frequently Asked Questions (FAQs)

### 1. Error in rule calcHMMrefFlat stating "Error in read.table(file = file, header = header, sep = sep, quote = quote,  : no lines available in input".

Please refer to the "SingleCellHMM_Run_combined_bam_HMM_features.log" created in the same directory as your Snakefile. You will likely see the following error:
```
cannot read: chr*_HMM.bed: No such file or directory
```
Please manually install [groHMM](https://bioconductor.org/packages/release/bioc/html/groHMM.html) and [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) and make sure you can load these packages in R. This error indicates that the groHMM R script did not finish running to generate groHMM bed files.