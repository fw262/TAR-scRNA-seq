# TAR-scRNA-seq (`from_STARsolo`)
Authors: Michael Wang (fw262@cornell.edu) & David McKellar (dwm269@cornell.edu)

## Outline
This workflow for TAR-scRNA-seq was written for datasets which have already been aligned with [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md). All output files are saved into the output directory (```sampleID_STARsolo/TAR```). See the main README file or the manuscript for details on this analysis.

## Required Software
This workflow requires the following packages listed below. Please ensure that tool can be called from the command line (i.e. the paths to each tool is in your path variable). We recommend using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to organize your packages. We also included a .yml (```scTAR_STARsolo.yml```) file which can be used to initialize a conda environment with all of the required dependencies.

### 1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)

### 2. [Drop-seq Computational Tools v2.3.0](https://github.com/broadinstitute/Drop-seq/releases)

Run the following command in your command line.
```
conda install -c bioconda dropseq_tools
```
### 3. [Picard Tools, version 2.18.29-0 or greater](https://broadinstitute.github.io/picard/)
```
conda install -c bioconda picard
```

### 4. [featureCounts](https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts)

featureCounts is a tool in the subreads package- we use it to tag the .bam file prior to UMI counting

<!-- ### num. [umi_tools](https://umi-tools.readthedocs.io/en/latest/index.html)
```
conda install -c bioconda -c conda-forge umi_tools==1.1.2
```
```
pip install umi_tools
``` -->
### 5. [R, version 3.6 or greater](https://www.r-project.org/)

Please also ensure that you have downloaded the following R packages. They will be used throughout the pipeline.
- [BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)
- [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [groHMM](https://www.bioconductor.org/packages/release/bioc/html/groHMM.html)
- [Seurat, version >= 4.0](https://satijalab.org/seurat/install.html)
- [data.table](https://github.com/Rdatatable/data.table)
- [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
- [stringr](https://cran.r-project.org/web/packages/stringr/readme/README.html)

### 6. [Samtools](http://www.htslib.org/)
```
conda install -c bioconda samtools
```
### 7. [GtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)
This tool is used to convert gtf annotation files to refFlat format.
```
conda install -c bioconda ucsc-gtftogenepred
```
### 8. [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
Please make sure this tool is available in your working environment.
```
conda install -c bioconda bedtools
```
### 9. [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

**Please also download the nt database.**

## Procedure

### 1. Align samples with STARsolo.
- Be sure that the aligned /bam files hav the `CB` and

### 2. Clone this repository.

Run the following command in your command line.
```
git clone https://github.com/fw262/TAR-scRNA-seq
```

### 3. Download required software listed above.

We recommend using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to manage packages, see above.

### 4. Align .fastq files with [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).

Place each STAR output into the same parent directory, so that it follows this file tree: (`DATADIR` in `STARsolo_config.yaml`).
```
DATADIR
├── sample_A
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Aligned.sortedByCoord.out.bam.bai
│   ├── Solo.out
│   │   ├── Gene
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
├── sample_B
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Aligned.sortedByCoord.out.bam.bai
│   ├── Solo.out
│   │   ├── Gene
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
├── sample_C
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Aligned.sortedByCoord.out.bam.bai
│   ├── Solo.out
│   │   ├── Gene
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
```

### 5. Edit the config.yaml file for your experiment.

Please change the variable names in the `config.yaml` as required for your analysis. This includes the following changes:
- **Samples**: STARsolo output directory names; each should be a directory in `DATADIR`
- **GENOME_FASTA**: path to the .fasta (genome sequence) file used to build STAR reference
- **GENOME_GTF**: path to the .gtf (gene annotations) file used to build STAR reference
- **DATADIR**: Path to where the STARsolo outputs (`.../DATADIR/sample_A`) are stored. Also note that outputs will be stored in each sample's directory, individually. See above for the expected file tree.
- **TMPDIR**: Directory to store temporary files.
- **PICARD**: Path to the picard tools .jar file.
- **DROPSEQ**: Path to the Dropseq tools folder.
- **GTFTOGENEPRED**: Path to gtfToGenePred tool. **NOTE** This step will only need to run once for each reference genome you use, and the resulting .refFlat file ("genes.refFlat") will be stored in the same directory as `GENOME_GTF`.
- **BLASTDB**: Directory in which the nt BLAST database is stored.
- **CORES**: Number of cores used in each step of the pipeline. To run multiple samples in parallel, please specify total number of cores in the snakemake command (i.e. "snakemake -j {total cores}").
- **MERGEBP**: Number of bases to merge in groHMM. Smaller numbers creates more TARs but takes longer to run. We recommend keeping the default value of 500.
- **THRESH**: Used to set TARs coverage threshold. This is sequence depth dependent. Default coverage threshold set at 1 in 10,000,000 uniquely aligned reads. For example, if there are 500,000,000 total aligned reads, TARs with at least 50 reads are kept when **THRESH** is set to 10,000,000. A higher **THRESH** value increases the number of TARs kept after filtering. We recommend keeping the default value of 10000000.

### 6. Run snakemake with the command ```snakemake -j [# total cores]```.

From the command line, `cd` into the directory ```.../TAR-SCRNA-seq/from_STARsolo```

Please ensure the Snakefile and STARsolo_config.yaml files as well as the scripts folder are in the directory where you intend to run the pipeline.

## Test datset
TODO- make toy dataset...

## Output files
All output files are stored in a directory inside of the STARsolo output  ```.../DATADIR/sample_A/TAR/```
- RefFlat format of TAR features with and without consideration of directionality stored in "**TAR_reads.bed.gz.withDir.genes.refFlat**" and "**TAR_reads.bed.gz.noDir.refFlat.refFlat**".
- Digital expression matrix for TAR features, with consideration of TAR directionality relative to annotated gene features, is stored in "**TAR_expression_matrix_withDir.txt.gz**". Additional rules within the **Snakefile** are provided to generate a matrix that does not consider strandedness- just uncomment them.
- A list of differentially expressed genes and uTARs in "**results_out/{sample}/{sample}\_diffMarkers.txt**".
- A list of differentially expressed uTARs and their labels based on BLASTn results in "**results_out/{sample}/{sample}\_TAR_diff_uTAR_Features_Labeled.txt**".
- Results of the BLASTn analysis for differentially expressed uTARs in "**TAR_blastResults.txt**".

#### See the README in the root directory of this repository for more details.
