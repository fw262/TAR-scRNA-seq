# TAR-scRNA-seq
Author: Michael Wang (fw262@cornell.edu)

Please cite the following paper when using this tool:

[Wang, Michael FZ, et al. "Uncovering transcriptional dark matter via gene annotation independent single-cell RNA sequencing analysis." _Nature Communications_ 12.1 (2021): 1-10.](https://www.nature.com/articles/s41467-021-22496-3)

## Abstract
Conventional scRNA-seq expression analyses rely on the availability of a high quality genome annotation. Yet, as we show here with scRNA-seq experiments and analyses spanning human, mouse, chicken, mole rat, lemur and sea urchin, genome annotations are often incomplete, in particular for organisms that are not routinely studied. To overcome this hurdle, we created a scRNA-seq analysis routine that recovers biologically relevant transcriptional activity beyond the scope of the best available genome annotation by performing scRNA-seq analysis on any region in the genome for which transcriptional products are detected. Our tool generates a single-cell expression matrix for all transcriptionally active regions (TARs), performs single-cell TAR expression analysis to identify biologically significant TARs, and then annotates TARs using gene homology analysis. This procedure uses single-cell expression analyses as a filter to direct annotation efforts to biologically significant transcripts and thereby uncovers biology to which scRNA-seq would otherwise be in the dark.

## Getting Started
We have written Snakemake workflows to perform this analysis. Select your workflow by that files which you start with. Current workflows support starting with either raw fastq files (*from_fastqs*) or if you have already aligned with cellranger count (*from_cellranger*). See the README files inside each directory for install and runtime directions.

## Frequently Asked Questions (FAQs)

### 1. Error in rule calcHMMrefFlat stating "Error in read.table(file = file, header = header, sep = sep, quote = quote,  : no lines available in input".

Please refer to the "SingleCellHMM_Run_combined_bam_HMM_features.log" created in the same directory as your Snakefile. You will likely see the following error:
```
cannot read: chr*_HMM.bed: No such file or directory
```
Please manually install [groHMM](https://bioconductor.org/packages/release/bioc/html/groHMM.html) and [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) and make sure you can load these packages in R. This error indicates that the groHMM R script did not finish running to generate groHMM bed files.
