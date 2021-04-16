#########################################################################
#    uTAR Snakemake pipeline
#    Copyright (C) 2020 Michael Wang
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##########################################################################

import pdb
########################################################################################################
# Configfile
########################################################################################################
configfile:'config.yaml'
########################################################################################################
# Variables and references
########################################################################################################
SAMPLEWDIR = config['SAMPWDIR']
MISMATCH = config['GLOBAL']['allowed_mismatch']
DATADIR = config['DATADIR']
GENOMEREF = config['GENOMEREF']
TMPDIR = config['TMPDIR']
gtffile=config['REFGTF']
MERGEBP=str(config['MERGEBP'])
THRESH=str(config['THRESH'])
#FULLSCRIPTPATH=str(config['FULLSCRIPTPATH'])
########################################################################################################
# Executables
########################################################################################################
TMPDIR = config['TMPDIR']
PICARD = config['PICARD']
DROPSEQ = config['DROPSEQ']
STAREXEC = config['STAREXEC']
CORES = config['CORES']
gtfToGenePred=config['GTFTOGENEPRED']

rule all:
############## original default call here
	input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_diffuTARMarkersLabeled.txt', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])

rule getMats:
	input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_getMats.txt', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])

#####################################################################################
# create sequence dictionary using picard tools
rule createGenomeDict:
        input: config['GENOMEREF']
        output: config['GENOMEREF']+'.dict'
        shell:
                """
                java -Dpicard.useLegacyParser=false -jar {PICARD} CreateSequenceDictionary \
                        -R {input} \
                        -O {output}
                """
###########################################################################

# create STAR index
rule generateStar:
        input:
              	fastaFile=config['GENOMEREF']
        output: directory('STAR_ind_noAnno')
        threads: CORES
        shell:
              	"""
                mkdir STAR_ind_noAnno
                \
                {STAREXEC}\
                        --runThreadN {CORES}\
                        --runMode genomeGenerate\
                        --genomeDir STAR_ind_noAnno\
                        --genomeFastaFiles {input.fastaFile}
                """
#####################################################################################

#####################################################################################
# convert GTF to REFFlat
# not run if REFFlat file exists

rule convertToRefFlat1:
        input:  {gtffile}
        output: "refFlat.refFlat"
        shell:
              	"""
                {gtfToGenePred} -genePredExt -geneNameAsName2 {input} refFlat.tmp
                paste <(cut -f 12 refFlat.tmp) <(cut -f 1-10 refFlat.tmp) > {output}
                rm refFlat.tmp
                """

########################################################################################################
#fastqc to sam
"""Create an empty bam file linking cell/UMI barcodes to reads"""
########################################################################################################
rule fastq_to_sam:
	input:
		r1=config['DATADIR']+'/{sample}_R1.fastq.gz',
		r2=config['DATADIR']+'/{sample}_R2.fastq.gz'
	output: temp('{path}/{sample}_unaligned.bam')
	threads: CORES
	shell:
		"""java -Dpicard.useLegacyParser=false -Djava.io.tmpdir={TMPDIR} -jar {PICARD} FastqToSam\
		-F1 {input.r1}\
		-F2 {input.r2}\
		-SM DS -O {output}"""

########################################################################################################
#Pre-alignment
########################################################################################################


rule stage1:
	input: '{path}/{sample}_unaligned.bam'
	output: temp('{path}/{sample}_tagged_unmapped.bam')
	params:
		BC_summary = '{sample}_CELL_barcode.txt',
		UMI_summary = '{sample}_UMI_barcode.txt',
		start_trim = '{sample}_start_trim.txt',
		polyA_trim = '{sample}_polyA_trim.txt',
		BC_range_1 = config['GLOBAL']['BC_range']['first'],
		BC_range_2 = config['GLOBAL']['BC_range']['last'],
		UMI_range_1 = config['GLOBAL']['UMI_range']['first'],
		UMI_range_2 = config['GLOBAL']['UMI_range']['last'],
		sample = '{sample}'
	threads: CORES
	shell:
		"""
		mkdir -p {SAMPLEWDIR}/{params.sample}
		{DROPSEQ}/TagBamWithReadSequenceExtended\
		SUMMARY={SAMPLEWDIR}/{params.sample}/{params.BC_summary}\
		BASE_RANGE={params.BC_range_1}-{params.BC_range_2}\
		BASE_QUALITY=10\
		BARCODED_READ=1\
		DISCARD_READ=false\
		TAG_NAME=XC\
		NUM_BASES_BELOW_QUALITY=1\
		INPUT={input}\
		OUTPUT=/dev/stdout COMPRESSION_LEVEL=0|\
		\
		{DROPSEQ}/TagBamWithReadSequenceExtended\
		SUMMARY={SAMPLEWDIR}/{params.sample}/{params.UMI_summary}\
		BASE_RANGE={params.UMI_range_1}-{params.UMI_range_2}\
		BASE_QUALITY=10\
		BARCODED_READ=1\
		DISCARD_READ=true\
		TAG_NAME=XM\
		NUM_BASES_BELOW_QUALITY=1\
		INPUT=/dev/stdin\
		OUTPUT=/dev/stdout COMPRESSION_LEVEL=0|\
		\
		{DROPSEQ}/FilterBam TAG_REJECT=XQ\
		INPUT=/dev/stdin\
		OUTPUT=/dev/stdout COMPRESSION_LEVEL=0|\
		\
		{DROPSEQ}/TrimStartingSequence\
		OUTPUT_SUMMARY={SAMPLEWDIR}/{params.sample}/{params.start_trim}\
		SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG\
		MISMATCHES=0\
		NUM_BASES=5\
		INPUT=/dev/stdin\
		OUTPUT=/dev/stdout COMPRESSION_LEVEL=0|\
		\
		{DROPSEQ}/PolyATrimmer\
		OUTPUT_SUMMARY={SAMPLEWDIR}/{params.sample}/{params.polyA_trim}\
		MISMATCHES=0\
		NUM_BASES=6\
		OUTPUT={output}\
		INPUT=/dev/stdin
		"""

########################################################################################################
# sam to fastq
########################################################################################################
rule sam_to_fastq:
	input: '{path}/{sample}_tagged_unmapped.bam'
	output: temp('{path}/{sample}_tagged_unmapped.fastq')
	shell:
		"""java -Dpicard.useLegacyParser=false -Xmx500m -jar -Djava.io.tmpdir={TMPDIR}	{PICARD} SamToFastq\
		-INPUT {input}\
		-FASTQ {output}"""

########################################################################################################
# Align the data with star
########################################################################################################
#--genomeLoad LoadAndKeep\   this part used to be here instead of "--genomeLoad NoSharedMemory\"

rule STAR_align:
	input:  fastq='{path}/{sample}_tagged_unmapped.fastq',
		STAR_ind='STAR_ind_noAnno'
	output: temp('{path}/{sample}_STAR_Aligned.out.sam')
	params:
		prefix = '{sample}_STAR_',
		mismatch = MISMATCH,
		sample = '{sample}'
	threads: CORES
	shell:"""{STAREXEC}\
			--genomeDir {input.STAR_ind}\
			--runThreadN {CORES}\
			--outFilterMismatchNmax {params.mismatch}\
			--readFilesIn {input.fastq}\
			--genomeLoad NoSharedMemory\
			--outFileNamePrefix {wildcards.path}/{params.prefix} || true"""

########################################################################################################
# Post-alignment
########################################################################################################
rule sort:
	input:
		samples = '{path}/{sample}_STAR_Aligned.out.sam'
	output: '{path}/{sample}_Aligned_sorted_2.bam'
	threads: CORES
	shell:
		"""java	-Dpicard.useLegacyParser=false -Djava.io.tmpdir={TMPDIR} -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar {PICARD} SortSam\
		-INPUT {input}\
		-OUTPUT {output}\
		-SORT_ORDER queryname\
		-TMP_DIR {TMPDIR}"""

rule sortCoord:
	input:
		samples = '{path}/{sample}_STAR_Aligned.out.sam'
	output: '{path}/{sample}_Aligned_coordSort.bam'
	threads: CORES
	shell:
		"""java	-Dpicard.useLegacyParser=false -Djava.io.tmpdir={TMPDIR} -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar {PICARD} SortSam\
		-INPUT {input}\
		-OUTPUT {output}\
		-SORT_ORDER coordinate\
		-TMP_DIR {TMPDIR}"""

rule indexBamFiles:
	input:
		samples='{path}/{sample}_Aligned_coordSort.bam'
	output: '{path}/{sample}_Aligned_coordSort.bam.bai'
	threads: CORES
	shell:
		"""samtools index {input}"""

rule mergeBamFiles:
	input: 
		samples=expand('{PIPELINE_MAJOR}/{sample}/{sample}_Aligned_sorted_2.bam',PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	output: 'combined_bam.bam'
	shell:
		"""samtools merge {output} {input}"""

# if combined file is too big, need to downsample
#rule# downSampleBam:
#	input: "combined_bam.bam"
#	output: "combined_bam.bam"
#	shell: """samtools view -s {downSamp} -b {input} > {output}
		

rule calcHMMbed:
	input: "combined_bam.bam"
	output: "TAR_reads.bed.gz"
	threads: CORES
	shell:
		"""
		bash scripts/SingleCellHMM_MW.bash {input} {CORES} {MERGEBP} {THRESH}
		"""

rule calcHMMrefFlat:
	input:	"TAR_reads.bed.gz"
	output: 'TAR_reads.bed.gz.noDir.refFlat.refFlat',
		'TAR_reads.bed.gz.withDir.refFlat.refFlat'
	shell:
		"""
		Rscript scripts/generate_refFlat_script_both.R refFlat.refFlat {input}
		"""

rule MergeBamAlignment:
	input:  unmapped = '{path}/{sample}_tagged_unmapped.bam', 
		mapped = '{path}/{sample}_Aligned_sorted_2.bam',
		dictFile = config['GENOMEREF']+'.dict'
	output: temp('{path}/{sample}_merged.bam')
        threads: CORES
        shell:
                """java -Dpicard.useLegacyParser=false -Djava.io.tmpdir={TMPDIR} -Xmx4000m -jar {PICARD} MergeBamAlignment\
                -REFERENCE_SEQUENCE {GENOMEREF}\
                -UNMAPPED_BAM {input.unmapped}\
              	-ALIGNED_BAM {input.mapped}\
                -INCLUDE_SECONDARY_ALIGNMENTS false\
                -PAIRED_RUN false\
                -OUTPUT {output}
		"""

rule stage3:
	input:  merged = '{path}/{sample}_merged.bam', reference = 'refFlat.refFlat'
        output: temp('{path}/{sample}_gene_exon_tagged.bam')
        threads: CORES
        shell:
                """
		{DROPSEQ}/TagReadWithGeneFunction\
                O={output}\
                I={input.merged}\
                ANNOTATIONS_FILE={input.reference}\
                #TAG=GE\
                CREATE_INDEX=true
                """
		
rule stage3_withDir:
	input:	merged = '{path}/{sample}_merged.bam',
		reference ='TAR_reads.bed.gz.withDir.refFlat.refFlat'
	output: temp('{path}/{sample}_TAR_tagged_withDir.bam')
	threads: CORES
        shell:
                """
                {DROPSEQ}/TagReadWithGeneFunction\
                O={output}\
                I={input.merged}\
                ANNOTATIONS_FILE={input.reference}\
                #TAG=GE\
                CREATE_INDEX=true
                """
rule stage3_noDir:
        input:  merged = '{path}/{sample}_merged.bam',           
	        reference = 'TAR_reads.bed.gz.noDir.refFlat.refFlat'
	output: temp('{path}/{sample}_TAR_tagged_noDir.bam')
	threads: CORES
        shell:
                """
                {DROPSEQ}/TagReadWithGeneFunction\
                O={output}\
                I={input.merged}\
                ANNOTATIONS_FILE={input.reference}\
                #TAG=GE\
                CREATE_INDEX=true
                """


########################################################################################################
# Extract expression for single species
########################################################################################################
rule extract_expression:
	input: 
		bam='{path}/{sample}_gene_exon_tagged.bam'
	output: dge='{path}/{sample}_gene_expression_matrix.txt.gz',
		cellSummary='{path}/{sample}_gene_dge.summary.txt'
	params:
		sample = '{sample}',
		numbarcodes = config["expectedCells"]
	shell:
		"""{DROPSEQ}/DigitalExpression\
		I={input.bam}\
		O={output.dge}\
		SUMMARY={output.cellSummary} \
		NUM_CORE_BARCODES={params.numbarcodes}"""

rule getCellsList:
	input:	'{path}/{sample}_gene_dge.summary.txt'
	output: '{path}/{sample}_cellList.txt'
	shell:
		"""
		sed '/^#/ d' < {input} > {input}.temp.txt
		tail -n +3 {input}.temp.txt > {input}.temp2.txt
		cut -f1 {input}.temp2.txt > {output}
		rm {input}.temp.txt
		rm {input}.temp2.txt
		"""

rule extract_HMM_expression_withDir:
	input: 
		bam='{path}/{sample}_TAR_tagged_withDir.bam',
		barcodes='{path}/{sample}_cellList.txt'
	output: '{path}/{sample}_TAR_expression_matrix_withDir.txt.gz'
	params:
		sample = '{sample}',
		numbarcodes = config["expectedCells"]
	shell:
		"""{DROPSEQ}/DigitalExpression\
		I={input.bam}\
		O={output}\
		SUMMARY={SAMPLEWDIR}/{params.sample}/{params.sample}_TAR_withDir_dge.summary.txt \
		CELL_BC_FILE={input.barcodes}"""
		
rule extract_HMM_expression_noDir:
	input: 
		bam='{path}/{sample}_TAR_tagged_noDir.bam',
		barcodes='{path}/{sample}_cellList.txt'
	output: '{path}/{sample}_TAR_expression_matrix_noDir.txt.gz'
	params:
		sample = '{sample}',
		numbarcodes = config["expectedCells"]
	shell:
		"""{DROPSEQ}/DigitalExpression\
		I={input.bam}\
		O={output}\
		SUMMARY={SAMPLEWDIR}/{params.sample}/{params.sample}_TAR_noDir_dge.summary.txt \
		CELL_BC_FILE={input.barcodes}"""


# generate differentially expressed genes and uTARs
rule getDiffMarkers:
	input:
		geneFile='{path}/{sample}_gene_expression_matrix.txt.gz',
		hmmFile='{path}/{sample}_TAR_expression_matrix_noDir.txt.gz'
	output: '{path}/{sample}_diffMarkers.txt'
	shell:
		"""
		Rscript scripts/analyzeExpressionMat.R {input.geneFile} {input.hmmFile}
		"""

# from diff uTARs, extract fasta region to blast
rule getDiffSeqsToBlast:
	input: 
		diffMarkers='{path}/{sample}_diffMarkers.txt',
		bamFile='{path}/{sample}_Aligned_coordSort.bam',
		bamIndex='{path}/{sample}_Aligned_coordSort.bam.bai'
	output:	fastaRegions='{path}/{sample}_seqsToBlast.txt',
		uTARDiffMarkers='{path}/{sample}_diffuTARMarkers.txt'
	shell:
		"""
		Rscript scripts/getFastasForBlast.R {input.diffMarkers} {input.bamFile}
		"""

# extract diff uTAR fastas
rule getDiffSeqsToBlastFa:
	input:
		fastaHeaders="{path}/{sample}_seqsToBlast.txt",
		fastaFile=config['GENOMEREF'],
		fastaIndex=config['GENOMEREF']+'.dict'
	output:
		"{path}/{sample}_seqsToBlast.fa"
	shell:
		"""
		samtools faidx -r {input.fastaHeaders} {input.fastaFile} > {output}
		"""

# run blast on fastas
rule ruleBlast:
	input:
		fastaFile="{path}/{sample}_seqsToBlast.fa",
		blastDB=config['BLASTDB']
	output:
		"{path}/{sample}_blastResults.txt"
	threads: CORES
	shell:
		"""
		blastn -db {input.blastDB}/nt -query {input.fastaFile} -out {output} -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 5 -num_threads {CORES}
		"""

# label uTARs based on best blast results
rule labelDiffuTARs:
        input:
              	diffMarkers='{path}/{sample}_diffuTARMarkers.txt',
                blastResult='{path}/{sample}_blastResults.txt'
        output: uTARSummary='{path}/{sample}_diffuTARMarkersLabeled.txt',
        shell:
                """
                Rscript scripts/examineBlastResults.R {input.diffMarkers} {input.blastResult}
                """

rule getMatsSteps:
	input: 	gene='{path}/{sample}_gene_expression_matrix.txt.gz',
		hmm1='{path}/{sample}_TAR_expression_matrix_withDir.txt.gz',
		hmm2='{path}/{sample}_TAR_expression_matrix_noDir.txt.gz',
		bamFile='{path}/{sample}_Aligned_coordSort.bam',
		baiFile='{path}/{sample}_Aligned_coordSort.bam.bai'

	output: '{path}/{sample}_getMats.txt'
	shell:
		"""echo "Expression matrices ready" > {output}"""
