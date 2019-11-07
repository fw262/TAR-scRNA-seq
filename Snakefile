#################################################
# Michael Wang			
# Drop-Seq Snakemake pipeline with de-novo generation of features integration			
# Adapted from the Dropseq pipeline by Hoohm:		
# https://github.com/Hoohm/dropseqpipe				
#################################################
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
MINCOV=str(config['MINCOV'])
########################################################################################################
# Executables
########################################################################################################
TMPDIR = config['TMPDIR']
MIXCR = config['MIXCR']
PICARD = config['PICARD']
DROPSEQ = config['DROPSEQ']
FASTQC = config['FASTQCEXEC']
STAREXEC = config['STAREXEC']
CORES = config['CORES']
gtfToGenePred=config['GTFTOGENEPRED']
SingleCellHMM=config['SingleCellHMM']
generate_refFlat_script=config['generate_refFlat_script']
Rscript=config['Rscript']

rule all:
############## original default call here
	input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_finish.txt', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_expression_matrix.txt.gz', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])

#####################################################################################
# create STAR index
rule generateStar:
        input:
              	fastaFile=config['GENOMEREF']
        output: directory('STAR_ind_noAnno')
        threads: 10
        shell:
              	"""
                mkdir STAR_ind_noAnno
                \
                {STAREXEC}\
                        --runThreadN 10\
                        --runMode genomeGenerate\
                        --genomeDir STAR_ind_noAnno\
                        --genomeFastaFiles {input.fastaFile}
                """
#####################################################################################

#####################################################################################
# convert GTF to REFFlat
rule convertToRefFlat1:
        input:  {gtffile}
        output: "refFlat.txt.gz"
        shell:
              	"""
                {gtfToGenePred} -genePredExt -geneNameAsName2 {input} refFlat.tmp
                paste <(cut -f 12 refFlat.tmp) <(cut -f 1-10 refFlat.tmp) > {output}
                rm refFlat.tmp
                """
				
rule unzipRefFlat:
	input: "refFlat.txt.gz"
	output: "refFlat.refFlat"
	shell:
		"""mv {input} {output}"""

#####################################################################################

#####################################################################################
# create sequence dictionary using picard tools
rule createGenomeDict:
	input: config['GENOMEREF']
	output: config['GENOMEREF']+'.dict'
	shell:
		"""
		java -jar {PICARD} CreateSequenceDictionary \
			R={input} \
			O={output}
		"""
#####################################################################################

########################################################################################################
#Create fastqc report, also create fastqc.pdf, with plots of key metrics
########################################################################################################
rule fastqc:
	input:
		r1='{DATADIR}/{sample}_R1.fastq.gz', 
		r2='{DATADIR}/{sample}_R2.fastq.gz'
	output: '{path}/{sample}_R1_fastqc.html'
	threads: 4
	shell:
		"{FASTQC} {input.r1} {input.r2} -t 2 -o {wildcards.path} --extract ; "
		# "Rscript R/fastqc.R {wildcards.path} "

########################################################################################################
#fastqc to sam for R1 only
"""Create an empty bam file linking cell/UMI barcodes to reads"""
########################################################################################################
rule fastq_to_sam_R1:
	input:
		r1='{DATADIR}/{sample}_R1.fastq.gz',
		fastqc='{path}/{sample}_R1_fastqc.html'
	output: '{path}/{sample}_unaligned_r1.bam'
	threads: 4
	shell:
		"""java -Djava.io.tmpdir={TMPDIR} -jar {PICARD} FastqToSam\
		F1={input.r1}\
		SM=DS O={output}"""

########################################################################################################
#fastqc to sam for R2 only
"""Create an empty bam file linking cell/UMI barcodes to reads"""
########################################################################################################
rule fastq_to_sam_R2:
	input:
		r2='{DATADIR}/{sample}_R2.fastq.gz',
		fastqc='{path}/{sample}_R1_fastqc.html'
	output: '{path}/{sample}_unaligned_r2.bam'
	threads: 4
	shell:
		"""java -Djava.io.tmpdir={TMPDIR} -jar {PICARD} FastqToSam\
		F1={input.r2}\
		SM=DS O={output}"""


########################################################################################################
#fastqc to sam
"""Create an empty bam file linking cell/UMI barcodes to reads"""
########################################################################################################
rule fastq_to_sam:
	input:
		r1=config['DATADIR']+'/{sample}_R1.fastq.gz',
		r2=config['DATADIR']+'/{sample}_R2.fastq.gz'
		#fastqc='{path}/{sample}_R1_fastqc.html'
	output: '{path}/{sample}_unaligned.bam'
	threads: 4
	shell:
		"""java -Djava.io.tmpdir={TMPDIR} -jar {PICARD} FastqToSam\
		F1={input.r1}\
		F2={input.r2}\
		SM=DS O={output}"""

########################################################################################################
#Pre-alignment
########################################################################################################


rule stage1:
	input: '{path}/{sample}_unaligned.bam'
	output: '{path}/{sample}_tagged_unmapped.bam'
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
	threads: 4
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
		OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 |\
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
		{DROPSEQ}/FilterBAM TAG_REJECT=XQ\
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
	output: '{path}/{sample}_tagged_unmapped.fastq'
	shell:
		"""java -Xmx500m -jar -Djava.io.tmpdir={TMPDIR}	{PICARD} SamToFastq\
		INPUT={input}\
		FASTQ={output}"""

########################################################################################################
# Align the data with star
########################################################################################################
#--genomeLoad LoadAndKeep\   this part used to be here instead of "--genomeLoad NoSharedMemory\"

rule STAR_align:
	input:  fastq='{path}/{sample}_tagged_unmapped.fastq',
		STAR_ind='STAR_ind_noAnno'
	output: '{path}/{sample}_STAR_Aligned.out.sam'
	params:
		prefix = '{sample}_STAR_',
		mismatch = MISMATCH,
		sample = '{sample}'
	threads: CORES
	shell:"""{STAREXEC}\
			--genomeDir {input.STAR_ind}\
			--runThreadN {CORES}\
			--outFilterMismatchNmax={params.mismatch}\
			--readFilesIn {input.fastq}\
			--genomeLoad NoSharedMemory\
			--outFileNamePrefix {wildcards.path}/{params.prefix} || true"""

rule mv_sam:
	input: '{path}/{sample}_STAR_Aligned.out.sam'
	output: '{path}/{sample}_Aligned.sam'
	shell:"""cp {input} {output};"""

########################################################################################################
# Post-alignment
########################################################################################################
rule sort:
	input:
		samples = '{path}/{sample}_Aligned.sam'
	output: temp('{path}/{sample}_Aligned_sorted.sam')
	threads: 4
	shell:
		"""java	-Djava.io.tmpdir={TMPDIR} -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar {PICARD} SortSam\
		INPUT={input}\
		OUTPUT={output}\
		SORT_ORDER=queryname\
		TMP_DIR={TMPDIR}"""

rule sortedSamToBam:
	input:
		samples='{path}/{sample}_Aligned_sorted.sam'
	output:temp('{path}/{sample}_Aligned_sorted.bam')
	threads:  4
	shell:
		"""samtools view -S -b -h -@ {threads} {input} > {output}"""
rule sort2:
	input:
		samples='{path}/{sample}_Aligned_sorted.bam'
	output:'{path}/{sample}_Aligned_sorted_2.bam'
	threads:  4
	shell:
		"""samtools sort {input} -o {output}"""
rule indexBamFiles:
	input:
		samples='{path}/{sample}_Aligned_sorted_2.bam'
	output:'{path}/{sample}_Aligned_sorted_2.bam.bai'
	threads: 4
	shell:
		"""samtools index {input}"""

rule mergeBamFiles:
	input: 
		samples=expand('{PIPELINE_MAJOR}/{sample}/{sample}_Aligned_sorted_2.bam',PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	output: 'combined_bam.bam'
	shell:
		"""samtools merge {output} {input}"""
		
rule calcHMMbed:
	input:
		samples='combined_bam.bam'
	output:
		out='combined_bam_HMM_features/combined_bam_merge500_'+str(config['MINCOV'])+'reads.bed.gz'
	run:
		shell("bash {SingleCellHMM} {input} {CORES} {MINCOV}")

rule calcHMMrefFlat:
	input:	'combined_bam_HMM_features/combined_bam_merge500_'+str(config['MINCOV'])+'reads.bed.gz'
	output: 'combined_bam_HMM_features/combined_bam_merge500_'+str(config['MINCOV'])+'reads.bed.gz.noDir.refFlat.refFlat',
			'combined_bam_HMM_features/combined_bam_merge500_'+str(config['MINCOV'])+'reads.bed.gz.withDir.refFlat.refFlat'
	shell:
		"""{Rscript} {generate_refFlat_script_both.R} refFlat.refFlat {input}"""


rule stage3:
	input:	unmapped = '{path}/{sample}_tagged_unmapped.bam',
			mapped = '{path}/{sample}_Aligned_sorted_2.bam',
			reference = 'refFlat.refFlat',
			dict=config['GENOMEREF']+'.dict'
	output: '{path}/{sample}_gene_exon_tagged.bam'
	threads: 4
	shell:
		"""java -Djava.io.tmpdir={TMPDIR} -Xmx4000m -jar {PICARD} MergeBamAlignment\
		REFERENCE_SEQUENCE={GENOMEREF}\
		UNMAPPED_BAM={input.unmapped}\
		ALIGNED_BAM={input.mapped}\
		INCLUDE_SECONDARY_ALIGNMENTS=false\
		PAIRED_RUN=false\
		OUTPUT=/dev/stdout COMPRESSION_LEVEL=0|\
		\
		{DROPSEQ}/TagReadWithGeneExon\
		O={output}\
		I=/dev/stdin\
		ANNOTATIONS_FILE={input.reference}\
		TAG=GE\
		CREATE_INDEX=true
		"""
		
rule stage3_withDir:
	input:	unmapped = '{path}/{sample}_tagged_unmapped.bam',
			mapped = '{path}/{sample}_Aligned_sorted_2.bam',
			reference='combined_bam_HMM_features/combined_bam_merge500_'+str(config['MINCOV'])+'reads.bed.gz.withDir.refFlat.refFlat'
	output: '{path}/{sample}_HMM_tagged_withDir.bam'
	threads: 4
	shell:
		"""java -Djava.io.tmpdir={TMPDIR} -Xmx4000m -jar {PICARD} MergeBamAlignment\
		REFERENCE_SEQUENCE={GENOMEREF}\
		UNMAPPED_BAM={input.unmapped}\
		ALIGNED_BAM={input.mapped}\
		INCLUDE_SECONDARY_ALIGNMENTS=false\
		PAIRED_RUN=false\
		OUTPUT=/dev/stdout COMPRESSION_LEVEL=0|\
		\
		{DROPSEQ}/TagReadWithGeneExon\
		O={output}\
		I=/dev/stdin\
		ANNOTATIONS_FILE={input.reference}\
		TAG=GE\
		CREATE_INDEX=true
		"""

rule stage3_noDir:
	input:	unmapped = '{path}/{sample}_tagged_unmapped.bam',
			mapped = '{path}/{sample}_Aligned_sorted_2.bam',
			reference='combined_bam_HMM_features/combined_bam_merge500_'+str(config['MINCOV'])+'reads.bed.gz.noDir.refFlat.refFlat'
	output: '{path}/{sample}_HMM_tagged_noDir.bam'
	threads: 4
	shell:
		"""java -Djava.io.tmpdir={TMPDIR} -Xmx4000m -jar {PICARD} MergeBamAlignment\
		REFERENCE_SEQUENCE={GENOMEREF}\
		UNMAPPED_BAM={input.unmapped}\
		ALIGNED_BAM={input.mapped}\
		INCLUDE_SECONDARY_ALIGNMENTS=false\
		PAIRED_RUN=false\
		OUTPUT=/dev/stdout COMPRESSION_LEVEL=0|\
		\
		{DROPSEQ}/TagReadWithGeneExon\
		O={output}\
		I=/dev/stdin\
		ANNOTATIONS_FILE={input.reference}\
		TAG=GE\
		CREATE_INDEX=true
		"""

########################################################################################################
# Extract expression for single species
########################################################################################################
rule extract_expression:
	input: 
		bam='{path}/{sample}_gene_exon_tagged.bam'
	output: '{path}/{sample}_expression_matrix.txt.gz'
	params:
		sample = '{sample}',
		numbarcodes = config["expectedCells"]
	shell:
		"""{DROPSEQ}/DigitalExpression\
		I={input.bam}\
		O={output}\
		SUMMARY={SAMPLEWDIR}/{params.sample}/{params.sample}_dge.summary.txt \
		NUM_CORE_BARCODES={params.numbarcodes}"""

rule extract_HMM_expression_withDir:
	input: 
		bam='{path}/{sample}_HMM_tagged_withDir.bam'
	output: '{path}/{sample}_expression_matrix_HMM_withDir.txt.gz'
	params:
		sample = '{sample}',
		numbarcodes = config["expectedCells"]
	shell:
		"""{DROPSEQ}/DigitalExpression\
		I={input.bam}\
		O={output}\
		SUMMARY={SAMPLEWDIR}/{params.sample}/{params.sample}_dge.summary.txt \
		NUM_CORE_BARCODES={params.numbarcodes}"""
		
rule extract_HMM_expression_noDir:
	input: 
		bam='{path}/{sample}_HMM_tagged_noDir.bam'
	output: '{path}/{sample}_expression_matrix_HMM_noDir.txt.gz'
	params:
		sample = '{sample}',
		numbarcodes = config["expectedCells"]
	shell:
		"""{DROPSEQ}/DigitalExpression\
		I={input.bam}\
		O={output}\
		SUMMARY={SAMPLEWDIR}/{params.sample}/{params.sample}_dge.summary.txt \
		NUM_CORE_BARCODES={params.numbarcodes}"""

		
rule dummyOut:
	input: gene='{path}/{sample}_expression_matrix.txt.gz',
	       hmm1='{path}/{sample}_expression_matrix_HMM_withDir.txt.gz',
		   hmm2='{path}/{sample}_expression_matrix_HMM_noDir.txt.gz'

	output: '{path}/{sample}_finish.txt'
	shell:
		"""echo "Expression matrices ready" > {output}"""