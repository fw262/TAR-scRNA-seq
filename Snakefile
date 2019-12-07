#################################################
# Michael Wang			#
# TAR scRNA-seq analysis			#
# last edited Dec. 6, 2019			#
#################################################
import pdb

########################################################################################################
# Configfile
########################################################################################################
# configfile:'config.B1-25.yaml'
configfile:'config.yaml'
########################################################################################################
# Variables and references
########################################################################################################
SAMPLEWDIR = config['SAMPWDIR']
MISMATCH = config['GLOBAL']['allowed_mismatch']
DATADIR = config['DATADIR']
GENOMEREF = config['GENOMEREF']
#REFFLAT = config['REFFLAT']
#REFFLAT_HMM = config['REFFLAT_HMM']
#REFFLAT_HMM_noDir = config['REFFLAT_HMM_noDir']
#REFFLAT_HMM_withDir = config['REFFLAT_HMM_withDir']
#METAREF = config['METAREF']
#RRNAINTERVALS=config['RRNAINTERVALS']
TMPDIR = config['TMPDIR']
gtffile=config['REFGTF']
MINCOV=str(config['MINCOV'])
MERGEBP=str(config['MERGEBP'])
########################################################################################################
# Executables
########################################################################################################
TMPDIR = config['TMPDIR']
#MIXCR = config['MIXCR']
PICARD = config['PICARD']
DROPSEQ = config['DROPSEQ']
FASTQC = config['FASTQCEXEC']
STAREXEC = config['STAREXEC']
CORES = config['CORES']
gtfToGenePred=config['GTFTOGENEPRED']

rule all:
############## original default call here
	input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_finish.txt', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_expression_matrix.txt.gz', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_Aligned_sorted_2.bam.bai', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_expression_matrix.txt.gz', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}.tagged.Bcell.numerated.cell.list', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples']) #### this works until the end
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}.allcells.descriptions.txt', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_STAR_Aligned.out.sam', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])
	#input: expand('{PIPELINE_MAJOR}/{sample}/{sample}_unaligned_r2.bam', PIPELINE_MAJOR=config['PIPELINE_MAJOR'], sample=config['Samples'])

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
		out='combined_bam_HMM_features/combined_bam_merge'+str(config['MERGEBP'])+'_'+str(config['MINCOV'])+'reads.bed.gz'
	run:
		shell("bash scripts/SingleCellHMM_2.bash {input} {CORES} {MINCOV} {MERGEBP}")

rule calcHMMrefFlat:
	input:	'combined_bam_HMM_features/combined_bam_merge'+str(config['MERGEBP'])+'_'+str(config['MINCOV'])+'reads.bed.gz'
	output: 'combined_bam_HMM_features/combined_bam_merge'+str(config['MERGEBP'])+'_'+str(config['MINCOV'])+'reads.bed.gz.noDir.refFlat.refFlat',
			'combined_bam_HMM_features/combined_bam_merge'+str(config['MERGEBP'])+'_'+str(config['MINCOV'])+'reads.bed.gz.withDir.refFlat.refFlat'
	shell:
		"""Rscript scripts/generate_refFlat_script_both.R refFlat.refFlat {input}"""

rule MergeBamAlignment:
	input:  unmapped = '{path}/{sample}_tagged_unmapped.bam', mapped = '{path}/{sample}_Aligned_sorted_2.bam'
        output: '{path}/{sample}_merged.bam'
        threads: 4
        shell:
                """java -Djava.io.tmpdir={TMPDIR} -Xmx4000m -jar {PICARD} MergeBamAlignment\
                REFERENCE_SEQUENCE={GENOMEREF}\
                UNMAPPED_BAM={input.unmapped}\
              	ALIGNED_BAM={input.mapped}\
                INCLUDE_SECONDARY_ALIGNMENTS=false\
                PAIRED_RUN=false\
                OUTPUT={output}
		"""

rule stage3:
	input:  merged = '{path}/{sample}_merged.bam', reference = 'refFlat.refFlat'
        output: '{path}/{sample}_gene_exon_tagged.bam'
        threads: 4
        shell:
                """
		{DROPSEQ}/TagReadWithGeneExon\
                O={output}\
                I={input.merged}\
                ANNOTATIONS_FILE={input.reference}\
                TAG=GE\
                CREATE_INDEX=true
                """
		
rule stage3_withDir:
	input:	merged = '{path}/{sample}_merged.bam',
		reference ='combined_bam_HMM_features/combined_bam_merge'+str(config['MERGEBP'])+'_'+str(config['MINCOV'])+'reads.bed.gz.withDir.refFlat.refFlat'
	output: '{path}/{sample}_HMM_tagged_withDir.bam'
	threads: 4
        shell:
                """
                {DROPSEQ}/TagReadWithGeneExon\
                O={output}\
                I={input.merged}\
                ANNOTATIONS_FILE={input.reference}\
                TAG=GE\
                CREATE_INDEX=true
                """
rule stage3_noDir:
        input:  merged = '{path}/{sample}_merged.bam',           
	        reference = 'combined_bam_HMM_features/combined_bam_merge'+str(config['MERGEBP'])+'_'+str(config['MINCOV'])+'reads.bed.gz.noDir.refFlat.refFlat'
	output: '{path}/{sample}_HMM_tagged_noDir.bam'
	threads: 4
        shell:
                """
                {DROPSEQ}/TagReadWithGeneExon\
                O={output}\
                I={input.merged}\
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


rule gunzip:
	input: '{path}/{sample}_expression_matrix.txt.gz'
	output: '{path}/{sample}_expression_matrix.txt'
	shell:
		"""gunzip -qf {input}"""

########################################################################################################
# Generate report
########################################################################################################
rule report:
	input:
		'{path}/{sample}_expression_matrix.txt'
	output:
		'{path}/{sample}.report.html'
	run:
		from snakemake.utils import report
		
		report("""
		A test dropseq pipeline
		""", output[0])


# Get Bcell related Cells and IG genes
########################################################################################################
rule Bcell_list:
	input: '{path}/{sample}_expression_matrix.txt'
	output: barcode = '{path}/{sample}.Bcell.barcodes.txt'
	shell: 
		"""
			 Rscript bin/Rscripts/DropSeq_analysis.barebones.justB.R {wildcards.sample} TRUE 100 ;
		"""

rule Bcell_raw_tags:
	input:	bam = '{path}/{sample}_tagged_unmapped.bam',
		barcode = '{path}/{sample}.Bcell.barcodes.txt'
	output: sam = temp('{path}/{sample}.tagged.Bcell.sam'),
		bam  = temp('{path}/{sample}.tagged.Bcell.bam'),
		fastq = '{path}/{sample}.tagged.Bcell.fastq',
		readid = '{path}/{sample}.tagged.Bcell.reads' 
	threads: 2
	shell:
		"""
			samtools view {input.bam} | LC_ALL=C grep -wFf {input.barcode} \
				 > {output.sam} ;
			awk '{{ print $1":"$12":"$14"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' {output.sam} \
				| samtools view -bS - > {output.bam} ;
			bamToFastq -i {output.bam} -fq {output.fastq} ;
			awk 'NR%4==1' {output.fastq} > {output.readid} ; 
		"""

rule allcell_raw_tags:
	input:  bam = '{path}/{sample}_tagged_unmapped.bam',
		celldesc = '{path}/{sample}.allcells.descriptions.txt'
	output: barcode = temp('{path}/{sample}.allcell.barcodes.txt'),
		sam = temp('{path}/{sample}.tagged.allcell.sam'),
		bam  = temp('{path}/{sample}.tagged.allcell.bam'),
		fastq = '{path}/{sample}.tagged.allcell.fastq',
		readid = '{path}/{sample}.tagged.allcell.reads'
	threads: 2
	shell:
		"""
			cut -f 1 {input.celldesc} > {output.barcode} ;
			samtools view {input.bam} | LC_ALL=C grep -wFf {output.barcode} \
				> {output.sam} ;
			awk '{{ print $1":"$12":"$14"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' {output.sam} \
				| samtools view -bS - > {output.bam} ;
			bamToFastq -i {output.bam} -fq {output.fastq} ;
			awk 'NR%4==1' {output.fastq} > {output.readid} ;
		"""



rule ig_process:
	input: '{path}/{sample}.tagged.Bcell.fastq'
	output:	vdj1 = '{path}/{sample}.tagged.Bcell.vdjca',
		rep = '{path}/{sample}.tagged.Bcell.reportTxt',
		rescue1 = '{path}/{sample}.tagged.Bcell.rescued1.vdjca',
		rescue2 = '{path}/{sample}.tagged.Bcell.rescued2.vdjca',
		clones = '{path}/{sample}.tagged.Bcell.clns',
		clnstxt = '{path}/{sample}.tagged.Bcell.clones.txt',
		trb = '{path}/{sample}.tagged.Bcell.clones.TRB.txt',
		igh = '{path}/{sample}.tagged.Bcell.clones.IGH.txt',
		alnAll = '{path}/{sample}.tagged.Bcell.alignmentsAll.txt',
		alns = '{path}/{sample}.tagged.Bcell.alignments.txt',
		viewalns = '{path}/{sample}.tagged.Bcell.aligned_view.txt'
	shell:
		"""
			{MIXCR} align -f -s hsa -p rna-seq -OallowPartialAlignments=true {input} {output.vdj1} > {output.rep} ;
			{MIXCR} assemblePartial -f {output.vdj1} {output.rescue1} ;
			{MIXCR} assemblePartial -f {output.rescue1} {output.rescue2} ;
			{MIXCR} assemble -f {output.rescue2} {output.clones} ;
			{MIXCR} exportClones -f {output.clones} {output.clnstxt} ;
			{MIXCR} exportClones -f -c TRB {output.clones} {output.trb} ;
			{MIXCR} exportClones -f -c IG {output.clones} {output.igh} ;
			{MIXCR} exportAlignments -f -c IG {output.rescue2} {output.alnAll} ;
			{MIXCR} exportAlignments -f -c IG -readID -vGene -vHitScore -dGene -dHitScore -jGene \
				-jHitScore -cGene -cHitScore -nMutations VRegion -aaMutations VRegion -nMutations JRegion -aaMutations JRegion -nMutations DRegion -aaMutations DRegion {output.rescue2} {output.alns} ;
			{MIXCR} exportAlignmentsPretty {output.rescue2} {output.viewalns} ;
		"""

rule ig_allcells_process:
	input: '{path}/{sample}.tagged.allcell.fastq'
	output: vdj1 = '{path}/{sample}.tagged.allcell.vdjca',
		rep = '{path}/{sample}.tagged.allcell.reportTxt',
		rescue1 = '{path}/{sample}.tagged.allcell.rescued1.vdjca',
		rescue2 = '{path}/{sample}.tagged.allcell.rescued2.vdjca',
		clones = '{path}/{sample}.tagged.allcell.clns',
		clnstxt = '{path}/{sample}.tagged.allcell.clones.txt',
		trb = '{path}/{sample}.tagged.allcell.clones.TRB.txt',
		igh = '{path}/{sample}.tagged.allcell.clones.IGH.txt',
		alnAll = '{path}/{sample}.tagged.allcell.alignmentsAll.txt',
		alns = '{path}/{sample}.tagged.allcell.alignments.txt',
		viewalns = '{path}/{sample}.tagged.allcell.aligned_view.txt'
	shell:
		"""
			{MIXCR} align -f -s hsa -p rna-seq -OallowPartialAlignments=true {input} {output.vdj1} > {output.rep} ;
			{MIXCR} assemblePartial -f {output.vdj1} {output.rescue1} ;
			{MIXCR} assemblePartial -f {output.rescue1} {output.rescue2} ;
			{MIXCR} assemble -f {output.rescue2} {output.clones} ;
			{MIXCR} exportClones -f {output.clones} {output.clnstxt} ;
			{MIXCR} exportClones -f -c TRB {output.clones} {output.trb} ;
			{MIXCR} exportClones -f -c IG {output.clones} {output.igh} ;
			{MIXCR} exportAlignments -f -c IG {output.rescue2} {output.alnAll} ;
			{MIXCR} exportAlignments -f -c IG -readID -vGene -vHitScore -dGene -dHitScore -jGene \
				-jHitScore -cGene -cHitScore -nMutations VRegion -aaMutations VRegion -nMutations JRegion -aaMutations JRegion -nMutations DRegion -aaMutations DRegion {output.rescue2} {output.alns} ;
			{MIXCR} exportAlignmentsPretty {output.rescue2} {output.viewalns} ;
		"""



rule ig_stats_Bcell:
	input: 	a='{path}/{sample}_expression_matrix.txt',
		b='{path}/{sample}.tagged.Bcell.alignments.txt'
	output: '{path}/{sample}.Bcell.IGstats.txt'
	shell: "Rscript bin/Rscripts/Mixcr_analysis.R {wildcards.sample} ; "

rule ig_stats_allcell:
	input:  a='{path}/{sample}_expression_matrix.txt',
		b='{path}/{sample}.tagged.allcell.alignments.txt'
	output: '{path}/{sample}.allcell.IGstats.txt'
	shell: "Rscript bin/Rscripts/Mixcr_analysis.allcell.R {wildcards.sample} ; "


rule colorTNSE:
	input: '{path}/{sample}.allcells.descriptions.txt'
	output: '{path}/{sample}.HCLC_color.eps'
	shell: "Rscript bin/Rscripts/Dropseq_tsne.bcell.hclcplot.R {wildcards.sample} ; "

rule enumerate_Bcell:
        input: '{path}/{sample}.tagged.Bcell.fastq'
        output: '{path}/{sample}.tagged.Bcell.numerated.cell.list'
        threads: 1
        shell:
                """
                   	awk 'NR % 4 ==1' {input} | cut -d':' -f 10 | nl > {output}
                """
