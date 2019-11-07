#UNCOMMENT BELOW BEFORE EXITING
# MAKE SURE CHROM IN REFERENCE AND HMM BED MATCH EACH OTHER
args=(commandArgs(TRUE))
inputRefFlat<-args[1]
HMMbedFile<-args[2]

#inputRefFlat<-"/workdir/References/Chick/GRCg6a/refFlat_ensembl.refFlat.filled.refFlat"
#HMMbedFile<-"/workdir/fw262/ShaoPei/chickenAll/all_Aligned_sorted_2_2019-07-01_11-29-22_MERGE500_MINREAD5/all_Aligned_sorted_2_merge500_5reads.bed"

#inputRefFlat<-"/workdir/fw262/references/equCab2_ucsc.refFlat"
#HMMbedFile<-"/workdir/fw262/ShaoPei/horseData/results_horse/8822_7858_62538_HGJT3BGX3_H2-Mock_TTCTGCCT/8822_7858_62538_HGJT3BGX3_H2-Mock_TTCTGCCT_Aligned_sorted_2_merge500_5reads.bed"
refName<-unlist(strsplit(inputRefFlat,"/"))[length(unlist(strsplit(inputRefFlat,"/")))]
# read in hg_38 ref genes
gene_ref <- read.delim(inputRefFlat, header=F, comment.char="#")
gene_ref_bare<-gene_ref[,c("V1","V3","V4","V5","V6")] ####### edit this based on file type
colnames(gene_ref_bare)<-c("gene","chr","direction","start","end")
gene_ref_bare$gene<-as.character(gene_ref_bare$gene)
gene_ref_bare$chr<-as.character(gene_ref_bare$chr)
gene_ref_bare$direction<-as.character(gene_ref_bare$direction)
gene_ref_bare$start<-as.numeric(gene_ref_bare$start)
gene_ref_bare$end<-as.numeric(gene_ref_bare$end)


# read in Shao-Pei's file # 500 bp blocks merge, at least 5 reads
input<-HMMbedFile
HMManno <- read.delim(input, header=FALSE)
HMManno_bare<-HMManno[,c("V1","V2","V3","V6")]
colnames(HMManno_bare)<-c("chr","start","end","direction")
HMManno_bare$chr<-as.character(HMManno_bare$chr)
HMManno_bare$direction<-as.character(HMManno_bare$direction)
HMManno_bare$start<-as.numeric(HMManno_bare$start)
HMManno_bare$end<-as.numeric(HMManno_bare$end)

outFile=paste0(input,".noDir.",refName)

########################################################
source("/workdir/fw262/ShaoPei/generate_refFlat_func.R")

df<-HMManno_bare
#HMManno_bare_sample<-df[sample(nrow(df),10),]
#HMManno$inAnno<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExist,exon_ref=gene_gtf_bare,gene_ref=gene_ref_bare)
HMManno$inGene<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene_noDir,gene_ref=gene_ref_bare)


########################### uncomment below to make refFlat file
# make dataframe into refFlat file format
HMManno$geneName<-paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,"_",HMManno$V6,"_",HMManno$V7,"_",HMManno$inGene,sep="")
#HMManno$name<-HMManno$geneName
HMManno$name<-paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,sep="")
HMManno$chrom<-HMManno$V1
HMManno$strand<-HMManno$V6
HMManno$txStart<-HMManno$V2
HMManno$txEnd<-HMManno$V3
HMManno$cdsStart<-HMManno$V2
HMManno$cdsEnd<-HMManno$V3
HMManno$exonCount<-1
HMManno$exonStarts<-HMManno$V2
HMManno$exonEnds<-HMManno$V3


HMMannoReady<-HMManno[,c("geneName","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")]

# export as GTF
write.table(HMMannoReady,outFile,sep="\t",row.names=F,col.names = F, quote=F)

###################################################################################
#### include direction below
###################################################################################
outFile=paste0(input,".withDir.",refName)

########################################################
source("/workdir/fw262/ShaoPei/generate_refFlat_func.R")

df<-HMManno_bare
#HMManno_bare_sample<-df[sample(nrow(df),10),]
#HMManno$inAnno<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExist,exon_ref=gene_gtf_bare,gene_ref=gene_ref_bare)
HMManno$inGene<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene2,gene_ref=gene_ref_bare)


########################### uncomment below to make refFlat file
# make dataframe into refFlat file format
HMManno$geneName<-paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,"_",HMManno$V6,"_",HMManno$V7,"_",HMManno$inGene,sep="")
#HMManno$name<-HMManno$geneName
HMManno$name<-paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,sep="")
HMManno$chrom<-HMManno$V1
HMManno$strand<-HMManno$V6
HMManno$txStart<-HMManno$V2
HMManno$txEnd<-HMManno$V3
HMManno$cdsStart<-HMManno$V2
HMManno$cdsEnd<-HMManno$V3
HMManno$exonCount<-1
HMManno$exonStarts<-HMManno$V2
HMManno$exonEnds<-HMManno$V3


HMMannoReady<-HMManno[,c("geneName","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")]

# export as GTF
write.table(HMMannoReady,outFile,sep="\t",row.names=F,col.names = F, quote=F)