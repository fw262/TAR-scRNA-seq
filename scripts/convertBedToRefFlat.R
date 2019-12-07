#UNCOMMENT BELOW BEFORE EXITING
# MAKE SURE CHROM IN REFERENCE AND HMM BED MATCH EACH OTHER
args=(commandArgs(TRUE))

input<-args[1]
HMManno <- read.delim(input, header=FALSE)
HMManno_bare<-HMManno[,c("V1","V2","V3","V6")]
colnames(HMManno_bare)<-c("chr","start","end","direction")
HMManno_bare$chr<-as.character(HMManno_bare$chr)
HMManno_bare$direction<-as.character(HMManno_bare$direction)
HMManno_bare$start<-as.numeric(HMManno_bare$start)
HMManno_bare$end<-as.numeric(HMManno_bare$end)

outFile=paste0(input,".refFlat")


########################### uncomment below to make refFlat file
# make dataframe into refFlat file format
HMManno$geneName<-paste(HMManno$V1,"_",HMManno$V2,"_",HMManno$V3,"_",HMManno$V6,"_",HMManno$V7,sep="")
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