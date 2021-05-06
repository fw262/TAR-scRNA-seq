# convert refFlat to gtf format
args=(commandArgs(TRUE))
inputRefFlat<-args[1]

TAR.refFlat <- read.delim(inputRefFlat, header=FALSE)

TAR.gtf<-data.frame(matrix(ncol = ncol(TAR.refFlat),nrow = nrow(TAR.refFlat)))
TAR.gtf$X1<-TAR.refFlat$V3 #seqname
TAR.gtf$X2<-"TAR_HMM" #source
TAR.gtf$X3<-"exon"
TAR.gtf$X4<-TAR.refFlat$V5
TAR.gtf$X5<-TAR.refFlat$V6
TAR.gtf$X6<-"."
TAR.gtf$X7<-TAR.refFlat$V4
TAR.gtf$X8<-"0"
TAR.gtf$X9<-paste0('gene_id "',TAR.refFlat$V1,'"; transcript_id "',TAR.refFlat$V1,'";')
TAR.gtf<-TAR.gtf[,1:9]

write.table(TAR.gtf,file=paste0(inputRefFlat,".gtf"),row.names=F,col.names = F,quote=F,sep="\t")
