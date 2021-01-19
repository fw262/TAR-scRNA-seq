## This script finds differentially expressed uTARs and returns a list of differentially expressed uTARs with genes
# saves 2 files (fasta region and diff uTAR markers)
args=(commandArgs(TRUE))
if(length(args)==2){
  inputDiffMarkersFile<-args[1]
  bamFile<-args[2]
} else {
  stop("Please specify differential markers and bamFile.")
}

print("Input arguments are good")
print("Loading required packages (data.table, dplyr, stringr).")

#library(Seurat, lib.loc = "/programs/R-3.5.0/library")
#library(Seurat) # please make sure Seurat is installed
if (!require('data.table')){
  install.packages("data.table")
}
library(data.table)

if (!require('dplyr')){
  install.packages("dplyr")
}
library(dplyr)

if (!require('stringr')){
  install.packages("stringr")
}
library(stringr)

print("Finished loaded packages (data.table, dplyr, stringr).")

diffMarkers <- read.delim(inputDiffMarkersFile, sep="\t")
print("Finished loading differentially expressed genes and uTARs.")

#### look at features and filter for peak
# extract only differentially expressed uTARs
diffMarkers$gene_or_uTAR<-as.character(diffMarkers$gene_or_uTAR)
numdash<-str_count(diffMarkers$gene_or_uTAR,"-")
outInd<-numdash>=5 & endsWith(diffMarkers$gene_or_uTAR,"-0")
inInd<-numdash>=5 & !endsWith(diffMarkers$gene_or_uTAR,"-0")
geneInd<-!(outInd|inInd)
uTARMarkers<-diffMarkers[outInd,]
colnames(uTARMarkers)<-c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","uTAR")

returnChromString<-function(input){
  out<-unlist(strsplit(input,"-",fixed = T))
  # check if neagtive strand or positive strand
  if(out[length(out)-2]==""){
    out2<-head(out,-6)
    out3<-paste(out2,collapse = "_")
    return(out3)
  } else {
    out2<-head(out,-5)
    out3<-paste(out2,collapse = "_")
    return(out3)
  }
}
returnStartString<-function(input){
  out<-unlist(strsplit(input,"-",fixed = T))
  # check if neagtive strand or positive strand
  if(out[length(out)-2]==""){ # negative strand case
    return(out[length(out)-5])
  } else {
    return(out[length(out)-4])
  }
}
returnEndString<-function(input){
  out<-unlist(strsplit(input,"-",fixed = T))
  # check if neagtive strand or positive strand
  if(out[length(out)-2]==""){ # negative strand case
    return(out[length(out)-4])
  } else {
    return(out[length(out)-3])
  }
}
returnDirString<-function(input){
  out<-unlist(strsplit(input,"-",fixed = T))
  # check if neagtive strand or positive strand
  if(out[length(out)-2]==""){ # negative strand case
    return("-")
  } else {
    return("+")
  }
}

uTARMarkers$chr<-as.character(sapply(as.character(uTARMarkers$uTAR),FUN=returnChromString))
uTARMarkers$start<-as.character(sapply(as.character(uTARMarkers$uTAR),FUN=returnStartString))
uTARMarkers$end<-as.character(sapply(as.character(uTARMarkers$uTAR),FUN=returnEndString))
uTARMarkers$dir<-as.character(sapply(as.character(uTARMarkers$uTAR),FUN=returnDirString))
uTARMarkers$fasta<-as.character(paste0(uTARMarkers$chr,":",uTARMarkers$start,"-",uTARMarkers$end))

# extract coverage
getCoverageOnly<-function(fasta,bamFile){
  interval<-strsplit(fasta,":")[[1]][2]
  interval<-as.numeric(strsplit(interval,"-")[[1]])
  covCom<-paste0("samtools depth -r ",fasta," ",bamFile," > temp",fasta)
  test<-system(covCom)
  #check if depth file is empty
  info = file.info(paste0("temp",fasta))
  if(info$size==0){
    system(paste0("rm temp",fasta))
    return(rep(1, diff(interval)+1))
  }
  temp <- read.delim(paste0("temp",fasta), header=FALSE)
  system(paste0("rm temp",fasta))
  coverage<-list(temp[,3])
  if(length(unlist(coverage))<diff(interval)*0.9){
    return(rep(1, diff(interval)+1))
  }
  return(unlist(coverage))
}
tissueMarkersCov<-lapply(X=uTARMarkers$fasta,FUN=getCoverageOnly,bamFile=bamFile)

#smooth coverage
smoothCovFunc<-function(input,span=0.5){
  x<-1:length(input)
  y<-as.numeric(input)
  smoothOut<-loess(y~x,span=span)
  return(predict(smoothOut))
}
smoothCov<-lapply(X=tissueMarkersCov, FUN=smoothCovFunc,span=0.25)

# get half width maximum of smoothed peak
findFHWM<-function(input){
  x<-1:length(input)
  y<-input
  xmax <- x[y==max(y)]
  
  if(max(y)> 100){
    x1<-rev(x[x < xmax])[which.min(rev(round(abs(y[x < xmax]-max(y)/2),digits=-1)))]
  } else {
    x1<-rev(x[x < xmax])[which.min(rev(round(abs(y[x < xmax]-max(y)/2))))]
  }
  
  x2 <- x[x > xmax][which.min((round(abs(y[x > xmax]-max(y)/2))))]
  
  #x1 <- rev(x[x < xmax])[(which(diff(sign(diff(abs(rev(y[x < xmax]-max(y)/2))*-1)))==-2)+1)[1]]
  #x2 <- x[x > xmax][(which(diff(sign(diff(abs((y[x > xmax]-max(y)/2))*-1)))==-2)+1)[1]]
  
  if(length(x1)==0){
    return(c(0,x2))
  } else if (length(x2)==0){
    return(c(x1,length(input)))
  }else{
    return(c(x1,x2))
  }
}
indicesToBlast<-lapply(X=smoothCov,FUN=findFHWM)

# add index shift to create fasta for peaks
shiftDF<-matrix(unlist(indicesToBlast), nrow=length(indicesToBlast), byrow=T)
tooShortInd<-abs(shiftDF[,2]-shiftDF[,1])<50
negativeInd<-(shiftDF[,2]-shiftDF[,1])<0
indTofix<-tooShortInd|negativeInd
shiftDF[indTofix,1]<-0
shiftDF[indTofix,2]<-(as.numeric(uTARMarkers$end[indTofix]) - as.numeric(uTARMarkers$start[indTofix]))

uTARMarkers<-cbind(uTARMarkers,data.frame(shiftDF))
uTARMarkers$fastaPeak<-paste0(uTARMarkers$chr,":",as.numeric(uTARMarkers$start)+uTARMarkers$X1,"-",as.numeric(uTARMarkers$start)+uTARMarkers$X2)
#indTooShort<-abs(uTARMarkers$X2 - uTARMarkers$X1)<50
#if peaks too short, just keep full uTAR
#uTARMarkers$fastaPeak[indTooShort]<-paste0(uTARMarkers$chr[indTooShort],":",as.numeric(uTARMarkers$start[indTooShort]),"-",as.numeric(uTARMarkers$end[indTooShort]))

if(dirname(inputDiffMarkersFile)!="."){
  sampleName<-tail(unlist(strsplit(dirname(inputDiffMarkersFile),"/")),n=1)
  output<-paste0(dirname(inputDiffMarkersFile),"/",sampleName,"_seqsToBlast.txt")
} else {
  output<-paste0("seqsToBlast.txt")
}

write.table(uTARMarkers$fastaPeak,file=output,quote = F, row.names = F,col.names = F) # save fasta format

if(dirname(inputDiffMarkersFile)!="."){
  sampleName<-tail(unlist(strsplit(dirname(inputDiffMarkersFile),"/")),n=1)
  output<-paste0(dirname(inputDiffMarkersFile),"/",sampleName,"_diffuTARMarkers.txt")
} else {
  output<-paste0("diffuTARMarkers.txt")
}
write.table(uTARMarkers,file=output,quote = F, row.names = F,col.names = T,sep="\t")
