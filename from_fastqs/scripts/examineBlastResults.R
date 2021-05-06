## This script finds differentially expressed uTARs and returns a list of differentially expressed uTARs with genes
args=(commandArgs(TRUE))
if(length(args)==2){
  inputDiffMarkersFile<-args[1]
  blastResultFile<-args[2]
} else {
  stop("Please specify differential markers and blast results.")
}

print("Input arguments are good")
print("Loading required packages (data.table, dplyr, stringr).")

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

diffMarkers <- read.delim(inputDiffMarkersFile,sep="\t")
blastResult <- read.delim(blastResultFile, header=FALSE)
print("Finished loading differentially expressed uTARs and fasta results.")

addInBlastResult<-function(blastResult,gonadUTarBed){
  getNFromList <- function(lst, n){
    sapply(lst, `[`, n)
  }
  # extract gene name in brackets
  temp<-getNFromList(strsplit(as.character(blastResult$V3),"(",fixed=TRUE),2)
  geneID<-getNFromList(strsplit(as.character(temp),")",fixed=TRUE),1)
  blastResult$geneID<-geneID
  
  # populate diff uTARs with blastResult
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  gonadUTarBed$blast<-NA
  gonadUTarBed$blastShort<-NA
  
  for (i in 1:nrow(gonadUTarBed)){
    blastTmp<-blastResult[blastResult$V1==gonadUTarBed$fastaPeak[i],]
    blastTmp<-blastTmp[order(blastTmp$V13,decreasing = T),]
    blastShorts<-as.character(blastTmp$geneID)
    blastShorts<-blastShorts[!is.na(blastShorts)]
    gonadUTarBed$blast[i]<-as.character(blastTmp$V3)[1]
    gonadUTarBed$blastShort[i]<-blastShorts[1]
  }
  
  gonadUTarBed$blastShort[is.na(gonadUTarBed$blastShort)]<-"undetermined"
  return(gonadUTarBed)
}
blastResult$V1<-as.character(blastResult$V1)
blastResult$V2<-as.character(blastResult$V2)
blastResult$V3<-as.character(blastResult$V3)
diffMarkers$uTAR<-as.character(diffMarkers$uTAR)
diffMarkers$dir<-as.character(diffMarkers$dir)
diffMarkers$fasta<-as.character(diffMarkers$fasta)
diffMarkers$fastaPeak<-as.character(diffMarkers$fastaPeak)

diffMarkers<-addInBlastResult(blastResult,diffMarkers)
diffMarkers$blast[is.na(diffMarkers$blast)]<-"undetermined"
print("Finished labeling based on best BLAST results.")

if(dirname(inputDiffMarkersFile)!="."){
  sampleName<-tail(unlist(strsplit(dirname(inputDiffMarkersFile),"/")),n=1)
  output<-paste0(dirname(inputDiffMarkersFile),"/",sampleName,"_diffuTARMarkersLabeled.txt")
} else {
  output<-paste0("diffuTARMarkersLabeled.txt")
}

write.table(diffMarkers,file=output,quote = F, row.names = F,col.names = T,sep="\t")
