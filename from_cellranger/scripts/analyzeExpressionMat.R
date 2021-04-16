## This script finds differentially expressed uTARs and returns a list of differentially expressed uTARs with genes
args=(commandArgs(TRUE))
if(length(args)==0){
  stop("Please specify gene and TAR expression matrices.")
} else if (length(args)==2){
  geneFile<-args[1] # first argument is gene File
  TARFile<-args[2] # second argument is TAR File
  nDiffFeat<-5
  cells<-2
  feats<-100
  nFeatMin<-200
  nFeatMax<-2500
  pcaDims<-10
  res_param<-0.2
} else if (length(args)==9){
  geneFile<-args[1] # first argument is gene File
  TARFile<-args[2] # second argument is TAR File
  nDiffFeat<-args[3]
  cells<-args[4]
  feats<-args[5]
  nFeatMin<-args[6]
  nFeatMax<-args[7]
  pcaDims<-args[8]
  res_param<-args[9]
} else {
  stop("Please specify all optional parameters.")
}

print("Input arguments are good")
print("Loading required packages (Seurat, data.table, dplyr, stringr).")

#library(Seurat, lib.loc = "/programs/R-3.5.0/library")
if (!require('Seurat')){
  install.packages('Seurat')
}
library(Seurat) # please make sure Seurat is installed

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

print("Finished loaded packages (Seurat, data.table, dplyr, stringr).")

createSeuratObj<-function(geneMat,hmmMat,nDiffFeat=5,cells=2,feats=100,nFeatMin=200,nFeatMax=2500,pcaDims=10,res_param=0.2){
  print("Loading in expression matrices.")
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  print("Creating gene expression matrix.")
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = cells, min.features = feats)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > nFeatMin & nFeature_RNA < nFeatMax)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM[,colnames(d4Gene)])
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  print("Finished creating Seurat object.")
  
  # identify aTARs and uTARs
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & !endsWith(rownames(sample_combined),"-0")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  sample_combined <- NormalizeData(sample_combined)
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  print("Finished normalizing and scaling data.")
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  print("Finished running PCA on uTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:pcaDims)
  sample_combined <- FindClusters(object=sample_combined,resolution=res_param)
  print("Finished running PCA on genes and cluster based on genes.")
  
  allMarkers <- FindAllMarkers(sample_combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
  
  numdash<-str_count(allMarkers$gene,"-")
  outInd<-numdash>=5 & endsWith(allMarkers$gene,"-0")
  inInd<-numdash>=5 & !endsWith(allMarkers$gene,"-0")
  geneInd<-!(outInd|inInd)
  
  markersUTAR<-allMarkers[outInd,]
  markersUTAR<-markersUTAR[order(markersUTAR$avg_logFC,decreasing = T),]
  markersUTAR <- markersUTAR %>% group_by(cluster) %>% top_n(n = nDiffFeat, wt = avg_logFC)
  
  markersGenes<-allMarkers[geneInd,]
  markersGenes<-markersGenes[order(markersGenes$avg_logFC,decreasing = T),]
  markersGenes <- markersGenes %>% group_by(cluster) %>% top_n(n = nDiffFeat, wt = avg_logFC)
  
  print("Finished finding differentially expressed genes and uTAR markers.")
    
  return(rbind(markersGenes,markersUTAR))
}

diffMarkers<-createSeuratObj(geneMat = geneFile,
                             hmmMat = TARFile,
                             nDiffFeat = nDiffFeat,
                             cells = cells,
                             feats = feats,
                             nFeatMin = nFeatMin,
                             nFeatMax = nFeatMax,
                             pcaDims = pcaDims,
                             res_param = res_param)

if(dirname(geneFile)!="."){
  sampleName<-tail(unlist(strsplit(dirname(geneFile),"/")),n=1)
  output<-paste0(dirname(geneFile),"/",sampleName,"_diffMarkers.txt")
} else {
  output<-paste0("diffMarkers.txt")
}

print(paste0("Output file stored in ",output,"."))
colnames(diffMarkers)<-c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene_or_uTAR")
write.table(diffMarkers,file=output,quote = F, row.names = F,col.names = T,sep="\t")
