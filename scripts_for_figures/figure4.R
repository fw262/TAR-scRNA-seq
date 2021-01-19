# Libraries
library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(pals);library(VennDiagram);library(viridis)
library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
#source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/functionsToAnalyse.R")
source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")

# spaceranger data folders
day14<-"/workdir/fw262/chickenSpatial/day14Spatial/outs"
day10<-"/workdir/fw262/chickenSpatial/day10Spatial/outs"
day7<-"/workdir/fw262/chickenSpatial/day7Spatial/outs"
day4<-"/workdir/fw262/chickenSpatial/day4Spatial/outs"

# load up gene markers
load(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure3/markersData_forPaper.RData")
View(diffMarkers$d4)
numDashes<-5
numdash<-str_count(diffMarkers$d4$gene,"-")
outInd<-numdash>=numDashes & endsWith(diffMarkers$d4$gene,"-0")
inInd<-numdash>=numDashes & endsWith(diffMarkers$d4$gene,"-1")
geneInd<-!(outInd|inInd)
geneOnly<-diffMarkers$d4$gene[geneInd]
inGene<-diffMarkers$d4$gene[inInd]
outGene<-diffMarkers$d4$gene[outInd]
d4Gene<-diffMarkers$d4[diffMarkers$d4$gene %in% geneOnly,]
top5D4Gene<-d4Gene %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)


# define mito genes
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
createSpatialSeu<-function(input,mito.genes=mito_genes,numDashes=5){
  output<-Load10X_Spatial(input)
  
  # calculate percent HMM out reads
  numdash<-str_count(rownames(output),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(output),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(output),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(output)[geneInd]
  inGene<-rownames(output)[inInd]
  outGene<-rownames(output)[outInd]
  
  out<-output[c(inGene,outGene),]
  output[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  output <- PercentageFeatureSet(object = output, features = mito.genes, col.name = "percent.mito")
  #VlnPlot(output,features=c("nCount_Spatial","percent.mito"),pt.size=0.1)
  output<-NormalizeData(output) %>% FindVariableFeatures() %>% ScaleData()
  return(output)
}

runDimReducFunction<-function(sample_combined, numDashes=5, res_param=1){
  # calculate percent HMM out reads
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=res_param)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=res_param)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  
  return(sample_combined)
}
checkIfExistGene2<-function(input,gene_ref){
  chrom<-input[[1]]
  startPos<-as.numeric(input[[2]])
  endPos<-as.numeric(input[[3]])
  direction<-input[[4]]
  gene_ref<-gene_ref # bed of entire genes
  #print(paste(chrom,startPos,endPos))
  
  # 1: either partially in or completely inside gene
  chrMatch_gene<-(chrom==gene_ref$chr)
  dirMatch_gene<-(direction==gene_ref$direction)
  qinrstartMatch<-(startPos>=gene_ref$start)&(startPos<=gene_ref$end)
  qinrendMatch<-(endPos<=gene_ref$end)&(endPos>=gene_ref$start)
  qinroutGeneAll_gene<-(chrMatch_gene*qinrstartMatch*qinrendMatch*dirMatch_gene) # query completely inside ref 
  rinqstartMatch<-(startPos<(gene_ref$start))
  rinqendMatch<-(endPos>(gene_ref$end))
  rinqoutGeneAll<-(chrMatch_gene*rinqstartMatch*rinqendMatch*dirMatch_gene) # ref completely inside query
  # also check partial overlap from beginning
  partialStart<-(chrMatch_gene*rinqstartMatch*qinrendMatch*dirMatch_gene)
  partialEnd<-(chrMatch_gene*qinrstartMatch*rinqendMatch*dirMatch_gene)
  if(sum(qinroutGeneAll_gene)|sum(rinqoutGeneAll)|sum(partialStart)|sum(partialEnd)){
    # g's are for gene names
    g1<-as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
    g2<-as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
    g3<-as.character(gene_ref$gene[as.logical(partialStart)])
    g4<-as.character(gene_ref$gene[as.logical(partialEnd)])
    gAll<-unique(c(g1,g2,g3,g4))
    # d's are for directions
    d1<-as.character(gene_ref$direction[as.logical(qinroutGeneAll_gene)])
    d2<-as.character(gene_ref$direction[as.logical(rinqoutGeneAll)])
    d3<-as.character(gene_ref$direction[as.logical(partialStart)])
    d4<-as.character(gene_ref$direction[as.logical(partialEnd)])
    dAll<-unique(c(d1,d2,d3,d4))
    return(paste(c(gAll,dAll,"1"),collapse="_"))
  }else{
    #return(paste(c(direction,"0"),collapse="_"))
    return(0)
  }
}

d4SpSeu<-createSpatialSeu(day4)
d7SpSeu<-createSpatialSeu(day7)
d10SpSeu<-createSpatialSeu(day10)
d14SpSeu<-createSpatialSeu(day14)

VlnPlot(d14SpSeu,features=c("nCount_Spatial","percent.mito","percent.outHMM"),pt.size=0)
d4Out<-SpatialFeaturePlot(d4SpSeu,features=c("percent.outHMM"))
d7Out<-SpatialFeaturePlot(d7SpSeu,features=c("percent.outHMM"))
d10Out<-SpatialFeaturePlot(d10SpSeu,features=c("percent.outHMM"))
d14Out<-SpatialFeaturePlot(d14SpSeu,features=c("percent.outHMM"))
plot_grid(d4Out,d7Out,d10Out,d14Out,nrow=2)
SpatialFeaturePlot(d14SpSeu,features=c("RUNX1T1"))
d4SpSeu<-runDimReducFunction(d4SpSeu)
d7SpSeu<-runDimReducFunction(d7SpSeu)
d10SpSeu<-runDimReducFunction(d10SpSeu)
d14SpSeu<-runDimReducFunction(d14SpSeu)

# merge datasets????
mergedSeu<-merge(d4SpSeu,y=c(d7SpSeu,d10SpSeu,d14SpSeu),add.cell.ids=c("D4","D7","D10","D14"))
# calculate percent HMM out reads
numDashes<-5
numdash<-str_count(rownames(mergedSeu),"-")
outInd<-numdash>=numDashes & endsWith(rownames(mergedSeu),"-0")
inInd<-numdash>=numDashes & endsWith(rownames(mergedSeu),"-1")
geneInd<-!(outInd|inInd)
geneOnly<-rownames(mergedSeu)[geneInd]
inGene<-rownames(mergedSeu)[inInd]
outGene<-rownames(mergedSeu)[outInd]

# list of differentially expressed uTARs from NON-spatial
nonSpUTar<-c("1-65683949-65711899-+-5762-0", #"SOX5"
"2-124878999-124895699---3773-0", #"RUNX1T1"
"21-1246149-1247899-+-169459-0", #"LOC419389"
"1-109375749-109386749-+-333313-0", #"SH3BGR"
"1-26785099-26808299-+-81344-0", #"LOC107052650"
"4-30658399-30689599---24422-0", #"ANAPC10"
"5-37573549-37577599---6394-0", #"CLEC14A"
"3-47571549-47666099-+-28210-0", #"SASH1"
"4-57304549-57388099---9473-0", #"LOC106038660"
"6-4454149-4522999---5590-0", #"NRG3"
"4-5778899-5816449---29933-0", #"DIAPH2"
"5-13398249-13417749-+-25439-0", #"CDKN1C"
"14-14579349-14597099-+-13435-0", #"LOC104913510"
"7-20808699-20938349-+-10524-0", #"KCNH7"
"13-9966349-10024349-+-3030-0", #"RPL26L1"
"1-125853499-125915949-+-4010-0", #"MID1"
"3-1536699-1613799---4759-0", #"LOC110120668"
"32-583299-584549---137241-0")#"LOC107050595"

uTARClust<-c(0,
             0,
             0,
             1,
             4,
             2,
             3,
             3,
             4,
             4,
             4,
             4,
             5,
             6,
             6,
             6,
             6,
             6)

nonSpUTarName<-c("SOX5_uTAR",
             "RUNX1T1_uTAR",
             "LOC419389_uTAR",
             "SH3BGR_uTAR",
             "LOC107052650_uTAR",
             "ANAPC10_uTAR",
             "CLEC14A_uTAR",
             "SASH1_uTAR",
             "LOC106038660_uTAR",
             "NRG3_uTAR",
             "DIAPH2_uTAR",
             "CDKN1C_uTAR",
             "LOC104913510_uTAR",
             "KCNH7_uTAR",
             "RPL26L1_uTAR",
             "MID1_uTAR",
             "LOC110120668_uTAR",
             "LOC107050595_uTAR")

uTARGeneList<-c("SOX5",
                "RUNX1T1",
                "LOC419389",
                "SH3BGR",
                "LOC107052650",
                "ANAPC10",
                "CLEC14A",
                "SASH1",
                "LOC106038660",
                "NRG3",
                "DIAPH2",
                "CDKN1C",
                "LOC104913510",
                "KCNH7",
                "RPL26L1",
                "MID1",
                "LOC110120668",
                "LOC107050595")
# non spatial uTAR df
nonSpatialDf<-data.frame(cbind(as.character(nonSpUTar),as.character(uTARGeneList)))
nonSpatialDf$chr<-getNFromList(strsplit(nonSpUTar,"-"),1)
nonSpatialDf$start<-as.numeric(getNFromList(strsplit(nonSpUTar,"-"),2))
nonSpatialDf$end<-as.numeric(getNFromList(strsplit(nonSpUTar,"-"),3))
nonSpatialDf$direction<-getNFromList(strsplit(nonSpUTar,"-"),4)
nonSpatialDf$direction[nonSpatialDf$direction!="+"]<-"-"
colnames(nonSpatialDf)[c(1,2)]<-c("uTar","gene")
nonSpatialDf$uTar<-as.character(nonSpatialDf$uTar)
nonSpatialDf$gene<-as.character(nonSpatialDf$gene)
str(nonSpatialDf)

# spatial uTAR df
spatialUTar<-data.frame(outGene)
spatialUTar$outGene<-as.character(spatialUTar$outGene)
spatialUTar$chr<-getNFromList(strsplit(as.character(spatialUTar$outGene),"-"),1)
spatialUTar$start<-as.numeric(getNFromList(strsplit(as.character(spatialUTar$outGene),"-"),2))
spatialUTar$end<-as.numeric(getNFromList(strsplit(as.character(spatialUTar$outGene),"-"),3))
spatialUTar$direction<-getNFromList(strsplit(as.character(spatialUTar$outGene),"-"),4)
spatialUTar$direction[spatialUTar$direction!="+"]<-"-"
HMManno_bare<-spatialUTar[,c("chr","start","end","direction")]
HMManno_bare$inGene<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene2,gene_ref=nonSpatialDf)
HMManno_bare$inGene_noDir<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene_noDir,gene_ref=nonSpatialDf)
HMManno_bare$uTAR<-spatialUTar$outGene

# filter for uTARs that overlap with gene annotations with direction
spaUTarDir<-HMManno_bare[HMManno_bare$inGene!=0,]

# myocyte expression and uTAR SH3BGR, genes TNNT2 TNNI1 TNNC1
SpatialFeaturePlot(d4SpSeu,features=c("TNNC1","TNNT2","TNNI1","1-109378099-109383399-+-156456-0"),pt.size.factor = 2,ncol = 4)
saveSpatialFunc(d4SpSeu,"TNNT2",maxScale = 4,dayText = "d4")
saveSpatialFunc(d4SpSeu,"1-109378099-109383399-+-156456-0",maxScale = 4,dayText = "d4")
SpatialFeaturePlot(d7SpSeu,features=c("TNNC1","TNNT2","TNNI1","1-109378099-109383399-+-156456-0"),pt.size.factor = 2,ncol = 4)
saveSpatialFunc(d7SpSeu,"TNNT2",maxScale = 5,dayText = "d7")
saveSpatialFunc(d7SpSeu,"1-109378099-109383399-+-156456-0",maxScale = 5,dayText = "d7")
SpatialFeaturePlot(d10SpSeu,features=c("TNNC1","TNNT2","TNNI1","1-109378099-109383399-+-156456-0"),pt.size.factor = 2,ncol = 4)
saveSpatialFunc(d10SpSeu,"TNNT2",maxScale = 4.5,dayText = "d10",ptsizefactor = 1.5)
saveSpatialFunc(d10SpSeu,"1-109378099-109383399-+-156456-0",maxScale = 4.5,dayText = "d10",ptsizefactor = 1.5)
SpatialFeaturePlot(d14SpSeu,features=c("TNNC1","TNNT2","TNNI1","1-109378099-109383399-+-156456-0"),pt.size.factor = 2,ncol = 4)
saveSpatialFunc(d14SpSeu,"TNNT2",maxScale = 4,dayText = "d14")
saveSpatialFunc(d14SpSeu,"1-109378099-109383399-+-156456-0",maxScale = 4,dayText = "d14")

## immature fibroblasts NRG3
#SpatialFeaturePlot(d4SpSeu,features=c("DCN","LPL","6-4506999-4507449---44-0"),pt.size.factor = 2)
#SpatialFeaturePlot(d7SpSeu,features=c("DCN","LPL","6-4506999-4507449---44-0"),pt.size.factor = 2)
#SpatialFeaturePlot(d10SpSeu,features=c("DCN","LPL","6-4506999-4507449---44-0"),pt.size.factor = 2)
#SpatialFeaturePlot(d14SpSeu,features=c("DCN","LPL","6-4506999-4507449---44-0"),pt.size.factor = 2)

# epicardial progenitor, uTAR RUNX1T1, gene CHODL VCAN CNMD COL1A1 HAPLN1
SpatialFeaturePlot(d4SpSeu,features=c("CHODL","COL1A1","CNMD","2-124889449-124890099---1997-0"),pt.size.factor = 2,ncol = 4)
saveSpatialFunc(d4SpSeu,"COL1A1",maxScale = 3.5,dayText = "d4")
saveSpatialFunc(d4SpSeu,"RUNX1T1",maxScale = 3.5,dayText = "d4")
saveSpatialFunc(d4SpSeu,"2-124889449-124890099---1997-0",maxScale = 3.5,dayText = "d4")
SpatialFeaturePlot(d7SpSeu,features=c("CHODL","COL1A1","CNMD","2-124889449-124890099---1997-0"),pt.size.factor = 2,ncol = 4)
saveSpatialFunc(d7SpSeu,"COL1A1",maxScale = 4,dayText = "d7")
saveSpatialFunc(d7SpSeu,"RUNX1T1",maxScale = 4,dayText = "d7")
saveSpatialFunc(d7SpSeu,"2-124889449-124890099---1997-0",maxScale = 4,dayText = "d7")
SpatialFeaturePlot(d10SpSeu,features=c("CHODL","COL1A1","CNMD","2-124889449-124890099---1997-0"),pt.size.factor = 2,ncol = 4)
saveSpatialFunc(d10SpSeu,"COL1A1",maxScale = 4.5,dayText = "d10",ptsizefactor = 1.5)
saveSpatialFunc(d10SpSeu,"RUNX1T1",maxScale = 4.5,dayText = "d10",ptsizefactor = 1.5)
saveSpatialFunc(d10SpSeu,"2-124889449-124890099---1997-0",maxScale = 4.5,dayText = "d10",ptsizefactor = 1.5)
SpatialFeaturePlot(d14SpSeu,features=c("CHODL","COL1A1","CNMD","2-124889449-124890099---1997-0"),pt.size.factor = 2,ncol = 4)
saveSpatialFunc(d14SpSeu,"COL1A1",maxScale = 5,dayText = "d14")
saveSpatialFunc(d14SpSeu,"RUNX1T1",maxScale = 5,dayText = "d14")
saveSpatialFunc(d14SpSeu,"2-124889449-124890099---1997-0",maxScale = 5,dayText = "d14")

# endothelial CLEC14A uTAR, uTAR 5-37574549-37576199---1225-0, gene POSTN CDH5 RARRES1 FHL1 SHE PECAM1
#SpatialFeaturePlot(d4SpSeu,features=c("PECAM1","RARRES1","5-37574549-37576199---1225-0"),pt.size.factor = 2)

# other, gene XIRP1 MYH7B POPDC2 TRIM55, uTAR KCNH7 7-20831849-20846849-+-8782-0
SpatialFeaturePlot(d4SpSeu,features=c("TRIM55","7-20831849-20846849-+-8782-0"),pt.size.factor = 2)
saveSpatialFunc(d4SpSeu,"TRIM55",maxScale = 3.75,dayText = "d4")
saveSpatialFunc(d4SpSeu,"7-20831849-20846849-+-8782-0",maxScale = 3.75,dayText = "d4")
SpatialFeaturePlot(d7SpSeu,features=c("TRIM55","7-20831849-20846849-+-8782-0"),pt.size.factor = 2)
saveSpatialFunc(d7SpSeu,"TRIM55",maxScale = 3.75,dayText = "d7")
saveSpatialFunc(d7SpSeu,"7-20831849-20846849-+-8782-0",maxScale = 3.75,dayText = "d7")
SpatialFeaturePlot(d10SpSeu,features=c("TRIM55","7-20831849-20846849-+-8782-0"),pt.size.factor = 2)
saveSpatialFunc(d10SpSeu,"TRIM55",maxScale = 2.5,dayText = "d10",ptsizefactor = 1.5)
saveSpatialFunc(d10SpSeu,"7-20831849-20846849-+-8782-0",maxScale = 2.5,dayText = "d10",ptsizefactor = 1.5)
SpatialFeaturePlot(d14SpSeu,features=c("TRIM55","7-20831849-20846849-+-8782-0"),pt.size.factor = 2)
saveSpatialFunc(d14SpSeu,"TRIM55",maxScale = 2.5,dayText = "d14")
saveSpatialFunc(d14SpSeu,"7-20831849-20846849-+-8782-0",maxScale = 2.5,dayText = "d14")

# HBA1,HBZ RBCs, uTAR 4-30676149-30679099---776-0 ANAPC10
SpatialFeaturePlot(d4SpSeu,features=c("HBZ","4-30676149-30679099---776-0"),pt.size.factor = 2)
saveSpatialFunc(d4SpSeu,"TRIM55",maxScale = 3.75,dayText = "d4")
saveSpatialFunc(d4SpSeu,"7-20831849-20846849-+-8782-0",maxScale = 3.75,dayText = "d4")
SpatialFeaturePlot(d7SpSeu,features=c("HBZ","4-30676149-30679099---776-0"),pt.size.factor = 2)
saveSpatialFunc(d7SpSeu,"TRIM55",maxScale = 3.75,dayText = "d7")
saveSpatialFunc(d7SpSeu,"7-20831849-20846849-+-8782-0",maxScale = 3.75,dayText = "d7")
SpatialFeaturePlot(d10SpSeu,features=c("HBZ","4-30676149-30679099---776-0"),pt.size.factor = 2)
saveSpatialFunc(d10SpSeu,"TRIM55",maxScale = 2.5,dayText = "d10",ptsizefactor = 1.5)
saveSpatialFunc(d10SpSeu,"7-20831849-20846849-+-8782-0",maxScale = 2.5,dayText = "d10",ptsizefactor = 1.5)
SpatialFeaturePlot(d14SpSeu,features=c("HBA1","4-30676149-30679099---776-0"),pt.size.factor = 2)
saveSpatialFunc(d14SpSeu,"TRIM55",maxScale = 2.5,dayText = "d14")
saveSpatialFunc(d14SpSeu,"7-20831849-20846849-+-8782-0",maxScale = 2.5,dayText = "d14")


# extract normalized data to perform pearson correlation test
d4SpSeuMat<-GetAssayData(object = d4SpSeu)
d7SpSeuMat<-GetAssayData(object = d7SpSeu)
d10SpSeuMat<-GetAssayData(object = d10SpSeu)
d14SpSeuMat<-GetAssayData(object = d14SpSeu)

cor.test(d4SpSeuMat["COL1A1",],d4SpSeuMat["1-109378099-109383399-+-156456-0",],method=c("pearson"))

# correlation tests
genesToTest<-unique(top5D4Gene$gene) # extract genes that are differentially expressed
genesToTest<-c(genesToTest,"TNNT2","COL1A1")
library("corrplot")
d4DfToTest<-d4SpSeuMat[c(genesToTest),]
geneCon<-data.frame(t(data.frame(d4DfToTest)))
colnames(geneCon)<-rownames(d4DfToTest)

#####################################################################
# add up expression in uTARs that correspond to same gene
for (i in unique(spaUTarDir$inGene)){
  uTARs<-spaUTarDir$uTAR[spaUTarDir$inGene==i]
  if(ncol(data.frame(d4SpSeuMat[uTARs,]))>1){
    temp<-data.frame(t(data.frame(d4SpSeuMat[uTARs,])))
    geneCon[,i]<-rowSums(temp)
  } else {
    temp<-data.frame(data.frame(d4SpSeuMat[uTARs,]))
    geneCon[,i]<-temp
  }
}

# pick spatial uTAR with the highest expression
for (i in unique(spaUTarDir$inGene)){
  uTARs<-spaUTarDir$uTAR[spaUTarDir$inGene==i]
  if(ncol(data.frame(d4SpSeuMat[uTARs,]))>1){
    l<-str_split(uTARs,"-") # get second last element in list to find expression
    out<-as.numeric(unlist(lapply(l, function(x) x[length(x) -1])))
    maxuTAR<-uTARs[which(out==max(out))]
    temp<-data.frame(data.frame(d4SpSeuMat[maxuTAR,]))
    geneCon[,i]<-temp
  } else {
    temp<-data.frame(data.frame(d4SpSeuMat[uTARs,]))
    geneCon[,i]<-temp
  }
}
#####################################################################

colnames(geneCon)<-str_replace_all(colnames(geneCon),pattern = "_\\+_1",replacement = "_uTAR")
colnames(geneCon)<-str_replace_all(colnames(geneCon),pattern = "_-_1",replacement = "_uTAR")
#keep only variables where mean is greater than 0
geneCon<-geneCon[,colMeans(geneCon)>0]
geneCon<-geneCon[,!startsWith(colnames(geneCon),"ENSGAL")]
# remove certain genes before comparison
geneCon<-geneCon[,!startsWith(colnames(geneCon),"LOC")]

M<-cor(geneCon,method = "pearson")
pvals<-cor.mtest(geneCon)$p

# corrplot version
head(round(M,2))
col<- colorRampPalette(viridis(n=3))(200)
corrplot(M,method="color",order = "hclust",type = "lower",p.mat = pvals, sig.level = 0.05, insig = "blank")

library(ggdendro)
# ggplot version
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  plot(hc)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
M <- reorder_cormat(M)
# Melt the correlation matrix
melted_M <- melt(M, na.rm = TRUE)

############### save dendrogram
dd <- as.dist((1-M)/2)
hc <- hclust(dd)
outFile<-paste0("/workdir/fw262/chickenSpatial/figures/dendrogram_noSum.pdf")
pdf(file=outFile, width=4.2, height=2.5, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F )
ggdendrogram(hc, rotate = FALSE, size = 2)
dev.off()
####################################


# define color
colorsForClust<-alphabet(n=length(unique(d4Gene$cluster)))
d4Gene$colors<-colorsForClust[as.numeric(levels(d4Gene$cluster))[d4Gene$cluster]+1]
colorDict<-d4Gene[rownames(d4DfToTest),"colors"]
names(colorDict)<-rownames(d4DfToTest)

temp<-unique(spaUTarDir$inGene)
temp<-str_replace_all(temp,pattern = "_\\+_1",replacement = "_uTAR")
temp<-str_replace_all(temp,pattern = "_-_1",replacement = "_uTAR")
colorDictuTAR<-colorsForClust[uTARClust+1]
names(colorDictuTAR)<-nonSpUTarName
colorDictAll<-c(colorDict,colorDictuTAR)

temp<-melted_M$Var1[melted_M$Var2==unique(melted_M$Var2)[1]]
colToPlot<-colorDictAll[as.character(temp)]
# Create a ggheatmap
#melted_M$value[melted_M$value>0.5]<-0.5
ggheatmap <- ggplot(melted_M, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  coord_fixed()+scale_fill_viridis()+#limits=c(min(melted_M$value),0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = colToPlot))+
  theme(axis.text.y = element_text(colour = colToPlot))


