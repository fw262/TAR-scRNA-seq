# Libraries
library(ggplot2); library(data.table); library(caTools);library(scales); library(dplyr)
library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
library(Seurat, lib.loc = "/programs/R-3.5.0/library")

load(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/allDataDiffParam.RData")   

addInBlastResult<-function(blastResult,gonadUTarBed){
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
      blastTmp<-blastResult[blastResult$V1==gonadUTarBed$fasta[i],]
      blastTmp<-blastTmp[order(blastTmp$V13,decreasing = T),]
      blastShorts<-as.character(blastTmp$geneID)
      blastShorts<-blastShorts[!is.na(blastShorts)]
      gonadUTarBed$blast[i]<-as.character(blastTmp$V3)[1]
      gonadUTarBed$blastShort[i]<-blastShorts[1]
   }
   
   gonadUTarBed$blastShort[is.na(gonadUTarBed$blastShort)]<-"undetermined"
   return(gonadUTarBed)
}

###########################################
#############################################
# look at smartseq data with TAR Antoine testes samples
###########################################
#############################################
# load up TAR Mat
TARMatTestes <- read.delim("/workdir/fw262/uTARAnalysisTool/TAR-scRNA-seq/smartseq2/mouseLemur/testesSS2/testesSS2TARs.mat.gz")
rownames(TARMatTestes)<-TARMatTestes[,1]
TARMat<-TARMatTestes[,-1]
colnames(TARMatTestes)<-paste0(getNFromList(strsplit(colnames(TARMatTestes),"[.]"),2),"-smartseq2")
TARMatTestes[is.na(TARMatTestes)] <- 0
TARMatTestes<-TARMatTestes[,order(colnames(TARMatTestes))]
TARMatTestes<-TARMatTestes[,!startsWith(colnames(TARMatTestes),"NA-smartseq2")]

# load up gene Mat 
geneMatTestes <- read.delim("/workdir/fw262/uTARAnalysisTool/TAR-scRNA-seq/smartseq2/mouseLemur/testesSS2/testesSS2genes.mat.gz")
rownames(geneMatTestes)<-geneMatTestes[,1]
geneMatTestes<-geneMatTestes[,-1]
colnames(geneMatTestes)<-paste0(getNFromList(strsplit(colnames(geneMatTestes),"[.]"),2),"-smartseq2")
geneMatTestes[is.na(geneMatTestes)] <- 0
geneMatTestes<-geneMatTestes[,order(colnames(geneMatTestes))]
geneMatTestes<-geneMatTestes[,!startsWith(colnames(geneMatTestes),"NA-smartseq2")]

sample_combined_mat<-rbind(geneMatTestes,TARMatTestes)
sample_combined<-CreateSeuratObject(counts = sample_combined_mat, min.cells=1)

numdash<-str_count(rownames(sample_combined),"-")
outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
inInd<-numdash>=5 & !endsWith(rownames(sample_combined),"-0")
geneInd<-!(outInd|inInd)
geneOnly<-rownames(sample_combined)[geneInd]
inGene<-rownames(sample_combined)[inInd]
outGene<-rownames(sample_combined)[outInd]
outputAll<-analysisOnFilSeurat(sample_combined)

# load up metadata of antoine spleen
metadataAll <- read.csv("/workdir/fw262/mouseLemur/LCA_SS2_10x_all_metadata_update.csv")
metadataAll$X<-paste(metadataAll$X,metadataAll$method,sep="")
metadata <- metadataAll[metadataAll$X %in% colnames(outputAll),]
metadataSS2<-metadataAll[metadataAll$method=="smartseq2",]

# add in cell ontology class
cellTypeLabel<-as.character(metadata$cell_ontology_class)
cellTypeLabel[is.na(cellTypeLabel)]<-"na"
names(cellTypeLabel)<-metadata$X

outputAll<-AddMetaData(
   object=outputAll,
   metadata=cellTypeLabel,
   col.name = "cellType")

outputAll<-getMoreStats(outputAll)

VlnPlot(outputAll,
        pt.size = 0.1,
        features=c("nFeature_RNA","nFeature_aTAR","nFeature_uTAR"),
        group.by = "cellType",
        cols=as.vector(alphabet(n=length(unique(Idents(outputAll))))))
VlnPlot(outputAll,
        pt.size=0.1,
        features=c("nCount_RNA","nCount_aTAR","nCount_uTAR"),
        group.by = "cellType",
        cols=as.vector(alphabet(n=length(unique(Idents(outputAll))))))
VlnPlot(outputAll,
        pt.size=0.1,
        features=c("percent.outHMM"),
        group.by = "cellType",
        cols=as.vector(alphabet(n=length(unique(Idents(outputAll))))))+NoLegend()

Idents(outputAll)<-outputAll$cellType

##### plotting umap and violin plot below
SS2figs<-createFig2(outputAll,60,addPoints=T)

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/SS2testes.pdf",
    width=(8/6)*4, height=((8/6)), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(SS2figs[[1]],SS2figs[[3]],SS2figs[[2]],SS2figs[[4]],nrow=1,align="b")
dev.off()


Idents(outputAll)<-outputAll$cellType
allMarkersSpleen <- FindAllMarkers(outputAll, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
geneMarkersSpleen<-allMarkersSpleen[allMarkersSpleen$gene %in% geneOnly,]
aTARMarkersSpleen<-allMarkersSpleen[allMarkersSpleen$gene %in% inGene,]
uTARMarkersSpleen<-allMarkersSpleen[allMarkersSpleen$gene %in% outGene,]
top3gene <- geneMarkersSpleen %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DoHeatmap(outputAll,features = unique(top3gene$gene))+NoLegend()
TARMarkersSpleen<-allMarkersSpleen[allMarkersSpleen$gene %in% c(inGene,outGene),]
top3TAR<-TARMarkersSpleen %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DoHeatmap(outputAll,features = unique(top3TAR$gene))+NoLegend()

top3uTAR<-uTARMarkersSpleen %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(outputAll,features = unique(top3uTAR$gene))+
   NoLegend()+
   coord_flip()+
   RotatedAxis()

###############################################
# go from diff uTARs to blast result
top3uTAR<-top3uTAR[!duplicated(top3uTAR$gene),]
test<-createOutGeneDf2(top3uTAR) ################### have this function ready
gonadUTarBed<-test
gonadUTarBed$fasta<-paste0(gonadUTarBed$chr,":",gonadUTarBed$start,"-",gonadUTarBed$end)
fastaFile<-"/workdir/fw262/mouseLemur/lemur_annotations/genome/ncbi-genomes-2019-10-23/GCF_000165445.2_Mmur_3.0_genomic.fna"
write.table(gonadUTarBed$fasta,file="uTARFasta.txt",quote = F, row.names = F,col.names = F) # save fasta format
extractFaCmd<-paste0("samtools faidx --region-file uTARFasta.txt ",fastaFile," > fullUTar.fa") # extract fasta sequence using faidx
system(extractFaCmd)
system("blastn -db /shared_data/genome_db/BLAST_NCBI/nt -query fullUTar.fa -out tt.txt -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 5 -num_threads 50")
blastResult <- read.delim("tt.txt", header=FALSE)

gonadUTarBed<-addInBlastResult(blastResult,gonadUTarBed)

DotPlot(outputAll, features = gonadUTarBed$TAR)+
   coord_flip()+RotatedAxis()+
   scale_x_discrete(breaks=gonadUTarBed$TAR,
                    labels=gonadUTarBed$blastShort)+
   ylab("cell types")+xlab("uTARs")

# load up coverage plots made in tesla server
load("/workdir/fw262/uTARAnalysisTool/TAR-scRNA-seq/smartseq2/mouseLemur/testesSS2/coverageForPaper.RData")
plot1<-plotDiffExpNolab("NC-033688.1-8631349-8633199-+-14-0",bamfile,"CCSER1",normalChrom = F,order = F) #
plot2<-plotDiffExpNolab("NC-033668.1-42790899-42792849---30-0",bamfile,"LOC108528520",normalChrom = F,order = F) #
plot3<-plotDiffExpNolab("NC-033661.1-52305649-52307099---31-0",bamfile,"LOC105867359",normalChrom = F,order = F) #
plot4<-plotDiffExpNolab("NC-033673.1-14169049-14169249-+-2-0",bamfile,"LOC112841550",normalChrom = F,order = F) #
plot5<-plotDiffExpNolab("NC-033683.1-6072949-6075199---38-0",bamfile,"TBXAS1",normalChrom = F,order = F) #

# make plots for paper
feat1<-featurePlotForPaper(outputAll,"NC-033688.1-8631349-8633199-+-14-0",T)+NoLegend()
feat2<-featurePlotForPaper(outputAll,"NC-033668.1-42790899-42792849---30-0",T)+NoLegend()
feat3<-featurePlotForPaper(outputAll,"NC-033661.1-52305649-52307099---31-0",T)+NoLegend()
feat4<-featurePlotForPaper(outputAll,"NC-033673.1-14169049-14169249-+-2-0",T)+NoLegend()
feat5<-featurePlotForPaper(outputAll,"NC-033683.1-6072949-6075199---38-0",T)+NoLegend()
gonadUTarBed$blastnResults<-gonadUTarBed$blastShort
rownames(gonadUTarBed)<-gonadUTarBed$TAR
testesDot<-plotDotPlotWithLabel(outputAll,gonadUTarBed,dotScale = 1.5)
testesUmap<-plotUmapColorUmap(outputAll,reducUse = "umap_gene")

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/testesSS2Umaps2.pdf",
    width=6, height=3/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(testesUmap,feat1,feat2,feat3,feat4,feat5,nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/testesSS2exCov2.pdf",
    width=5, height=2/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(plot1,plot2,plot3,plot4,plot5,nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/testesSS2Dot2.pdf",
    width=1.2, height=1+(2/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(testesDot+NoLegend(),nrow=1,align="hv")
dev.off()


analysisOnPlasma<-function(sample_combined,res_param=0.2){
   library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
   library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
   library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
   #source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
   library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
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
   
   sample_combined <- RunPCA(object = sample_combined, features = outGene)
   sample_combined <- FindNeighbors(object=sample_combined, dims=1:10)
   sample_combined <- FindClusters(object=sample_combined, resolution=res_param)
   sample_combined <- RunUMAP(sample_combined,dims=1:10,check_duplicates = F)
   return(sample_combined)
}

##############################################################################
##############################################################################
### look at higher resolution spleen TAR data
geneSpleen<-system("ls -d /workdir/fw262/mouseLemur/dge_final/*ANTOINE_SPLEEN*gene*",intern = T)
tarSpleen<-system("ls -d /fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/newSpleenData/MLCA_ANTOINE_SPLEEN_S10_MW_merged_TARTagged.HMM.txt.gz",intern = T)

seuratSpleen<-generateSeuratForMouseLemur(hmmMat=tarSpleen,geneMat=geneSpleen)

# load annotations
mouseLemurMeta10X <- read.csv("/workdir/fw262/mouseLemur/LCA_SS2_10x_all_metadata_update.csv")
library(readxl)
sameNamesMeta <- data.frame(read_excel("/workdir/fw262/mouseLemur/sameNamesMeta.xlsx",col_names = FALSE))
rownames(sameNamesMeta)<-sameNamesMeta[,1]
sameNamesMeta[is.na(sameNamesMeta)]<-""
# add in conversion to meta file
mouseLemurMeta10X<-mouseLemurMeta10X[mouseLemurMeta10X$method!="smartseq2",]
mouseLemurMeta10X$mwSample<-sameNamesMeta[as.character(mouseLemurMeta10X$channel),2]

mouseLemurMeta10X<-mouseLemurMeta10X[startsWith(as.character(mouseLemurMeta10X$X),"Antoine_Spleen_10X"),]
# get cells list
spleenCells<-getNFromList(strsplit(as.character(mouseLemurMeta10X$X),"_"),4)
seuratSpleen<-subset(seuratSpleen,cells = spleenCells)

cellAnno<-as.character(mouseLemurMeta10X$cell_ontology_class_update)
names(cellAnno)<-as.character(mouseLemurMeta10X$X)

# add cell type metadata
seuratSpleen<-AddMetaData(
   object=seuratSpleen,
   metadata=cellAnno,
   col.name = "cellType")

seuratSpleen <- NormalizeData(seuratSpleen)
seuratSpleen <- ScaleData(seuratSpleen, features = rownames(seuratSpleen))
seuratSpleen <- analysisOnFilSeurat(seuratSpleen)

Idents(seuratSpleen)<-seuratSpleen$cellType
geneP<-DimPlot(seuratSpleen,reduction="umap_gene",cols=as.vector(alphabet(n=length(unique(Idents(seuratSpleen))))),pt.size = 0.5)+NoLegend()+coord_fixed()+theme(legend.position = "right")
uTARP<-DimPlot(seuratSpleen,reduction="umap_uTAR",cols=as.vector(alphabet(n=length(unique(Idents(seuratSpleen))))),pt.size = 0.5)+NoLegend()#+coord_fixed()+NoLegend()
plot_grid(geneP,uTARP)

seuratSpleenFigs<-createFig2(seuratSpleen,20,addPoints=T)
pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/mouseLemurSpleen.pdf",
    width=(8/6)*4, height=((8/6)), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(seuratSpleenFigs[[1]],seuratSpleenFigs[[2]],seuratSpleenFigs[[4]],nrow=1,align="b",rel_widths = c(1,1,2))
dev.off()
pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/mouseLemurSpleenLabels.pdf",
    width=(8/6)*4, height=((8/6)*4), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(geneP)
dev.off()

# subset for plasma cells
spleenPlasma<-subset(seuratSpleen,idents = "plasma cell")
plasmaMat<-GetAssayData(object = spleenPlasma, slot = "counts")
spleenPlasma<-CreateSeuratObject(plasmaMat)
spleenPlasma<-analysisOnFilSeurat(spleenPlasma)
Idents(spleenPlasma)<-spleenPlasma$uTAR_snn_res.0.2
geneP<-umapPlainLegend(spleenPlasma,"umap_gene")
uTARP<-umapPlainLegend(spleenPlasma,"umap_uTAR")
spleenPlasmaMarkers <- FindAllMarkers(spleenPlasma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

numdash<-str_count(rownames(spleenPlasma),"-")
outInd<-numdash>=5 & endsWith(rownames(spleenPlasma),"-0")
inInd<-numdash>=5 & !endsWith(rownames(spleenPlasma),"-0")
geneInd<-!(outInd|inInd)
geneOnly<-rownames(spleenPlasma)[geneInd]
inGene<-rownames(spleenPlasma)[inInd]
outGene<-rownames(spleenPlasma)[outInd]

spermatoMarkersOut<-spleenPlasmaMarkers[spleenPlasmaMarkers$gene %in% outGene,]

top3uTAR<-spermatoMarkersOut #%>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#FeaturePlot(spleenPlasma,features="NC-033672.1-80355499-80356749---297323-0")
DotPlot(spleenPlasma,features = unique(top3uTAR$gene))+
   #NoLegend()+
   coord_flip()+
   RotatedAxis()

top3uTAR<-top3uTAR[!duplicated(top3uTAR$gene),]
test<-createOutGeneDf2(top3uTAR) 
gonadUTarBed<-test
gonadUTarBed$fasta<-paste0(gonadUTarBed$chr,":",gonadUTarBed$start,"-",gonadUTarBed$end)
write.table(gonadUTarBed$fasta,file="uTARFasta.txt",quote = F, row.names = F,col.names = F) # save fasta format
fastaFile<-"/workdir/fw262/mouseLemur/lemur_annotations/genome/ncbi-genomes-2019-10-23/GCF_000165445.2_Mmur_3.0_genomic.fna"
extractFaCmd<-paste0("samtools faidx --region-file uTARFasta.txt ",fastaFile," > fullUTar.fa") # extract fasta sequence using faidx
system(extractFaCmd)
system("blastn -db /shared_data/genome_db/BLAST_NCBI/nt -query fullUTar.fa -out tt.txt -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 5 -num_threads 50")
blastResult <- read.delim("tt.txt", header=FALSE)

gonadUTarBed<-addInBlastResult(blastResult,gonadUTarBed)
gonadUTarBed["NC-033663.1-32390449-32393449---6954-0","blastShort"]<-"IGKV2-28" # high bitscore and e-value, Macaca mulatta immunoglobulin locus, # VK2 in rhesus
gonadUTarBed["NC-033663.1-32334449-32339049---5676-0","blastShort"]<-"IGKV2-28" # high bitscore and e-value, Macaca mulatta immunoglobulin locus, # VK2 in rhesus
gonadUTarBed["NC-033669.1-91222449-91227099-+-74229-0","blastShort"]<-"red fox ncRNA LOC112927970" # high bitscore and e-value, Macaca mulatta immunoglobulin locus, # VK2 in rhesus

featuresToPlot<-gonadUTarBed[gonadUTarBed$blastShort!="undetermined",]
featuresToPlot<-featuresToPlot[order(featuresToPlot$avg_logFC,decreasing = T),]
# label featuresToPlot based on igBLAST results

vioP<-VlnPlot(spleenPlasma,features=c("NC-033663.1-32390449-32393449---6954-0",
                                      "NC-033663.1-32334449-32339049---5676-0"))+theme(plot.title = element_blank())+NoLegend()
feaP<-FeaturePlot(spleenPlasma,order = T,reduction = "umap_uTAR",features=c("NC-033663.1-32390449-32393449---6954-0",
                                          "NC-033663.1-32334449-32339049---5676-0"))#+theme(plot.title = element_blank())+NoLegend()

featuresToPlot$blastnResults<-featuresToPlot$blastShort
rownames(featuresToPlot)<-featuresToPlot$TAR
featuresToPlot<-featuresToPlot[order(featuresToPlot$pct.2),]
featuresToPlot$blastnResults<-featuresToPlot$blastShort
rownames(featuresToPlot)<-featuresToPlot$TAR
plasmaDot<-plotDotPlotWithLabel(spleenPlasma,featuresToPlot[1:4,],dotScale = 1.5)

featurePlotForPaperPlasma<-function(input,featureIn,orderPoints=T){
   out<-FeaturePlot(input,
                    pt.size=0.1,
                    features = featureIn,
                    reduction = "umap_uTAR",
                    cols=c("light grey","dark red"),
                    order=orderPoints)+
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks=element_blank(),
            axis.title = element_blank(),
            plot.title = element_blank(),
            #legend.text=element_text(size=12),
            #legend.position="bottom",
            aspect.ratio=1)+
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
   out<-addSmallLegendDot(out,spaceLegend =0.5,textSize = 8)
   return(out)
}

spleenfeat1<-featurePlotForPaperPlasma(spleenPlasma,featuresToPlot$TAR[1],T)+NoLegend()
spleenfeat2<-featurePlotForPaperPlasma(spleenPlasma,featuresToPlot$TAR[2],T)+NoLegend()
spleenfeat3<-featurePlotForPaperPlasma(spleenPlasma,featuresToPlot$TAR[3],T)+NoLegend()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/plasmaUmaps.pdf",
    width=5, height=(3/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(geneP,uTARP,spleenfeat1,spleenfeat2,spleenfeat3,nrow=1,align="b")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/plasmaDot.pdf",
    width=1.5, height=1+(2/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(plasmaDot,nrow=1,align="hv")
dev.off()

# load up coverage
load("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/coverageForPaperSpleen.RData")
pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/SpleenexCov2.pdf",
    width=3, height=2/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(plot1,plot2,plot3,nrow=1,align="hv")
dev.off()



