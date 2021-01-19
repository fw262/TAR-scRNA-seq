# Libraries
library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(pals);library(VennDiagram)
library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
#source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/functionsToAnalyse.R")
source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr")
#scanorama<-import('scanorama')
setwd("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei")

# humans data
############ 
filePath38<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg38_forPaper/hg38_pbmc_major/pbmc4k/"
filePath16<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg16/hg16_pbmc_major/pbmc4k/"
geneMat16<-paste0(filePath16,"pbmc4k_expression_matrix.txt.gz")
noDirMat16<-paste0(filePath16,"pbmc4k_expression_matrix_HMM_noDir_merge500_30reads.txt.gz")
geneMat38<-paste0(filePath38,"pbmc4k_expression_matrix.txt.gz")
noDirMat38<-paste0(filePath38,"pbmc4k_expression_matrix_HMM_noDir.txt.gz")

hg38<-createSeurathg38(geneMat38,noDirMat38)
hg16<-createSeuarthg16(geneMat16,noDirMat16,colnames(hg38))
############
# other human data
hg17gene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg17/hg17_pbmc_major/pbmc4k/pbmc4k_expression_matrix.txt.gz")
hg17TAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg17/hg17_pbmc_major/pbmc4k/pbmc4k_expression_matrix_HMM_noDir_merge500_31reads.txt.gz")
hg17<-createSeuarthg16(hg17gene,hg17TAR,colnames(hg38))
hg18gene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg18/hg18_pbmc_major/pbmc4k/pbmc4k_expression_matrix.txt.gz")
hg18TAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg18/hg18_pbmc_major/pbmc4k/pbmc4k_expression_matrix_HMM_noDir_merge500_30reads.txt.gz")
hg18<-createSeuarthg16(hg18gene,hg18TAR,colnames(hg38))
hg19gene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg19/hg19_pbmc_major/pbmc4k/pbmc4k_expression_matrix.txt.gz")
hg19TAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg19/hg19_pbmc_major/pbmc4k/pbmc4k_expression_matrix_HMM_31reads.txt.gz")
hg19<-createSeuarthg16(hg19gene,hg19TAR,colnames(hg38))

############ 
# chicken data
############
d4gene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99532_HGH27BGXB_1_D4_S1/99532_HGH27BGXB_1_D4_S1_expression_matrix.txt.gz")
d4TAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99532_HGH27BGXB_1_D4_S1/99532_HGH27BGXB_1_D4_S1_expression_matrix_HMM_noDir.txt.gz")
d4<-createSeuratChicken(d4gene,d4TAR)
d7RVgene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99534_HGH27BGXB_3_RV_S3/99534_HGH27BGXB_3_RV_S3_expression_matrix.txt.gz")
d7RVTAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99534_HGH27BGXB_3_RV_S3/99534_HGH27BGXB_3_RV_S3_expression_matrix_HMM_noDir.txt.gz")
d7RV<-createSeuratChicken(d7RVgene,d7RVTAR)
d10RVgene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days10_v2/results_chicken/72255_HGGFCBGX5_RV_S2_all/72255_HGGFCBGX5_RV_S2_all_expression_matrix.txt.gz")
d10RVTAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days10_v2/results_chicken/72255_HGGFCBGX5_RV_S2_all/72255_HGGFCBGX5_RV_S2_all_expression_matrix_HMM_noDir.txt.gz")
d10RV<-createSeuratChicken(d10RVgene,d10RVTAR)
d14RVgene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days14_v2/results_chicken/75019_HWYLJBGX5_RV_S2_all/75019_HWYLJBGX5_RV_S2_all_expression_matrix.txt.gz")
d14RVTAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days14_v2/results_chicken/75019_HWYLJBGX5_RV_S2_all/75019_HWYLJBGX5_RV_S2_all_expression_matrix_HMM_noDir.txt.gz")
d14RV<-createSeuratChicken(d14RVgene,d14RVTAR)
############
# mole rat
moleRatGene<-("/workdir/fw262/moleRat/results_moleRat/nmr_1.1/nmr_1.1_expression_matrix.txt.gz")
moleRatTAR<-("/workdir/fw262/moleRat/results_moleRat/nmr_1.1/nmr_1.1_expression_matrix_HMM_noDir.txt.gz")
moleRat<-createSeuratMoleRat(moleRatGene,moleRatTAR)
############
# sea urchin
urchinGene<-("/workdir/fw262/seaUrchin/results_seaUrchin_fromFastq/D1/D1_expression_matrix.txt.gz")
urchinTAR<-("/workdir/fw262/seaUrchin/results_seaUrchin_fromFastq/D1/D1_expression_matrix_HMM_noDir.txt.gz")
urchin<-createSeuratUrchin(urchinGene,urchinTAR)
############
# mouse lemur testis
testisGene<-("/workdir/fw262/mouseLemur/testis/results_testis/MLCA_ANTOINE_TESTIS_S7/MLCA_ANTOINE_TESTIS_S7_expression_matrix.txt.gz")
testisTAR<-("/workdir/fw262/mouseLemur/testis/results_testis/MLCA_ANTOINE_TESTIS_S7/MLCA_ANTOINE_TESTIS_S7_expression_matrix_HMM_noDir_26reads.txt.gz")
testis<-createSeuratTestis(testisGene,testisTAR)
############
# mouse lemur lungs
lungsCD31Gene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/olga/results_lungs/CD31_gene_exon_tagged.txt.gz")
lungsCD31TAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/olga/results_lungs/CD31_HMM_tagged_noDir_86reads.txt.gz")
lungsCD31<-createSeuratLung(lungsCD31Gene,lungsCD31TAR)
############
# mouse lemur lungs
lungsEpcamGene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/olga/results_lungs/epcam_gene_exon_tagged.txt.gz")
lungsEpcamTAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/olga/results_lungs/epcam_HMM_tagged_noDir_86reads.txt.gz")
lungsEpcam<-createSeuratLung(lungsEpcamGene,lungsEpcamTAR)
############

# generate list of seurats
# TIPS: run analysisOnFilSeurat after filtering and creating seurat object for each dataset
allSeurats<-list(hg38,
                 hg16,
                 d4,
                 d10RV,
                 d14RV,
                 moleRat,
                 urchin,
                 testis)
# run same analysis after generating seurat objects
allSeurats<-lapply(X=allSeurats,FUN=analysisOnFilSeurat)
names(allSeurats)<-c("hg38",
                     "hg16",
                     "d4",
                     "d10RV",
                     "d14RV",
                     "moleRat",
                     "urchin",
                     "testis")
allSeurats$d7RV<-analysisOnFilSeurat(d7RV)
allSeurats$hg17<-analysisOnFilSeurat(hg17)
allSeurats$hg18<-analysisOnFilSeurat(hg18)
allSeurats$hg19<-analysisOnFilSeurat(hg19)
allSeurats$lungsCD31<-analysisOnFilSeurat(lungsCD31)
allSeurats$lungsEpcam<-analysisOnFilSeurat(lungsEpcam)
allSeurats$hg38_4000<-analysisOnFilSeuratWithVarGene(hg38)
allSeurats$hg16_4000<-analysisOnFilSeuratWithVarGene(hg16)
save(allSeurats,file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/allDataDiffParam.RData")

load(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/allDataDiffParam.RData")

# set identities
for(i in allSeurats){
  Idents(i)<-i$RNA_snn_res.0.2
}
  

# plot figure 2 elements
load("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/mouseAtlasSpleen.RData")
load("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/mouseAtlasKidney.RData")

allSeurats$mouseAtlasSpleen<-mergedSeurat
allSeurats$mouseAtlasKidney<-seuratKidney
allSeuratsFig2<-lapply(X=allSeurats,FUN=createFig2)
names(allSeuratsFig2)<-names(allSeurats)
#allSeuratsFig2$mouseAtlas<-createFig2(mergedSeurat)
#allSeuratsFig2$hg38_4000<-createFig2(allSeurats$hg38_4000)
#allSeuratsFig2$hg16_4000<-createFig2(allSeurats$hg16_4000)

# calculate silhouette values
# make sure Idents(allSeurats) is set to cell type
silVals<-lapply(X=allSeurats,FUN=calcSilfromUMAPMeanFinal)#,tt=0)
meanSilVal<-data.frame(matrix(unlist(silVals), nrow=length(silVals), byrow=T))
colnames(meanSilVal)<-c("genes","aTARs","uTARs")
rownames(meanSilVal)<-names(allSeurats)
meanSilVal$experiment<-rownames(meanSilVal)#<-factor(rownames(meanSilVal),levels = c("hg38","hg16","chicken d4","chicken d7","chicken d10","chicken d14","lemur lungs","lemur testis","mole rat","sea urchin"))
#meanSilVal["mouseAtlas",]<-c(silValMouse,"mouseAtlas")
meanSilVal<-meanSilVal[c("mouseAtlasSpleen",
                         "mouseAtlasKidney",
                         "hg16_4000",
                         "hg38_4000",
                         "d4",
                         "d7RV",
                         "d10RV",
                         "d14RV",
                         "lungsEpcam",
                         "moleRat",
                         "urchin"),]
dataToPlot<-melt(meanSilVal,id.vars = "experiment")
dataToPlot$value<-as.numeric(dataToPlot$value)
dataToPlot$experiment<-factor(dataToPlot$experiment,
                              levels = c("mouseAtlasSpleen",
                                         "mouseAtlasKidney",
                                         "hg16_4000",
                                         "hg38_4000",
                                         "d4",
                                         "d7RV",
                                         "d10RV",
                                         "d14RV",
                                         "lungsEpcam",
                                         "moleRat",
                                         "urchin"))

silValBar<-ggplot(dataToPlot,aes(x=experiment,y=value,fill=variable))+
  geom_bar(position="dodge",stat="identity")+
  xlab("sample")+ylab("silhouette coefficient")+labs(fill="features")+
  scale_fill_manual(values=c("#0000B2","#FF0000","#640000"))
silValBar<-adjustThemeGG(silValBar)
silValBar<-addSmallLegend(silValBar)

############################################3
##############################################
##### look at pseudobulk counts
# load up true pseudo bulk values

allSeuratsBulkPca<-lapply(X=allSeurats[c("mouseAtlasSpleen",
                                         "mouseAtlasKidney",
                                         "hg16_4000",
                                         "hg38_4000",
                                         "d4",
                                         "d7RV",
                                         "d10RV",
                                         "d14RV",
                                         "lungsEpcam",
                                         "moleRat",
                                         "urchin")],
                          FUN=drawScatterTrueBulk2,pca="pca_uTAR",PCs = 5,numDashes = 5,minBulk = 0,returnData = T)
names(allSeuratsBulkPca)<-c("mouseAtlasSpleen",
                            "mouseAtlasKidney",
                            "hg16_4000",
                            "hg38_4000",
                            "d4",
                            "d7RV",
                            "d10RV",
                            "d14RV",
                            "lungsEpcam",
                            "moleRat",
                            "urchin")
#allSeuratsBulkPca$mouseAtlas<-drawScatterTrueBulk2(mergedSeurat,pca="pca_uTAR",PCs = 5,numDashes = 5,minBulk = 0,returnData = T)
bulkPCAComp<-do.call(rbind, unname(allSeuratsBulkPca))

# make scatter plot
stdDevP<-ggplot(bulkPCAComp,aes(x=sumExpUTar,y=sumPCLoadUTar))+
  geom_bin2d(bins=100)+scale_fill_continuous(type="viridis")+theme_bw()+scale_y_log10()+scale_x_log10()+
  labs(title="2D histogram of uTARs",y="uTAR PC loading",x="uTAR pseudo-bulk read coverage")
stdDevP<-stdDevP+
  theme(plot.title = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.text.x=element_text(size = 8,face="bold"),
        axis.text.y=element_text(size=8,face="bold"))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
stdDevP<-stdDevP+geom_hline(yintercept = 0.5)+geom_vline(xintercept =10000)
stdDevP<-addSmallLegend2(stdDevP,spaceLegend =0.5)#+theme(aspect.ratio = 1)

venDia<-drawVennDiagramFromData(bulkPCAComp,pcThresh=0.5, bulkThresh=10000)

# calculate R2 value
allScatter.lm = lm(sumPCLoadUTar ~ sumExpUTar, data=bulkPCAComp)
summary(allScatter.lm)$r.squared

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/silCoeff4.pdf",
    width=3, height=2, paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(silValBar+NoLegend(),nrow=1)
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/scatterPlot4.pdf",
    width=3, height=2, paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(stdDevP,nrow=1)
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/venDia4.pdf",
    width=2, height=2, paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(venDia,nrow=1)
dev.off()

# plot top of figure 2
pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/mouseAtlasSpleen.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$mouseAtlasSpleen[[1]],allSeuratsFig2$mouseAtlasSpleen[[2]],allSeuratsFig2$mouseAtlasSpleen[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/mouseAtlasKidney.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$mouseAtlasKidney[[1]],allSeuratsFig2$mouseAtlasKidney[[2]],allSeuratsFig2$mouseAtlasKidney[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()


pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/hg38_4000.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$hg38_4000[[1]],allSeuratsFig2$hg38_4000[[2]],allSeuratsFig2$hg38_4000[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/hg16_4000.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$hg16_4000[[1]],allSeuratsFig2$hg16_4000[[2]],allSeuratsFig2$hg16_4000[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()


pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/d4.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$d4[[1]],allSeuratsFig2$d4[[2]],allSeuratsFig2$d4[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/d14RV.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$d14RV[[1]],allSeuratsFig2$d14RV[[2]],allSeuratsFig2$d14RV[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/lungsEpcam.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$lungsEpcam[[1]],allSeuratsFig2$lungsEpcam[[2]],allSeuratsFig2$lungsEpcam[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/moleRat.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$moleRat[[1]],allSeuratsFig2$moleRat[[2]],allSeuratsFig2$moleRat[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/urchin.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$urchin[[1]],allSeuratsFig2$urchin[[2]],allSeuratsFig2$urchin[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

# plot supplement
pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/hg16.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$hg16[[1]],allSeuratsFig2$hg16[[2]],allSeuratsFig2$hg16[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/d7RV.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$d7RV[[1]],allSeuratsFig2$d7RV[[2]],allSeuratsFig2$d7RV[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/d10RV.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(allSeuratsFig2$d10RV[[1]],allSeuratsFig2$d10RV[[2]],allSeuratsFig2$d10RV[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()




#### look at differential gene expression and classify cells
diffMarkers<-lapply(X=allSeurats,FUN=FindDiffExp)

# filter for uTAR markers
diffMarkersUtar<-lapply(X=diffMarkers,FUN=FindDiffExpUTAR)

day7bam<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99534_HGH27BGXB_3_RV_S3/99534_HGH27BGXB_3_RV_S3_Aligned_sorted_2.bam"
day10bam<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days10_v2/results_chicken/72255_HGGFCBGX5_RV_S2_all/72255_HGGFCBGX5_RV_S2_all_Aligned_sorted_2.bam"
fastaFile<-"/fs/cbsuvlaminck2/workdir/References/Chick/GRCg6a/Gallus_gallus.GRCg6a.dna.toplevel.fa"


runBlastOnDiffMarkers<-function(markers,fastaFile,bamFile,normalChromT=T){
  test<-createOutGeneDf3(markers)
  test$dot<-"."
  test<-test[order(test$avg_logFC,decreasing = T),]
  #test<-test[!duplicated(test$gene),]
  cellTypeMarkersBed<-test[,c("chr","start","end","TAR","dot","strand","avg_logFC","pct.1","pct.2","p_val_adj","cluster","fasta")]
  system("mkdir -p temp")
  cellTypeMarkersCov<-lapply(FUN=getCoverageOnly,X=cellTypeMarkersBed[,4],bamFile=bamFile,normalChrom=normalChromT)
  
  #smooth coverage
  smoothCov<-lapply(X=cellTypeMarkersCov, FUN=smoothCovFunc,span=0.25)
  # get half width maximum of smoothed peak
  indicesToBlast<-lapply(X=smoothCov,FUN=findFHWM)
  
  # add index shift to create fasta for peaks
  cellTypeMarkersBed<-cbind(cellTypeMarkersBed,data.frame(matrix(unlist(indicesToBlast), nrow=length(indicesToBlast), byrow=T)))
  cellTypeMarkersBed$fastaPeak<-paste0(cellTypeMarkersBed$chr,":",cellTypeMarkersBed$start+cellTypeMarkersBed$X1,"-",cellTypeMarkersBed$start+cellTypeMarkersBed$X2)
  cellTypeMarkersBed$length<-cellTypeMarkersBed$X2-cellTypeMarkersBed$X1
  
  write.table(cellTypeMarkersBed$fastaPeak,file="cellTypeseqToBlast.txt",quote = F, row.names = F,col.names = F) # save fasta format
  extractFaCmd<-paste0("samtools faidx -r cellTypeseqToBlast.txt ",fastaFile," > cellTypeseqToBlast.fa") # extract fasta sequence using faidx
  system(extractFaCmd)
  
  #blast command
  system("blastn -db /shared_data/genome_db/BLAST_NCBI/nt -query cellTypeseqToBlast.fa -out tt.txt -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 5 -num_threads 10")
  
  return()

}

