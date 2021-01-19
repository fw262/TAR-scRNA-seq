# Libraries
library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(pals);library(VennDiagram)
library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
s#ource("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/functionsToAnalyse.R")
source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr")
#scanorama<-import('scanorama')
setwd("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei")

############ 
# chicken data
############
d4gene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99532_HGH27BGXB_1_D4_S1/99532_HGH27BGXB_1_D4_S1_expression_matrix.txt.gz")
d4TAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99532_HGH27BGXB_1_D4_S1/99532_HGH27BGXB_1_D4_S1_expression_matrix_HMM_withDir.txt.gz")
d4<-createSeuratChicken(d4gene,d4TAR)
d7RVgene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99534_HGH27BGXB_3_RV_S3/99534_HGH27BGXB_3_RV_S3_expression_matrix.txt.gz")
d7RVTAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99534_HGH27BGXB_3_RV_S3/99534_HGH27BGXB_3_RV_S3_expression_matrix_HMM_withDir.txt.gz")
d7RV<-createSeuratChicken(d7RVgene,d7RVTAR)
d10RVgene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days10_v2/results_chicken/72255_HGGFCBGX5_RV_S2_all/72255_HGGFCBGX5_RV_S2_all_expression_matrix.txt.gz")
d10RVTAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days10_v2/results_chicken/72255_HGGFCBGX5_RV_S2_all/72255_HGGFCBGX5_RV_S2_all_expression_matrix_HMM_withDir.txt.gz")
d10RV<-createSeuratChicken(d10RVgene,d10RVTAR)
d14RVgene<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days14_v2/results_chicken/75019_HWYLJBGX5_RV_S2_all/75019_HWYLJBGX5_RV_S2_all_expression_matrix.txt.gz")
d14RVTAR<-("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days14_v2/results_chicken/75019_HWYLJBGX5_RV_S2_all/75019_HWYLJBGX5_RV_S2_all_expression_matrix_HMM_withDir.txt.gz")
d14RV<-createSeuratChicken(d14RVgene,d14RVTAR)

d4withDir   <-analysisOnFilSeurat(d4)
d7RVwithDir <-analysisOnFilSeurat(d7RV)
d10RVwithDir<-analysisOnFilSeurat(d10RV)
d14RVwithDir<-analysisOnFilSeurat(d14RV)

d4withDir_fig<-createFig2(d4withDir)
d7RVwithDir_fig<-createFig2(d7RVwithDir)
d10RVwithDir_fig<-createFig2(d10RVwithDir)
d14RVwithDir_fig<-createFig2(d14RVwithDir)

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/d4withDir.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(d4withDir_fig[[1]],d4withDir_fig[[2]],d4withDir_fig[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/d7withDir.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(d7RVwithDir_fig[[1]],d7RVwithDir_fig[[2]],d7RVwithDir_fig[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/d10withDir.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(d10RVwithDir_fig[[1]],d10RVwithDir_fig[[2]],d10RVwithDir_fig[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/d14withDir.pdf",
    width=(8/7), height=((8/7)*2.5), paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(d14RVwithDir_fig[[1]],d14RVwithDir_fig[[2]],d14RVwithDir_fig[[4]],ncol=1,align="vh",rel_heights = c(1,1,0.5))
dev.off()

# calculate silhouette values
silValsDir<-lapply(X=list(d4withDir,
                       d7RVwithDir,
                       d10RVwithDir,
                       d14RVwithDir),FUN=calcSilfromUMAPMean)
meanSilVal<-data.frame(matrix(unlist(silValsDir), nrow=length(silValsDir), byrow=T))
colnames(meanSilVal)<-c("genes","aTARs","uTARs")
rownames(meanSilVal)<-c("d4","d7","d10","d14")
meanSilVal$experiment<-rownames(meanSilVal)#<-factor(rownames(meanSilVal),levels = c("hg38","hg16","chicken d4","chicken d7","chicken d10","chicken d14","lemur lungs","lemur testis","mole rat","sea urchin"))
dataToPlot<-melt(meanSilVal,id.vars = "experiment")
dataToPlot$value<-as.numeric(dataToPlot$value)
dataToPlot$experiment<-factor(dataToPlot$experiment,
                              levels = c("d4",
                                         "d7",
                                         "d10",
                                         "d14"))

silValBar<-ggplot(dataToPlot,aes(x=experiment,y=value,fill=variable))+
  geom_bar(position="dodge",stat="identity")+
  xlab("sample")+ylab("silhouette coefficient")+labs(fill="features")+
  scale_fill_manual(values=c("#0000B2","#FF0000","#640000"))
silValBar<-adjustThemeGG(silValBar)
silValBar<-addSmallLegend(silValBar)

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/silCoeff2_withDir.pdf",
    width=2.5, height=2, paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(silValBar,nrow=1)
dev.off()

# find differential marker
d4withDirMarkers<-FindDiffExp(d4withDir)
d4DirUtar<-FindDiffExpUTAR(d4withDirMarkers)
d4DirUtar<-d4DirUtar[!endsWith(rownames(d4DirUtar),"-01"),]#remove duplicates

d4DirUtar$blastnResults[rownames(d4DirUtar)=="3-59669599-59692249-+-5940-0"]<-"CENPW (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-65683949-65711899-+-5762-0"]<-"SOX5"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-52651049-52688749---6170-0"]<-"LOC107053589"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="24-3293949-3349649---5242-0"]<-"LOC101749178"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="2-102630649-102657549---9144-0"]<-"GATA6 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="12-2836549-2872299-+-36062-0"]<-"DAG1 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="2-124878999-124895699---3773-0"]<-"RUNX1"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="21-1246149-1247899-+-169459-0"]<-"LOC419389"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-109375749-109386749-+-333313-0"]<-"SH3BGR"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-34406949-34434949---64324-0"]<-"HMGA2 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-34313499-34346099---15667-0"]<-"HMGA2 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="2-65517999-65537899---63139-0"]<-"LYRM4 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-26785099-26808299-+-81344-0"]<-"LOC107052650"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="4-30658399-30689599---24422-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="2-98071749-98154149-+-10705-0"]<-"MOG"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="27-5429599-5434549---3614-0"]<-"SLC4A1 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="5-37573549-37577599---6394-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="2-39103949-39146049---165236-0"]<-"RBMS3 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="3-47571549-47666099-+-28210-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="2-53292399-53421599-+-35029-0"]<-"CNTNAP2"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="4-57304549-57388099---9473-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="3-77228699-77246349---3383-0"]<-"TBX18 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="6-4454149-4522999---5590-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="4-5778899-5816449---29933-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="5-13398249-13417749-+-25439-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="20-4666849-4698049-+-11175-0"]<-"MAFB (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="14-14579349-14597099-+-13435-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="7-20808699-20938349-+-10524-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="13-9966349-10024349-+-3030-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="13-9912199-9965749-+-11435-0"]<-"STC2 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-125853499-125915949-+-4010-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="3-1536699-1613799---4759-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="3-59941549-59977949-+-5296-0"]<-"HEY2 (-)"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="32-583299-584549---137241-0"]<-"unknown"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-65683949-65711899-+-5762-0"]<-"SOX5"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="2-124878999-124895699---3773-0"]<-"RUNX1T1"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="21-1246149-1247899-+-169459-0"]<-"LOC419389"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-109375749-109386749-+-333313-0"]<-"SH3BGR"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-26785099-26808299-+-81344-0"]<-"LOC107052650"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="4-30658399-30689599---24422-0"]<-"ANAPC10"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="5-37573549-37577599---6394-0"]<-"CLEC14A"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="3-47571549-47666099-+-28210-0"]<-"SASH1"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="4-57304549-57388099---9473-0"]<-"LOC106038660"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="6-4454149-4522999---5590-0"]<-"NRG3"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="4-5778899-5816449---29933-0"]<-"DIAPH2"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="5-13398249-13417749-+-25439-0"]<-"CDKN1C"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="14-14579349-14597099-+-13435-0"]<-"LOC104913510"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="7-20808699-20938349-+-10524-0"]<-"KCNH7"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="13-9966349-10024349-+-3030-0"]<-"RPL26L1"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="1-125853499-125915949-+-4010-0"]<-"MID1"
d4DirUtar$blastnResults[rownames(d4DirUtar)=="3-1536699-1613799---4759-0"]<-"LOC110120668"

cluster_labels<-c("epicardial progenitor",
                  "myocytes",
                  "red blood cells",
                  "endothelial",
                  "immature fibroblasts",
                  "immune",
                  "other")
names(cluster_labels) <- levels(d4withDir)
d4withDir<-RenameIdents(d4withDir,cluster_labels)

cluster_labels<-c("fibroblasts",
                  "vascular endothelial",
                  "myocytes",
                  "endothelial",
                  "smooth muscle cells",
                  "epicardial progenitor",
                  "erythrocytes")
names(cluster_labels) <- levels(d14RVwithDir)
d14RVwithDir<-RenameIdents(d14RVwithDir,cluster_labels)

d4DotDir<-plotDotPlotWithLabelSupp(d4withDir,d4DirUtar)
d7DotDir<-plotDotPlotWithLabelSupp(d7RVwithDir,d4DirUtar)
d10DotDir<-plotDotPlotWithLabelSupp(d10RVwithDir,d4DirUtar)
d14DotDir<-plotDotPlotWithLabelSupp(d14RVwithDir,d4DirUtar)

plotDotPlotWithLabelSupp<-function(seuratOb,markers,dotScale=4){
  dotPlotOutF<-DotPlot(seuratOb, features = rownames(markers),dot.scale = dotScale)+
    coord_flip()+RotatedAxis()+
    scale_x_discrete(breaks=rownames(markers),
                     labels=markers$blastnResults)+
    labs(x="uTAR features")+
    theme(legend.position="right",
          legend.direction = "vertical",
          legend.text=element_text(size=6),
          legend.title = element_text(size=6),
          legend.spacing.x=unit(1,"mm"),
          legend.spacing.y=unit(1,"mm"),
          axis.text.y=element_text(size=6),
          axis.text.x=element_text(size=6,
                                   hjust=1,
                                   angle=30,
                                   face = "bold",
                                   colour=as.vector(alphabet(n=length(unique(Idents(seuratOb)))))),
          axis.title.y=element_blank(),#element_text(size=10, face="bold"),
          axis.title.x=element_blank())
  dotPlotOutF<-addSmallLegendDot(dotPlotOutF,spaceLegend =0.5,textSize = 8)+ #+theme(aspect.ratio = 1)
    theme(legend.title = element_blank())
  
  return(dotPlotOutF)
  
}

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/withDirDot.pdf",
    width=4, height=3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(d4DotDir,d14DotDir,align="hv",nrow=1)
dev.off()

#### plot sense and antisense

d4anti<-d4DirUtar[endsWith(d4DirUtar$blastnResults,"(-)"),]
d4anti<-d4anti[!duplicated(d4anti$blastnResults),]

d4sense<-d4anti
rownames(d4sense)<-getNFromList(str_split(d4sense$blastnResults," "),1)
d4sense$blastnResults<-rownames(d4sense)

d4antiD<-plotDotPlotWithLabelSupp(d4withDir,d4anti)
d4senseD<-plotDotPlotWithLabelSupp(d4withDir,d4sense)

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/senseAntisense.pdf",
    width=4, height=3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats = F)
plot_grid(d4antiD,d4senseD,align="hv",nrow=1)
dev.off()
