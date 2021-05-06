library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library("pals")
library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
#source("/workdir/fw262/ShaoPei/allOrganisms/functionsToAnalyse.R")
source("/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")

########## loading paper ready data
readsData <- read.delim("/workdir/fw262/ShaoPei/paperScripts_wksp/figure1/assemblyInfo.txt")
rownames(readsData)<-readsData$assembly
readsData<-readsData[!(rownames(readsData) %in% "hg19"),] # remove 19
readsData<-readsData[!(rownames(readsData) %in% "horse"),] # remove horse

readsData$assembly<-factor(readsData$assembly,levels=c("hg38","hg16","mouse","chicken","mouse_lemur","naked_mole_rat","purple_sea_urchin"))
readsData$total.bases<-readsData$total.non.N.bases+readsData$repeating.bases
transcriptsP<-ggplot(readsData)+
              geom_bar(aes(x=factor(readsData$assembly),
                           y=readsData$transcripts/readsData$total.bases),
                       fill=as.vector(c("blue","#399cbd","#009600","#960000","#696969","black","#b08f26")),
                       stat="identity",
                       width=0.75)+
              NoLegend()+
              ylab("relative # of\nannotated transcripts")+
              theme_classic()+
              theme(axis.title.x=element_blank(),
                    axis.title.y=element_text(size=8,face="bold"),
                    axis.text.y=element_text(size=8,angle = 90,hjust = 0.5))+
              theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
transcriptsP<-adjustThemeGG(transcriptsP)
transcriptsP<-addSmallLegend(transcriptsP)+NoLegend()

readsData_m<-melt(readsData[,1:5],id.vars="assembly")
readsData_m<-readsData_m[grep("bases",readsData_m$variable),]
readsData_m$variable<-factor(readsData_m$variable,levels=c("transcripts","chromosomes","repeating.bases","total.non.N.bases"))
basesP<-ggplot(readsData)+
  geom_bar(aes(x=factor(readsData$assembly),y=readsData$total.bases),
           fill=as.vector(c("blue","#399cbd","#009600","#960000","#696969","black","#b08f26")),
           stat="identity",width=0.75)+
  NoLegend()+
  ylab("assembly length (bps)")+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=8,face="bold"),
        axis.text.y=element_text(size=8,angle = 90,hjust = 0.5))
basesP<-adjustThemeGG(basesP)
basesP<-addSmallLegend(basesP)+NoLegend()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

plot_grid(basesP,transcriptsP,ncol=2,nrow=1)

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/assemblyStats.pdf",
    width=3.5, height=2.2, paper="special", bg="white",
    fonts="Helvetica", pointsize=6,useDingbats=F)
plot_grid(basesP,transcriptsP,ncol=2,nrow=1)
dev.off()

# load up reads outside all organisms data, matching figure 2

hg38uTARper<-allSeurats$hg38_4000$percent.outHMM
hg16uTARper<-allSeurats$hg16_4000$percent.outHMM
mouseuTARper<-allSeurats$mouseAtlasSpleen$percent.outHMM
chickenuTARper<-c(allSeurats$d4$percent.outHMM,
                  allSeurats$d7RV$percent.outHMM,
                  allSeurats$d10RV$percent.outHMM,
                  allSeurats$d14RV$percent.outHMM)
lemuruTARper<-allSeurats$lungsEpcam$percent.outHMM
moleuTARper<-allSeurats$moleRat$percent.outHMM
urchinuTARper<-allSeurats$urchin$percent.outHMM
allUTARs<-rbind(data.frame(cbind("hg38",as.vector(hg38uTARper))),
      data.frame(cbind("hg16",as.vector(hg16uTARper))),
      data.frame(cbind("mouse",as.vector(mouseuTARper))),
      data.frame(cbind("chicken",as.vector(chickenuTARper))),           
      data.frame(cbind("mouse lemur",as.vector(lemuruTARper))),
      data.frame(cbind("mole rat",as.vector(moleuTARper))),
      data.frame(cbind("sea urchin",as.vector(urchinuTARper))))
allUTARs$X2<-as.numeric(levels(allUTARs$X2))[allUTARs$X2]

readsP<-ggplot(allUTARs,aes(x=X1,y=X2,fill=X1,color=X1))+
  geom_violin()+
  #facet_wrap(facets = readsToPlot$X2,ncol = length(unique(readsToPlot$X2)))+
  stat_summary(fun.data=mean_sdl,geom="pointrange",size=0.2,color="black")+
  ylab("reads outside of\nannotations (%)")+xlab("sample")+
  scale_fill_manual(values=c("blue","#399cbd","#009600","#960000","#696969","black","#b08f26"))+
  scale_color_manual(values=c("blue","#399cbd","#009600","#960000","#696969","black","#b08f26"))
readsP<-adjustThemeGG(readsP)
readsP<-addSmallLegend(readsP)+NoLegend()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/assemblyStats_readsOnly2.pdf",
    width=4/3, height=2, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats=F)
plot_grid(readsP,nrow=1,align="h")
dev.off()


pdf(file="/workdir/fw262/ShaoPei/paperFigs/assemblyStats_v2.pdf",
    width=4, height=2, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats=F)
plot_grid(basesP,transcriptsP,readsP,nrow=1,align="h")
dev.off()
