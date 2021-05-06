# Libraries
library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot)
library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
#source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/functionsToAnalyse.R")
source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")

# load up seurat objects
load(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/allDataDiffParam.RData")
# load up markers data
load(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure3/markersData_forPaper.RData")

#### fine bam files to look at coverage
hg16bam<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg16/hg16_pbmc_major/pbmc4k/pbmc4k_Aligned_sorted_2.bam"
hg38bam<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/humanPBMCHg38_forPaper/hg38_pbmc_major/pbmc4k/pbmc4k_Aligned_sorted_2.bam "
day4bam<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days4_7_v3/results_chicken/99532_HGH27BGXB_1_D4_S1/99532_HGH27BGXB_1_D4_S1_Aligned_sorted_2.bam"
day14bam<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/chickenAll/days14_v2/results_chicken/75019_HWYLJBGX5_RV_S2_all/75019_HWYLJBGX5_RV_S2_all_Aligned_sorted_2.bam"
lungbam<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/olga/results_lungs/epcam_HMM_tagged_noDir.bam"
testisbam<-"/workdir/fw262/mouseLemur/testis/results_testis/MLCA_ANTOINE_TESTIS_S7/MLCA_ANTOINE_TESTIS_S7_Aligned_sorted_2.bam"
molebam<-"/workdir/fw262/moleRat/results_moleRat/nmr_1.1/nmr_1.1_Aligned_sorted_2.bam"
urchinbam<-"/workdir/fw262/seaUrchin/results_seaUrchin_fromFastq/D1/D1_Aligned_sorted_2.bam"

lemurfasta<-"/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/olga/lemur_annotations/genome/ncbi-genomes-2019-10-23/GCF_000165445.2_Mmur_3.0_genomic.fna"
molefasta<-"/workdir/fw262/moleRat/genome/GCF_000247695.1_HetGla_female_1.0_genomic.fna"
urchinfasta<-"/workdir/fw262/seaUrchin/genome/genome_assemblies/ncbi-genomes-2019-11-21/GCF_000002235.4_Spur_4.2_genomic.fna"
############################################


# relabel cluster day 4
cluster_labels<-c("epicardial progenitor",
                  "myocytes",
                  "red blood cells",
                  "endothelial cells",
                  "immature fibroblasts",
                  "immune cells",
                  "other")
cluster_labels<-c("0",
                  "1",
                  "2",
                  "3",
                  "4",
                  "5",
                  "6")

names(cluster_labels) <- levels(allSeurats$d4)
allSeurats$d4<-RenameIdents(allSeurats$d4,cluster_labels)

day4MarkersOut<-diffMarkersUtar$d4
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="1-65683949-65711899-+-5762-0"]<-"SOX5"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="2-124878999-124895699---3773-0"]<-"RUNX1T1"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="21-1246149-1247899-+-169459-0"]<-"LOC419389"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="1-109375749-109386749-+-333313-0"]<-"SH3BGR"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="1-26785099-26808299-+-81344-0"]<-"LOC107052650"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="4-30658399-30689599---24422-0"]<-"ANAPC10"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="5-37573549-37577599---6394-0"]<-"CLEC14A"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="3-47571549-47666099-+-28210-0"]<-"SASH1"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="4-57304549-57388099---9473-0"]<-"LOC106038660"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="6-4454149-4522999---5590-0"]<-"NRG3"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="4-5778899-5816449---29933-0"]<-"DIAPH2"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="5-13398249-13417749-+-25439-0"]<-"CDKN1C"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="14-14579349-14597099-+-13435-0"]<-"LOC104913510"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="7-20808699-20938349-+-10524-0"]<-"KCNH7"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="13-9966349-10024349-+-3030-0"]<-"RPL26L1"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="1-125853499-125915949-+-4010-0"]<-"MID1"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="3-1536699-1613799---4759-0"]<-"LOC110120668"
day4MarkersOut$blastnResults[rownames(day4MarkersOut)=="32-583299-584549---137241-0"]<-"LOC107050595"
day4MarkersOut<-day4MarkersOut[!endsWith(rownames(day4MarkersOut),"-01"),]

d4Dot<-plotDotPlotWithLabel(allSeurats$d4,day4MarkersOut)

day4FeatOut<-findSiguTARmarkers(allSeurats$d4,day4MarkersOut,shift = 0,minpc.1 = 0.0,fcThresh = 0)
day4exp4<-plotDiffExpNolab(allSeurats$d4,"2-124878999-124895699---3773-0",day4bam,"RUNX1T1",normalChrom = T) # searched for
day4exp5<-plotDiffExpNolab(allSeurats$d4,"5-37573549-37577599---6394-0",day4bam,"CLEC14A",normalChrom = T) # searched for
day4exp6<-plotDiffExpNolab(allSeurats$d4,"1-65683949-65711899-+-5762-0",day4bam,"SOX5",normalChrom = T) # searched for
day4exp7<-plotDiffExpNolab(allSeurats$d4,"21-1246149-1247899-+-169459-0",day4bam,"LOC419389",normalChrom = T) # searched for
day4exp8<-plotDiffExpNolab(allSeurats$d4,"6-4454149-4522999---5590-0",day4bam,"NRG3",normalChrom = T) # searched for
day4exp9<-plotDiffExpNolab(allSeurats$d4,"13-9966349-10024349-+-3030-0",day4bam,"RPL26L1",normalChrom = T) # searched for

day4Umap<-plotUmapColorUmap(allSeurats$d4,reducUse = "umap_gene")

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/day4Umaps2.pdf",
    width=6, height=3/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(day4Umap,day4exp4[[2]],day4exp5[[2]],day4exp6[[2]],day4exp7[[2]],day4exp8[[2]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/day4exCov2.pdf",
    width=5, height=2/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(day4exp4[[1]],day4exp5[[1]],day4exp6[[1]],day4exp7[[1]],day4exp8[[1]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/day4Dot2.pdf",
    width=1.2, height=1+(2/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(d4Dot+NoLegend(),nrow=1,align="hv")
dev.off()
####################################################
# day 14 
d14Dot<-plotDotPlotWithLabel(allSeurats$d14RV,day4MarkersOut)
day14Umap<-plotUmapColorUmap(allSeurats$d14RV,reducUse = "umap_gene")
day14exp4<-plotDiffExpNolab(allSeurats$d14RV,"2-124878999-124895699---3773-0",day14bam,"RUNX1T1",normalChrom = T) # searched for
day14exp5<-plotDiffExpNolab(allSeurats$d14RV,"5-37573549-37577599---6394-0",day14bam,"CLEC14A",normalChrom = T) # searched for
day14exp6<-plotDiffExpNolab(allSeurats$d14RV,"1-65683949-65711899-+-5762-0",day14bam,"SOX5",normalChrom = T) # searched for
day14exp7<-plotDiffExpNolab(allSeurats$d14RV,"21-1246149-1247899-+-169459-0",day14bam,"LOC419389",normalChrom = T) # searched for
day14exp8<-plotDiffExpNolab(allSeurats$d14RV,"6-4454149-4522999---5590-0",day14bam,"NRG3",normalChrom = T) # searched for
day14exp9<-plotDiffExpNolab(allSeurats$d14RV,"13-9966349-10024349-+-3030-0",day14bam,"RPL26L1",normalChrom = T) # searched for

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/day14Umaps2.pdf",
    width=6, height=3/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(day14Umap,day14exp4[[2]],day14exp5[[2]],day14exp6[[2]],day14exp7[[2]],day14exp8[[2]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/day14exCov2.pdf",
    width=5, height=2/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(day14exp4[[1]],day14exp5[[1]],day14exp6[[1]],day14exp7[[1]],day14exp8[[1]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/day14Dot2.pdf",
    width=1.2, height=1+(2/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(d14Dot+NoLegend(),nrow=1,align="hv")
dev.off()
####################################################
# mole rat
moleRatMarkers<-findSiguTARmarkers2(allSeurats$moleRat,diffMarkersUtar$moleRat,shift = 1,fcThresh = 0.5,minpc.1 = 0.5)

#fcThresh = 0.5,minpc.1 = 0.5, positive only
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624733.1-17257599-17263349-+-240134-0"]<-"DUSP1"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624817.1-3730249-3739699---53834-0"]<-"LOC105860288"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624775.1-8146199-8147899---54887-0"]<-"Fth1"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624747.1-24035499-24038099---135965-0"]<-"LOC101722152"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624775.1-1545749-1551199-+-84462-0"]<-"LOC102005716"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624824.1-634899-641249-+-1065080-0"]<-"LOC104847932"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004625219.1-799-1199---191586-0"]<-"Rpl23a"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624743.1-6144399-6144949---378141-0"]<-"LOC101699547"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624788.1-5697449-5699299---126669-0"]<-"Rpl24"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624778.1-13208149-13222649---23668-0"]<-"Plxdc1"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624793.1-13075449-13089049---12919-0"]<-"LOC110343856"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624733.1-17202749-17211199---66206-0"]<-"Tsku"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624777.1-49-2399---43827-0"]<-"Mad2l1"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624849.1-4971349-4974149-+-5586-0"]<-"NATD1"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624771.1-9176449-9180799---9603-0"]<-"CUNH1orf105"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624740.1-23019849-23026049-+-6817-0"]<-"TRG"
moleRatMarkers$blastnResults[moleRatMarkers$gene=="NW-004624826.1-6519149-6520599-+-12931-0"]<-"XCL1"
moleRatMarkers<-moleRatMarkers[!endsWith(rownames(moleRatMarkers),"-01"),]

moleDot<-plotDotPlotWithLabel(allSeurats$moleRat,moleRatMarkers,dotScale = 1.25)

moleexp1<-plotDiffExpNolab(allSeurats$moleRat,"NW-004624733.1-17257599-17263349-+-240134-0",molebam,"DUSP1",normalChrom = F) # searched for
moleexp2<-plotDiffExpNolab(allSeurats$moleRat,"NW-004624775.1-8146199-8147899---54887-0",molebam,"Fth1",normalChrom = F) # searched for
moleexp3<-plotDiffExpNolab(allSeurats$moleRat,"NW-004624778.1-13208149-13222649---23668-0",molebam,"Plxdc1",normalChrom = F) # searched for
moleexp4<-plotDiffExpNolab(allSeurats$moleRat,"NW-004624849.1-4971349-4974149-+-5586-0",molebam,"NATD1",normalChrom = F) # searched for
moleexp5<-plotDiffExpNolab(allSeurats$moleRat,"NW-004624740.1-23019849-23026049-+-6817-0",molebam,"TRG",normalChrom = F) # searched for
moleUmap<-plotUmapColorUmap(allSeurats$moleRat,reducUse = "umap_gene")

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/moleUmaps2.pdf",
    width=6, height=3/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(moleUmap,moleexp1[[2]],moleexp2[[2]],moleexp3[[2]],moleexp4[[2]],moleexp5[[2]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/moleexCov2.pdf",
    width=5, height=2/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(moleexp1[[1]],moleexp2[[1]],moleexp3[[1]],moleexp4[[1]],moleexp5[[1]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/moleDot2.pdf",
    width=1.2, height=1+(2/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(moleDot+NoLegend(),nrow=1,align="hv")
dev.off()

###########################################
# mouse lemur epcam
#find lungs lemur example
############## differential marker analysis
lungMarkers<-findSiguTARmarkers2(allSeurats$lungsEpcam,diffMarkersUtar$lungsEpcam,shift = 1,fcThresh = 0.75,minpc.1 = 0.5)

# annotate uTARs fcThresh = 0.75,minpc.1 = 0.5, positive =T
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033669.1-3366899-3370949-")]<-"H3F3C"
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033678.1-40045849-40046349-")]<-"LOC105876721"
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033676.1-46824199-46826449-")]<-"RPLP1"
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033674.1-27742249-27745449-")]<-"RPS20"
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033664.1-8257299-8265799-")]<-"MS4A8"
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033683.1-1184549-1186249-")]<-"BST2"
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033687.1-7523249-7524299-")]<-"LOC105866554"
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033672.1-69686699-69689399-")]<-"TGFBRAP1"
lungMarkers$blastnResults[startsWith(lungMarkers$gene,"NC-033662.1-16611599-16612799-")]<-"IGHA"
lungMarkers<-lungMarkers[!endsWith(rownames(lungMarkers),"-01"),]

lungDot<-plotDotPlotWithLabel(allSeurats$lungsEpcam,lungMarkers,dotScale = 1.5)

lungexp1<-plotDiffExpNolab(allSeurats$lungsEpcam,"NC-033669.1-3366899-3370949-+-637924-0",lungbam,"H3F3C",normalChrom = F,order = F) # searched for
lungexp2<-plotDiffExpNolab(allSeurats$lungsEpcam,"NC-033676.1-46824199-46826449---143403-0",lungbam,"RPLP1",normalChrom = F,order = F) # searched for
lungexp3<-plotDiffExpNolab(allSeurats$lungsEpcam,"NC-033674.1-27742249-27745449---926795-0",lungbam,"RPS20",normalChrom = F,order = F) # searched for
lungexp4<-plotDiffExpNolab(allSeurats$lungsEpcam,"NC-033664.1-8257299-8265799---1410537-0",lungbam,"MS4A8",normalChrom = F,order = F) # searched for
lungexp5<-plotDiffExpNolab(allSeurats$lungsEpcam,"NC-033683.1-1184549-1186249---107026-0",lungbam,"BST2",normalChrom = F,order = F) # searched for
lungUmap<-plotUmapColorUmap(allSeurats$lungsEpcam,reducUse = "umap_gene")

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/lungUmaps2.pdf",
    width=6, height=3/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(lungUmap,lungexp1[[2]],lungexp2[[2]],lungexp3[[2]],lungexp4[[2]],lungexp5[[2]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/lungexCov2.pdf",
    width=5, height=2/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(lungexp1[[1]],lungexp2[[1]],lungexp3[[1]],lungexp4[[1]],lungexp5[[1]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/lungDot2.pdf",
    width=1.2, height=1+(2/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(lungDot+NoLegend(),nrow=1,align="hv")
dev.off()

##################################################
#find sea urchin example
############## differential marker analysis

cluster_labels<-c("endoderm",
                  "aboral ECD1",
                  "aboral ECD2",
                  "pigment",
                  "Veg2+",
                  "neuronal",
                  "LOC576614+",
                  "skeleton")


urchinFeatOut<-findSiguTARmarkers2(allSeurats$urchin,diffMarkersUtar$urchin,shift = 1,fcThresh = 2,minpc.1 = 0.5)
urchinFeatGene<-findSigGenemarkers2(allSeurats$urchin,diffMarkers$urchin,shift = 1,fcThresh = 0.5,minpc.1 = 0.25)

urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011997178.1-133249-136749")]<-"LOC764887"
urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011989268.1-19699-23999")]<-"LOC100893340"
urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011991856.1-522999-526249")]<-"LOC115926763"
urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011982291.1-29599-43249")]<-"LOC577764"
urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011969894.1-25849-31999")]<-"LOC105437727"
urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011980210.1-2449-4299")]<-"LOC115926624"
urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011992260.1-24399-26449")]<-"LOC105440873"
urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011993144.1-339299-341349")]<-"LOC115920005"
urchinFeatOut$blastnResults[startsWith(urchinFeatOut$gene,"NW-011980580.1-1275899-1283599")]<-"LOC105439981"

urchinDot<-plotDotPlotWithLabel(allSeurats$urchin,urchinFeatOut,dotScale = 1.5)

urchinexp1<-plotDiffExpNolab(allSeurats$urchin,"NW-011997178.1-133249-136749-+-47678-0",urchinbam,"LOC764887",normalChrom = F,order = F) # searched for
urchinexp2<-plotDiffExpNolab(allSeurats$urchin,"NW-011989268.1-19699-23999---36027-0",urchinbam,"LOC100893340",normalChrom = F,order = F) # searched for
urchinexp3<-plotDiffExpNolab(allSeurats$urchin,"NW-011991856.1-522999-526249---137297-0",urchinbam,"LOC115926763",normalChrom = F,order = F) # searched for
urchinexp4<-plotDiffExpNolab(allSeurats$urchin,"NW-011980210.1-2449-4299---83336-0",urchinbam,"LOC115926624",normalChrom = F,order = F) # searched for
urchinexp5<-plotDiffExpNolab(allSeurats$urchin,"NW-011993144.1-339299-341349---391388-0",urchinbam,"LOC115920005",normalChrom = F,order = F) # searched for
urchinUmap<-plotUmapColorUmap(allSeurats$urchin,reducUse = "umap_gene")

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/urchinUmaps2.pdf",
    width=6, height=3/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(urchinUmap,urchinexp1[[2]],urchinexp2[[2]],urchinexp3[[2]],urchinexp4[[2]],urchinexp5[[2]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/urchinexCov2.pdf",
    width=5, height=2/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(urchinexp1[[1]],urchinexp2[[1]],urchinexp3[[1]],urchinexp4[[1]],urchinexp5[[1]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/urchinDot2.pdf",
    width=1.2, height=1+(2/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(urchinDot+NoLegend(),nrow=1,align="hv")
dev.off()
##################################################
# mouse atlas
############### setting the statistics
load("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/mouseAtlasSpleen.RData")
load(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperScripts_wksp/figure2/mouseAtlasSpleenMarkers.RData")
mouseFeatOut<-cellTypeMarkersBed
mouseFeatOut<-mouseFeatOut[!duplicated(mouseFeatOut$TAR),]
rownames(mouseFeatOut)<-mouseFeatOut$TAR
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr1:21325481-21325744"]<-"GTH1" #new
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr17:21474944-21475613"]<-"unknown" #
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr1:88567488-88567557"]<-"PRPF8"
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr16:19267406-19267865"]<-"SNX29" #
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr16:19274835-19275494"]<-"ATF7IP2" # Igl complex
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr1:21320030-21320758"]<-"Phf20" #pseudogene
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chrX:94909038-94915832"]<-"ABITRAM" #pseudogene
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr16:19327491-19328415"]<-"PDXDC1" #high confidence
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr3:97992116-97994425"]<-"unknown" #pseudogene
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr2:150301377-150301556"]<-"unknown" #pseudogene
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr14:51173988-51174239"]<-"unknown" #pseudogene
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr14:69824433-69825050"]<-"unknown" #pseudogene
mouseFeatOut$blastnResults[mouseFeatOut$fastaPeak=="chr7:67813500-67813562"]<-"unknown" #pseudogene
mouseBam<-"/workdir/fw262/mouseAtlas/results_MouseAtlas_fromFastq/10X_P7_6/10X_P7_6_HMM_tagged_noDir.bam"
orderOfFeat<-order(mouseFeatOut$cluster) # order features
mouseFeatOut<-mouseFeatOut[orderOfFeat,]
mouseDot<-plotDotPlotWithLabel(mergedSeurat,mouseFeatOut,dotScale = 2.5)

mouseexp1<-plotDiffExpNolab(mergedSeurat,"chr1-21324999-21327499---792-0",mouseBam,"GTH1",normalChrom = T,order = F) # searched for
mouseexp2<-plotDiffExpNolab(mergedSeurat,"chr1-88567149-88567849-+-2781-0",mouseBam,"PRPF8",normalChrom = T,order = F) # searched for
mouseexp3<-plotDiffExpNolab(mergedSeurat,"chr16-19263949-19268649-+-559-0",mouseBam,"SNX29",normalChrom = T,order = F) # searched for
mouseexp4<-plotDiffExpNolab(mergedSeurat,"chr16-19270249-19276499-+-466-0",mouseBam,"ATF7IP2",normalChrom = T,order = F) # searched for
mouseexp5<-plotDiffExpNolab(mergedSeurat,"chr1-21318099-21324249---377-0",mouseBam,"Phf20",normalChrom = T,order = F) # searched for
mouseexp6<-plotDiffExpNolab(mergedSeurat,"chr16-19325399-19332999-+-326-0",mouseBam,"PDXDC1",normalChrom = T,order = F)
mouseexp7<-plotDiffExpNolab(mergedSeurat,"chrX-94858199-94918049-+-20098-0",mouseBam,"ABITRAM",normalChrom = T,order = F)

cluster_labels<-c("0",
                  "1",
                  "2",
                  "3",
                  "4")
names(cluster_labels) <- levels(mergedSeurat)
mergedSeurat<-RenameIdents(mergedSeurat,cluster_labels)

mouseUmap<-plotUmapColorUmap(mergedSeurat,reducUse = "umap_gene")

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/mouseUmaps3.pdf",
    width=6, height=3/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(mouseUmap,mouseexp1[[2]],mouseexp2[[2]],mouseexp3[[2]],mouseexp4[[2]],mouseexp5[[2]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/mouseexCov2.pdf",
    width=5, height=2/3, paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(mouseexp1[[1]],mouseexp2[[1]],mouseexp3[[1]],mouseexp4[[1]],mouseexp5[[1]],nrow=1,align="hv")
dev.off()

pdf(file="/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/paperFigs/mouseDot3.pdf",
    width=1.2, height=1+(2/3), paper="special", bg="white",
    fonts="Helvetica", pointsize=6, useDingbats = F)
plot_grid(mouseDot+NoLegend(),nrow=1,align="hv")
dev.off()
