readsDataFunc<-function(mouseMuscle_seurat){
  df<-data.frame(GetAssayData(mouseMuscle_seurat, slot = "counts"))
  outGene<-rownames(df)[endsWith(rownames(df),"-0")]
  inGene<-rownames(df)[endsWith(rownames(df),"-1")]
  genes<-rownames(df)[!endsWith(rownames(df),"-0")&!endsWith(rownames(df),"-1")]
  readsData<-data.frame(rbind(colSums(df[genes,]),colSums(df[inGene,]),colSums(df[outGene,])))
  rownames(readsData)<-c("genes","HMMinGene","HMMoutGene")
  readsData<-data.frame(t(readsData))
  readsData<-readsData[order(readsData$HMMoutGene),]
  readsData$cell<-rownames(readsData)
  readsData$outpercent = readsData$HMMoutGene / (readsData$HMMinGene + readsData$HMMoutGene)
  readsData$genesPercent=readsData$genes/(readsData$HMMinGene + readsData$HMMoutGene + readsData$genes)
  readsData$cell2 = factor(readsData$cell, levels = as.character(readsData[order(readsData$outpercent),]$cell))
  return(readsData)
}

# first input is the non seurat object
# second input is the seurat object
readsDataFuncQuick<-function(mouseBrain,mouseBrain_seurat){
  df<-mouseBrain[,colnames(mouseBrain_seurat)]
  rownames(df)<-gsub("_","-",rownames(df))
  outGene<-rownames(mouseBrain_seurat)[endsWith(rownames(mouseBrain_seurat),"-0")]
  inGene<-rownames(mouseBrain_seurat)[endsWith(rownames(mouseBrain_seurat),"-1")]
  genes<-rownames(mouseBrain_seurat)[!endsWith(rownames(mouseBrain_seurat),"-0")&!endsWith(rownames(mouseBrain_seurat),"-1")]
  readsData<-data.frame(rbind(colSums(df[genes,]),colSums(df[inGene,]),colSums(df[outGene,])))
  rownames(readsData)<-c("genes","HMMinGene","HMMoutGene")
  readsData<-data.frame(t(readsData))
  readsData<-readsData[order(readsData$HMMoutGene),]
  readsData$cell<-rownames(readsData)
  readsData$outpercent = readsData$HMMoutGene / (readsData$HMMinGene + readsData$HMMoutGene)
  readsData$genesPercent=readsData$genes/(readsData$HMMinGene + readsData$HMMoutGene + readsData$genes)
  readsData$cell2 = factor(readsData$cell, levels = as.character(readsData[order(readsData$outpercent),]$cell))
  return(readsData)
}

convertHMMDirHorse<-function(Sample_HMM,refFlat){
  rownames(Sample_HMM)<-paste0("chr",rownames(Sample_HMM))
  rownames(Sample_HMM)<-substr(rownames(Sample_HMM),1,nchar(rownames(Sample_HMM))-2)
  
  x<-vector(mode="character", length=nrow(Sample_HMM))
  for(i in 1:nrow(Sample_HMM)){
    x[i]<-as.character(refFlat$V1[rownames(Sample_HMM)[i]==refFlat$V2])
  }
  rownames(Sample_HMM)<-x
  return(Sample_HMM)
  
}

convertHMMDir<-function(Sample_HMM,refFlat){
  #rownames(Sample_HMM)<-paste0("chr",rownames(Sample_HMM))
  rownames(Sample_HMM)<-substr(rownames(Sample_HMM),1,nchar(rownames(Sample_HMM))-2)
  
  x<-vector(mode="character", length=nrow(Sample_HMM))
  for(i in 1:nrow(Sample_HMM)){
    x[i]<-as.character(refFlat$V1[rownames(Sample_HMM)[i]==refFlat$V2])
  }
  rownames(Sample_HMM)<-x
  return(Sample_HMM)
  
}

convertHMMDirHorse<-function(Sample_HMM,refFlat){
  rownames(Sample_HMM)<-paste0("chr",rownames(Sample_HMM))
  rownames(Sample_HMM)<-substr(rownames(Sample_HMM),1,nchar(rownames(Sample_HMM))-2)
  
  x<-vector(mode="character", length=nrow(Sample_HMM))
  for(i in 1:nrow(Sample_HMM)){
    x[i]<-as.character(refFlat$V1[rownames(Sample_HMM)[i]==refFlat$V2])
  }
  rownames(Sample_HMM)<-x
  return(Sample_HMM)
  
}

loadExpMat<-function(file,projName){
  Sample_HMM=fread(paste0(file),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sueratOb<-CreateSeuratObject(counts = Sample_HMM, project = projName,min.cells = 2, min.features = 100)
  #x1<-minNFeature_RNA
  #y1<-maxNFeature_RNA
  #out<-subset(sueratOb, subset = nFeature_RNA > x1 & nFeature_RNA < y1)
  return(sueratOb)
}

loadExpMatWithHMM<-function(HMMfile,genefile,projName){
  Sample_HMM=fread(HMMfile,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  Sample_genes=fread(genefile,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_genes <- as.data.frame(Sample_genes) # force to data frame
  rownames(Sample_genes) <- Sample_genes$GENE # make the name of rows GENEs
  Sample_genes <- Sample_genes[,-1] # take out first column
  #mito.genes<-c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")
  #rownames(Sample_genes)[(rownames(Sample_genes) %in% mito.genes)]<-paste0("MT-",rownames(Sample_genes)[rownames(Sample_genes) %in% mito.genes]) #incorporate MT genes
  
  
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  cells<-intersect(colnames(Sample_HMM),colnames(Sample_genes))

  combinedData<-rbind(Sample_genes[,cells],Sample_HMM[,cells])
  seuratOb<-CreateSeuratObject(counts = combinedData, project = projName,min.cells = 2, min.features = 100)
  return(seuratOb)
}

normVarScaleData<-function(seuratOb,numFeat=2000){
  seuratOb <- NormalizeData(seuratOb)
  seuratOb <- FindVariableFeatures(seuratOb, selection.method = "vst", nfeatures=numFeat)
  seuratOb <- ScaleData(seuratOb, features = rownames(seuratOb))
  return(seuratOb)
}

normVarScaleDataNoVarGene<-function(seuratOb){
  seuratOb <- NormalizeData(seuratOb)
  seuratOb <- ScaleData(seuratOb, features = rownames(seuratOb))
  return(seuratOb)
}


createOutGeneDf<-function(OutMarkSig){
  OutMarkSig$gene<-as.character(rownames(OutMarkSig))
  OutMarkSig$chr<-as.character(sapply(strsplit(OutMarkSig$gene,"-"),"[",1))
  OutMarkSig$start<-as.numeric(sapply(strsplit(OutMarkSig$gene,"-"),"[",2))
  OutMarkSig$stop<-as.numeric(sapply(strsplit(OutMarkSig$gene,"-"),"[",3))
  OutMarkSig$fasta<-paste0(OutMarkSig$chr,":",OutMarkSig$start,"-",OutMarkSig$stop)
  OutMarkSig$fasta<-substr(OutMarkSig$gene,1,nchar(OutMarkSig$gene)-5)
  i <- 2
  j <- 2
  OutMarkSig$fasta<-sapply(strsplit(OutMarkSig$fasta, "-"), function(x) {
    g <- seq_along(x)
    g[g < i] <- i
    g[g > j + 1] <- j+1
    paste(tapply(x, g, paste, collapse = "-"), collapse = ":")
  })
  
  i <- 1
  j <- 1
  OutMarkSig$fasta<-sapply(strsplit(OutMarkSig$fasta, "-"), function(x) {
    g <- seq_along(x)
    g[g < i] <- i
    g[g > j + 1] <- j+1
    paste(tapply(x, g, paste, collapse = "-"), collapse = "_")
  })
  
  
  OutMarkSig$featureLength<-OutMarkSig$stop-OutMarkSig$start+1
  OutMarkSig$totalReads<-as.numeric(unlist(lapply(strsplit(OutMarkSig$gene,"-"),function(x) x[length(x)-1])))
  OutMarkSig$readsPerLen<-as.numeric(OutMarkSig$totalReads)/OutMarkSig$featureLength
  OutMarkSig$blastnResults<-NA
  
  return(OutMarkSig)
}

createOutGeneDf2<-function(OutMarkSig){
  OutMarkSig$TAR<-OutMarkSig$gene
  OutMarkSig$chr<-paste(as.character(sapply(strsplit(OutMarkSig$TAR,"-"),"[",1)),
                        as.character(sapply(strsplit(OutMarkSig$TAR,"-"),"[",2)),
                        sep="_")
  OutMarkSig$start<-as.numeric(sapply(strsplit(OutMarkSig$TAR,"-"),"[",3))
  OutMarkSig$end<-as.numeric(sapply(strsplit(OutMarkSig$TAR,"-"),"[",4))
  OutMarkSig$fasta<-paste0(OutMarkSig$chr,":",OutMarkSig$start,"-",OutMarkSig$end)
  #OutMarkSig$fasta<-substr(OutMarkSig$TAR,1,nchar(OutMarkSig$TAR)-5)
  
  OutMarkSig$strand<-as.character(sapply(strsplit(OutMarkSig$TAR,"-"),"[",5))
  OutMarkSig$strand[OutMarkSig$strand==""]<-"-"
  
  i <- 2
  j <- 2
  OutMarkSig$fasta<-sapply(strsplit(OutMarkSig$fasta, "-"), function(x) {
    g <- seq_along(x)
    g[g < i] <- i
    g[g > j + 1] <- j+1
    paste(tapply(x, g, paste, collapse = "-"), collapse = ":")
  })
  
  i <- 1
  j <- 1
  OutMarkSig$fasta<-sapply(strsplit(OutMarkSig$fasta, "-"), function(x) {
    g <- seq_along(x)
    g[g < i] <- i
    g[g > j + 1] <- j+1
    paste(tapply(x, g, paste, collapse = "-"), collapse = "_")
  })
  
  
  OutMarkSig$featureLength<-OutMarkSig$end-OutMarkSig$start+1
  OutMarkSig$totalReads<-as.numeric(unlist(lapply(strsplit(OutMarkSig$TAR,"-"),function(x) x[length(x)-1])))
  OutMarkSig$readsPerLen<-as.numeric(OutMarkSig$totalReads)/OutMarkSig$featureLength
  #OutMarkSig$blastnResults<-NA
  
  return(OutMarkSig)
}

createOutGeneDf4<-function(OutMarkSig){
  OutMarkSig$TAR<-OutMarkSig$gene
  OutMarkSig$chr<-paste(as.character(sapply(strsplit(OutMarkSig$TAR,"-"),"[",1)),
                        as.character(sapply(strsplit(OutMarkSig$TAR,"-"),"[",2)),
                        sep="_")
  OutMarkSig$start<-as.numeric(sapply(strsplit(OutMarkSig$TAR,"-"),"[",3))
  OutMarkSig$end<-as.numeric(sapply(strsplit(OutMarkSig$TAR,"-"),"[",4))
  OutMarkSig$fasta<-paste0(OutMarkSig$chr,":",OutMarkSig$start,"-",OutMarkSig$end)
  #OutMarkSig$fasta<-substr(OutMarkSig$TAR,1,nchar(OutMarkSig$TAR)-5)
  
  OutMarkSig$strand<-as.character(sapply(strsplit(OutMarkSig$TAR,"-"),"[",5))
  OutMarkSig$strand[OutMarkSig$strand==""]<-"-"
  
  OutMarkSig$featureLength<-OutMarkSig$end-OutMarkSig$start+1
  OutMarkSig$totalReads<-as.numeric(unlist(lapply(strsplit(OutMarkSig$TAR,"-"),function(x) x[length(x)-1])))
  OutMarkSig$readsPerLen<-as.numeric(OutMarkSig$totalReads)/OutMarkSig$featureLength
  #OutMarkSig$blastnResults<-NA
  
  return(OutMarkSig)
}


plotListOfCov<-function(outMarkerDf){
  out<-list()
  for (i in 1:nrow(outMarkerDf)){
    
    eval(substitute(df<-data.frame(as.numeric(outMarkerDf$coverage[[i]]))    ,list(i=i)))

    eval(substitute(p<-ggplot(df,aes(x=1:nrow(df),y=df[,1]),width=1)+geom_area(stat="identity",fill="blue")    ,list(i=i)))

    eval(substitute(print(i)    ,list(i=i)))

    eval(substitute( out[[i]]<-p    ,list(i=i)))

  }
  plotOut<-plot_grid(out,ncol=1,align="hv")
  return(plotOut)
}

plotListOfCov2<-function(input,colorInput="black"){
  df<-data.frame(as.numeric(input))
  out<-ggplot(df,aes(x=1:nrow(df),y=df[,1]),wdith=1)+
    geom_area(stat="identity",fill=colorInput)+ theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8))+
    scale_y_continuous(limits=c(0,max(df[,1])),
                       breaks = c(0,max(df[])))+
    scale_x_continuous(limits=c(0,nrow(df)),
                       breaks = c(0,nrow(df)))
  
  return(out)
}

plotUmapColorUmap<-function(input,reducUse="umap"){
  umapGenes<-DimPlot(input,
                     reduction=reducUse,
                     label=T,
                     pt.size=0.01,
                     label.size = 0,
                     cols=as.vector(alphabet(n=length(unique(Idents(input))))))+
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank(),
          #legend.text=element_text(size=12),
          #legend.position="bottom",
          aspect.ratio=1)+
    NoLegend()+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  return(umapGenes)
}


plotListOfCovNoLabel<-function(input,colorInput="dark red"){
  df<-data.frame(as.numeric(input))
  out<-ggplot(df,aes(x=1:nrow(df),y=df[,1]),wdith=1)+
    geom_area(stat="identity",fill=colorInput)+ theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=15,face = "bold"),
          axis.text.y=element_text(size=15,face = "bold"))+
    scale_y_continuous(limits=c(0,max(df[,1])),
                       breaks = c(0,max(df[])))+
    scale_x_continuous(limits=c(0,nrow(df)),
                       breaks = c(0,nrow(df)))
  
  return(out)
}



plotListOfCov1.5<-function(input,colorInput="dark red"){
  df<-data.frame(as.numeric(input))
  out<-ggplot(df,aes(x=1:nrow(df),y=df[,1]),wdith=1)+
    geom_area(stat="identity",fill=colorInput)+ theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=15,face = "bold"),
          axis.text.y=element_text(size=15,face = "bold"))

  return(out)
}



generateClusters<-function(chicken_seurat){
  all_Seurat_chicken<-chicken_seurat
  varHMMoutGene<-VariableFeatures(object=all_Seurat_chicken)[endsWith(VariableFeatures(object=all_Seurat_chicken),"-0")]
  varHMMoutGene<-varHMMoutGene[!grepl("random",varHMMoutGene)]
  varHMMinGene<-VariableFeatures(object=all_Seurat_chicken)[endsWith(VariableFeatures(object=all_Seurat_chicken),"-1")]
  varHMMinGene<-varHMMinGene[!grepl("random",varHMMinGene)]
  varGenes<-VariableFeatures(object=all_Seurat_chicken)[!(endsWith(VariableFeatures(object=all_Seurat_chicken),"-0"))&!(endsWith(VariableFeatures(object=all_Seurat_chicken),"-1"))]
  all_Seurat_chicken_Genes<-RunPCA(object = all_Seurat_chicken, features = varGenes)
  all_Seurat_chicken_Genes<-FindNeighbors(object=all_Seurat_chicken_Genes,dims=1:10)
  all_Seurat_chicken_Genes<-FindClusters(object=all_Seurat_chicken_Genes,resolution=0.5)
  all_Seurat_chicken_inGene<-RunPCA(object = all_Seurat_chicken, features = varHMMinGene)
  Idents(object=all_Seurat_chicken_inGene)<-Idents(object=all_Seurat_chicken_Genes) #assign same cluster identity as all_Seurat_chicken_Genes
  all_Seurat_chicken_outGene<-RunPCA(object = all_Seurat_chicken, features = varHMMoutGene)
  Idents(object=all_Seurat_chicken_outGene)<-Idents(object=all_Seurat_chicken_Genes) #assign same cluster identity as all_Seurat_chicken_Genes
  all_Seurat_chicken_HMM<-RunPCA(object = all_Seurat_chicken, features = c(varHMMinGene,varHMMoutGene))
  Idents(object=all_Seurat_chicken_HMM)<-Idents(object=all_Seurat_chicken_Genes) #assign same cluster identity as all_Seurat_chicken_Genes
  
  all_Seurat_chicken_Genes<-RunUMAP(all_Seurat_chicken_Genes, dims=1:10,do.fast=T,check_duplicates=F)
  all_Seurat_chicken_inGene<-RunUMAP(all_Seurat_chicken_inGene, dims=1:10,do.fast=T,check_duplicates=F)
  all_Seurat_chicken_outGene<-RunUMAP(all_Seurat_chicken_outGene, dims=1:10,do.fast=T,check_duplicates=F)
  all_Seurat_chicken_HMM<-RunUMAP(all_Seurat_chicken_HMM, dims=1:10,do.fast=T,check_duplicates=F)
  
  chickenCluster_HMMout<-RunPCA(object = all_Seurat_chicken, features = varHMMoutGene)
  chickenCluster_HMMout<-FindNeighbors(object=chickenCluster_HMMout,dims=1:10)
  chickenCluster_HMMout<-FindClusters(object=chickenCluster_HMMout,resolution=0.5)
  chickenCluster_HMMout<-RunUMAP(chickenCluster_HMMout, dims=1:10,do.fast=T,check_duplicates=F)
  
  output<-c(all_Seurat_chicken_Genes,all_Seurat_chicken_HMM,all_Seurat_chicken_outGene,chickenCluster_HMMout)
  return(output)
  
}

doDimPlots<-function(cluster){
  p1<-DimPlot(cluster[[1]], pt.size=0.3,reduction = "umap",no.axes=T)+labs(title="gene exp")+theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())#,legend.position="none")
  p2<-DimPlot(cluster[[2]], pt.size=0.3,reduction = "umap",no.axes=T)+labs(title="HMM features")+NoLegend()+theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none")
  p3<-DimPlot(cluster[[3]],pt.size=0.3, reduction = "umap",no.axes=T)+labs(title="HMM out")+NoLegend()+theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none")
  p4<-DimPlot(cluster[[4]],pt.size=0.3, reduction = "umap",no.axes=T)+labs(title="HMM out cluster")+NoLegend()+theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none")
  mouseMuscle_cluster<-plot_grid(p1,p2,p3,p4,ncol=2, nrow=2)
  
  return(mouseMuscle_cluster)
}


######## run cken2 first, then use those inputs to generate equivalent output of combined without Scanorama
combineNoScanorama<-function(LVdata,RVdata,dayNumber){
  LVdatadf<-as.data.frame(t(data.frame(GetAssayData(object = LVdata[[1]],slot="counts"))))
  RVdatadf<-as.data.frame(t(data.frame(GetAssayData(object = RVdata[[1]],slot="counts"))))
  rownames(LVdatadf)<-paste0(rownames(LVdatadf),paste0("_",dayNumber,"LV"))
  rownames(RVdatadf)<-paste0(rownames(RVdatadf),paste0("_",dayNumber,"RV"))
  Matt<-rbind.fill(LVdatadf,RVdatadf)
  Matt[is.na(Matt)]<-0
  rownames(Matt)<-c(rownames(LVdatadf),rownames(RVdatadf))
  output<-CreateSeuratObject(counts=t(Matt))
  output<-NormalizeData(output)
  outputL<-generateDatafromScanorama(output)
  outputL[[1]][["sample"]]<-sapply(strsplit(colnames(outputL[[1]]),"_"),"[",2)
  outputL[[2]][["sample"]]<-sapply(strsplit(colnames(outputL[[2]]),"_"),"[",2)
  outputL[[3]][["sample"]]<-sapply(strsplit(colnames(outputL[[3]]),"_"),"[",2)
  outputL[[4]]<-getReadsInfo(outputL[[1]])
  return(outputL)
  
}


cken<-function(geneMat,noDirMat,withDirMat){
  
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  d4Gene<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  d4Gene[["percent.mt"]] <- PercentageFeatureSet(object = d4Gene, pattern = "^MT-")
  VlnPlot(d4Gene,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
  d4Gene <- subset(x = d4Gene, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(d4Gene)
  
  d4noDir<-loadExpMat(noDirMat,"d4noDir")
  d4withDir<-loadExpMat(withDirMat,"d4withDir")
  d4noDir <- subset(x = d4noDir, cells = d4Cells)
  d4withDir <- subset(x = d4withDir, cells=d4Cells)
  
  d4Gene<-normVarScaleData(d4Gene)
  d4Gene <- RunPCA(object = d4Gene, features = VariableFeatures(d4Gene))
  d4Gene <- FindNeighbors(object=d4Gene,dims=1:10)
  d4Gene <- FindClusters(object=d4Gene,resolution=0.2)
  
  d4noDir<-normVarScaleData(d4noDir)
  d4noDiroutGene<-rownames(d4noDir)[endsWith(rownames(d4noDir),"-0")]
  d4noDirinGene<-rownames(d4noDir)[endsWith(rownames(d4noDir),"-1")]
  d4noDirIn <- RunPCA(object = d4noDir, features = d4noDirinGene)
  d4noDirOut <- RunPCA(object = d4noDir, features = d4noDiroutGene)
  Idents(d4noDirIn)<-Idents(d4Gene)
  Idents(d4noDirOut)<-Idents(d4Gene)
  
  d4withDir<-normVarScaleData(d4withDir)
  d4withDiroutGene<-rownames(d4withDir)[endsWith(rownames(d4withDir),"-0")]
  d4withDirinGene<-rownames(d4withDir)[endsWith(rownames(d4withDir),"-1")]
  d4withDirIn <- RunPCA(object = d4withDir, features = d4withDirinGene)
  d4withDirOut <- RunPCA(object = d4withDir, features = d4withDiroutGene)
  Idents(d4withDirIn)<-Idents(d4Gene)
  Idents(d4withDirOut)<-Idents(d4Gene)
  
  d4Gene <- RunUMAP(d4Gene,dims=1:10,check_duplicates = F)
  d4noDirIn <- RunUMAP(d4noDirIn,dims=1:10,check_duplicates = F)
  d4noDirOut <- RunUMAP(d4noDirOut,dims=1:10,check_duplicates = F)
  d4withDirIn <- RunUMAP(d4withDirIn,dims=1:10,check_duplicates = F)
  d4withDirOut <- RunUMAP(d4withDirOut,dims=1:10,check_duplicates = F)
  
  d4Gene <- RunTSNE(d4Gene,dims=1:10,check_duplicates = F)
  d4noDirIn <- RunTSNE(d4noDirIn,dims=1:10,check_duplicates = F)
  d4noDirOut <- RunTSNE(d4noDirOut,dims=1:10,check_duplicates = F)
  d4withDirIn <- RunTSNE(d4withDirIn,dims=1:10,check_duplicates = F)
  d4withDirOut <- RunTSNE(d4withDirOut,dims=1:10,check_duplicates = F)
  
  return(c(d4Gene,d4noDirIn,d4noDirOut,d4withDirIn,d4withDirOut))
}



generateGeneWithHMMout<-function(geneMat,hmmMat,importantOut){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  geneExpMat<-GetAssayData(geneOnly,slot="counts")
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  rownames(Sample_HMM)<-gsub("_","-",rownames(Sample_HMM))
  Sample_HMM<-Sample_HMM[importantOut,d4Cells]
  
  out<-rbind(data.frame(geneExpMat),Sample_HMM)
}
  

# generates data that has values of gene expression + hmm features
# useful for featurePlot of gene expression of hmm features
generateDataforChicken2<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  ####//////// setting the statistics
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  #geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  TARSeu<-sample_combined[c(inGene,outGene),]
  aTARSeu<-sample_combined[c(inGene),]
  uTARSeu<-sample_combined[c(outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(TARSeu,features = outGene)
  sample_combined[["nFeature_RNA"]]<-geneOnly[["nFeature_RNA"]]
  sample_combined[["nCount_RNA"]]<-geneOnly[["nCount_RNA"]]
  sample_combined[["nFeature_aTAR"]]<-aTARSeu[["nFeature_RNA"]]
  sample_combined[["nFeature_uTAR"]]<-uTARSeu[["nFeature_RNA"]]
  sample_combined[["nCount_aTAR"]]<-aTARSeu[["nCount_RNA"]]
  sample_combined[['nCount_uTAR']]<-uTARSeu[["nCount_RNA"]]
  ############################################
  
  # cluster on gene expression only
  sample_combined<-normVarScaleData(sample_combined,1)
  geneOnly_seurat <- RunPCA(object = sample_combined, features = rownames(geneOnly))
  geneOnly_seurat <- FindNeighbors(object=geneOnly_seurat,dims=1:10)
  geneOnly_seurat <- FindClusters(object=geneOnly_seurat,resolution=0.2)
  
  # run pca on in and out hmm features
  #outGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_0")]
  #inGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_1")]
  ########################## The above used to be uncommented
  
  outGene_seurat <- RunPCA(object = sample_combined, features = gsub("_","-",outGene))
  inGene_seurat <- RunPCA(object = sample_combined, features = gsub("_","-",inGene))
  Idents(outGene_seurat)<-Idents(geneOnly_seurat)
  Idents(inGene_seurat)<-Idents(geneOnly_seurat)
  
  geneOnly_seurat <- RunUMAP(geneOnly_seurat,dims=1:10,check_duplicates = F)
  inGene_seurat <- RunUMAP(inGene_seurat,dims=1:10,check_duplicates = F)
  outGene_seurat <- RunUMAP(outGene_seurat,dims=1:10,check_duplicates = F)

  #geneOnly_seurat <- RunTSNE(geneOnly_seurat,dims=1:10,check_duplicates = F)
  #inGene_seurat <- RunTSNE(inGene_seurat,dims=1:10,check_duplicates = F)
  #outGene_seurat <- RunTSNE(outGene_seurat,dims=1:10,check_duplicates = F)

  df<-sample_combined_mat[,colnames(sample_combined)]
  #rownames(df)<-gsub("_","-",rownames(df))
  genes<-rownames(geneOnly)
  readsData<-data.frame(rbind(colSums(df[genes,]),colSums(df[inGene,]),colSums(df[outGene,])))
  rownames(readsData)<-c("genes","HMMinGene","HMMoutGene")
  readsData<-data.frame(t(readsData))
  readsData<-readsData[order(readsData$HMMoutGene),]
  readsData$cell<-rownames(readsData)
  readsData$outpercent = readsData$HMMoutGene / (readsData$HMMinGene + readsData$HMMoutGene)
  readsData$genesPercent=readsData$genes/(readsData$HMMinGene + readsData$HMMoutGene + readsData$genes)
  readsData$cell2 = factor(readsData$cell, levels = as.character(readsData[order(readsData$outpercent),]$cell))
  
  return(c(geneOnly_seurat,inGene_seurat,outGene_seurat,list(readsData)))
}

##### get more read stats function from a created seurat object
getMoreStats<-function(input,numDashes){
  ####//////// setting the statistics
  numdash<-str_count(rownames(input),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(input),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(input),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(input)[geneInd]
  inGene<-rownames(input)[inInd]
  outGene<-rownames(input)[outInd]
  
  TARSeu<-input[c(inGene,outGene),]
  aTARSeu<-input[c(inGene),]
  uTARSeu<-input[c(outGene),]
  geneSeu<-input[geneOnly,]
  input[["percent.outHMM"]]<-PercentageFeatureSet(TARSeu,features = outGene)
  input[["nFeature_RNA"]]<-geneSeu[["nFeature_RNA"]]
  input[["nCount_RNA"]]<-geneSeu[["nCount_RNA"]]
  input[["nFeature_aTAR"]]<-aTARSeu[["nFeature_RNA"]]
  input[["nFeature_uTAR"]]<-uTARSeu[["nFeature_RNA"]]
  input[["nCount_aTAR"]]<-aTARSeu[["nCount_RNA"]]
  input[['nCount_uTAR']]<-uTARSeu[["nCount_RNA"]]
  ############################################
  return(input)
  
}

getpercentOut<-function(input,numDashes=5){
  ####//////// setting the statistics
  numdash<-str_count(rownames(input),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(input),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(input),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(input)[geneInd]
  inGene<-rownames(input)[inInd]
  outGene<-rownames(input)[outInd]
  
  TARSeu<-input[c(inGene,outGene),]
  aTARSeu<-input[c(inGene),]
  uTARSeu<-input[c(outGene),]
  geneSeu<-input[geneOnly,]
  input[["percent.outHMM"]]<-PercentageFeatureSet(TARSeu,features = outGene)

  return(input)
}

# generates data that has values of gene expression + hmm features
# useful for featurePlot of gene expression of hmm features
generateDataforLemur<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  #geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  # cluster on gene expression only
  sample_combined<-normVarScaleData(sample_combined)
  geneOnly_seurat <- RunPCA(object = sample_combined, features = rownames(geneOnly))
  geneOnly_seurat <- FindNeighbors(object=geneOnly_seurat,dims=1:10)
  geneOnly_seurat <- FindClusters(object=geneOnly_seurat,resolution=0.2)
  
  # run pca on in and out hmm features
  #outGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_0")]
  #inGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_1")]
  
  numdash<-str_count(rownames(Sample_HMM),"_")
  outInd<-numdash>=5 & endsWith(rownames(Sample_HMM),"_0")
  inInd<-numdash>=5 & endsWith(rownames(Sample_HMM),"_1")
  inGene<-rownames(Sample_HMM)[inInd]
  outGene<-rownames(Sample_HMM)[outInd]
  
  
  
  outGene_seurat <- RunPCA(object = sample_combined, features = gsub("_","-",outGene))
  #outGene_seurat <- FindNeighbors(object=outGene_seurat,dims=1:10)
  #outGene_seurat <- FindClusters(object=outGene_seurat,resolution=0.2)
  
  inGene_seurat <- RunPCA(object = sample_combined, features = gsub("_","-",inGene))
  #inGene_seurat <- FindNeighbors(object=inGene_seurat,dims=1:10)
  #inGene_seurat <- FindClusters(object=inGene_seurat,resolution=0.2)
  
  Idents(outGene_seurat)<-Idents(geneOnly_seurat)
  Idents(inGene_seurat)<-Idents(geneOnly_seurat)
  
  geneOnly_seurat <- RunUMAP(geneOnly_seurat,dims=1:10,check_duplicates = F)
  inGene_seurat <- RunUMAP(inGene_seurat,dims=1:10,check_duplicates = F)
  outGene_seurat <- RunUMAP(outGene_seurat,dims=1:10,check_duplicates = F)
  
  #geneOnly_seurat <- RunTSNE(geneOnly_seurat,dims=1:10,check_duplicates = F)
  #inGene_seurat <- RunTSNE(inGene_seurat,dims=1:10,check_duplicates = F)
  #outGene_seurat <- RunTSNE(outGene_seurat,dims=1:10,check_duplicates = F)
  
  df<-sample_combined_mat[,colnames(sample_combined)]
  #rownames(df)<-gsub("_","-",rownames(df))
  genes<-rownames(geneOnly)
  readsData<-data.frame(rbind(colSums(df[genes,]),colSums(df[inGene,]),colSums(df[outGene,])))
  rownames(readsData)<-c("genes","HMMinGene","HMMoutGene")
  readsData<-data.frame(t(readsData))
  readsData<-readsData[order(readsData$HMMoutGene),]
  readsData$cell<-rownames(readsData)
  readsData$outpercent = readsData$HMMoutGene / (readsData$HMMinGene + readsData$HMMoutGene)
  readsData$genesPercent=readsData$genes/(readsData$HMMinGene + readsData$HMMoutGene + readsData$genes)
  readsData$cell2 = factor(readsData$cell, levels = as.character(readsData[order(readsData$outpercent),]$cell))
  
  return(c(geneOnly_seurat,inGene_seurat,outGene_seurat,list(readsData)))
}

generateDataforLemur_objectToMerge<-function(geneMat,hmmMat,project_name){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 200 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = project_name, min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  
  return(sample_combined)
}

loadGene_TARs<-function(geneMat,hmmMat,project_name){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = project_name, min.cells = 2)

  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  
  return(list(sample_combined,sample_combined_mat))
}

loadGene_TARs_mouseBrain<-function(geneMat,hmmMat,project_name){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined_mat<-sample_combined_mat[,d4Cells]
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  # cluster on gene expression only
  sample_combined<-normVarScaleData(sample_combined,1)
  geneOnly_seurat <- RunPCA(object = sample_combined, features = geneOnly)
  geneOnly_seurat <- FindNeighbors(object=geneOnly_seurat,dims=1:10)
  geneOnly_seurat <- FindClusters(object=geneOnly_seurat,resolution=0.2)
  
  # run pca on in and out hmm features
  #outGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_0")]
  #inGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_1")]
  ########################## The above used to be uncommented
  
  outGene_seurat <- RunPCA(object = sample_combined, features = outGene)
  inGene_seurat <- RunPCA(object = sample_combined, features = inGene)
  Idents(outGene_seurat)<-Idents(geneOnly_seurat)
  Idents(inGene_seurat)<-Idents(geneOnly_seurat)
  
  geneOnly_seurat <- RunUMAP(geneOnly_seurat,dims=1:10,check_duplicates = F)
  inGene_seurat <- RunUMAP(inGene_seurat,dims=1:10,check_duplicates = F)
  outGene_seurat <- RunUMAP(outGene_seurat,dims=1:10,check_duplicates = F)
  
  return(c(geneOnly_seurat,inGene_seurat,outGene_seurat,sample_combined_mat))
}



generateDataforLemur_objectToMerge2<-function(geneMat1,geneMat2,hmmMat1,hmmMat2,project_name){
  
  # load gene expression matrix and find valid cells
  d4Gene1=fread(paste0(geneMat1),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene1 <- as.data.frame(d4Gene1) # force to data frame
  rownames(d4Gene1) <- d4Gene1$GENE # make the name of rows GENEs
  d4Gene1 <- d4Gene1[,-1] # take out first column
  d4Gene2=fread(paste0(geneMat2),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene2 <- as.data.frame(d4Gene2) # force to data frame
  rownames(d4Gene2) <- d4Gene2$GENE # make the name of rows GENEs
  d4Gene2 <- d4Gene2[,-1] # take out first column
  
  test<-merge(d4Gene1, d4Gene2, by = "row.names", all = TRUE)
  d4Gene<-as.matrix(test[-1])
  rownames(d4Gene) <- test[,1]
  d4Gene[is.na(d4Gene)] <- 0
  colnames(d4Gene)<-substr(colnames(d4Gene),1,nchar(colnames(d4Gene))-2)
  
  
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 200 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = project_name, min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  
  return(sample_combined)
}


generateDataforSeaUrchin<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 400 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  
  return(sample_combined)
}


generateDataforLemur_objectToMerge_lightFilter<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  #geneOnly <- subset(x = geneOnly, nFeature_RNA > 200 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  
  return(sample_combined)
}


# input is a seurat object
getReadsInfo<-function(input){
  df<-data.frame(GetAssayData(object = input, slot = "counts"))
  
  numdash<-str_count(rownames(df),"-")
  outInd<-numdash>=5 & endsWith(rownames(df),"-0")
  inInd<-numdash>=5 & endsWith(rownames(df),"-1")
  geneInd<-!(outInd|inInd)
  
  geneOnly<-rownames(df)[geneInd]
  inGene<-rownames(df)[inInd]
  outGene<-rownames(df)[outInd]
  
  
  readsData<-data.frame(rbind(colSums(df[geneOnly,]),colSums(df[inGene,]),colSums(df[outGene,])))
  rownames(readsData)<-c("genes","HMMinGene","HMMoutGene")
  readsData<-data.frame(t(readsData))
  readsData<-readsData[order(readsData$HMMoutGene),]
  readsData$cell<-rownames(readsData)
  readsData$outpercent = readsData$HMMoutGene / (readsData$HMMinGene + readsData$HMMoutGene)
  readsData$genesPercent=readsData$genes/(readsData$HMMinGene + readsData$HMMoutGene + readsData$genes)
  readsData$cell2 = factor(readsData$cell, levels = as.character(readsData[order(readsData$outpercent),]$cell))
  
  return(readsData)
}

# input is a seurat object
getReadsInfoHuman<-function(input){
  df<-data.frame(GetAssayData(object = input, slot = "counts"))
  
  numdash<-str_count(rownames(df),"-")
  outInd<-numdash>=3 & endsWith(rownames(df),"-0")
  inInd<-numdash>=3 & endsWith(rownames(df),"-1")
  geneInd<-!(outInd|inInd)
  
  geneOnly<-rownames(df)[geneInd]
  inGene<-rownames(df)[inInd]
  outGene<-rownames(df)[outInd]
  
  
  readsData<-data.frame(rbind(colSums(df[geneOnly,]),colSums(df[inGene,]),colSums(df[outGene,])))
  rownames(readsData)<-c("genes","HMMinGene","HMMoutGene")
  readsData<-data.frame(t(readsData))
  readsData<-readsData[order(readsData$HMMoutGene),]
  readsData$cell<-rownames(readsData)
  readsData$outpercent = readsData$HMMoutGene / (readsData$HMMinGene + readsData$HMMoutGene)
  readsData$genesPercent=readsData$genes/(readsData$HMMinGene + readsData$HMMoutGene + readsData$genes)
  readsData$cell2 = factor(readsData$cell, levels = as.character(readsData[order(readsData$outpercent),]$cell))
  
  return(readsData)
}




# gener all the output based on input of scanorama combined dataset
generateDatafromScanorama<-function(sample_combined){
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  # cluster on gene expression only
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  geneOnly_seurat <- RunPCA(object = sample_combined, features = geneOnly)
  geneOnly_seurat <- FindNeighbors(object=geneOnly_seurat,dims=1:10)
  geneOnly_seurat <- FindClusters(object=geneOnly_seurat,resolution=0.2)
  
  # run pca on in and out hmm features
  outGene_seurat <- RunPCA(object = sample_combined, features = outGene)
  outGene_seurat <- FindNeighbors(object=outGene_seurat,dims=1:10)
  outGene_seurat <- FindClusters(object=outGene_seurat,resolution=0.2)
  
  inGene_seurat <- RunPCA(object = sample_combined, features = inGene)
  inGene_seurat <- FindNeighbors(object=inGene_seurat,dims=1:10)
  inGene_seurat <- FindClusters(object=inGene_seurat,resolution=0.2)
  
  Idents(outGene_seurat)<-Idents(geneOnly_seurat)
  Idents(inGene_seurat)<-Idents(geneOnly_seurat)
  
  geneOnly_seurat <- RunUMAP(geneOnly_seurat,dims=1:10,check_duplicates = F)
  inGene_seurat <- RunUMAP(inGene_seurat,dims=1:10,check_duplicates = F)
  outGene_seurat <- RunUMAP(outGene_seurat,dims=1:10,check_duplicates = F)
  
  ##geneOnly_seurat <- RunTSNE(geneOnly_seurat,dims=1:10,check_duplicates = F)
  ##inGene_seurat <- RunTSNE(inGene_seurat,dims=1:10,check_duplicates = F)
  ##outGene_seurat <- RunTSNE(outGene_seurat,dims=1:10,check_duplicates = F)
  #df<-data.frame(GetAssayData(sample_combined@assays$RNA,slot="counts"))
  ##genes<-rownames(geneOnly)
  #readsData<-data.frame(rbind(colSums(df[geneOnly,]),colSums(df[inGene,]),colSums(df[outGene,])))
  #rownames(readsData)<-c("genes","HMMinGene","HMMoutGene")
  #readsData<-data.frame(t(readsData))
  #readsData<-readsData[order(readsData$HMMoutGene),]
  #readsData$cell<-rownames(readsData)
  #readsData$outpercent = readsData$HMMoutGene / (readsData$HMMinGene + readsData$HMMoutGene)
  #readsData$genesPercent=readsData$genes/(readsData$HMMinGene + readsData$HMMoutGene + readsData$genes)
  #readsData$cell2 = factor(readsData$cell, levels = as.character(readsData[order(readsData$outpercent),]$cell))
  
  return(c(geneOnly_seurat,inGene_seurat,outGene_seurat))
}


# generates data that has values of gene expression + hmm features
# useful for featurePlot of gene expression of hmm features
generateDataforHumans<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  #geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  # cluster on gene expression only
  sample_combined<-normVarScaleData(sample_combined)
  geneOnly_seurat <- RunPCA(object = sample_combined, features = rownames(geneOnly))
  geneOnly_seurat <- FindNeighbors(object=geneOnly_seurat,dims=1:10)
  geneOnly_seurat <- FindClusters(object=geneOnly_seurat,resolution=0.2)
  
  # run pca on in and out hmm features
  outGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_0")]
  inGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_1")]
  
  outGene_seurat <- RunPCA(object = sample_combined, features = gsub("_","-",outGene))
  outGene_seurat <- FindNeighbors(object=outGene_seurat,dims=1:10)
  outGene_seurat <- FindClusters(object=outGene_seurat,resolution=0.2)
  
  inGene_seurat <- RunPCA(object = sample_combined, features = gsub("_","-",inGene))
  inGene_seurat <- FindNeighbors(object=inGene_seurat,dims=1:10)
  inGene_seurat <- FindClusters(object=inGene_seurat,resolution=0.2)
  
  Idents(outGene_seurat)<-Idents(geneOnly_seurat)
  Idents(inGene_seurat)<-Idents(geneOnly_seurat)
  
  geneOnly_seurat <- RunUMAP(geneOnly_seurat,dims=1:10,check_duplicates = F)
  inGene_seurat <- RunUMAP(inGene_seurat,dims=1:10,check_duplicates = F)
  outGene_seurat <- RunUMAP(outGene_seurat,dims=1:10,check_duplicates = F)
  
  #geneOnly_seurat <- RunTSNE(geneOnly_seurat,dims=1:10,check_duplicates = F)
  #inGene_seurat <- RunTSNE(inGene_seurat,dims=1:10,check_duplicates = F)
  #outGene_seurat <- RunTSNE(outGene_seurat,dims=1:10,check_duplicates = F)
  
  df<-sample_combined_mat[,colnames(sample_combined)]
  #rownames(df)<-gsub("_","-",rownames(df))
  genes<-rownames(geneOnly)
  readsData<-data.frame(rbind(colSums(df[genes,]),colSums(df[inGene,]),colSums(df[outGene,])))
  rownames(readsData)<-c("genes","HMMinGene","HMMoutGene")
  readsData<-data.frame(t(readsData))
  readsData<-readsData[order(readsData$HMMoutGene),]
  readsData$cell<-rownames(readsData)
  readsData$outpercent = readsData$HMMoutGene / (readsData$HMMinGene + readsData$HMMoutGene)
  readsData$genesPercent=readsData$genes/(readsData$HMMinGene + readsData$HMMoutGene + readsData$genes)
  readsData$cell2 = factor(readsData$cell, levels = as.character(readsData[order(readsData$outpercent),]$cell))
  
  return(c(geneOnly_seurat,inGene_seurat,outGene_seurat,list(readsData)))
}


# generates data that has values of gene expression + hmm features
# useful for featurePlot of gene expression of hmm features
generateDataHMM<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  #subset here
  outGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_0")]
  Sample_HMM<-Sample_HMM[outGene,d4Cells]
  Sample_HMM<-CreateSeuratObject(counts = Sample_HMM, project = "sample_HMM")

  # cluster on gene expression only
  geneOnly<-normVarScaleData(geneOnly)
  geneOnly_seurat <- RunPCA(object = geneOnly, features = rownames(geneOnly))
  geneOnly_seurat <- FindNeighbors(object=geneOnly_seurat,dims=1:10)
  geneOnly_seurat <- FindClusters(object=geneOnly_seurat,resolution=0.2)
  Sample_HMM<-normVarScaleData(Sample_HMM)
  Sample_all <- RunPCA(object = Sample_HMM, features = rownames(Sample_HMM))
  Sample_all <- FindNeighbors(object=Sample_all,dims=1:10)
  Sample_all <- FindClusters(object=Sample_all,resolution=0.2)
  Idents(Sample_all)<-Idents(geneOnly_seurat)
  
  Sample_variable <- RunPCA(object = Sample_HMM, features = VariableFeatures(Sample_HMM))
  Sample_variable <- FindNeighbors(object=Sample_variable,dims=1:10)
  Sample_variable <- FindClusters(object=Sample_variable,resolution=0.2)
  Idents(Sample_variable)<-Idents(geneOnly_seurat)
  
  Sample_all <- RunUMAP(Sample_all,dims=1:10)
  Sample_variable <- RunUMAP(Sample_variable,dims=1:10)
  
  return(c(Sample_all,Sample_variable))
}

generateDataGeneOnly<-function(geneMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)

  out<-GetAssayData(object = geneOnly, slot = "counts")

  return(data.frame(t(out)))
}



# input is list of seurat objects with all features listed in each seurat object
# look at only variable HMM features and variable gene features
calculateSilRawVariable<-function(input){
  day4GeneMat<-GetAssayData(object = input[[1]])#, slot = "scale.data")
  day4Gene<-t(as.matrix(day4GeneMat[VariableFeatures(input[[1]]),])) %>% dist# %>% hclust(method="ward.D2")# %>% as.dendrogram
  
  hmmMat<-GetAssayData(object = input[[2]])#, slot = "scale.data")
  outGene<-VariableFeatures(input[[2]])[endsWith(VariableFeatures(input[[2]]),"-0")]
  inGene<-VariableFeatures(input[[2]])[endsWith(VariableFeatures(input[[2]]),"-1")]

  day4out<-t(as.matrix(hmmMat[outGene,])) %>% dist #%>% hclust(method="ward.D2")# %>% as.dendrogram
  day4in<-t(as.matrix(hmmMat[inGene,])) %>% dist #%>% hclust(method="ward.D2")# %>% as.dendrogram
  
  geneSil<-silhouette(x=as.numeric(Idents(input[[1]])),dist=as.matrix(day4Gene))
  inSil<-silhouette(x=as.numeric(Idents(input[[2]])),dist=as.matrix(day4in))
  outSil<-silhouette(x=as.numeric(Idents(input[[2]])),dist=as.matrix(day4out))
  
  return(list(geneSil,inSil,outSil))
}

# look at only variable HMM features and variable gene features
calcSilfromUMAP<-function(input){
  
  umapGene<-data.frame(input[[1]]@reductions$umap@cell.embeddings)
  umapIn<-data.frame(input[[2]]@reductions$umap@cell.embeddings)
  umapOut<-data.frame(input[[3]]@reductions$umap@cell.embeddings)
  
  geneSil<-silhouette(x=as.numeric(Idents(input[[1]])),dist=dist(umapGene))
  inSil<-silhouette(x=as.numeric(Idents(input[[2]])),dist=dist(umapIn))
  outSil<-silhouette(x=as.numeric(Idents(input[[3]])),dist=dist(umapOut))
  
  out<-cbind(geneSil[,3],inSil[,3],outSil[,3])
  
  return(out)
}

# used for when everything is contained in 1 seurat object
calcSilfromUMAP2<-function(input){
  
  umapGene<-data.frame(input@reductions$umap_gene@cell.embeddings)
  umapIn<-data.frame(input@reductions$umap_aTAR@cell.embeddings)
  umapOut<-data.frame(input@reductions$umap_uTAR@cell.embeddings)
  
  labels<-input[["RNA_snn_res.0.2"]][,1]
  labels<-as.numeric(levels(labels))[labels]
  geneSil<-silhouette(x=labels,dist=dist(umapGene))
  inSil<-silhouette(x=labels,dist=dist(umapIn))
  outSil<-silhouette(x=labels,dist=dist(umapOut))
  
  out<-cbind(geneSil[,3],inSil[,3],outSil[,3])
  
  return(out)
}

calcSilfromUMAPMean<-function(input){
  
  umapGene<-data.frame(input@reductions$umap_gene@cell.embeddings)
  umapIn<-data.frame(input@reductions$umap_aTAR@cell.embeddings)
  umapOut<-data.frame(input@reductions$umap_uTAR@cell.embeddings)
  
  labels<-input[["RNA_snn_res.0.2"]][,1]
  labels<-as.numeric(levels(labels))[labels]
  geneSil<-silhouette(x=labels,dist=dist(umapGene))
  inSil<-silhouette(x=labels,dist=dist(umapIn))
  outSil<-silhouette(x=labels,dist=dist(umapOut))
  
  out<-cbind(geneSil[,3],inSil[,3],outSil[,3])
  
  return(colMeans(out))
}

calcSilfromUMAPMeanFinal<-function(input){
  
  umapGene<-data.frame(input@reductions$umap_gene@cell.embeddings)
  umapIn<-data.frame(input@reductions$umap_aTAR@cell.embeddings)
  umapOut<-data.frame(input@reductions$umap_uTAR@cell.embeddings)
  
  labels<-unique(Idents(input))
  nums<-1:length(unique(Idents(input)))
  names(nums)<-labels
  tt<-nums[Idents(input)]
  names(tt)<-names(Idents(input))
  geneSil<-silhouette(x=tt,dist=dist(umapGene))
  inSil<-silhouette(x=tt,dist=dist(umapIn))
  outSil<-silhouette(x=tt,dist=dist(umapOut))
  
  out<-cbind(geneSil[,3],inSil[,3],outSil[,3])
  
  return(colMeans(out))
}

calcSilfromUMAPMeanFinalLem<-function(input){
  print("finished 1 in list")
  umapGene<-data.frame(input@reductions$umap_gene@cell.embeddings)
  umapIn<-data.frame(input@reductions$umap_aTAR@cell.embeddings)
  umapOut<-data.frame(input@reductions$umap_uTAR@cell.embeddings)
  umapnonVarGenes<-data.frame(input@reductions$umap_nonVarGenes@cell.embeddings)
  
  labels<-unique(Idents(input))
  nums<-1:length(unique(Idents(input)))
  names(nums)<-labels
  tt<-nums[Idents(input)]
  names(tt)<-names(Idents(input))
  geneSil<-silhouette(x=tt,dist=dist(umapGene))
  inSil<-silhouette(x=tt,dist=dist(umapIn))
  outSil<-silhouette(x=tt,dist=dist(umapOut))
  nonVarGenesSil<-silhouette(x=tt,dist=dist(umapnonVarGenes))
  
  out<-cbind(geneSil[,3],inSil[,3],outSil[,3],nonVarGenesSil[,3])
  
  return(colMeans(out))
}


calcSilfromUMAPMean2<-function(input){
  input<-input[[1]]
  umapGene<-data.frame(input@reductions$umap_gene@cell.embeddings)
  umapRand<-data.frame(input@reductions$umap_uTARrand@cell.embeddings)
  umapOut<-data.frame(input@reductions$umap_uTAR@cell.embeddings)
  
  labels<-unique(Idents(input))
  nums<-1:length(unique(Idents(input)))
  names(nums)<-labels
  tt<-nums[Idents(input)]
  names(tt)<-names(Idents(input))
  geneSil<-silhouette(x=tt,dist=dist(umapGene))
  randSil<-silhouette(x=tt,dist=dist(umapRand))
  outSil<-silhouette(x=tt,dist=dist(umapOut))
  
  out<-cbind(geneSil[,3],outSil[,3],randSil[,3])
  
  return(colMeans(out))
}
calcSilfromUMAPMean3<-function(input){
  input<-input[[1]]
  umapGene<-data.frame(input@reductions$umap_gene@cell.embeddings)
  umapCont<-data.frame(input@reductions$umap_nonVarGene@cell.embeddings)
  umapOut<-data.frame(input@reductions$umap_uTAR@cell.embeddings)
  
  labels<-unique(Idents(input))
  nums<-1:length(unique(Idents(input)))
  names(nums)<-labels
  tt<-nums[Idents(input)]
  names(tt)<-names(Idents(input))
  geneSil<-silhouette(x=tt,dist=dist(umapGene))
  contSil<-silhouette(x=tt,dist=dist(umapCont))
  outSil<-silhouette(x=tt,dist=dist(umapOut))
  
  out<-cbind(geneSil[,3],outSil[,3],contSil[,3])
  
  return(colMeans(out))
}


calcSilfromGraph<-function(input){
  
  geneSil<-silhouette(x=as.numeric(Idents(input[[1]])),dist=dist(input[[1]]@graphs$RNA_nn))
  inSil<-silhouette(x=as.numeric(Idents(input[[2]])),dist=dist(input[[2]]@graphs$RNA_nn))
  outSil<-silhouette(x=as.numeric(Idents(input[[3]])),dist=dist(input[[3]]@graphs$RNA_nn))
  
  out<-cbind(geneSil[,3],inSil[,3],outSil[,3])
  
  return(out)
}


calcSilfromUMAPSingle<-function(input){
  umapGene<-data.frame(input@reductions$umap@cell.embeddings)
  geneSil<-silhouette(x=as.numeric(Idents(input)),dist=dist(umapGene))
  return(geneSil)
}


generateHMMMat<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  #subset here
  outGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_0")]
  Sample_HMM<-Sample_HMM[outGene,d4Cells]

  out<-data.frame(t(Sample_HMM))
  return(out)
}

generateHMMMatIn<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  #subset here
  outGene<-rownames(Sample_HMM)[endsWith(rownames(Sample_HMM),"_1")]
  Sample_HMM<-Sample_HMM[outGene,d4Cells]
  
  out<-data.frame(t(Sample_HMM))
  return(out)
}


extractRNA_chicken <- function(seurat.object, sample_name){
  return(t(as.matrix(GetAssayData(seurat.object, slot = "counts")))[colnames(seurat.object)[seurat.object$orig.ident == sample_name],])
}

scanoramaFunc<-function(listOfMat){
  
  day4hmmMat<-listOfMat[[1]]
  day7RVhmmMat<-listOfMat[[1]]
  day7LVhmmMat<-listOfMat[[1]]
  day10RVhmmMat<-listOfMat[[1]]
  day10LVhmmMat<-listOfMat[[1]]
  day14RVhmmMat<-listOfMat[[1]]
  day14LVhmmMat<-listOfMat[[1]]
  
  rownames(day4hmmMat)<-paste0(rownames(day4hmmMat),"_day4")
  rownames(day7RVhmmMat)<-paste0(rownames(day7RVhmmMat),"_day7RV")
  rownames(day7LVhmmMat)<-paste0(rownames(day7LVhmmMat),"_day7LV")
  rownames(day10RVhmmMat)<-paste0(rownames(day10RVhmmMat),"_day10RV")
  rownames(day10LVhmmMat)<-paste0(rownames(day10LVhmmMat),"_day10LV")
  rownames(day14RVhmmMat)<-paste0(rownames(day14RVhmmMat),"_day14RV")
  rownames(day14LVhmmMat)<-paste0(rownames(day14LVhmmMat),"_day14LV")
  
  allDaysMat<-rbind.fill(day4hmmMat,day7RVhmmMat,day7LVhmmMat,day10RVhmmMat,day10LVhmmMat,day14RVhmmMat,day14LVhmmMat)
  allDaysMat[is.na(allDaysMat)]<-0
  
  rownames(allDaysMat)<-c(rownames(day4hmmMat),
                          rownames(day7RVhmmMat),
                          rownames(day7LVhmmMat),
                          rownames(day10RVhmmMat),
                          rownames(day10LVhmmMat),
                          rownames(day14RVhmmMat),
                          rownames(day14LVhmmMat))
  
  datasets<-list(as.matrix(allDaysMat[rownames(day4hmmMat),]),
                 as.matrix(allDaysMat[rownames(day7RVhmmMat),]),
                 as.matrix(allDaysMat[rownames(day7LVhmmMat),]),
                 as.matrix(allDaysMat[rownames(day10RVhmmMat),]),
                 as.matrix(allDaysMat[rownames(day10LVhmmMat),]),
                 as.matrix(allDaysMat[rownames(day14RVhmmMat),]),
                 as.matrix(allDaysMat[rownames(day14LVhmmMat),]))
  
  genes_list<-list(colnames(allDaysMat),
                   colnames(allDaysMat),
                   colnames(allDaysMat),
                   colnames(allDaysMat),
                   colnames(allDaysMat),
                   colnames(allDaysMat),
                   colnames(allDaysMat))
  
  integrated.corrected.data <- scanorama$correct(datasets, genes_list,
                                                 return_dimred=TRUE, return_dense=TRUE,verbose=T)
  
  ### combining corrected_scanorama and non scanorama data matrix
  corrected_scanorama <- t(rbind(integrated.corrected.data[[2]][[1]], integrated.corrected.data[[2]][[2]], integrated.corrected.data[[2]][[3]], integrated.corrected.data[[2]][[4]], integrated.corrected.data[[2]][[5]], integrated.corrected.data[[2]][[6]], integrated.corrected.data[[2]][[7]]))
  original_data<-t(allDaysMat)
  rownames(original_data)<-gsub("_","-",rownames(original_data))
  colnames(corrected_scanorama)<-colnames(original_data)
  rownames(corrected_scanorama)<-rownames(original_data)
  chicken<-CreateSeuratObject(counts = original_data, project = "original_data")
  chicken[["sample"]]<-sapply(strsplit(colnames(chicken),"_"),"[",2)
  
  cluster_idents<-c(paste0("day4_",Idents(day4Data[[1]])),
                    paste0("day7RV_",Idents(day7RVData[[1]])),
                    paste0("day7LV_",Idents(day7LVData[[1]])),
                    paste0("day10RV_",Idents(day10RVData[[1]])),
                    paste0("day10LV_",Idents(day10LVData[[1]])),
                    paste0("day14RV_",Idents(day14RVData[[1]])),
                    paste0("day14LV_",Idents(day14LVData[[1]])))
  chicken[["cluster_ids"]]<-cluster_idents
  
  ### from madhav
  new.assay <- new(
    Class = 'Assay',
    counts =  new(Class = "dgCMatrix"),
    data = corrected_scanorama,
    scale.data = matrix(),
    var.features = vector(),
    meta.features = data.frame(row.names = rownames(x = corrected_scanorama)),
    misc = NULL
  )
  chicken[["scanorama"]] <- new.assay
  DefaultAssay(chicken) <- "scanorama"
  
  # run standard scRNA-seq analysis steps
  chicken<-FindVariableFeatures(chicken, selection.method = "vst", nfeatures=2000)
  chicken<-ScaleData(chicken, features = rownames(chicken))
  chicken<-RunPCA(object = chicken, features = rownames(chicken))
  chicken<-RunUMAP(chicken,dims=1:10,check_duplicates = F)
  DimPlot(chicken,reduction = "umap")
  
  return(chicken)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

ggsaveMW<-function(pp,outFile,wid,hei){
  ggsave(filename = outFile,
         plot=pp,
         width=wid,
         height=hei,
         units="in")
}

plot_data_column = function (data, column) {
  ggplot(data, aes_string(x = column)) +
    geom_histogram(fill = "lightgreen") +
    xlab(column)
}

umapPlain<-function(data){
  out<-DimPlot(data,
               reduction="umap",
               label=F,
               pt.size = 0.01,
               cols=as.vector(alphabet(n=length(unique(Idents(data))))))+
    #NoLegend()+
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank())+
          #aspect.ratio=1)+
    NoLegend()+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  return(out)
}

umapPlainLegend<-function(data,reduc="umap_gene",group="ident",noLegend=T){
  out<-DimPlot(data,
               reduction=reduc,
               label=F,
               pt.size = 0.001,
               group.by = group,
               cols=as.vector(alphabet(n=length(unique(Idents(data))))))+
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank(),
          aspect.ratio=1)+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  if(noLegend){
    out<-out+NoLegend()
  }
  return(out)
}


umapPlainGG<-function(test){
  umapCoords<-data.frame(test@reductions$umap@cell.embeddings)
  umapCoords$colr<-Idents(test)
  
  out<-ggplot(umapCoords,
              aes(x=UMAP_1,
                  y=UMAP_2,
                  color=colr))+
    geom_point(size=0.01)+
    scale_color_manual(values=as.vector(alphabet(n=length(unique(Idents(test))))))+
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank(),
          aspect.ratio=1)+
    NoLegend()+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  return(out)
}



makeNewScanoramaAssay<-function(seuratObj,scanorama_data){
  new.assay <- new(
    Class = 'Assay',
    counts =  new(Class = "dgCMatrix"),
    data = scanorama_data,
    scale.data = matrix(),
    var.features = vector(),
    meta.features = data.frame(row.names = rownames(x = scanorama_data)),
    misc = NULL
  )
  
  seuratObj[["scanorama"]] <- new.assay
  DefaultAssay(seuratObj) <- "scanorama"
  
  return(seuratObj)
}

makeNewRawAssay<-function(seuratObj,raw.data){
  new.assay <- new(
    Class = 'Assay',
    counts =  new(Class = "dgCMatrix"),
    data = raw.data,
    scale.data = matrix(),
    var.features = vector(),
    meta.features = data.frame(row.names = rownames(x = raw.data)),
    misc = NULL
  )
  
  seuratObj[["raw.data"]] <- new.assay

  return(seuratObj)
}

calculatePCLoadingForChicken<-function(input){
  # define genes, aTARs, uTARs
  numdash<-str_count(rownames(input),"-")
  outInd<-numdash>=5 & endsWith(rownames(input),"-0")
  inInd<-numdash>=5 & endsWith(rownames(input),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(input)[geneInd]
  inGene<-rownames(input)[inInd]
  outGene<-rownames(input)[outInd]
  
  # run PCA on all TARs
  uTARs_seurat<-RunPCA(object=input,features=c(inGene,outGene))
  
  pcLoadings<-data.frame(uTARs_seurat[["pca"]]@feature.loadings)
  pcLoadingsAb<-abs(pcLoadings)
  
  # update aTARs and uTARs, some have 0 variance and not in pcLoadings
  numdash<-str_count(rownames(pcLoadingsAb),"-")
  outInd<-numdash>=5 & endsWith(rownames(pcLoadingsAb),"-0")
  inInd<-numdash>=5 & endsWith(rownames(pcLoadingsAb),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(pcLoadingsAb)[geneInd]
  inGene<-rownames(pcLoadingsAb)[inInd]
  outGene<-rownames(pcLoadingsAb)[outInd]
  
  
  loadings<-data.frame(t(rbind(data.frame(colSums(pcLoadingsAb[inGene,]),colSums(pcLoadingsAb[outGene,])))))
  #loadings$variable<-c("aTARs","uTARs")
  uTARLoading<-loadings[2,]/(colSums(loadings))
  jackStrawGG<-data.frame(as.numeric(uTARs_seurat[["pca"]]@stdev))
  output<-cbind(jackStrawGG,t(loadings))
  colnames(output)<-c("std","aTAR_PCload","uTAR_PCload")
  
  return(t(output))
  
  
  loadingsGG<-melt(loadings)
  
  
  
  
  
  
  colnames(loadingsGG)<-c("Var1","Var2","Value")
  pcLoadingsTotal<-ggplot(loadingsGG,aes(x=Var2,y=Value,fill=Var1))+
    geom_bar(stat="identity")+
    scale_x_discrete(labels=as.character(c(10,20,30,40,50)),
                     breaks=c("PC_10","PC_20","PC_30","PC_40","PC_50"))+
    theme_bw()+
    theme(legend.position="top",
          axis.text.y=element_text(size=14,face="bold"),
          axis.text.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=14,face="bold"),
          axis.title.x=element_text(size=14,face="bold"))+
    ylab("total abs. PC loading values")+
    xlab("PC")+
    labs(fill="features")
  
  jackStrawGG$Var2<-as.factor(1:50)
  colnames(jackStrawGG)<-c("Value","Var2")
  jackStrawP<-ggplot(jackStrawGG,aes(x=Var2,y=Value))+
    geom_bar(stat="identity")+
    ylab("standard deviation")+
    xlab("PC")+
    scale_x_discrete(labels=as.character(c(10,20,30,40,50)),
                     breaks=c(10,20,30,40,50))+
    theme_bw()+
    theme(axis.text.y=element_text(size=14,face="bold"),
          axis.text.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=14,face="bold"),
          axis.title.x=element_text(size=14,face="bold"))
  
  plot_grid(pcLoadings,jackStrawP,nrow=2,align="v")
  
}

addSmallLegend <- function(myPlot, legPosY=1, pointSize = 2, textSize = 6, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize,face="bold"), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          #legend.position = c(0.5,legPosY),
          #legend.direction="horizontal",
          legend.background = element_blank())+
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
}
addSmallLegend2 <- function(myPlot, legPosY=1, pointSize = 2, textSize = 6, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize,face="bold"), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          #legend.position = c(0.5,legPosY),
          #legend.direction="horizontal",
          legend.background = element_blank())+
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
}

addSmallLegend3 <- function(myPlot, legPosY=1, pointSize = 2, textSize = 6, spaceLegend = 0.5) {
  myPlot +
    #guides(shape = guide_legend(override.aes = list(size = pointSize)),
    #       color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize,face="bold"), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          #legend.position = c(0.5,legPosY),
          #legend.direction="horizontal",
          legend.background = element_blank())+
    guides(color=guide_legend(override.aes=list(fill=NA)))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
}

addSmallLegendDot <- function(myPlot, pointSize = 2, textSize = 6, spaceLegend = 0.5) {
  myPlot +
    #guides(shape = guide_legend(override.aes = list(size = pointSize)),
    #       color = guide_legend(override.aes = list(size = pointSize))
    #       ) +
    theme(legend.title = element_text(size = textSize,face="bold"), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          #legend.position = c(0.5,legPosY),
          #legend.direction="horizontal",
          legend.background = element_blank())+
    #guides(color=guide_legend(override.aes=list(fill=NA)))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
}


adjustThemeGG<-function(myPlot){
  myPlot+
    theme(plot.title = element_text(size=8,face="bold"),
          axis.title.x = element_text(size=8,face="bold"),
          axis.title.y = element_text(size=8,face="bold"),
          axis.text.x=element_text(size = 6,face="bold"),
          axis.text.y=element_text(size=6,face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
}

adjustThemeGGSize<-function(myPlot,title.size=8,x.title.size=8,y.title.size=8,text.x=6,text.y=6){
  myPlot+
    theme(plot.title = element_text(size=title.size,face="bold"),
          axis.title.x = element_text(size=x.title.size,face="bold"),
          axis.title.y = element_text(size=y.title.size,face="bold"),
          axis.text.x=element_text(size = text.x,face="bold"),
          axis.text.y=element_text(size=text.y,face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
}


adjustThemeGGDot<-function(myPlot){
  myPlot+
    theme(plot.title = element_text(size=8,face="bold"),
          axis.title.x = element_blank(),#element_text(size=8,face="bold"),
          axis.title.y = element_blank(),#element_text(size=8,face="bold"),
          axis.text.x=element_text(size = 6,face="bold"),
          axis.text.y=element_text(size=6,face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
}

adjustThemeGGumap<-function(myPlot,titleSize=8){
  myPlot+
    theme(plot.title = element_text(size=titleSize,face="bold"),
          axis.title.x = element_blank(),#element_text(size=8,face="bold"),
          axis.title.y = element_blank(),#element_text(size=8,face="bold"),
          axis.text.x=element_blank(),#element_text(size=8,face="bold"),
          axis.text.y=element_blank(),#element_text(size=8,face="bold"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
}


# TARs is a 3 column vector without direction for now
checkInOrOutGene_noDir<-function(HMManno_bare,refFlat){
  refName<-unlist(strsplit(refFlat,"/"))[length(unlist(strsplit(refFlat,"/")))]
  # read in hg_38 ref genes
  gene_ref <- read.delim(refFlat, header=F, comment.char="#")
  gene_ref_bare<-gene_ref[,c("V1","V3","V4","V5","V6")] ####### edit this based on file type
  colnames(gene_ref_bare)<-c("gene","chr","direction","start","end")
  gene_ref_bare$gene<-as.character(gene_ref_bare$gene)
  gene_ref_bare$chr<-as.character(gene_ref_bare$chr)
  gene_ref_bare$direction<-as.character(gene_ref_bare$direction)
  gene_ref_bare$start<-as.numeric(gene_ref_bare$start)
  gene_ref_bare$end<-as.numeric(gene_ref_bare$end)
  
  colnames(HMManno_bare)<-c("chr","start","end")
  HMManno_bare$chr<-as.character(levels(HMManno_bare$chr))[HMManno_bare$chr]
  HMManno_bare$start<-as.numeric(levels(HMManno_bare$start))[HMManno_bare$start]
  HMManno_bare$end<-as.numeric(levels(HMManno_bare$end))[HMManno_bare$end]
  HMManno_bare[is.na(HMManno_bare)]<-0
  
  ########################################################
  checkIfExistGene_noDir<-function(input,gene_ref){
    chrom<-input[[1]]
    startPos<-as.numeric(input[[2]])
    endPos<-as.numeric(input[[3]])
    gene_ref<-gene_ref # bed of entire genes
    
    # 1: either partially in or completely inside gene
    chrMatch_gene<-(chrom==gene_ref$chr)
    qinrstartMatch<-(startPos>=gene_ref$start)&(startPos<=gene_ref$end)
    qinrendMatch<-(endPos<=gene_ref$end)&(endPos>=gene_ref$start)
    qinroutGeneAll_gene<-(chrMatch_gene*qinrstartMatch*qinrendMatch) # query completely inside ref 
    rinqstartMatch<-(startPos<(gene_ref$start))
    rinqendMatch<-(endPos>(gene_ref$end))
    rinqoutGeneAll<-(chrMatch_gene*rinqstartMatch*rinqendMatch) # ref completely inside query
    # also check partial overlap from beginning
    partialStart<-(chrMatch_gene*rinqstartMatch*qinrendMatch)
    partialEnd<-(chrMatch_gene*qinrstartMatch*rinqendMatch)
    if(sum(qinroutGeneAll_gene)|sum(rinqoutGeneAll)|sum(partialStart)|sum(partialEnd)){
      # g's are for gene names
      g1<-as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
      g2<-as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
      g3<-as.character(gene_ref$gene[as.logical(partialStart)])
      g4<-as.character(gene_ref$gene[as.logical(partialEnd)])
      gAll<-unique(c(g1,g2,g3,g4))
      return(paste(c(gAll,"1"),collapse="_"))
    }else{
      return(0)
    }
  }
    
  df<-HMManno_bare
  #HMManno_bare_sample<-df[sample(nrow(df),10),]
  #HMManno$inAnno<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExist,exon_ref=gene_gtf_bare,gene_ref=gene_ref_bare)
  library(parallel)
  num_cores<-detectCores()
  num_coresUse<-floor(num_cores/5)
  clust<-makeCluster(num_coresUse)
  HMManno_bare$inGene<-parApply(cl = clust,X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene_noDir,gene_ref=gene_ref_bare)
  stopCluster(clust)
  
  ########################### uncomment below to make refFlat file
  # make dataframe into refFlat file format
  HMManno_bare$geneName<-paste0(HMManno_bare$chr,"_",HMManno_bare$start,"_",HMManno_bare$end,"_",HMManno_bare$inGene,sep="")
  return(HMManno_bare$geneName)
}


# TARs is a 3 column vector without direction for now
checkInOrOutGene_noDir_v2<-function(HMManno_bare,refFlat){
  refName<-unlist(strsplit(refFlat,"/"))[length(unlist(strsplit(refFlat,"/")))]
  # read in hg_38 ref genes
  gene_ref <- read.delim(refFlat, header=F, comment.char="#")
  gene_ref_bare<-gene_ref[,c("V1","V3","V4","V5","V6")] ####### edit this based on file type
  colnames(gene_ref_bare)<-c("gene","chr","direction","start","end")
  gene_ref_bare$gene<-as.character(gene_ref_bare$gene)
  gene_ref_bare$chr<-as.character(gene_ref_bare$chr)
  gene_ref_bare$direction<-as.character(gene_ref_bare$direction)
  gene_ref_bare$start<-as.numeric(gene_ref_bare$start)
  gene_ref_bare$end<-as.numeric(gene_ref_bare$end)
  
  colnames(HMManno_bare)<-c("chr","start","end")
  HMManno_bare$chr<-as.character(levels(HMManno_bare$chr))[HMManno_bare$chr]
  HMManno_bare$start<-as.numeric(levels(HMManno_bare$start))[HMManno_bare$start]
  HMManno_bare$end<-as.numeric(levels(HMManno_bare$end))[HMManno_bare$end]
  HMManno_bare[is.na(HMManno_bare)]<-0
  
  ########################################################
  checkIfExistGene_noDir<-function(input,gene_ref){
    chrom<-input[[1]]
    startPos<-as.numeric(input[[2]])
    endPos<-as.numeric(input[[3]])
    gene_ref<-gene_ref # bed of entire genes
    
    # 1: either partially in or completely inside gene
    chrMatch_gene<-(chrom==gene_ref$chr)
    qinrstartMatch<-(startPos>=gene_ref$start)&(startPos<=gene_ref$end)
    qinrendMatch<-(endPos<=gene_ref$end)&(endPos>=gene_ref$start)
    qinroutGeneAll_gene<-(chrMatch_gene*qinrstartMatch*qinrendMatch) # query completely inside ref 
    rinqstartMatch<-(startPos<(gene_ref$start))
    rinqendMatch<-(endPos>(gene_ref$end))
    rinqoutGeneAll<-(chrMatch_gene*rinqstartMatch*rinqendMatch) # ref completely inside query
    # also check partial overlap from beginning
    partialStart<-(chrMatch_gene*rinqstartMatch*qinrendMatch)
    partialEnd<-(chrMatch_gene*qinrstartMatch*rinqendMatch)
    if(sum(qinroutGeneAll_gene)|sum(rinqoutGeneAll)|sum(partialStart)|sum(partialEnd)){
      # g's are for gene names
      g1<-as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
      g2<-as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
      g3<-as.character(gene_ref$gene[as.logical(partialStart)])
      g4<-as.character(gene_ref$gene[as.logical(partialEnd)])
      gAll<-unique(c(g1,g2,g3,g4))
      return(paste(c(gAll,"1"),collapse="_"))
    }else{
      return(0)
    }
  }
  
  df<-HMManno_bare
  #HMManno_bare_sample<-df[sample(nrow(df),10),]
  #HMManno$inAnno<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExist,exon_ref=gene_gtf_bare,gene_ref=gene_ref_bare)
  library(parallel)
  num_cores<-detectCores()
  num_coresUse<-floor(num_cores/5)
  clust<-makeCluster(num_coresUse)
  HMManno_bare$inGene<-parApply(cl = clust,X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene_noDir,gene_ref=gene_ref_bare)
  stopCluster(clust)
  
  ########################### uncomment below to make refFlat file
  # make dataframe into refFlat file format
  #HMManno_bare$geneName<-paste0(HMManno_bare$chr,"_",HMManno_bare$start,"_",HMManno_bare$end,"_",HMManno_bare$inGene,sep="")
  return(HMManno_bare$inGene)
}




calcPercentOutType<-function(helloInput,n){
  numdash<-str_count(rownames(helloInput),"-")
  outInd<-numdash>=n & endsWith(rownames(helloInput),"-0")
  inInd<-numdash>=n & endsWith(rownames(helloInput),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(helloInput)[geneInd]
  inGene<-rownames(helloInput)[inInd]
  outGene<-rownames(helloInput)[outInd]
  
  out<-helloInput[c(inGene,outGene),]
  helloInput[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  return(helloInput)

}

plotVlnPlotType<-function(input,yscale){
  p<-VlnPlot(input,
          features="percent.outHMM",
          pt.size=0,
          cols=as.vector(alphabet(n=length(unique(Idents(input))))),
          fill=NULL)+#as.vector(alphabet(n=length(unique(Idents(input))))))+
    ylim(yscale[1],yscale[2])+
    NoLegend()+
    ggtitle("reads outside annotations %")+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_text(size=8,face="bold"),
          title = element_text(size=8,face="bold"))
  return(p)
}


plotVlnPlotTypeGG<-function(input,yscale=c(0,20),addPoints=F){
  dataToPlot<-data.frame(cbind(input$percent.outHMM,Idents(input)))
  colnames(dataToPlot)<-c("percent.out","cellType")
  dataToPlot$cellType<-factor(dataToPlot$cellType,levels = 1:length(Idents(input)))
  
  colours=as.vector(alphabet(n=length(unique(Idents(input)))))
  
  for(i in 1:length(unique(dataToPlot$cellType))){
    dataToPlot$color[dataToPlot$cellType==i]<-colours[i]
  }
  
  if(addPoints){
    readsP<-ggplot(dataToPlot,aes(x=dataToPlot$cellType,y=dataToPlot$percent.out,fill=dataToPlot$cellType,color=dataToPlot$cellType))+
      geom_point(size=0.00001)+
      coord_cartesian(ylim=c(yscale[1],yscale[2]))+
      NoLegend()+
      geom_violin()+
      scale_fill_manual(values=colours)+
      scale_color_manual(values=colours)
  } else {
    readsP<-ggplot(dataToPlot,aes(x=dataToPlot$cellType,y=dataToPlot$percent.out,fill=dataToPlot$cellType,color=dataToPlot$cellType))+
      geom_violin()+
      coord_cartesian(ylim=c(yscale[1],yscale[2]))+
      NoLegend()+
      scale_fill_manual(values=colours)+
      scale_color_manual(values=colours)
  }
  
  readsP<-adjustThemeGG(readsP)
  readsP<-addSmallLegend(readsP)+NoLegend()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          #axis.text.y=element_blank(),
          axis.ticks.x = element_blank())+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  
  return(readsP)
}



plotTotalOut<-function(input,yscale){
  dataToPlot<-data.frame(cbind(input$nCount_aTAR,input$nCount_uTAR,Idents(input)))
  dataGG<-data.frame()
  for(i in 1: length(unique(dataToPlot[,3]))){
    dataGG<-rbind(dataGG,c(colSums(dataToPlot[dataToPlot$X3==i,c(1,2)]),i))
  }
  dataGG$ident<-as.character(dataGG[,3])
  dataGG$out<-(dataGG[,2]/(dataGG[,2]+dataGG[,1]))*100
  
  readsP<-ggplot(dataGG,aes(x=ident,
                    y=out,
                    col=ident,
                    fill=ident,
                    width=0.75))+
    geom_bar(stat="identity")+
    ylim(yscale[1],yscale[2])+
    scale_fill_manual(values=as.vector(alphabet(n=nrow(dataGG))))+
    scale_color_manual(values=as.vector(alphabet(n=nrow(dataGG))))
  
  readsP<-adjustThemeGG(readsP)
  readsP<-addSmallLegend(readsP)+NoLegend()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
    theme(axis.text.x = element_blank())
  
  return(readsP)
}

getNFromList <- function(lst, n){
  sapply(lst, `[`, n)
}

drawVennDiagram<-function(outGeneSeurat,numDashes,topPer){
  
  numdash<-str_count(rownames(outGeneSeurat),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(outGeneSeurat)[geneInd]
  inGene<-rownames(outGeneSeurat)[inInd]
  outGene<-rownames(outGeneSeurat)[outInd]
  
  # look at ven diagram
  rawMat<-as.matrix(GetAssayData(object = outGeneSeurat[outGene,], slot = "counts"))
  
  # look at total expression of each tar
  sumExpUTar<-rowSums(rawMat[outGene,])
  # look at PC loading of just uTARs, top 10 PCs
  pcLoadings<-data.frame(outGeneSeurat[["pca"]]@feature.loadings)
  pcLoadingsAb<-abs(pcLoadings)
  sumPCLoadUTar<-sort(rowSums(pcLoadingsAb[,1:10]),decreasing=T)
  
  # look for common uTARs
  commonTARs<-intersect(names(sumExpUTar),names(sumPCLoadUTar))
  sumExpUTar<-sumExpUTar[commonTARs]
  sumPCLoadUTar<-sumPCLoadUTar[commonTARs]
  
  sumExpUTar<-sort(sumExpUTar,decreasing=T)
  sumExpUTar<-sumExpUTar[1:(length(sumExpUTar)*topPer)]# filter for top 10% uTARs of this metric
  namesExp<-names(sumExpUTar)
  
  sumPCLoadUTar<-sort(sumPCLoadUTar,decreasing=T)
  sumPCLoadUTar<-sumPCLoadUTar[1:(length(sumPCLoadUTar)*topPer)]# filter for top 10% uTARs of this metric
  namesPC<-names(sumPCLoadUTar)
  
  # draw venn diagram make sure to library(VennDiagram)
  VennDay4<-draw.pairwise.venn(area1=length(namesExp),
                               area2=length(namesPC),
                               cross.area = length(intersect(namesExp,namesPC)),
                               category=c("top 10% uTARs\ntotal reads",
                                          "top 10% uTARs\nloading, 10PCs"),
                               fill=c("#0000B2","#640000"),
                               fontface="bold",
                               cat.fontface = "bold",
                               cex=1,
                               cat.cex = 1,
                               ind=F,
                               aspect.ratio=1)
  #out<-grid.draw(VennDay4)
  return(VennDay4)
}

drawVennDiagram2<-function(outGeneSeurat,pca,rawMat,numDashes=5,PCs,topPer){
  # look at ven diagram
  numdash<-str_count(rownames(outGeneSeurat),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(outGeneSeurat)[geneInd]
  inGene<-rownames(outGeneSeurat)[inInd]
  outGene<-rownames(outGeneSeurat)[outInd]
  
  # look at total expression of each tar
  sumExpUTar<-rowSums(rawMat[outGene,])
  # look at PC loading of just uTARs, top 10 PCs
  eigenVec<-data.frame(outGeneSeurat[[pca]]@feature.loadings)
  pcLoadings<-sweep(eigenVec,2,outGeneSeurat[[pca]]@stdev,"*")
  sumPCLoadUTar<-sort(rowSums(abs(pcLoadings[,1:PCs])),decreasing=T)
  
  # look for common uTARs
  commonTARs<-intersect(names(sumExpUTar),names(sumPCLoadUTar))
  sumExpUTar<-sumExpUTar[commonTARs]
  sumPCLoadUTar<-sumPCLoadUTar[commonTARs]
  
  sumExpUTar<-sort(sumExpUTar,decreasing=T)
  sumExpUTar<-sumExpUTar[1:(length(sumExpUTar)*topPer)]# filter for top 10% uTARs of this metric
  namesExp<-names(sumExpUTar)
  
  sumPCLoadUTar<-sort(sumPCLoadUTar,decreasing=T)
  sumPCLoadUTar<-sumPCLoadUTar[1:(length(sumPCLoadUTar)*topPer)]# filter for top 10% uTARs of this metric
  namesPC<-names(sumPCLoadUTar)
  
  # draw venn diagram make sure to library(VennDiagram)
  VennDay4<-draw.pairwise.venn(area1=length(namesExp),
                               area2=length(namesPC),
                               cross.area = length(intersect(namesExp,namesPC)),
                               category=c("uTARs total reads",
                                          "uTARs loadings"),
                               fill=c("#0000B2","#640000"),
                               fontface="bold",
                               cat.fontface = "bold",
                               cex=1,
                               cat.cex = 1,
                               ind=F)
  #out<-grid.draw(VennDay4)
  return(VennDay4)
}

drawVennDiagramTrueBulk<-function(outGeneSeurat,pca,rawMat,numDashes=5,PCs,minPCLoad=0.5,minBulkCov=10000){
  # look at ven diagram
  numdash<-str_count(rownames(outGeneSeurat),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(outGeneSeurat)[geneInd]
  inGene<-rownames(outGeneSeurat)[inInd]
  outGene<-rownames(outGeneSeurat)[outInd]
  
  # look at total expression of each tar
  sumExpUTar<-rowSums(rawMat)
  # look at PC loading of just uTARs, top 10 PCs
  eigenVec<-data.frame(outGeneSeurat[[pca]]@feature.loadings)
  pcLoadings<-sweep(eigenVec,2,outGeneSeurat[[pca]]@stdev,"*")
  sumPCLoadUTar<-sort(rowSums(abs(pcLoadings[,1:PCs])),decreasing=T)
  
  # look for common uTARs
  commonTARs<-intersect(names(sumExpUTar),names(sumPCLoadUTar))
  sumExpUTar<-sumExpUTar[commonTARs]
  sumPCLoadUTar<-sumPCLoadUTar[commonTARs]
  
  sumExpUTar<-sort(sumExpUTar,decreasing=T)
  sumExpUTar<-sumExpUTar[sumExpUTar > minBulkCov ]# filter for min PC loading top 10% uTARs of this metric
  namesExp<-names(sumExpUTar)
  
  sumPCLoadUTar<-sort(sumPCLoadUTar,decreasing=T)
  sumPCLoadUTar<-sumPCLoadUTar[sumPCLoadUTar > minPCLoad]# filter for top 10% uTARs of this metric
  namesPC<-names(sumPCLoadUTar)
  
  # draw venn diagram make sure to library(VennDiagram)
  VennDay4<-draw.pairwise.venn(area1=length(namesExp),
                               area2=length(namesPC),
                               cross.area = length(intersect(namesExp,namesPC)),
                               category=c("uTARs total reads",
                                          "uTARs loadings"),
                               fill=c("#0000B2","#640000"),
                               fontface="bold",
                               cat.fontface = "bold",
                               cex=1,
                               cat.cex = 1,
                               ind=F)
  #out<-grid.draw(VennDay4)
  return(VennDay4)
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

checkIfExistGene_noDir<-function(input,gene_ref){
  chrom<-input[[1]]
  startPos<-as.numeric(input[[2]])
  endPos<-as.numeric(input[[3]])
  gene_ref<-gene_ref # bed of entire genes
  
  # 1: either partially in or completely inside gene
  chrMatch_gene<-(chrom==gene_ref$chr)
  qinrstartMatch<-(startPos>=gene_ref$start)&(startPos<=gene_ref$end)
  qinrendMatch<-(endPos<=gene_ref$end)&(endPos>=gene_ref$start)
  qinroutGeneAll_gene<-(chrMatch_gene*qinrstartMatch*qinrendMatch) # query completely inside ref 
  rinqstartMatch<-(startPos<(gene_ref$start))
  rinqendMatch<-(endPos>(gene_ref$end))
  rinqoutGeneAll<-(chrMatch_gene*rinqstartMatch*rinqendMatch) # ref completely inside query
  # also check partial overlap from beginning
  partialStart<-(chrMatch_gene*rinqstartMatch*qinrendMatch)
  partialEnd<-(chrMatch_gene*qinrstartMatch*rinqendMatch)
  if(sum(qinroutGeneAll_gene)|sum(rinqoutGeneAll)|sum(partialStart)|sum(partialEnd)){
    # g's are for gene names
    g1<-as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
    g2<-as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
    g3<-as.character(gene_ref$gene[as.logical(partialStart)])
    g4<-as.character(gene_ref$gene[as.logical(partialEnd)])
    gAll<-unique(c(g1,g2,g3,g4))
    return(paste(c(gAll,"1"),collapse="_"))
  }else{
    return(0)
  }
}

drawVennDiagramFromData<-function(dataInput,pcThresh=0.5, bulkThresh=10000){
  
  sumExpUTar<-dataInput[,1]
  names(sumExpUTar)<-rownames(dataInput)
  sumExpUTar<-sort(sumExpUTar,decreasing=T)
  sumExpUTar<-sumExpUTar[sumExpUTar > bulkThresh ]# filter for min PC loading top 10% uTARs of this metric
  namesExp<-names(sumExpUTar)
  
  sumPCLoadUTar<-dataInput[,2]
  names(sumPCLoadUTar)<-rownames(dataInput)
  sumPCLoadUTar<-sort(sumPCLoadUTar,decreasing=T)
  sumPCLoadUTar<-sumPCLoadUTar[sumPCLoadUTar > pcThresh]# filter for top 10% uTARs of this metric
  namesPC<-names(sumPCLoadUTar)
  
  
  # draw venn diagram make sure to library(VennDiagram)
  VennDay4<-draw.pairwise.venn(area1=length(namesExp),
                               area2=length(namesPC),
                               cross.area = length(intersect(namesExp,namesPC)),
                               category=c("uTARs total reads",
                                          "uTARs loadings"),
                               fill=c("#0000B2","#640000"),
                               fontface="bold",
                               cat.fontface = "bold",
                               cex=1,
                               cat.cex = 1,
                               ind=F)
  #out<-grid.draw(VennDay4)
  return(VennDay4)
  
}

drawVennDiagramFromDataNoTitle<-function(dataInput,pcThresh=0.5, bulkThresh=10000){
  
  sumExpUTar<-dataInput[,1]
  names(sumExpUTar)<-rownames(dataInput)
  sumExpUTar<-sort(sumExpUTar,decreasing=T)
  sumExpUTar<-sumExpUTar[sumExpUTar > bulkThresh ]# filter for min PC loading top 10% uTARs of this metric
  namesExp<-names(sumExpUTar)
  
  sumPCLoadUTar<-dataInput[,2]
  names(sumPCLoadUTar)<-rownames(dataInput)
  sumPCLoadUTar<-sort(sumPCLoadUTar,decreasing=T)
  sumPCLoadUTar<-sumPCLoadUTar[sumPCLoadUTar > pcThresh]# filter for top 10% uTARs of this metric
  namesPC<-names(sumPCLoadUTar)
  
  
  # draw venn diagram make sure to library(VennDiagram)
  VennDay4<-draw.pairwise.venn(area1=length(namesExp),
                               area2=length(namesPC),
                               cross.area = length(intersect(namesExp,namesPC)),
                               fill=c("#0000B2","#640000"),
                               fontface="bold",
                               cat.fontface = "bold",
                               cex=1,
                               cat.cex = 1,
                               ind=F)
  #out<-grid.draw(VennDay4)
  return(VennDay4)
  
}

drawScatter<-function(outGeneSeurat,pca,rawMat,numDashes=5,PCs=10,minBulk=100){
  # look at ven diagram
  numdash<-str_count(rownames(outGeneSeurat),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(outGeneSeurat)[geneInd]
  inGene<-rownames(outGeneSeurat)[inInd]
  outGene<-rownames(outGeneSeurat)[outInd]
  
  # look at total expression of each tar
  sumExpUTar<-rowSums(rawMat[outGene,])
  # look at PC loading of just uTARs, top 10 PCs
  eigenVec<-data.frame(outGeneSeurat[[pca]]@feature.loadings)
  pcLoadings<-sweep(eigenVec,2,outGeneSeurat[[pca]]@stdev,"*")
  sumPCLoadUTar<-sort(rowSums(abs(pcLoadings[,1:PCs])),decreasing=T)
  
  # look for common uTARs
  commonTARs<-intersect(names(sumExpUTar),names(sumPCLoadUTar))
  sumExpUTar<-sumExpUTar[commonTARs]
  sumPCLoadUTar<-sumPCLoadUTar[commonTARs]
  
  toPlot<-data.frame(cbind(sumExpUTar,sumPCLoadUTar))
  toPlotTest<-toPlot[toPlot$sumExpUTar>minBulk,]
  stdDevP<-ggplot(toPlotTest,aes(x=sumExpUTar,y=sumPCLoadUTar))+
    geom_bin2d(bins=100)+scale_fill_continuous(type="viridis")+theme_bw()+scale_y_log10()+scale_x_log10()+
    labs(title="2D histogram of uTARs",y="uTAR PC loading",x="uTAR pseudo-bulk read coverage")
  stdDevP<-stdDevP+
      theme(plot.title = element_text(size=10,face="bold"),
            axis.title.x = element_text(size=10,face="bold"),
            axis.title.y = element_text(size=10,face="bold"),
            axis.text.x=element_text(size = 8,face="bold"),
            axis.text.y=element_text(size=8,face="bold"))+
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  #out<-grid.draw(VennDay4)
  return(stdDevP)
}

drawScatterTrueBulk<-function(outGeneSeurat,pca,rawMat,numDashes=5,PCs=10,minBulk=100,returnData=F){
  # look at ven diagram
  numdash<-str_count(rownames(outGeneSeurat),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(outGeneSeurat)[geneInd]
  inGene<-rownames(outGeneSeurat)[inInd]
  outGene<-rownames(outGeneSeurat)[outInd]
  
  # look at total expression of each tar
  sumExpUTar<-rowSums(rawMat)
  # look at PC loading of just uTARs, top 10 PCs
  eigenVec<-data.frame(outGeneSeurat[[pca]]@feature.loadings)
  pcLoadings<-sweep(eigenVec,2,outGeneSeurat[[pca]]@stdev,"*")
  sumPCLoadUTar<-sort(rowSums(abs(pcLoadings[,1:PCs])),decreasing=T)
  
  # look for common uTARs
  commonTARs<-intersect(names(sumExpUTar),names(sumPCLoadUTar))
  sumExpUTar<-sumExpUTar[commonTARs]
  sumPCLoadUTar<-sumPCLoadUTar[commonTARs]
  
  toPlot<-data.frame(cbind(sumExpUTar,sumPCLoadUTar))
  toPlotTest<-toPlot[toPlot$sumExpUTar>minBulk,]
  if(returnData==T){
    return(toPlotTest)
  }
  stdDevP<-ggplot(toPlotTest,aes(x=sumExpUTar,y=sumPCLoadUTar))+
    geom_bin2d(bins=100)+scale_fill_continuous(type="viridis")+theme_bw()+scale_y_log10()+scale_x_log10()+
    labs(title="2D histogram of uTARs",y="uTAR PC loading",x="uTAR pseudo-bulk read coverage")
  stdDevP<-stdDevP+
    theme(plot.title = element_text(size=10,face="bold"),
          axis.title.x = element_text(size=10,face="bold"),
          axis.title.y = element_text(size=10,face="bold"),
          axis.text.x=element_text(size = 8,face="bold"),
          axis.text.y=element_text(size=8,face="bold"))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  #out<-grid.draw(VennDay4)
  return(stdDevP)
}

# taken from TAR expression matrix
drawScatterTrueBulk2<-function(outGeneSeurat,pca,numDashes=5,PCs=10,minBulk=100,returnData=F){
  # look at ven diagram
  numdash<-str_count(rownames(outGeneSeurat),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(outGeneSeurat),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(outGeneSeurat)[geneInd]
  inGene<-rownames(outGeneSeurat)[inInd]
  outGene<-rownames(outGeneSeurat)[outInd]
  
  # look at total expression of each tar
  bulkCounts<-as.numeric(unlist(lapply(str_split(outGene,"-"),function(x) x[length(x)-1])))
  rawMat<-data.frame(bulkCounts)
  rownames(rawMat)<-outGene
  
  sumExpUTar<-rowSums(rawMat)
  # look at PC loading of just uTARs, top 10 PCs
  eigenVec<-data.frame(outGeneSeurat[[pca]]@feature.loadings)
  pcLoadings<-sweep(eigenVec,2,outGeneSeurat[[pca]]@stdev,"*")
  sumPCLoadUTar<-sort(rowSums(abs(pcLoadings[,1:PCs])),decreasing=T)
  
  # look for common uTARs
  commonTARs<-intersect(names(sumExpUTar),names(sumPCLoadUTar))
  sumExpUTar<-sumExpUTar[commonTARs]
  sumPCLoadUTar<-sumPCLoadUTar[commonTARs]
  
  toPlot<-data.frame(cbind(sumExpUTar,sumPCLoadUTar))
  toPlotTest<-toPlot[toPlot$sumExpUTar>minBulk,]
  if(returnData==T){
    return(toPlotTest)
  }
  stdDevP<-ggplot(toPlotTest,aes(x=sumExpUTar,y=sumPCLoadUTar))+
    geom_bin2d(bins=100)+scale_fill_continuous(type="viridis")+theme_bw()+scale_y_log10()+scale_x_log10()+
    labs(title="2D histogram of uTARs",y="uTAR PC loading",x="uTAR pseudo-bulk read coverage")
  stdDevP<-stdDevP+
    theme(plot.title = element_text(size=10,face="bold"),
          axis.title.x = element_text(size=10,face="bold"),
          axis.title.y = element_text(size=10,face="bold"),
          axis.text.x=element_text(size = 8,face="bold"),
          axis.text.y=element_text(size=8,face="bold"))+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  #out<-grid.draw(VennDay4)
  return(stdDevP)
}



plotSingleSilBar<-function(meanSilVal,yscale){
  silValBar<-ggplot(meanSilVal,aes(x=Var2,y=value,fill=c("#0000B2","#640000","#FF0000")))+
    ylim(yscale[1],yscale[2])+
    geom_bar(position="dodge",stat="identity",width = 0.5)+
    xlab("sample")+ylab("silhouette coefficient")+labs(fill="features")+
    scale_fill_manual(values=c("#0000B2","#FF0000","#640000"))
  #scale_fill_brewer(palette = "Set1")+
  #theme(axis.title.x = element_text(size=10,face="bold"),
  #      axis.title.y = element_text(size=10,face="bold"),
  #      axis.text.x=element_text(size = 8,angle=45,hjust=1,face="bold"),
  #      axis.text.y=element_text(size=8,face="bold"),
  #      legend.position = "top",
  #      legend.text = element_text(size=8,face="bold"),
  #      legend.title=element_blank(),
  #      panel.background = element_blank(),
  #      axis.line = element_line(colour = "black"))+
  silValBar<-adjustThemeGG(silValBar)
  silValBar<-addSmallLegend(silValBar)+NoLegend()+theme(axis.title.x=element_blank(),
                                                        axis.text.x=element_text(size=8),
                                                        axis.title.y=element_blank())+
    geom_hline(yintercept=0)
  return(silValBar)
}

# from 3 column data, each column is a set of silhouette values
plotSingleSilBar2<-function(silVals,yscale){
  dataToPlot<-data.frame(cbind(c("genes","aTARs","uTARs"),as.numeric(silVals)))
  colnames(dataToPlot)<-c("feature","silhouette")
  dataToPlot$silhouette<-as.numeric(as.character(dataToPlot$silhouette))
  dataToPlot$feature<-factor(dataToPlot$feature,levels=c("genes","aTARs","uTARs"))
  
  
  silValBar<-ggplot(data=dataToPlot,aes(x=dataToPlot$feature,
                                        y=dataToPlot$silhouette,
                                        fill=dataToPlot$feature))+
    ylim(yscale[1],yscale[2])+
    geom_bar(stat="identity")+
    xlab("sample")+ylab("silhouette coefficient")+
    scale_fill_manual(values=c("#0000B2","#FF0000","#640000"))
  
  silValBar<-adjustThemeGG(silValBar)
  silValBar<-addSmallLegend(silValBar)+NoLegend()+theme(axis.title.x=element_blank(),
                                                        axis.text.x=element_text(size=8),
                                                        axis.title.y=element_blank())+
    geom_hline(yintercept=0)
  return(silValBar)
}


createNewAssay<-function(object, assay_name, set_assasy_as_default = FALSE, data_matrix, slot){
  counts = new(Class = "dgCMatrix") 
  data = matrix()
  scale.data = matrix()
  var.features = vector()
  meta.features = data.frame(row.names = rownames(x = data_matrix))
  if (slot == "counts"){
      counts <- as(data_matrix, "dgCMatrix")
  }
  if (slot == "data"){
      data <- data_matrix
  }
  if (slot == "scale.data"){
      scale.data <- data_matrix
  }
  new.assay <- new(Class = 'Assay', counts = counts, data = data, scale.data = scale.data, misc = NULL)
  object[[assay_name]] <- new.assay
  if (set_assasy_as_default){
    DefaultAssay(object) <- assay_name
  }
  return(object)
}


humanRunPCAtars<-function(input, numDashes){
  numdash<-str_count(rownames(input),"-")
  outInd<-numdash>=numDashes & endsWith(rownames(input),"-0")
  inInd<-numdash>=numDashes & endsWith(rownames(input),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(input)[geneInd]
  inGene<-rownames(input)[inInd]
  outGene<-rownames(input)[outInd]
  
  input<-normVarScaleData(input,1)
  
  geneOnly_seurat <- RunPCA(object = input, features = geneOnly)
  geneOnly_seurat <- FindNeighbors(object=input,dims=1:10)
  geneOnly_seurat <- FindClusters(object=input,resolution=0.2)
  
  outGene_seurat <- RunPCA(object = input, features = outGene)
  outGene_seurat <- FindNeighbors(object=outGene_seurat,dims=1:10)
  outGene_seurat <- FindClusters(object=outGene_seurat,resolution=0.2)
  inGene_seurat <- RunPCA(object = input, features = inGene)
  inGene_seurat <- FindNeighbors(object=inGene_seurat,dims=1:10)
  inGene_seurat <- FindClusters(object=inGene_seurat,resolution=0.2)
  Idents(outGene_seurat)<-Idents(geneOnly_seurat)
  Idents(inGene_seurat)<-Idents(geneOnly_seurat)
  
  geneOnly_seurat <- RunUMAP(geneOnly_seurat,dims=1:10,check_duplicates = F)
  inGene_seurat <- RunUMAP(inGene_seurat,dims=1:10,check_duplicates = F)
  outGene_seurat <- RunUMAP(outGene_seurat,dims=1:10,check_duplicates = F)
  
  return(list(geneOnly_seurat,inGene_seurat,outGene_seurat))
  
}

featurePlotForPaper<-function(input,featureIn,orderPoints=T){
  out<-FeaturePlot(input,
                   pt.size=0.1,
                   features = featureIn ,
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

featurePlotForPaperNoAsp<-function(input,featureIn){
  out<-FeaturePlot(input,
                   pt.size=0.01,
                   features = featureIn ,
                   cols=c("light grey","dark red"),
                   order=T)+
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank(),
          plot.title = element_blank()
          #legend.text=element_text(size=12),
          #legend.position="bottom"
          )+
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  out<-addSmallLegendDot(out,spaceLegend =0.5,textSize = 8)
  return(out)
}


generateDataforMouseBrain<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 500 & nFeature_RNA < 6000)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=3 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=3 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  return(sample_combined)
}

  
generateDataforAllOrgSameParam<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 200 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  # normalize, variable feature, scale
  sample_combined <- NormalizeData(sample_combined)
  sample_combined <- FindVariableFeatures(sample_combined, selection.method = "mean.var.plot",dispersion.cutoff=c(1,Inf),mean.cutoff=c(0,Inf))
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene[outGene %in% VariableFeatures(sample_combined)])
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene[inGene %in% VariableFeatures(sample_combined)])
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=0.2)
  print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly[geneOnly %in% VariableFeatures(sample_combined)]) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  return(sample_combined)
}

generateSeuratForMouseAtlas<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)

  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  return(sample_combined)
}

generateSeuratForMouseLemur<-function(geneMat,hmmMat){

  print(geneMat)
  print(hmmMat)
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T,fill = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  d4Gene<-d4Gene[!duplicated(d4Gene$GENE),] # remove duplicate rows to account for trachea data
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  
  d4Gene <- d4Gene[,-1] # take out first column
  
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T,fill = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  Sample_HMM<-Sample_HMM[!duplicated(Sample_HMM$GENE),] # remove duplicate rows to account for trachea data
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, min.cells=1, project = "sample_combined")
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  return(sample_combined)
}

generateSeuratForMouseLemurGeneUTAR<-function(geneMat,hmmMat){
  
  print(geneMat)
  print(hmmMat)
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T,fill = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  d4Gene<-d4Gene[!duplicated(d4Gene$GENE),] # remove duplicate rows to account for trachea data
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  
  d4Gene <- d4Gene[,-1] # take out first column
  
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T,fill = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  Sample_HMM<-Sample_HMM[!duplicated(Sample_HMM$GENE),] # remove duplicate rows to account for trachea data
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, min.cells=1, project = "sample_combined")
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  TARSeu<-sample_combined[c(inGene,outGene),]
  aTARSeu<-sample_combined[c(inGene),]
  uTARSeu<-sample_combined[c(outGene),]
  geneSeu<-sample_combined[c(geneOnly),]
  sample_combined<-sample_combined[c(geneOnly,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(TARSeu,features = outGene)
  sample_combined[["nFeature_RNA"]]<-geneSeu[["nFeature_RNA"]]
  sample_combined[["nCount_RNA"]]<-geneSeu[["nCount_RNA"]]
  sample_combined[["nFeature_aTAR"]]<-aTARSeu[["nFeature_RNA"]]
  sample_combined[["nFeature_uTAR"]]<-uTARSeu[["nFeature_RNA"]]
  sample_combined[["nCount_aTAR"]]<-aTARSeu[["nCount_RNA"]]
  sample_combined[['nCount_uTAR']]<-uTARSeu[["nCount_RNA"]]
  
  return(sample_combined)
}

generateSeuratForMouseLemur2<-function(geneMat,hmmMat){
  
  print(geneMat)
  print(hmmMat)
  # load gene expression matrix and find valid cells
  if (grepl("TRACHEA",geneMat)){
    d4Gene=read.delim(geneMat)
  } else {
    d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T,fill = T)
  }
  
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  d4Gene<-d4Gene[!duplicated(d4Gene$GENE),] # remove duplicate rows to account for trachea data
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  
  d4Gene <- d4Gene[,-1] # take out first column
  
  if (grepl("TRACHEA",hmmMat)){
    Sample_HMM=read.delim(hmmMat)
  } else {
    Sample_HMM=fread(paste0(hmmMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T,fill = T)
  }
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  Sample_HMM<-Sample_HMM[!duplicated(Sample_HMM$GENE),] # remove duplicate rows to account for trachea data
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, min.cells=1, project = "sample_combined")
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  return(sample_combined)
}

generateMatForMouseLemur<-function(hmmMat){
  library(data.table);library(stringr)
  if (grepl("TRACHEA",hmmMat)){
    Sample_HMM=read.delim(hmmMat)
  } else {
    Sample_HMM=fread(paste0(hmmMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T,fill = T)
  }
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  Sample_HMM<-Sample_HMM[!duplicated(Sample_HMM$GENE),] # remove duplicate rows to account for trachea data
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  sample_combined <- Sample_HMM[,-1] # take out first column
  
  numdash<-str_count(rownames(sample_combined),"_")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"_0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"_1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  return(sample_combined[outGene,])
}



generateSeuratForMouseLemurTAR<-function(hmmMat){
  
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T,fill = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  sample_combined_mat<-Sample_HMM
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, min.cells = 1,project = "sample_combined")#,fill = T)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  return(sample_combined[outGene,])
}

generateSeuratForMouseLemurTAR2<-function(hmmMat){
  
  Sample_HMM=read.delim(hmmMat)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  sample_combined_mat<-Sample_HMM
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, min.cells = 1,project = "sample_combined")#,fill = T)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  return(sample_combined[outGene,])
}

generateSeuratForMouseLemurTAR3<-function(hmmMat){
  
  Sample_HMM=read.delim(hmmMat)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  sample_combined_mat<-Sample_HMM
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, min.cells = 1,project = "sample_combined")#,fill = T)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  return(sample_combined)
}

generateSeuratForMouseLemurGene<-function(hmmMat){
  
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  
  sample_combined_mat<-Sample_HMM
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, min.cells = 1,project = "sample_combined")
  
  return(sample_combined)
}


calcUtarStats<-function(input){
  numdash<-str_count(rownames(input),"-")
  outInd<-numdash>=5 & endsWith(rownames(input),"-0")
  inInd<-numdash>=5 & endsWith(rownames(input),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(input)[geneInd]
  inGene<-rownames(input)[inInd]
  outGene<-rownames(input)[outInd]
  
  mat<-GetAssayData(object = input, slot = "counts")
  
  return(c(length(inGene),length(outGene),sum(mat[inGene,]),sum(mat[outGene,])))
  
}


generateDataforAllOrgSameParamAllTARs<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 400 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  # normalize, variable feature, scale
  sample_combined <- NormalizeData(sample_combined)
  sample_combined <- FindVariableFeatures(sample_combined, selection.method = "mean.var.plot",dispersion.cutoff=c(0.1,Inf),mean.cutoff=c(0,Inf))
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=0.2)
  print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)

  return(sample_combined)
}

runAnalysisOnMouse<-function(sample_combined){
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  # normalize, variable feature, scale
  sample_combined <- NormalizeData(sample_combined)
  #sample_combined <- FindVariableFeatures(sample_combined, selection.method = "mean.var.plot",dispersion.cutoff=c(0.1,Inf),mean.cutoff=c(0,Inf))
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=0.2)
  print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  
  return(sample_combined)
}

runAnalysisOnLemur<-function(sample_combined){
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  # normalize, variable feature, scale
  sample_combined <- NormalizeData(sample_combined)
  #sample_combined <- FindVariableFeatures(sample_combined, selection.method = "mean.var.plot",dispersion.cutoff=c(0.1,Inf),mean.cutoff=c(0,Inf))
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  #sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  #sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  #sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene)
  #sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  #sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=0.2)
  #print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10)
  #sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  
  return(sample_combined)
}

runAnalysisOnLemur_geneUTar<-function(sample_combined){
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  sample_combined<-sample_combined[c(geneOnly,outGene),]
  # normalize, variable feature, scale
  sample_combined <- NormalizeData(sample_combined)
  #sample_combined <- FindVariableFeatures(sample_combined, selection.method = "mean.var.plot",dispersion.cutoff=c(0.1,Inf),mean.cutoff=c(0,Inf))
  sample_combined <- ScaleData(sample_combined, features = c(geneOnly,outGene))
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  #sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  #sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  #print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  #sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene)
  #sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  #sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=0.2)
  #print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  #sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  #sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  #print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:50)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:50)
  #sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  
  return(sample_combined)
}

runAnalysisOnLemurGene<-function(sample_combined){
  # normalize, variable feature, scale
  sample_combined <- NormalizeData(sample_combined)
  #sample_combined <- FindVariableFeatures(sample_combined, selection.method = "mean.var.plot",dispersion.cutoff=c(0.1,Inf),mean.cutoff=c(0,Inf))
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10)
  #sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  
  return(sample_combined)
}

runAnalysisOnLemurTAR<-function(sample_combined){
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  # normalize, variable feature, scale
  sample_combined <- NormalizeData(sample_combined)
  #sample_combined <- FindVariableFeatures(sample_combined, selection.method = "mean.var.plot",dispersion.cutoff=c(0.1,Inf),mean.cutoff=c(0,Inf))
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10)
  return(sample_combined)
}


convertWinstonChr<-function(input){
  
}

createSeuratUrchin<-function(geneMat,hmmMat){

  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 400 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  return(sample_combined)
}

createSeuratTestis<-function(geneMat,hmmMat){
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 400 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  return(sample_combined)
}

createSeuratMoleRat<-function(geneMat,hmmMat){
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 400 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  return(sample_combined)
}

createSeuratLung<-function(geneMat,hmmMat){

  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly <- subset(x = geneOnly, nFeature_RNA > 200 & nFeature_RNA < 2500)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM[,colnames(d4Gene)])
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "project_name", min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  return(sample_combined)
}

createSeurathg38<-function(geneMat,hmmMat){
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 3, min.features = 200)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(geneOnly, pattern = "^MT-")
  geneOnly <- subset(geneOnly, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  return(sample_combined)
}

createSeuarthg16<-function(geneMat,hmmMat,cellsFromHg38){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 3, min.features = 200)
  #geneOnly[["percent.mt"]] <- PercentageFeatureSet(geneOnly, pattern = "^MT-")
  geneOnly <- subset(geneOnly, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)# & percent.mt < 5)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2)
  sample_combined<-subset(x = sample_combined, cells=cellsFromHg38)
  return(sample_combined)
}

createSeuratChicken<-function(geneMat,hmmMat){
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  return(sample_combined)
}



generateDataforHumansHg38<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 3, min.features = 200)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(geneOnly, pattern = "^MT-")
  geneOnly <- subset(geneOnly, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  sample_combined <- NormalizeData(sample_combined)
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=0.2)
  print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  return(sample_combined)
  
}

createFig2<-function(inputSeurat,maxVal=20,addPoints=F){
  geneUmap<-umapPlainLegend(inputSeurat,"umap_gene")
  uTARUmap<-umapPlainLegend(inputSeurat,"umap_uTAR")
  aTARUmap<-umapPlainLegend(inputSeurat,"umap_aTAR")
  vlnP<-plotVlnPlotTypeGG(inputSeurat,c(0,maxVal),addPoints=addPoints)
  return(list(geneUmap,uTARUmap,aTARUmap,vlnP))
}

analysisOnFilSeurat<-function(sample_combined,res_param=0.2){
  library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
  library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
  library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
  #source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
  library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  sample_combined <- NormalizeData(sample_combined)
  sample_combined <- ScaleData(sample_combined, features = rownames(sample_combined))
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=res_param)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=res_param)
  print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=res_param)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  return(sample_combined)
}

analysisOnFilSeuratWithVarGene<-function(sample_combined,res_param=0.2){
  library(stringr); library(Seurat, lib.loc = "/programs/R-3.5.0/library")
  #library(ggplot2); library(data.table); library(caTools);library(scales); library(Seurat, lib.loc = "/programs/R-3.5.0/library"); library(dplyr)
  #library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
  #library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
  #source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
  #library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  normVarScaleData<-function(seuratOb,numFeat=2000){
    seuratOb <- NormalizeData(seuratOb)
    seuratOb <- FindVariableFeatures(seuratOb, selection.method = "vst", nfeatures=numFeat)
    seuratOb <- ScaleData(seuratOb, features = rownames(seuratOb))
    return(seuratOb)
  }
  
  sample_combined <- normVarScaleData(sample_combined,10000)

  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, npcs = min(50,ncol(sample_combined)), reduction.name="pca_uTAR",features = outGene[outGene %in% VariableFeatures(sample_combined)])
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=res_param)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, npcs = min(50,ncol(sample_combined)), reduction.name="pca_aTAR",features = inGene[inGene %in% VariableFeatures(sample_combined)])
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=res_param)
  print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, npcs = min(50,ncol(sample_combined)), features = geneOnly[geneOnly %in% VariableFeatures(sample_combined)]) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=res_param)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunPCA(object = sample_combined, npcs = min(50,ncol(sample_combined)), reduction.name="pca_nonVarGenes", features = geneOnly[!(geneOnly %in% VariableFeatures(sample_combined))]) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,reduction="pca_nonVarGenes",graph.name="nonVarGenes_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="nonVarGenes_snn",resolution=res_param)
  print("Finished running PCA on NON variable genes and clustered based on NON variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_nonVarGenes",reduction="pca_nonVarGenes",dims=1:10,check_duplicates = F)
  return(sample_combined)
}



generateDataforChicken3<-function(geneMat,hmmMat){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 2, min.features = 100)
  geneOnly[["percent.mt"]] <- PercentageFeatureSet(object = geneOnly, pattern = "^MT-")
  geneOnly <- subset(x = geneOnly, subset = percent.mt < 20 & nFeature_RNA > 200)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2, min.features = 100)
  sample_combined<-subset(x = sample_combined, cells=d4Cells)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  sample_combined<-normVarScaleData(sample_combined,1)
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=0.2)
  print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  return(sample_combined)
  
}


generateDataforHumansHg16<-function(geneMat,hmmMat,cellsFromHg38){
  
  # load gene expression matrix and find valid cells
  d4Gene=fread(paste0(geneMat),sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  d4Gene <- as.data.frame(d4Gene) # force to data frame
  rownames(d4Gene) <- d4Gene$GENE # make the name of rows GENEs
  d4Gene <- d4Gene[,-1] # take out first column
  # incorporate mito info
  #mito.genes<-c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
  #rownames(d4Gene)[(rownames(d4Gene) %in% mito.genes)]<-paste0("MT-",rownames(d4Gene)[rownames(d4Gene) %in% mito.genes]) #incorporate MT genes
  geneOnly<-CreateSeuratObject(counts = d4Gene, project = "d4Gene",min.cells = 3, min.features = 200)
  #geneOnly[["percent.mt"]] <- PercentageFeatureSet(geneOnly, pattern = "^MT-")
  geneOnly <- subset(geneOnly, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)# & percent.mt < 5)
  d4Cells<-colnames(geneOnly)
  
  # load hmm matrix and combine with gene expression matrix and subset for valid cells
  Sample_HMM=fread(hmmMat,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
  rownames(Sample_HMM)<-tolower(rownames(Sample_HMM))
  Sample_HMM <- as.data.frame(Sample_HMM) # force to data frame
  rownames(Sample_HMM) <- Sample_HMM$GENE # make the name of rows GENEs
  Sample_HMM <- Sample_HMM[,-1] # take out first column
  sample_combined_mat<-rbind(d4Gene,Sample_HMM)
  sample_combined<-CreateSeuratObject(counts = sample_combined_mat, project = "sample_combined",min.cells = 2)
  sample_combined<-subset(x = sample_combined, cells=cellsFromHg38)
  
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  out<-sample_combined[c(inGene,outGene),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(out,features = outGene)
  
  sample_combined<-normVarScaleData(sample_combined,1)
  
  print(paste(length(VariableFeatures(sample_combined))," variable features in total, genes and TARs included."))
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_uTAR",features = outGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_uTAR", graph.name="uTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="uTAR_snn",resolution=0.2)
  print("Finished running PCA on variable uTARs and clustered based on variable uTARs.")
  sample_combined <- RunPCA(object = sample_combined, reduction.name="pca_aTAR",features = inGene)
  sample_combined <- FindNeighbors(object=sample_combined, reduction="pca_aTAR", graph.name="aTAR_snn",dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,graph.name="aTAR_snn",resolution=0.2)
  print("Finished running PCA on variable aTARs and clustered based on variable aTARs.")
  sample_combined <- RunPCA(object = sample_combined, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  sample_combined <- FindNeighbors(object=sample_combined,dims=1:10)
  sample_combined <- FindClusters(object=sample_combined,resolution=0.2)
  print("Finished running PCA on variable genes and clustered based on variable genes.")
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_gene",reduction="pca",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:10,check_duplicates = F)
  sample_combined <- RunUMAP(sample_combined,reduction.name="umap_aTAR",reduction="pca_aTAR",dims=1:10,check_duplicates = F)
  return(sample_combined)
  
}






plotDiffExp<-function(inputSeurat,uTAR,bamFile,titleForPlot="uTAR",normalChrom=T){
  featureP<-featurePlotForPaperNoAsp(inputSeurat,uTAR)
  
  if (normalChrom){
    chr<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    fasta<-paste0(chr,":",start,"-",stop)
  } else {
    chr1<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    chr2<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",4))
    fasta<-paste0(chr1,"_",chr2,":",start,"-",stop)
  }
  covCom<-paste0("samtools depth -r ",fasta," ",bamFile," > temp.txt")
  test<-system(covCom)
  temp <- read.delim("temp.txt", header=FALSE)
  system("rm temp.txt")
  coverage<-list(temp[,3])
  
  allCovPlot<-plotListOfCov2(unlist(coverage))
  allCovPlot<-allCovPlot+labs(title=titleForPlot,
                              x="position along uTAR",
                              y="coverage")+theme(axis.title.y=element_text(angle=90))
  allCovPlot<-adjustThemeGGSize(allCovPlot,10,10,10,8,8)
  
  return(plot_grid(allCovPlot,featureP,ncol=1,align = "v",axis="lr",rel_heights = c(1,4)))
  
}


plotDiffExpNolab<-function(inputSeurat,uTAR,bamFile,titleForPlot="uTAR",normalChrom=T,order=T){
  featureP<-featurePlotForPaper(inputSeurat,uTAR,order)+NoLegend()
  
  if (normalChrom){
    chr<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    fasta<-paste0(chr,":",start,"-",stop)
  } else {
    chr1<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    chr2<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",4))
    fasta<-paste0(chr1,"_",chr2,":",start,"-",stop)
  }
  covCom<-paste0("samtools depth -r ",fasta," ",bamFile," > temp.txt")
  test<-system(covCom)
  temp <- read.delim("temp.txt", header=FALSE)
  system("rm temp.txt")
  coverage<-list(temp[,3])
  
  allCovPlot<-plotListOfCov2(unlist(coverage))
  allCovPlot<-allCovPlot+labs(title=titleForPlot)+
    theme(plot.title = element_text(size=8,face = "bold"),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6))

  #out<-plot_grid(allCovPlot,featureP,ncol=1,align = "v",axis="lr",rel_heights = c(2,3))
  return(list(allCovPlot,featureP))
}

plotDiffExpNolabGene<-function(inputSeurat,uTAR,bamFile,titleForPlot="uTAR",normalChrom=T,order=T){
  featureP<-featurePlotForPaper(inputSeurat,titleForPlot,order)+NoLegend()
  
  if (normalChrom){
    chr<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    fasta<-paste0(chr,":",start,"-",stop)
  } else {
    chr1<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    chr2<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",4))
    fasta<-paste0(chr1,"_",chr2,":",start,"-",stop)
  }
  covCom<-paste0("samtools depth -r ",fasta," ",bamFile," > temp.txt")
  test<-system(covCom)
  temp <- read.delim("temp.txt", header=FALSE)
  system("rm temp.txt")
  coverage<-list(temp[,3])
  
  allCovPlot<-plotListOfCov2(unlist(coverage))
  allCovPlot<-allCovPlot+labs(title=titleForPlot)+
    theme(plot.title = element_text(size=8,face = "bold"),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6))

  #out<-plot_grid(allCovPlot,featureP,ncol=1,align = "v",axis="lr",rel_heights = c(2,3))
  return(list(allCovPlot,featureP))
  
}


findSiguTARmarkers<-function(seuratOb,moleRatMarkers,shift=0,fcThresh=0.5,minpc.1=0.5){
  numdash<-str_count(rownames(moleRatMarkers),"-")
  outInd<-numdash>=5 & endsWith(rownames(moleRatMarkers),"-0")
  inInd<-numdash>=5 & endsWith(rownames(moleRatMarkers),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(moleRatMarkers)[geneInd]
  inGene<-rownames(moleRatMarkers)[inInd]
  outGene<-rownames(moleRatMarkers)[outInd]
  
  #filter for uTARs
  markers<-moleRatMarkers[outGene,]
  #filter by FC
  markers<-markers[markers$avg_logFC>=fcThresh,]
  #filter by minimum group1 marker
  markers<-markers[markers$pct.1>=minpc.1,]
  
  dotPlotOutF<-plotDotPlotWithLabel(seuratOb,markers)
  
  #dotPlotOutF<-DotPlot(seuratOb, features = rownames(markers),dot.scale = 4)+
  #  coord_flip()+RotatedAxis()+
  #  labs(x="uTAR features")+
  #  theme(legend.position="right",
  #        legend.direction = "vertical",
  #        legend.text=element_text(size=4),
  #        legend.title = element_text(size=4),
  #        legend.spacing.x=unit(1,"mm"),
  #        legend.spacing.y=unit(1,"mm"),
  #        axis.text.y=element_blank(),#element_blank(size=4),
  #        axis.text.x=element_text(size=8,
  #                                 angle = 0,
  #                                 colour=as.vector(alphabet(n=length(unique(Idents(seuratOb)))))),
  #        axis.title.y=element_text(size=10, face="bold"),
  #        axis.title.x=element_blank())
  #dotPlotOutF<-addSmallLegendDot(dotPlotOutF,spaceLegend =0.5,textSize = 8)+ #+theme(aspect.ratio = 1)
  #  theme(legend.title = element_blank())
  
  start<-as.numeric(sapply(strsplit(rownames(markers),"-"),"[",2+shift))
  stop<-as.numeric(sapply(strsplit(rownames(markers),"-"),"[",3+shift))
  
  return(list(dotPlotOutF,markers,stop-start))
}

findSigGenemarkers<-function(seuratOb,moleRatMarkers,shift=0,fcThresh=0.5,minpc.1=0.5){
  numdash<-str_count(rownames(moleRatMarkers),"-")
  outInd<-numdash>=5 & endsWith(rownames(moleRatMarkers),"-0")
  inInd<-numdash>=5 & endsWith(rownames(moleRatMarkers),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(moleRatMarkers)[geneInd]
  inGene<-rownames(moleRatMarkers)[inInd]
  outGene<-rownames(moleRatMarkers)[outInd]
  
  #filter for uTARs
  markers<-moleRatMarkers[geneOnly,]
  #filter by FC
  markers<-markers[markers$avg_logFC>=fcThresh,]
  #filter by minimum group1 marker
  markers<-markers[markers$pct.1>=minpc.1,]
  
  dotPlotOutF<-plotDotPlotWithLabel(seuratOb,markers)

  start<-as.numeric(sapply(strsplit(rownames(markers),"-"),"[",2+shift))
  stop<-as.numeric(sapply(strsplit(rownames(markers),"-"),"[",3+shift))
  
  return(list(dotPlotOutF,markers,stop-start))
}

  
findSiguTARmarkers2<-function(seuratOb,moleRatMarkers,shift=0,fcThresh=0.5,minpc.1=0.5){
  numdash<-str_count(rownames(moleRatMarkers),"-")
  outInd<-numdash>=5 & endsWith(rownames(moleRatMarkers),"-0")
  inInd<-numdash>=5 & endsWith(rownames(moleRatMarkers),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(moleRatMarkers)[geneInd]
  inGene<-rownames(moleRatMarkers)[inInd]
  outGene<-rownames(moleRatMarkers)[outInd]
  
  #filter for uTARs
  markers<-moleRatMarkers[outGene,]
  #filter by FC
  markers<-markers[markers$avg_logFC>=fcThresh,]
  #filter by minimum group1 marker
  markers<-markers[markers$pct.1>=minpc.1,]
  
  return(markers)
}

findSigGenemarkers2<-function(seuratOb,moleRatMarkers,shift=0,fcThresh=0.5,minpc.1=0.5){
  numdash<-str_count(rownames(moleRatMarkers),"-")
  outInd<-numdash>=5 & endsWith(rownames(moleRatMarkers),"-0")
  inInd<-numdash>=5 & endsWith(rownames(moleRatMarkers),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(moleRatMarkers)[geneInd]
  inGene<-rownames(moleRatMarkers)[inInd]
  outGene<-rownames(moleRatMarkers)[outInd]
  #filter for uTARs
  markers<-moleRatMarkers[geneOnly,]
  #filter by FC
  markers<-markers[markers$avg_logFC>=fcThresh,]
  #filter by minimum group1 marker
  markers<-markers[markers$pct.1>=minpc.1,]
  return(markers)
}

plotDotPlotWithLabel<-function(seuratOb,markers,dotScale=1.5){
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
          axis.text.y=element_text(size=4),
          axis.text.x=element_text(size=6,
                                   hjust=0.5,
                                   angle=0,
                                   face = "bold",
                                   colour=as.vector(alphabet(n=length(unique(Idents(seuratOb)))))),
          axis.title.y=element_blank(),#element_text(size=10, face="bold"),
          axis.title.x=element_blank())
  dotPlotOutF<-addSmallLegendDot(dotPlotOutF,spaceLegend =0.5,textSize = 8)+ #+theme(aspect.ratio = 1)
    theme(legend.title = element_blank())
  
  return(dotPlotOutF)
  
}

plotDotPlotWithLabelLemur<-function(seuratOb,metaData,markers,dotScale=1.5){
  dotPlotOutF<-DotPlot(seuratOb, 
                       features = markers$gene,
                       dot.scale = dotScale,
                       cols = c("white","brown"))+
    coord_flip()+RotatedAxis()+
    scale_x_discrete(breaks=markers$gene,
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
                                   angle=45,
                                   face = "bold",
                                   colour=as.vector(alphabet(n=length(levels(seuratOb$cellType))))),
          axis.title.y=element_blank(),#element_text(size=10, face="bold"),
          axis.title.x=element_blank())
  dotPlotOutF<-addSmallLegendDot(dotPlotOutF,spaceLegend =0.5,textSize = 6)+ #+theme(aspect.ratio = 1)
    theme(legend.title = element_blank())
  
  return(dotPlotOutF)
}

plotDotPlotWithLabelLemurFlip<-function(seuratOb,metaData,markers,dotScale=1.5){
  dotPlotOutF<-DotPlot(seuratOb, 
                       features = markers$gene,
                       dot.scale = dotScale,
                       cols = c("white","brown"))+
    scale_x_discrete(breaks=markers$gene,
                     labels=markers$blastnResults)+
    labs(x="uTAR features")+
    theme(legend.position="right",
          legend.direction = "vertical",
          legend.text=element_text(size=6),
          legend.title = element_text(size=6),
          legend.spacing.x=unit(1,"mm"),
          legend.spacing.y=unit(1,"mm"),
          axis.text.y=element_text(size=6,
                                   face = "bold",
                                   colour=as.vector(alphabet(n=length(levels(seuratOb$cellType))))),
          axis.text.x=element_text(size=6,
                                   hjust=1,
                                   angle=45),
          axis.title.y=element_blank(),#element_text(size=10, face="bold"),
          axis.title.x=element_blank())
  dotPlotOutF<-addSmallLegendDot(dotPlotOutF,spaceLegend =0.5,textSize = 6)+ #+theme(aspect.ratio = 1)
    theme(legend.title = element_blank())
  
  return(dotPlotOutF)
}


plotDotPlotWithLabelSupp<-function(seuratOb,markers,dotScale=2){
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

getSeqCov<-function(uTAR,bamFile,fastaFile,titleForPlot="uTAR",normalChrom=T){
  if (normalChrom){
    chr<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    fasta<-paste0(chr,":",start,"-",stop)
  } else {
    chr1<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    chr2<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",4))
    fasta<-paste0(chr1,"_",chr2,":",start,"-",stop)
  }
  covCom<-paste0("samtools depth -r ",fasta," ",bamFile," > temp.txt")
  test<-system(covCom)
  temp <- read.delim("temp.txt", header=FALSE)
  system("rm temp.txt")
  coverage<-list(temp[,3])
  
  seqCom<-paste0("samtools faidx -r ",fasta," ",fastaFile," > temp.txt")
  test<-system(seqCom)
  temp <- read.delim("temp.txt", header=FALSE)
  system("rm temp.txt")
  coverage<-list(temp[,3])
  
  
  allCovPlot<-plotListOfCov2(unlist(coverage))
  allCovPlot<-allCovPlot+labs(title=titleForPlot,
                              x="position along uTAR",
                              y="coverage")+theme(axis.title.y=element_text(angle=90))
  allCovPlot<-adjustThemeGGSize(allCovPlot,10,10,10,8,8)
  
  return(plot_grid(allCovPlot,featureP,ncol=1,align = "vl",rel_heights = c(1,4)))
  
}

calcBulk<-function(bulkFile,SeuObj){
  TARBulkCounts <- read.delim(bulkFile, header=FALSE)
  TARBulkCounts[,1]<-str_replace_all(TARBulkCounts[,1],"_","-")
  numdash<-str_count(TARBulkCounts[,1],"-")
  outInd<-numdash>=5 & endsWith(TARBulkCounts[,1],"-0")
  inInd<-numdash>=5 & endsWith(TARBulkCounts[,1],"-1")
  geneInd<-!(outInd|inInd)
  geneOnlyBulk<-TARBulkCounts[,1][geneInd]
  inGeneBulk<-TARBulkCounts[,1][inInd]
  outGeneBulk<-TARBulkCounts[,1][outInd]
  rownames(TARBulkCounts)<-TARBulkCounts[,1]
  uTARBulkCounts<-TARBulkCounts[outGeneBulk,]
  uTARBulkCountsLemur<-data.frame(uTARBulkCounts[,-1],row.names = outGeneBulk)
  testisScatter<-drawScatterTrueBulk(SeuObj,pca="pca_uTAR",rawMat = uTARBulkCountsLemur,PCs = 5,numDashes = 5,minBulk = 0,returnData = T)
  
  return(testisScatter)
}

FindDiffExp<-function(input){
  mouseMarkers <- FindAllMarkers(input, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
  
  numdash<-str_count(mouseMarkers$gene,"-")
  outInd<-numdash>=5 & endsWith(mouseMarkers$gene,"-0")
  inInd<-numdash>=5 & endsWith(mouseMarkers$gene,"-1")
  geneInd<-!(outInd|inInd)
  markersUTAR<-mouseMarkers[outInd,]

  return(mouseMarkers)
  
}

FindDiffExpUTAR<-function(mouseMarkers){

  numdash<-str_count(mouseMarkers$gene,"-")
  outInd<-numdash>=5 & endsWith(mouseMarkers$gene,"-0")
  inInd<-numdash>=5 & endsWith(mouseMarkers$gene,"-1")
  geneInd<-!(outInd|inInd)
  markersUTAR<-mouseMarkers[outInd,]
  
  return(markersUTAR)
  
}

checkIfExistGene_noDir_num<-function(input,gene_ref){
  chrom<-input[[2]]
  startPos<-input[[4]]
  endPos<-input[[5]]
  gene_ref<-gene_ref # bed of entire genes
  
  # 1: either partially in or completely inside gene
  chrMatch_gene<-(chrom==gene_ref$chr)
  qinrstartMatch<-(startPos>=gene_ref$start)&(startPos<=gene_ref$end)
  qinrendMatch<-(endPos<=gene_ref$end)&(endPos>=gene_ref$start)
  qinroutGeneAll_gene<-(chrMatch_gene*qinrstartMatch*qinrendMatch) # query completely inside ref 
  rinqstartMatch<-(startPos<(gene_ref$start))
  rinqendMatch<-(endPos>(gene_ref$end))
  rinqoutGeneAll<-(chrMatch_gene*rinqstartMatch*rinqendMatch) # ref completely inside query
  # also check partial overlap from beginning
  partialStart<-(chrMatch_gene*rinqstartMatch*qinrendMatch)
  partialEnd<-(chrMatch_gene*qinrstartMatch*rinqendMatch)
  if(sum(qinroutGeneAll_gene)|sum(rinqoutGeneAll)|sum(partialStart)|sum(partialEnd)){
    # g's are for gene names
    g1<-as.character(gene_ref$gene[as.logical(qinroutGeneAll_gene)])
    g2<-as.character(gene_ref$gene[as.logical(rinqoutGeneAll)])
    g3<-as.character(gene_ref$gene[as.logical(partialStart)])
    g4<-as.character(gene_ref$gene[as.logical(partialEnd)])
    gAll<-unique(c(g1,g2,g3,g4))
    return(1)
  }else{
    return(0)
  }
}

calcPerCellStats<-function(sample_combined){
  ####//////// setting the statistics
  numdash<-str_count(rownames(sample_combined),"-")
  outInd<-numdash>=5 & endsWith(rownames(sample_combined),"-0")
  inInd<-numdash>=5 & endsWith(rownames(sample_combined),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(sample_combined)[geneInd]
  inGene<-rownames(sample_combined)[inInd]
  outGene<-rownames(sample_combined)[outInd]
  
  TARSeu<-sample_combined[c(inGene,outGene),]
  aTARSeu<-sample_combined[c(inGene),]
  uTARSeu<-sample_combined[c(outGene),]
  geneSeu<-sample_combined[c(geneOnly),]
  sample_combined[["percent.outHMM"]]<-PercentageFeatureSet(TARSeu,features = outGene)
  sample_combined[["nFeature_RNA"]]<-geneSeu[["nFeature_RNA"]]
  sample_combined[["nCount_RNA"]]<-geneSeu[["nCount_RNA"]]
  sample_combined[["nFeature_aTAR"]]<-aTARSeu[["nFeature_RNA"]]
  sample_combined[["nFeature_uTAR"]]<-uTARSeu[["nFeature_RNA"]]
  sample_combined[["nCount_aTAR"]]<-aTARSeu[["nCount_RNA"]]
  sample_combined[['nCount_uTAR']]<-uTARSeu[["nCount_RNA"]]
  
  return(sample_combined)
}

plot3D<-function(mergedSeurat,embedding="umap_uTAR_3D"){
  umap_1 <- mergedSeurat[[embedding]]@cell.embeddings[,1]
  umap_2 <- mergedSeurat[[embedding]]@cell.embeddings[,2]
  umap_3 <- mergedSeurat[[embedding]]@cell.embeddings[,3]
  cellIdents<-mergedSeurat$tissue
  
  dataToPlot<-data.frame(cbind(as.character(cellIdents),umap_1,umap_2,umap_3))
  plot_ly(data = dataToPlot, 
          x = ~umap_1, y = ~umap_2, z = ~umap_3, 
          color = ~V1,
          type = "scatter3d", 
          mode = "markers")
          
}

removeAxisLabel<-function(input){
  input<-input+
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()
          )
  return(input)
}

removeAxisLabelBox<-function(input){
  input<-input+
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          #panel.grid = element_blank(),
          axis.ticks = element_blank()
          #axis.line = element_blank()
    )
  return(input)
}


getGeneNameFromTARs<-function(tarNames){
  i <- 6
  j <- 6
  tarsMod<-sapply(strsplit(tarNames, "_"), function(x) {
    g <- seq_along(x)
    g[g < i] <- i
    g[g > j + 1] <- j+1
    paste(tapply(x, g, paste, collapse = "_"), collapse = "!MW!")
  })
  ListOfGenes<-getNFromList(strsplit(tarsMod,"!MW!"),2)
  ListOfGenes<-unique(unlist(strsplit(ListOfGenes,"_")))
  ListOfGenes<-ListOfGenes[!(ListOfGenes=="0"|ListOfGenes=="1")]
  return(ListOfGenes)
}

checkTARinWinston<-function(input,gene_ref_bare){
  outGene<-input
  outGene<-str_replace_all(outGene,pattern = "_",replacement = "-")
  outGene<-outGene[startsWith(outGene,"NC-0336")]
  length(outGene)
  
  HMManno_bare<-data.frame(nrow=length(outGene))
  chromCode<-paste(getNFromList(strsplit(outGene,"-"),1),getNFromList(strsplit(outGene,"-"),2),sep = "_")
  
  chrVar<-chromCode
  startVar<-getNFromList(strsplit(outGene,"-"),3)
  endVar<-getNFromList(strsplit(outGene,"-"),4)
  dirVar<-getNFromList(strsplit(outGene,"-"),5)
  dirVar[dirVar==""]<-"-"
  HMManno_bare<-data.frame(cbind(chrVar,startVar,endVar,dirVar))
  colnames(HMManno_bare)<-c("chr","start","end","direction")
  HMManno_bare$chr<-as.character(HMManno_bare$chr)
  HMManno_bare$direction<-as.character(HMManno_bare$direction)
  HMManno_bare$start<-as.numeric(as.character(HMManno_bare$start))
  HMManno_bare$end<-as.numeric(as.character(HMManno_bare$end))
  HMManno_bare$chr[is.na(HMManno_bare$chr)]<-"unannoChr"
  HMManno_bare$TAR<-outGene
  HMManno_bare$noDir<-apply(X=HMManno_bare,MARGIN=1,FUN=checkIfExistGene_noDir,gene_ref=gene_ref_bare)
  #subset outGene that match with winston annotations
  outGeneWinston<-data.frame(outGene[HMManno_bare$noDir!=0])
  HMManno_bareWinston<-HMManno_bare[HMManno_bare$noDir!=0,]
  winsTransInUTAR<-unique(unlist(strsplit(as.character(HMManno_bareWinston$noDir),"; _")))
  winsTransInUTAR<-winsTransInUTAR[!winsTransInUTAR=="1"]
  
  return(winsTransInUTAR)
}

getTARsFilterlogFC<-function(allTARsInput,thresFC){
  out<-list()
  for(logFCThres in thresFC){
    temp<-gsub(pattern = "-",x = allTARsInput$gene[allTARsInput$avg_logFC>logFCThres],replacement = "_")
    tarNames<-gsub(pattern="___",x=temp,replacement = "_-_")
    genesInDiffTar<-getGeneNameFromTARs(tarNames)
    out<-c(out,list(genesInDiffTar))
  }
  return(out)
}

getTARsFilterPVal<-function(allTARsInput,thresFC){
  out<-list()
  for(logFCThres in thresFC){
    temp<-gsub(pattern = "-",x = allTARsInput$gene[allTARsInput$p_val<logFCThres],replacement = "_")
    tarNames<-gsub(pattern="___",x=temp,replacement = "_-_")
    genesInDiffTar<-getGeneNameFromTARs(tarNames)
    out<-c(out,list(genesInDiffTar))
  }
  return(out)
}


getTARsFilterlogFCTarName<-function(allTARsInput,thresFC){
  out<-list()
  for(logFCThres in thresFC){
    temp<-gsub(pattern = "-",x = allTARsInput$gene[allTARsInput$avg_logFC>logFCThres],replacement = "_")
    tarNames<-gsub(pattern="___",x=temp,replacement = "_-_")
    out<-c(out,list(unique(tarNames)))
  }
  return(out)
}

getTARsFilterPValTarName<-function(allTARsInput,thresFC){
  out<-list()
  for(logFCThres in thresFC){
    temp<-gsub(pattern = "-",x = allTARsInput$gene[allTARsInput$p_val<logFCThres],replacement = "_")
    tarNames<-gsub(pattern="___",x=temp,replacement = "_-_")
    out<-c(out,list(unique(tarNames)))
  }
  return(out)
}


sumLengthFunc<-function(genesInAllTar,features){
  return(sum(features%in%genesInAllTar))
}

cellTypeReadsFunc<-function(mergedSeuratMeta,column="percent.outHMM"){
  # calculate median percent out for each 
  cellTypes<-unique(mergedSeuratMeta$cellType)
  cellTypesOutPer<-vector()
  for(cells in cellTypes){
    cellTypesOutPer<-c(cellTypesOutPer,median(mergedSeuratMeta[,column][mergedSeuratMeta$cellType==cells]))
  }
  cellTypeInfo<-data.frame(cbind(cellTypes,cellTypesOutPer))
  cellTypeInfo$cellTypesOutPer<-as.numeric(levels(cellTypeInfo$cellTypesOutPer))[cellTypeInfo$cellTypesOutPer]
  cellTypeInfo <- cellTypeInfo[order(-cellTypeInfo$cellTypesOutPer),]
  cellTypeInfo$cellTypes<-as.character(cellTypeInfo$cellTypes)
  
  for(i in 1:nrow(cellTypeInfo)){
    cellTypeInfo$colour[i]<-unique(mergedSeuratMeta[mergedSeuratMeta$cellType==cellTypeInfo$cellTypes[i],"colour"])
  }
  
  mergedSeuratMeta$cellType<-factor(mergedSeuratMeta$cellType,
                                    levels=cellTypeInfo$cellTypes)
  mergedSeuratMeta$colour<-factor(mergedSeuratMeta$colour,levels=cellTypeInfo$colour)
  library(scales)
  percentUTarCellType<-ggplot(mergedSeuratMeta)+
    geom_boxplot(aes(x=cellType,
                     y=mergedSeuratMeta[,column]+1,
                     fill=cellType),
                 outlier.size = 0.1,
                 lwd=0.25)+
    scale_fill_manual(values = levels(mergedSeuratMeta$colour))+
    theme_bw()+
    NoLegend()+
    ylab(paste0(column))+
    #pseudo_log_trans(sigma = 1, base = exp(1))+
    scale_y_log10()+
    theme(axis.text.x=element_text(colour = levels(mergedSeuratMeta$colour),size=4,angle=45,hjust=1),
          axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(size=8,face="bold"),
          axis.title.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.minor.y = element_blank())
  
  return(percentUTarCellType)
}

cellTypeReadsFuncSample<-function(mergedSeuratMeta,column="percent.outHMM"){
  # calculate median percent out for each 
  channel<-unique(as.character(mergedSeuratMeta$channel))
  channelOutPer<-vector()
  for(chann in channel){
    channelOutPer<-c(channelOutPer,median(mergedSeuratMeta[,column][mergedSeuratMeta$channel==chann]))
  }
  channelInfo<-data.frame(cbind(channel,channelOutPer))
  channelInfo$channelOutPer<-as.numeric(levels(channelInfo$channelOutPer))[channelInfo$channelOutPer]
  channelInfo <- channelInfo[order(-channelInfo$channelOutPer),]
  channelInfo$channel<-as.character(channelInfo$channel)
  
  #for(i in 1:nrow(cellTypeInfo)){
  #  cellTypeInfo$colour[i]<-unique(mergedSeuratMeta[mergedSeuratMeta$cellType==cellTypeInfo$cellTypes[i],"colour"])
  #}
  #
  mergedSeuratMeta$channel<-factor(mergedSeuratMeta$channel,
                                    levels=channelInfo$channel)
  #mergedSeuratMeta$colour<-factor(mergedSeuratMeta$colour,levels=cellTypeInfo$colour)
  
  percentUTarCellType<-ggplot(mergedSeuratMeta)+
    geom_boxplot(aes(x=mergedSeuratMeta$channel,
                     y=mergedSeuratMeta[,column],
                     fill=mergedSeuratMeta$channel),
                 outlier.size = 0.1,
                 lwd=0.25)+
    #scale_fill_manual(values = levels(mergedSeuratMeta$colour))+
    theme_bw()+
    NoLegend()+
    ylab(column)+
    scale_y_log10()+
    theme(axis.text.x=element_text(size=4,angle=45,hjust=1),
          axis.text.y=element_text(colour = "black"),
          axis.title.y = element_text(size=8,face="bold"),
          axis.title.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.minor.y = element_blank())
  
  return(percentUTarCellType)
}


runIndividualAnalysis<-function(ind1){
  # Libraries
  library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
  library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
  library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
  source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
  library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
  
  # load up metafile
  mouseLemurMeta10X <- read.csv("/workdir/fw262/mouseLemur/tabula-microcebus-tenx-alltissues.csv")
  # load up conversion from channel to sample same in /mnt/ibm_sm/michaelw/TARAnalysis/pipeline_major
  library(readxl)
  sameNamesMeta <- data.frame(read_excel("/workdir/fw262/mouseLemur/sameNamesMeta.xlsx",col_names = FALSE))
  rownames(sameNamesMeta)<-sameNamesMeta[,1]
  sameNamesMeta[is.na(sameNamesMeta)]<-""
  # add in conversion to meta file
  mouseLemurMeta10X$mwSample<-sameNamesMeta[as.character(mouseLemurMeta10X$channel),2]
  #####################################################
  
  # look at plasma cells
  geneCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*gene*")
  tarCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*TARs*")
  geneSpleen<-system(geneCmd,intern = T)
  tarSpleen<-system(tarCmd,intern = T)
  
  ### this is important here
  seuratSpleen<-generateSeuratForMouseLemurGeneUTAR(hmmMat=tarSpleen,geneMat=geneSpleen)
  seuratSpleen <- NormalizeData(seuratSpleen)
  seuratSpleen <- ScaleData(seuratSpleen, features = rownames(seuratSpleen))
  
  mouseLemurMeta10X<-mouseLemurMeta10X[!duplicated(mouseLemurMeta10X$X),]
  cellAnno<-mouseLemurMeta10X[,c("compartment","cell_ontology_class","tissue","mwSample")]
  rownames(cellAnno)<-mouseLemurMeta10X$X
  cellAnno<-cellAnno[cellAnno$mwSample==ind1,]
  
  # create metadata of cell types and tissue
  cellTypes<-as.character(cellAnno$cell_ontology_class)
  cellTypes[cellTypes==""]<-"doublet"
  names(cellTypes)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=cellTypes,
    col.name = "cellType")
  
  # create metadata of cell types and tissue
  compartment<-as.character(cellAnno$compartment)
  compartment[compartment==""]<-"doublet"
  names(compartment)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=compartment,
    col.name = "compartment")
  
  # create metadata of cell types and tissue
  tissue<-as.character(cellAnno$tissue)
  tissue[tissue==""]<-"doublet"
  names(tissue)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=tissue,
    col.name = "tissue")
  
  numdash<-str_count(rownames(seuratSpleen),"-")
  outInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-0")
  inInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(seuratSpleen)[geneInd]
  inGene<-rownames(seuratSpleen)[inInd]
  outGene<-rownames(seuratSpleen)[outInd]
  
  # find differential genes based on cell types
  Idents(seuratSpleen)<-seuratSpleen$cellType
  cellTypeMarkers <- FindAllMarkers(seuratSpleen, features = outGene, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
  top5<-cellTypeMarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
  top5<-top5[!duplicated(top5$gene),]
  
  # filter for differential genes in all cells
  #load("/local/workdir/fw262/mouseLemur/cellTypeMarkersoutGene.RData")
  #top5<-cellTypeMarkersoutGene %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
  #top5<-top5[(top5$cluster %in% unique(seuratSpleen$cellType)),]
  #top5<-top5[!duplicated(top5$gene),]
  
  #seuratSpleen <- RunPCA(object = seuratSpleen, reduction.name="pca_uTAR_gene",features=c(geneOnly,outGene))
  seuratSpleen <- RunPCA(object = seuratSpleen, reduction.name="pca_uTAR",features = outGene)
  seuratSpleen <- RunPCA(object = seuratSpleen, features = geneOnly) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  seuratSpleen <- RunUMAP(seuratSpleen,reduction.name="umap_gene",reduction="pca",dims=1:50)
  seuratSpleen <- RunUMAP(seuratSpleen,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:50)
  #seuratSpleen <- RunUMAP(seuratSpleen,reduction.name="umap_uTAR_gene",reduction="pca_uTAR_gene",dims=1:50)

  ######## randomly assign read counts
  test<-data.frame(GetAssayData(object = seuratSpleen, slot = "counts"))
  test2<-test
  for(i in 1:ncol(test)){
    set.seed(i)
    test2[outGene,i]<-sample(test[outGene,i])
  }
  seurRand<-CreateAssayObject(test2)
  seuratSpleen@assays$randomUTAR<-seurRand
  DefaultAssay(seuratSpleen)<-"randomUTAR"
  
  seuratSpleen <- NormalizeData(seuratSpleen,assay = "randomUTAR")
  seuratSpleen <- ScaleData(seuratSpleen,features = rownames(seurRand),assay = "randomUTAR")
  seuratSpleen <- RunPCA(object = seuratSpleen,assay = "randomUTAR",reduction.name="pca_uTARrand",features = outGene)
  seuratSpleen <- RunUMAP(seuratSpleen,assay = "randomUTAR",reduction.name="umap_uTARrand",reduction="pca_uTARrand",dims=1:50)
  

  return(c(seuratSpleen,cellTypeMarkers))
}

runIndividualAnalysis2<-function(ind1,nfeat=5000){
  # Libraries
  library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
  library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
  library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
  source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
  library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
  
  # load up metafile
  mouseLemurMeta10X <- read.csv("/workdir/fw262/mouseLemur/tabula-microcebus-tenx-alltissues.csv")
  # load up conversion from channel to sample same in /mnt/ibm_sm/michaelw/TARAnalysis/pipeline_major
  library(readxl)
  sameNamesMeta <- data.frame(read_excel("/workdir/fw262/mouseLemur/sameNamesMeta.xlsx",col_names = FALSE))
  rownames(sameNamesMeta)<-sameNamesMeta[,1]
  sameNamesMeta[is.na(sameNamesMeta)]<-""
  # add in conversion to meta file
  mouseLemurMeta10X$mwSample<-sameNamesMeta[as.character(mouseLemurMeta10X$channel),2]
  #####################################################
  
  # look at plasma cells
  geneCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*gene*")
  tarCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*TARs*")
  geneSpleen<-system(geneCmd,intern = T)
  tarSpleen<-system(tarCmd,intern = T)
  
  ### this is important here
  seuratSpleen<-generateSeuratForMouseLemurTAR3(hmmMat=tarSpleen)
  #seuratSpleen<-FindVariableFeatures(seuratSpleen,nfeatures=nfeat)
  seuratSpleen <- NormalizeData(seuratSpleen)
  seuratSpleen <- ScaleData(seuratSpleen, features = rownames(seuratSpleen))
  
  mouseLemurMeta10X<-mouseLemurMeta10X[!duplicated(mouseLemurMeta10X$X),]
  cellAnno<-mouseLemurMeta10X[,c("compartment","cell_ontology_class","tissue","mwSample")]
  rownames(cellAnno)<-mouseLemurMeta10X$X
  cellAnno<-cellAnno[cellAnno$mwSample==ind1,]
  
  # create metadata of cell types and tissue
  cellTypes<-as.character(cellAnno$cell_ontology_class)
  cellTypes[cellTypes==""]<-"doublet"
  names(cellTypes)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=cellTypes,
    col.name = "cellType")
  
  # create metadata of cell types and tissue
  compartment<-as.character(cellAnno$compartment)
  compartment[compartment==""]<-"doublet"
  names(compartment)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=compartment,
    col.name = "compartment")
  
  # create metadata of cell types and tissue
  tissue<-as.character(cellAnno$tissue)
  tissue[tissue==""]<-"doublet"
  names(tissue)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=tissue,
    col.name = "tissue")
  
  numdash<-str_count(rownames(seuratSpleen),"-")
  outInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-0")
  inInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(seuratSpleen)[geneInd]
  inGene<-rownames(seuratSpleen)[inInd]
  outGene<-rownames(seuratSpleen)[outInd]
  
  # find differential genes based on cell types
  Idents(seuratSpleen)<-seuratSpleen$cellType
  #cellTypeMarkersUTAR <- FindAllMarkers(seuratSpleen, features = outGene, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
  #cellTypeMarkersATAR <- FindAllMarkers(seuratSpleen, features = inGene, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
  #cellTypeMarkersGENE <- FindAllMarkers(seuratSpleen, features = geneOnly, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
  #cellTypeMarkers<-rbind(cellTypeMarkersUTAR,cellTypeMarkersATAR,cellTypeMarkersGENE)
  cellTypeMarkers <- FindAllMarkers(seuratSpleen, features = rownames(seuratSpleen), only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
  cellTypeMarkers$mwSample<-ind1
  diffUTars<-cellTypeMarkers$gene[endsWith(cellTypeMarkers$gene,"-0")]
  seuratSpleen[["percent.ouHMM_0.5FC"]]<-PercentageFeatureSet(seuratSpleen,features = diffUTars)
  # filter for differential genes in all cells
  #load("/local/workdir/fw262/mouseLemur/cellTypeMarkersoutGene.RData")
  #top5<-cellTypeMarkersoutGene %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
  #top5<-top5[(top5$cluster %in% unique(seuratSpleen$cellType)),]
  #top5<-top5[!duplicated(top5$gene),]
  
  #seuratSpleen <- RunPCA(object = seuratSpleen, reduction.name="pca_uTAR_gene",features=c(geneOnly,outGene))
  #seuratSpleen <- RunPCA(object = seuratSpleen, reduction.name="pca_uTAR",features = outGene[outGene %in% VariableFeatures(seuratSpleen)])
  #seuratSpleen <- RunPCA(object = seuratSpleen, features = geneOnly[(geneOnly %in% VariableFeatures(seuratSpleen))]) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  #seuratSpleen <- RunPCA(object = seuratSpleen, reduction.name="pca_nonVarGene",features = geneOnly[!(geneOnly %in% VariableFeatures(seuratSpleen))]) # run gene pca, clustering last so Idents gets filled with gene exp clustering
  #
  #seuratSpleen <- RunUMAP(seuratSpleen,reduction.name="umap_gene",reduction="pca",dims=1:50)
  #seuratSpleen <- RunUMAP(seuratSpleen,reduction.name="umap_uTAR",reduction="pca_uTAR",dims=1:50)
  #seuratSpleen <- RunUMAP(seuratSpleen,reduction.name="umap_nonVarGene",reduction="pca_nonVarGene",dims=1:50)
  #seuratSpleen <- RunUMAP(seuratSpleen,reduction.name="umap_uTAR_gene",reduction="pca_uTAR_gene",dims=1:50)
  
  ######## randomly assign read counts
  #test<-data.frame(GetAssayData(object = seuratSpleen, slot = "counts"))
  #test2<-test
  #for(i in 1:ncol(test)){
  #  set.seed(i)
  #  test2[outGene,i]<-sample(test[outGene,i])
  #}
  #seurRand<-CreateAssayObject(test2)
  #seuratSpleen@assays$randomUTAR<-seurRand
  #DefaultAssay(seuratSpleen)<-"randomUTAR"
  #
  #seuratSpleen <- NormalizeData(seuratSpleen,assay = "randomUTAR")
  #seuratSpleen <- ScaleData(seuratSpleen,features = rownames(seurRand),assay = "randomUTAR")
  #seuratSpleen <- RunPCA(object = seuratSpleen,assay = "randomUTAR",reduction.name="pca_uTARrand",features = outGene)
  #seuratSpleen <- RunUMAP(seuratSpleen,assay = "randomUTAR",reduction.name="umap_uTARrand",reduction="pca_uTARrand",dims=1:50)
  
  return(cellTypeMarkers)
  
  #return(c(seuratSpleen,list(cellTypeMarkers)))
}

runIndividualAnalysis3<-function(ind1,nfeat=5000,logfcVal=0.5){
  # Libraries
  library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
  library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
  library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
  source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
  library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
  
  # load up metafile
  mouseLemurMeta10X <- read.csv("/workdir/fw262/mouseLemur/tabula-microcebus-tenx-alltissues.csv")
  # load up conversion from channel to sample same in /mnt/ibm_sm/michaelw/TARAnalysis/pipeline_major
  library(readxl)
  sameNamesMeta <- data.frame(read_excel("/workdir/fw262/mouseLemur/sameNamesMeta.xlsx",col_names = FALSE))
  rownames(sameNamesMeta)<-sameNamesMeta[,1]
  sameNamesMeta[is.na(sameNamesMeta)]<-""
  # add in conversion to meta file
  mouseLemurMeta10X$mwSample<-sameNamesMeta[as.character(mouseLemurMeta10X$channel),2]
  #####################################################
  
  # look at plasma cells
  geneCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*gene*")
  tarCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*TARs*")
  geneSpleen<-system(geneCmd,intern = T)
  tarSpleen<-system(tarCmd,intern = T)
  
  ### this is important here
  seuratSpleen<-generateSeuratForMouseLemurTAR3(hmmMat=tarSpleen)
  #seuratSpleen<-FindVariableFeatures(seuratSpleen,nfeatures=nfeat)
  seuratSpleen <- NormalizeData(seuratSpleen)
  seuratSpleen <- ScaleData(seuratSpleen, features = rownames(seuratSpleen))
  
  mouseLemurMeta10X<-mouseLemurMeta10X[!duplicated(mouseLemurMeta10X$X),]
  cellAnno<-mouseLemurMeta10X[,c("compartment","cell_ontology_class","tissue","mwSample")]
  rownames(cellAnno)<-mouseLemurMeta10X$X
  cellAnno<-cellAnno[cellAnno$mwSample==ind1,]
  
  # create metadata of cell types and tissue
  cellTypes<-as.character(cellAnno$cell_ontology_class)
  cellTypes[cellTypes==""]<-"doublet"
  names(cellTypes)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=cellTypes,
    col.name = "cellType")
  
  # create metadata of cell types and tissue
  compartment<-as.character(cellAnno$compartment)
  compartment[compartment==""]<-"doublet"
  names(compartment)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=compartment,
    col.name = "compartment")
  
  # create metadata of cell types and tissue
  tissue<-as.character(cellAnno$tissue)
  tissue[tissue==""]<-"doublet"
  names(tissue)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=tissue,
    col.name = "tissue")
  
  numdash<-str_count(rownames(seuratSpleen),"-")
  outInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-0")
  inInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(seuratSpleen)[geneInd]
  inGene<-rownames(seuratSpleen)[inInd]
  outGene<-rownames(seuratSpleen)[outInd]
  
  # find differential genes based on cell types
  Idents(seuratSpleen)<-seuratSpleen$cellType
  cellTypeMarkers <- FindAllMarkers(seuratSpleen, features = rownames(seuratSpleen), only.pos = T, min.pct = 0.25, logfc.threshold = logfcVal)
  cellTypeMarkers$mwSample<-ind1
  diffUTars<-unique(cellTypeMarkers$gene[endsWith(cellTypeMarkers$gene,"-0")])
  seuratSpleen[["percent.outHMM_0.5FC"]]<-PercentageFeatureSet(seuratSpleen,features = diffUTars)
  metaOut<-seuratSpleen@meta.data
  metaOut$mwSample<-ind1
  
  return(list(metaOut,cellTypeMarkers))
}

# run diff analysis on cell type per compartment
runIndividualAnalysis4<-function(ind1,nfeat=5000,logfcVal=0.5){
  # Libraries
  library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
  library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
  library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
  source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
  library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
  
  # load up metafile
  mouseLemurMeta10X <- read.csv("/workdir/fw262/mouseLemur/tabula-microcebus-tenx-alltissues.csv")
  # load up conversion from channel to sample same in /mnt/ibm_sm/michaelw/TARAnalysis/pipeline_major
  library(readxl)
  sameNamesMeta <- data.frame(read_excel("/workdir/fw262/mouseLemur/sameNamesMeta.xlsx",col_names = FALSE))
  rownames(sameNamesMeta)<-sameNamesMeta[,1]
  sameNamesMeta[is.na(sameNamesMeta)]<-""
  # add in conversion to meta file
  mouseLemurMeta10X$mwSample<-sameNamesMeta[as.character(mouseLemurMeta10X$channel),2]
  #####################################################
  
  # look at plasma cells
  geneCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*gene*")
  tarCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*TARs*")
  geneSpleen<-system(geneCmd,intern = T)
  tarSpleen<-system(tarCmd,intern = T)
  
  ### this is important here
  seuratSpleen<-generateSeuratForMouseLemurTAR3(hmmMat=tarSpleen)
  #seuratSpleen<-FindVariableFeatures(seuratSpleen,nfeatures=nfeat)
  seuratSpleen <- NormalizeData(seuratSpleen)
  seuratSpleen <- ScaleData(seuratSpleen, features = rownames(seuratSpleen))
  
  mouseLemurMeta10X<-mouseLemurMeta10X[!duplicated(mouseLemurMeta10X$X),]
  cellAnno<-mouseLemurMeta10X[,c("compartment","cell_ontology_class","tissue","mwSample")]
  rownames(cellAnno)<-mouseLemurMeta10X$X
  cellAnno<-cellAnno[cellAnno$mwSample==ind1,]
  
  # create metadata of cell types and tissue
  cellTypes<-as.character(cellAnno$cell_ontology_class)
  cellTypes[cellTypes==""]<-"doublet"
  names(cellTypes)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=cellTypes,
    col.name = "cellType")
  
  # create metadata of cell types and tissue
  compartment<-as.character(cellAnno$compartment)
  compartment[compartment==""]<-"doublet"
  names(compartment)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=compartment,
    col.name = "compartment")
  
  # create metadata of cell types and tissue
  tissue<-as.character(cellAnno$tissue)
  tissue[tissue==""]<-"doublet"
  names(tissue)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=tissue,
    col.name = "tissue")
  
  # add in cell type for each compartment
  typeCompart<-paste(cellTypes,compartment,sep="-")
  names(typeCompart)<-names(cellTypes)
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=typeCompart,
    col.name = "typeCompart")
  
  numdash<-str_count(rownames(seuratSpleen),"-")
  outInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-0")
  inInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(seuratSpleen)[geneInd]
  inGene<-rownames(seuratSpleen)[inInd]
  outGene<-rownames(seuratSpleen)[outInd]
  
  # find differential genes based on cell types
  Idents(seuratSpleen)<-seuratSpleen$typeCompart
  cellTypeMarkers <- FindAllMarkers(seuratSpleen, features = rownames(seuratSpleen), only.pos = T, min.pct = 0.25, logfc.threshold = logfcVal)
  cellTypeMarkers$mwSample<-ind1
  diffUTars<-unique(cellTypeMarkers$gene[endsWith(cellTypeMarkers$gene,"-0")])
  seuratSpleen[["percent.outHMM_0.5FC"]]<-PercentageFeatureSet(seuratSpleen,features = diffUTars)
  metaOut<-seuratSpleen@meta.data
  metaOut$mwSample<-ind1
  
  return(list(metaOut,cellTypeMarkers,seuratSpleen))
}



runIndividualAnalysis_getSeu<-function(ind1,logfcVal=0.5){
  # Libraries
  library(ggplot2); library(data.table); library(caTools);library(scales); library("Seurat"); library(dplyr)
  library(ontologyPlot);library(ontologyIndex);library(dendextend);library(cowplot);library(VennDiagram)
  library("BCellMA");library(ggplot2); library(data.table); library(caTools);library(plyr);library(stringr);library(data.table); library(caTools); library("gridExtra");library(ggpubr)
  source("/fs/cbsuvlaminck2/workdir/fw262/ShaoPei/allOrganisms/analyzeAllOrgFunc.R")
  library(reticulate); library("cluster"); library(qpcR); library(gtools);library("rowr");library("pals")
  
  # load up metafile
  mouseLemurMeta10X <- read.csv("/workdir/fw262/mouseLemur/tabula-microcebus-tenx-alltissues.csv")
  # load up conversion from channel to sample same in /mnt/ibm_sm/michaelw/TARAnalysis/pipeline_major
  library(readxl)
  sameNamesMeta <- data.frame(read_excel("/workdir/fw262/mouseLemur/sameNamesMeta.xlsx",col_names = FALSE))
  rownames(sameNamesMeta)<-sameNamesMeta[,1]
  sameNamesMeta[is.na(sameNamesMeta)]<-""
  # add in conversion to meta file
  mouseLemurMeta10X$mwSample<-sameNamesMeta[as.character(mouseLemurMeta10X$channel),2]
  #####################################################
  
  # look at plasma cells
  geneCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*gene*")
  tarCmd<-paste0("ls -d /workdir/fw262/mouseLemur/dge_final/*",ind1,"*TARs*")
  geneSpleen<-system(geneCmd,intern = T)
  tarSpleen<-system(tarCmd,intern = T)
  
  ### this is important here
  seuratSpleen<-generateSeuratForMouseLemurTAR3(hmmMat=tarSpleen)
  #seuratSpleen<-FindVariableFeatures(seuratSpleen,nfeatures=nfeat)
  seuratSpleen <- NormalizeData(seuratSpleen)
  seuratSpleen <- ScaleData(seuratSpleen, features = rownames(seuratSpleen))
  
  mouseLemurMeta10X<-mouseLemurMeta10X[!duplicated(mouseLemurMeta10X$X),]
  cellAnno<-mouseLemurMeta10X[,c("compartment","cell_ontology_class","tissue","mwSample")]
  rownames(cellAnno)<-mouseLemurMeta10X$X
  cellAnno<-cellAnno[cellAnno$mwSample==ind1,]
  
  # create metadata of cell types and tissue
  cellTypes<-as.character(cellAnno$cell_ontology_class)
  cellTypes[cellTypes==""]<-"doublet"
  names(cellTypes)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=cellTypes,
    col.name = "cellType")
  
  # create metadata of cell types and tissue
  compartment<-as.character(cellAnno$compartment)
  compartment[compartment==""]<-"doublet"
  names(compartment)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=compartment,
    col.name = "compartment")
  
  # create metadata of cell types and tissue
  tissue<-as.character(cellAnno$tissue)
  tissue[tissue==""]<-"doublet"
  names(tissue)<-colnames(seuratSpleen)
  # add cell type metadata
  seuratSpleen<-AddMetaData(
    object=seuratSpleen,
    metadata=tissue,
    col.name = "tissue")
  
  numdash<-str_count(rownames(seuratSpleen),"-")
  outInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-0")
  inInd<-numdash>=5 & endsWith(rownames(seuratSpleen),"-1")
  geneInd<-!(outInd|inInd)
  geneOnly<-rownames(seuratSpleen)[geneInd]
  inGene<-rownames(seuratSpleen)[inInd]
  outGene<-rownames(seuratSpleen)[outInd]
  
  # find differential genes based on cell types
  Idents(seuratSpleen)<-seuratSpleen$cellType
  cellTypeMarkers <- FindAllMarkers(seuratSpleen, features = rownames(seuratSpleen), only.pos = T, min.pct = 0.25, logfc.threshold = logfcVal)
  cellTypeMarkers$mwSample<-ind1
  #diffUTars<-unique(cellTypeMarkers$gene[endsWith(cellTypeMarkers$gene,"-0")])
  #seuratSpleen[["percent.outHMM_0.5FC"]]<-PercentageFeatureSet(seuratSpleen,features = diffUTars)
  #metaOut<-seuratSpleen@meta.data
  #metaOut$mwSample<-ind1
  
  return(list(seuratSpleen,cellTypeMarkers))
}




#### functions here
getNFromList <- function(lst, n){
  sapply(lst, `[`, n)
}

smoothCovFunc<-function(input,span=0.5){
  x<-1:length(input)
  y<-as.numeric(input)
  smoothOut<-loess(y~x,span=span)
  return(predict(smoothOut))
}

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

argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

test <- function(xin,yin,w, span) {
  peaks <- argmax(x=xin, y=yin, w=w, span=span)
  
  plot(xin, yin, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""))
  lines(xin, peaks$y.hat,  lwd=2) #$
  y.min <- min(yin)
  sapply(peaks$i, function(i) lines(c(xin[i],xin[i]), c(y.min, peaks$y.hat[i]),
                                    col="Red", lty=2))
  points(xin[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
}


plotListOfCov2<-function(input,colorInput="black"){
  df<-data.frame(as.numeric(input))
  out<-ggplot(df,aes(x=1:nrow(df),y=df[,1]),wdith=1)+
    geom_area(stat="identity",fill=colorInput)+ theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8))+
    scale_y_continuous(limits=c(0,max(df[,1])),
                       breaks = c(0,max(df[])))+
    scale_x_continuous(limits=c(0,nrow(df)),
                       breaks = c(0,nrow(df)))
  
  return(out)
}


plotDiffExpNolab<-function(input,uTAR,bamFile,titleForPlot="uTAR",normalChrom=T,order=T){
  library(ggplot2)
  plotListOfCov2<-function(input,colorInput="black"){
    df<-data.frame(as.numeric(input))
    out<-ggplot(df,aes(x=1:nrow(df),y=df[,1]),wdith=1)+
      geom_area(stat="identity",fill=colorInput)+ theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8))+
      scale_y_continuous(limits=c(0,max(df[,1])),
                         breaks = c(0,max(df[])))+
      scale_x_continuous(limits=c(0,nrow(df)),
                         breaks = c(0,nrow(df)))
    
    return(out)
  }
  
  if (normalChrom){
    chr<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    fasta<-paste0(chr,":",start,"-",stop)
  } else {
    chr1<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    chr2<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",4))
    fasta<-paste0(chr1,"_",chr2,":",start,"-",stop)
  }
  covCom<-paste0("samtools depth -r ",fasta," ",bamFile," > temp",uTAR,".txt")
  test<-system(covCom)
  readTempFile<-paste0("temp",uTAR,".txt")
  temp <- read.delim(readTempFile, header=FALSE)
  rmCmd<-paste0("rm ",readTempFile)
  system(rmCmd)
  coverage<-list(temp[,3])
  
  allCovPlot<-plotListOfCov2(unlist(coverage))
  allCovPlot<-allCovPlot+labs(title=titleForPlot)+
    theme(plot.title = element_text(size=8,face = "bold"),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6))
  
  #out<-plot_grid(allCovPlot,featureP,ncol=1,align = "v",axis="lr",rel_heights = c(2,3))
  return(allCovPlot)
}

getCoverageOnly<-function(uTAR,bamFile,titleForPlot="uTAR",normalChrom=T,order=T){
  print(uTAR)

  if (normalChrom){
    chr<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    fasta<-paste0(chr,":",start,"-",stop)
  } else {
    chr1<-as.character(sapply(strsplit(uTAR,"-"),"[",1))
    chr2<-as.character(sapply(strsplit(uTAR,"-"),"[",2))
    start<-as.character(sapply(strsplit(uTAR,"-"),"[",3))
    stop<-as.character(sapply(strsplit(uTAR,"-"),"[",4))
    fasta<-paste0(chr1,"_",chr2,":",start,"-",stop)
  }
  covCom<-paste0("samtools depth -r ",fasta," ",bamFile," > temp/temp",uTAR,".txt")
  test<-system(covCom)
  readTempFile<-paste0("temp/temp",uTAR,".txt")
  temp <- read.delim(readTempFile, header=FALSE)
  #rmCmd<-paste0("rm ",readTempFile)
  #system(rmCmd)
  coverage<-list(temp[,3])
  
  return(unlist(coverage))
  
}

saveSpatialFunc<-function(input,feature,maxScale,dayText,ptsizefactor = 2){
  out<-SpatialPlot(input,features=feature,pt.size.factor = ptsizefactor)+
    #scale_fill_gradient(low="grey",high="brown",limits=c(0,maxScale))
    #scale_fill_brewer(palette = "RdBu")
    scale_fill_viridis(option = "inferno",limits=c(0,maxScale))
    #scale_color_gradient(limits=c(0,maxScale))
  
  outFile<-paste0("/workdir/fw262/chickenSpatial/figures/",feature,"_",dayText,".pdf")
  pdf(file=outFile, width=2, height=2, paper="special", bg="white",
      fonts="Helvetica", pointsize=6, useDingbats = F, )
  plot(out)
  dev.off()
  return()
}
