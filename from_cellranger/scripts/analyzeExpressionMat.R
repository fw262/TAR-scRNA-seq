## This script finds differentially expressed uTARs and returns a list of differentially expressed uTARs with genes
args = (commandArgs(TRUE))
if (length(args) == 0) {
  stop("Please specify gene and TAR expression matrices.")
} else if (length(args) == 2) {
  geneFile <- args[1] # first argument is gene File
  TARFile <- args[2] # second argument is TAR File
  nDiffFeat <- 2000
  cells <- 2
  feats <- 100
  nFeatMin <- 200
  nFeatMax <- 2500
  pcaDims <- 20
  res_param <- 0.8
} else if (length(args) == 9) {
  geneFile <- args[1] # first argument is gene File
  TARFile <- args[2] # second argument is TAR File
  nDiffFeat <- args[3]
  cells <- args[4]
  feats <- args[5]
  nFeatMin <- args[6]
  nFeatMax <- args[7]
  pcaDims <- args[8]
  res_param <- args[9]
} else {
  stop("Please specify all parameters.")
}

cat("Input arguments are good\n")
cat("Loading required packages (Seurat, data.table, dplyr, stringr).\n")

if (!require('Seurat')) {
  install.packages('Seurat')
}
library(Seurat,  verbose=F) # Please use Seurat >= v4.0

if (!require('data.table')) {
  install.packages("data.table")
}
library(data.table,  verbose=F)

if (!require('dplyr')) {
  install.packages("dplyr")
}
library(dplyr,  verbose=F)

if (!require('stringr')) {
  install.packages("stringr")
}
library(stringr, verbose=F)

cat("Finished loading packages (Seurat, data.table, dplyr, stringr).\n")

FindDiffExprFeatures <-
  function(
    geneFile,
    hmmFile,
    nDiffFeat = 5,
    cells = 2,
    feats = 100,
    nFeatMin = 200,
    nFeatMax = 2500,
    pcaDims = 20,
    res_param = 0.8
  ){
    cat("Loading in expression matrix...")
    geneMat <- Seurat::Read10X(geneFile)
    cat("Done.\n")
    
    # load HMM matrix and combine with gene expression matrix and subset for valid cells ####
    cat("Loading in HMM matrix...")
    
    #TODO- clean this up...
    hmmMat <- fread(
      hmmFile,
      sep = "\t",
      header = T,
      stringsAsFactors = F,
      showProgress = F
    ) 
    rownames(hmmMat) <- tolower(rownames(hmmMat))
    hmmMat <- as.data.frame(hmmMat) # force to data frame
    rownames(hmmMat) <- hmmMat$GENE # set the rownames as GENEs
    hmmMat <- hmmMat[, -1] # take out first column
    hmmMat <- as.sparse(hmmMat) #force to sparse to match with gene matrix
    cat("Done.\n")
    
    combined.mat <- rbind(geneMat, hmmMat[,colnames(geneMat)])
    combined.seu <- 
      CreateSeuratObject(counts = combined.mat) 
    cat("Finished creating Seurat object.\n")
    
    # identify aTARs and uTARs
    numdash <- str_count(rownames(combined.seu), "-")
    outInd <- numdash >= 5 & endsWith(rownames(combined.seu), "-0")
    inInd <- numdash >= 5 & !endsWith(rownames(combined.seu), "-0")
    geneInd <- !(outInd | inInd)
    geneOnly.features <- rownames(combined.seu)[geneInd]
    inGene <- rownames(combined.seu)[inInd]
    outGene <- rownames(combined.seu)[outInd]
    
    out <- combined.seu[c(inGene, outGene), ]
    combined.seu[["percent.outHMM"]] <-
      PercentageFeatureSet(
        out, 
        features = outGene
      )
    
    # Preprocess data ####
    cat("Preprocessing data...\n") 
    combined.seu <- 
      NormalizeData(combined.seu)%>%
      ScaleData(features = rownames(combined.seu)) %>%
      RunPCA(
        reduction.name = "pca_uTAR",
        reduction.key="PCuTAR_",
        features = outGene,
        verbose=F
      ) %>%
      RunPCA(# run gene pca, clustering last so Idents gets filled with gene exp clustering
        features = geneOnly.features,
        verbose=F
      )  %>% FindNeighbors(
        dims = 1:pcaDims
      ) %>% 
      FindClusters(
        resolution = res_param
      )
    cat("Done preprocessing.\n")
    
    # DGE ####
    cat("Running DGE...\n")
    allMarkers <-
      FindAllMarkers(
        combined.seu,
        only.pos = T,
        min.pct = 0.25,
        logfc.threshold = 0.25
      )
    cat("Done with DGE\n")
    
    numdash <- str_count(allMarkers$gene, "-")
    outInd <- numdash >= 5 & endsWith(allMarkers$gene, "-0")
    inInd <- numdash >= 5 & !endsWith(allMarkers$gene, "-0")
    geneInd <- !(outInd | inInd)
    
    markersUTAR <- allMarkers[outInd, ]
    markersUTAR <- markersUTAR[order(markersUTAR$avg_log2FC, decreasing = T), ]
    markersUTAR <- markersUTAR %>% group_by(cluster) %>% top_n(n = nDiffFeat, wt = avg_log2FC)
    
    markersGenes <- allMarkers[geneInd, ]
    markersGenes <- markersGenes[order(markersGenes$avg_log2FC, decreasing = T), ]
    markersGenes <- markersGenes %>% group_by(cluster) %>% top_n(n = nDiffFeat, wt = avg_log2FC)
    
    cat("Finished finding differentially expressed genes and uTAR markers.\n")
    
    return(rbind(markersGenes, markersUTAR))
  }

diffFeatures <- FindDiffExprFeatures(
  geneFile = geneFile,
  hmmFile = TARFile,
  nDiffFeat = nDiffFeat,
  cells = cells,
  feats = feats,
  nFeatMin = nFeatMin,
  nFeatMax = nFeatMax,
  pcaDims = pcaDims,
  res_param = res_param
)

if(dirname(TARFile) != "."){
  sampleName <- tail(unlist(strsplit(dirname(TARFile), "/")), n = 1)
  output <-
    paste0(dirname(TARFile), "/", "diffFeatures.txt")
}else{
  output <- paste0("diffFeatures.txt")
}

colnames(diffFeatures) <-
  c(
    "p_val",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "cluster",
    "gene_or_uTAR"
  )

write.table(
  diffFeatures,
  file = output,
  quote = F,
  row.names = F,
  col.names = T,
  sep = "\t"
)
cat(paste0("Output file stored in ", output, ".\n"))
