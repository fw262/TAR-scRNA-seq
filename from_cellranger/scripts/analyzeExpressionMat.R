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
  n.pcs <- 50
  res_param <- 0.8
} else if (length(args) == 9) {
  geneFile <- args[1] # first argument is gene File
  TARFile <- args[2] # second argument is TAR File
  nDiffFeat <- args[3]
  cells <- args[4]
  feats <- args[5]
  nFeatMin <- args[6]
  nFeatMax <- args[7]
  n.pcs <- args[8]
  res_param <- args[9]
} else {
  stop("Please specify all parameters.")
}

cat("Input arguments are good\n")
cat("Loading required packages (Matrix, Seurat, data.table, dplyr, stringr).\n")

if (!require('R.utils', quietly=T)) {
  install.packages('R.utils', repo="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(R.utils, verbose=F))

# if (!require('Matrix')) {
#   install.packages('Matrix', repo="https://cloud.r-project.org")
# }
suppressPackageStartupMessages(library(Matrix,  verbose=F))

if (!require('Seurat', quietly=T)) {
  install.packages('Seurat', repo="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(Seurat,  verbose=F)) # Please use Seurat >= v4.0

if (!require('data.table', quietly=T)) {
  install.packages("data.table", repo="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(data.table,  verbose=F))

if (!require('dplyr', quietly=T)) {
  install.packages("dplyr", repo="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(dplyr,  verbose=F))

if (!require('stringr', quietly=T)) {
  install.packages("stringr", repo="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(stringr, verbose=F))

cat("Finished loading packages (Seurat, data.table, dplyr, stringr).\n")

# Borrowed from the Marioni Lab, DropletUtils package (https://rdrr.io/github/MarioniLab/DropletUtils/src/R/write10xCounts.R)
#   (Had R version issues getting it to work as a dependency)
#' @importFrom utils write.table
#' @importFrom Matrix writeMM
#' @importFrom R.utils gzip
.write_sparse <- function(path, x, barcodes, gene.id, gene.symbol, gene.type) {
  # dir.create(path, showWarnings=FALSE)
  gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)

  gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
  mhandle <- file.path(path, "matrix.mtx")
  bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
  fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
  on.exit({
    close(bhandle)
    close(fhandle)
  })

  writeMM(x, file=mhandle)
  write(barcodes, file=bhandle)
  write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

  # Annoyingly, writeMM doesn't take connection objects.
  gzip(mhandle)

  return(NULL)
}

# Calculate the number of PCs that contain some proportion (95%) of the variance
.find_npcs <- function(seu, var.toal=0.95, reduction="pca"){
  if(is.null(seu@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }

  tmp.var <- (seu@reductions[[reduction]]@stdev)^2
  var.cut <- var.toal*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }

  return(n.pcs)
}

FindDiffExprFeatures <-
  function(
    geneFile,
    TARFile,
    nDiffFeat = 5,
    cells = 2,
    feats = 100,
    nFeatMin = 200,
    nFeatMax = 2500,
    n.pcs = 50,
    res_param = 0.8
  ){
    cat("Loading in expression matrix... ")
    geneMat <- Seurat::Read10X(geneFile)
    cat("Done.\n")

    # load TAR matrix and combine with gene expression matrix and subset for valid cells ####
    cat("Loading in TAR matrix... ")
    #TODO- clean this up...
    TARMat <- fread(
      TARFile,
      sep = "\t",
      header = T,
      stringsAsFactors = F,
      showProgress = F
    )
    rownames(TARMat) <- tolower(rownames(TARMat))
    TARMat <- as.data.frame(TARMat) # force to data frame
    rownames(TARMat) <- TARMat$GENE # set the rownames as GENEs
    TARMat <- TARMat[, -1] # take out first column
    TARMat <- as.sparse(TARMat) #force to sparse to match with gene matrix

    cat("Done.\n")
    cat("Loaded in ", ncol(TARMat), "cells and ", nrow(TARMat), "TARs.\n")

    cat("Saving TAR matrix in MTX format...")
    newMatDirName <- paste0(
      stringr::str_remove(TARFile,"TAR_expression_matrix_withDir.txt.gz"),
      "TAR_feature_bc_matrix"
    )
    if(!dir.exists(newMatDirName)){
      dir.create(newMatDirName)
    }
    if(!file.exists(paste0(newMatDirName,"/matrix.mtx.gz"))){
      .write_sparse(
        path=newMatDirName,
        x=TARMat,
        barcodes=colnames(TARMat),
        gene.id=rownames(TARMat),
        gene.symbol=rownames(TARMat),
        gene.type="TAR"
      )
    }
    cat("Done.\n")

    combined.mat <- rbind(geneMat, TARMat[,colnames(geneMat)])
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
    combined.seu[["percent.outTAR"]] <-
      PercentageFeatureSet(
        out,
        features = outGene
      )

    # Preprocess data ####
    cat("Preprocessing data... \n")
    combined.seu <-
      NormalizeData(combined.seu, verbose=F) %>%
      ScaleData(features = rownames(combined.seu),verbose=F)

    combined.seu <- RunPCA(
      combined.seu,
      reduction.name = "pca_uTAR",
      reduction.key="PCuTAR_",
      features = outGene,
      n.pcs = n.pcs,
      nfeatures.print = 0,
      verbose=F
    )

    # run gene pca, clustering last so Idents gets filled with gene exp clustering
    combined.seu <- RunPCA(
      combined.seu,
      reduction.name = "pca_genes",
      reduction.key="PCgenes_",
      features = geneOnly.features,
      n.pcs = n.pcs,
      nfeatures.print = 0,
      verbose=F
    )

    # Find number of PCs that account for 95% of variance
    npcs_significant <- .find_npcs(
      combined.seu,
      var.toal=0.95,
      reduction="pca_genes"
    )

    combined.seu <- FindNeighbors(
      combined.seu,
      reduction="pca_genes",
      dims=1:npcs_significant
    )

    combined.seu <- FindClusters(
      combined.seu,
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
    markersUTAR <- markersUTAR %>%
      group_by(cluster) %>%
      top_n(n = nDiffFeat, wt = avg_log2FC)

    markersGenes <- allMarkers[geneInd, ]
    markersGenes <- markersGenes[order(markersGenes$avg_log2FC, decreasing = T), ]
    markersGenes <- markersGenes %>%
      group_by(cluster) %>%
      top_n(n = nDiffFeat, wt = avg_log2FC)

    cat("Finished finding differentially expressed genes and uTAR markers.\n")

    return(rbind(markersGenes, markersUTAR))
  }

diffFeatures <- FindDiffExprFeatures(
  geneFile = geneFile,
  TARFile = TARFile,
  nDiffFeat = nDiffFeat,
  cells = cells,
  feats = feats,
  nFeatMin = nFeatMin,
  nFeatMax = nFeatMax,
  n.pcs = n.pcs,
  res_param = res_param
)

if(dirname(TARFile) != "."){
  sampleName <- tail(unlist(strsplit(dirname(TARFile), "/")), n = 1)
  output <- paste0(dirname(TARFile), "/", "diff_Expressed_Features.txt")
}else{
  output <- paste0("diff_Expressed_Features.txt")
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
