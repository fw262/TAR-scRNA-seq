## This script finds differentially expressed uTARs and returns a list of differentially expressed uTARs with genes
# saves 2 files (fasta region and diff uTAR markers)
args = (commandArgs(TRUE))
if (length(args) == 2) {
  inputDiffFeaturesFile <- args[1]
  bamFile <- args[2]
} else {
  stop("Please specify differential markers and bamFile.")
}

cat("Input arguments are good\n")
cat("Loading required packages (data.table, dplyr, stringr).\n")

if (!require('data.table')) {
  install.packages("data.table")
}
suppressPackageStartupMessages(library(data.table))

if (!require('dplyr')) {
  install.packages("dplyr")
}
suppressPackageStartupMessages(library(dplyr))

if (!require('stringr')) {
  install.packages("stringr")
}
suppressPackageStartupMessages(library(stringr))

cat("Finished loaded packages (data.table, dplyr, stringr).\n")

diffFeatures <- read.delim(inputDiffFeaturesFile, sep = "\t")
cat("Finished loading differentially expressed genes and uTARs.\n")

#### look at features and filter for peak
# extract only differentially expressed uTARs
diffFeatures$gene_or_uTAR <- as.character(diffFeatures$gene_or_uTAR)
numdash <- str_count(diffFeatures$gene_or_uTAR, "-")
outInd <- numdash >= 5 & endsWith(diffFeatures$gene_or_uTAR, "-0")
inInd <- numdash >= 5 & !endsWith(diffFeatures$gene_or_uTAR, "-0")
geneInd <- !(outInd | inInd)
uTARFeatures <- diffFeatures[outInd, ]
colnames(uTARFeatures) <-
  c("p_val",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "cluster",
    "uTAR")

returnChromString <- function(input) {
  out <- unlist(strsplit(input, "-", fixed = T))
  # check if negative strand or positive strand
  if (out[length(out) - 2] == "") {
    out2 <- head(out, -6)
    out3 <- paste(out2, collapse = "_")
    return(out3)
  } else {
    out2 <- head(out, -5)
    out3 <- paste(out2, collapse = "_")
    return(out3)
  }
}
returnStartString <- function(input) {
  out <- unlist(strsplit(input, "-", fixed = T))
  # check if negative strand or positive strand
  if (out[length(out) - 2] == "") {
    # negative strand case
    return(out[length(out) - 5])
  } else {
    return(out[length(out) - 4])
  }
}
returnEndString <- function(input) {
  out <- unlist(strsplit(input, "-", fixed = T))
  # check if negative strand or positive strand
  if (out[length(out) - 2] == "") {
    # negative strand case
    return(out[length(out) - 4])
  } else {
    return(out[length(out) - 3])
  }
}
returnDirString <- function(input) {
  out <- unlist(strsplit(input, "-", fixed = T))
  # check if negative strand or positive strand
  if (out[length(out) - 2] == "") {
    # negative strand case
    return("-")
  } else {
    return("+")
  }
}

uTARFeatures$chr <-
  as.character(sapply(as.character(uTARFeatures$uTAR), FUN = returnChromString))
uTARFeatures$start <-
  as.character(sapply(as.character(uTARFeatures$uTAR), FUN = returnStartString))
uTARFeatures$end <-
  as.character(sapply(as.character(uTARFeatures$uTAR), FUN = returnEndString))
uTARFeatures$dir <-
  as.character(sapply(as.character(uTARFeatures$uTAR), FUN = returnDirString))
uTARFeatures$fasta <-
  as.character(paste0(uTARFeatures$chr, ":", uTARFeatures$start, "-", uTARFeatures$end))

# extract coverage
getCoverageOnly <- function(fasta, bamFile) {
  interval <- strsplit(fasta, ":")[[1]][2]
  interval <- as.numeric(strsplit(interval, "-")[[1]])
  covCom <-
    paste0("samtools depth -r ", fasta, " ", bamFile, " > temp", fasta)
  test <- system(covCom)
  #check if depth file is empty
  info = file.info(paste0("temp", fasta))
  if (info$size == 0) {
    system(paste0("rm temp", fasta))
    return(rep(1, diff(interval) + 1))
  }
  temp <- read.delim(paste0("temp", fasta), header = FALSE)
  system(paste0("rm temp", fasta))
  coverage <- list(temp[, 3])
  if (length(unlist(coverage)) < diff(interval) * 0.9) {
    return(rep(1, diff(interval) + 1))
  }
  return(unlist(coverage))
}
tissueFeaturesCov <-
  lapply(X = uTARFeatures$fasta,
         FUN = getCoverageOnly,
         bamFile = bamFile)

#smooth coverage
smoothCovFunc <- function(input, span = 0.5) {
  x <- 1:length(input)
  y <- as.numeric(input)
  smoothOut <- loess(y ~ x, span = span)
  return(predict(smoothOut))
}
smoothCov <- lapply(X = tissueFeaturesCov, FUN = smoothCovFunc, span = 0.25)

# get half width maximum of smoothed peak
findFHWM <- function(input) {
  x <- 1:length(input)
  y <- input
  xmax <- x[y == max(y)]

  if (max(y) > 100) {
    x1 <-
      rev(x[x < xmax])[which.min(rev(round(abs(y[x < xmax] - max(
        y
      ) / 2), digits = -1)))]
  } else {
    x1 <-
      rev(x[x < xmax])[which.min(rev(round(abs(y[x < xmax] - max(
        y
      ) / 2))))]
  }

  x2 <- x[x > xmax][which.min((round(abs(y[x > xmax] - max(
    y
  ) / 2))))]

  #x1 <- rev(x[x < xmax])[(which(diff(sign(diff(abs(rev(y[x < xmax]-max(y)/2))*-1)))==-2)+1)[1]]
  #x2 <- x[x > xmax][(which(diff(sign(diff(abs((y[x > xmax]-max(y)/2))*-1)))==-2)+1)[1]]

  if (length(x1) == 0) {
    return(c(0, x2))
  } else if (length(x2) == 0) {
    return(c(x1, length(input)))
  } else{
    return(c(x1, x2))
  }
}
indicesToBlast <- lapply(X = smoothCov, FUN = findFHWM)

# add index shift to create fasta for peaks
shiftDF <-
  matrix(
    unlist(indicesToBlast),
    nrow = length(indicesToBlast),
    byrow = T
  )
tooShortInd <- abs(shiftDF[, 2] - shiftDF[, 1]) < 50
negativeInd <- (shiftDF[, 2] - shiftDF[, 1]) < 0
indTofix <- tooShortInd | negativeInd
shiftDF[indTofix, 1] <- 0
shiftDF[indTofix, 2] <-
  (as.numeric(uTARFeatures$end[indTofix]) - as.numeric(uTARFeatures$start[indTofix]))

uTARFeatures <- cbind(uTARFeatures, data.frame(shiftDF))
uTARFeatures$fastaPeak <-
  paste0(
    uTARFeatures$chr,
    ":",
    as.numeric(uTARFeatures$start) + uTARFeatures$X1,
    "-",
    as.numeric(uTARFeatures$start) + uTARFeatures$X2
  )
#indTooShort<-abs(uTARFeatures$X2 - uTARFeatures$X1)<50
#if peaks too short, just keep full uTAR
#uTARFeatures$fastaPeak[indTooShort]<-paste0(uTARFeatures$chr[indTooShort],":",as.numeric(uTARFeatures$start[indTooShort]),"-",as.numeric(uTARFeatures$end[indTooShort]))

if (dirname(inputDiffFeaturesFile) != ".") {
  sampleName <- tail(unlist(strsplit(dirname(inputDiffFeaturesFile), "/")), n = 1)
  output <- paste0(dirname(inputDiffFeaturesFile), "/", sampleName,"_seqsToBlast.txt")
} else {
  output <- paste0("seqsToBlast.txt")
}

write.table(
  uTARFeatures$fastaPeak,
  file = output,
  quote = F,
  row.names = F,
  col.names = F
) # save fasta format

if (dirname(inputDiffFeaturesFile) != ".") {
  sampleName <- tail(unlist(strsplit(dirname(inputDiffFeaturesFile), "/")), n = 1)
  output <- paste0(dirname(inputDiffFeaturesFile), "/", "TAR_diff_uTAR_Features.txt")
} else {
  output <- paste0("TAR_diff_uTAR_Features.txt")
}
write.table(
  uTARFeatures,
  file = output,
  quote = F,
  row.names = F,
  col.names = T,
  sep = "\t"
)
