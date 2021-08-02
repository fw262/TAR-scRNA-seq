# R --vanilla --slave --args $(pwd) PREFIX.bed < SingleCellHMM.R

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input_f = args[2] #PREFIX.bed.gz #05-007B1_gene_exon_tagged.REF_chr22_split.bed.gz
if (!requireNamespace("BiocManager", quietly = TRUE)){
  #TODO: check global options for default repo, and set if it isn't found
  # options()

  install.packages("BiocManager", repos = "http://cran.us.r-project.org") # DWM- explicitly set repo here to avoid install error
}

if (!require('rtracklayer')){
  BiocManager::install("rtracklayer")
}
if (!require('groHMM')){
  BiocManager::install("groHMM")
}

suppressPackageStartupMessages(library(groHMM, quietly = TRUE))
suppressPackageStartupMessages(library(rtracklayer, quietly = TRUE))

print("installed/loaded packages (rtracklayer & groHMM)")

GRangeTobed <- function(gr, f_name){
  df <- data.frame(seqnames=seqnames(gr),
    starts=start(gr)-1,
    ends=end(gr),
    names=c(rep(".", length(gr))),
    scores=c(rep(".", length(gr))),
    strands=strand(gr))

  write.table(df, file=f_name, quote=F, sep="\t", row.names=F, col.names=F)
}

detectTranscripts_AllEM <- function(reads=NULL, Fp=NULL, Fm=NULL, LtProbA=-5,
    LtProbB=-200, UTS=5, size=50, threshold=0.1, debug=TRUE, ...) {

    stopifnot(!is.null(reads)|(!is.null(Fp) & !is.null(Fm)))

    ## Setup/Window Analysis/Casting.
    epsilon <- 0.001

    ## Allow equilavent form of Fp and Fm to be spcified in the function
    ## automatically.
    if(is.null(Fp) & is.null(Fm)) {
     Fp <- windowAnalysis(reads=reads, strand="+", windowSize=size, ...)
     Fm <- windowAnalysis(reads=reads, strand="-", windowSize=size, ...)
    }

    nFp <- NROW(Fp)
    nFm <- NROW(Fm)
    CHRp <- as.character(names(Fp))
    CHRm <- as.character(names(Fm))

    size <- as.integer(size)
    ANS <- NULL

    ## Set up initial HMM variables.
    HMM <- list()
    HMM$nstates <- as.integer(2)
    HMM$ePrDist <- c("dgamma", "dgamma")
    ## CGD: 3-3-13: Still legacy. Switch to integrating gamma between read
    ## and read+1

    HMM$iProb <- as.double(log(c(1.0,0.0)))
                                    ## Non-transcribed,  transcribed.
    HMM$ePrVars <- as.list(data.frame(c(UTS, 1/UTS, -1), c(0.5, 10, -1)))
    HMM$tProb <- as.list(data.frame(c(log(1-exp(LtProbA)), LtProbA),
        c(LtProbB, log(1-exp(LtProbB))) ))

    ## Cast counts to a real, and combine +/- strand into one list variable.
    ##  Treat like separate training sequences (they really are).
    FT <- list()    # MHC; 7/2/2014, bioconductor complains F thinking as False
    for(i in 1:nFp) FT[[i]]<- as.double(Fp[[i]]+1)
    ## CGD: 3-3-13: Still legacy.  Switch to integrating gamma between read and
    ## read+1

    for(i in 1:nFm) FT[[i+nFp]] <- as.double(Fm[[i]]+1)
    ## CGD: 3-3-13: Still legacy.  Switch to integrating gamma between read and
    ## read+1

    ## In case the above command copies, rather than points ...
    ## free unused memory.
    remove(Fp)
    remove(Fm)

    ## Run EM algorithm.
    BWem <- .Call("RBaumWelchEM", HMM$nstates, FT, as.integer(1),
                HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb,
                as.double(threshold), c(TRUE, TRUE), c(TRUE, TRUE),
                as.integer(1), TRUE, PACKAGE="groHMM")
               # Update Transitions, Emissions.

    ## Translate these into transcript positions.
    for(i in seq_along(CHRp)) {
        ans <- .Call("getTranscriptPositions", as.double(BWem[[3]][[i]]),
                    as.double(0.5), size, PACKAGE="groHMM")
        Nrep <- NROW(ans$Start)
        # ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep),
        #           chromStart =ans$Start, chromEnd =ans$End,
        #   name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("+",
        #           Nrep)))
        ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep),
            start=ans$Start, end =ans$End, strand =rep("+", Nrep)))
    }

    for(i in seq_along(CHRm)) {
        ans <- .Call("getTranscriptPositions", as.double(BWem[[3]][[i+nFp]]),
                    as.double(0.5), size, PACKAGE="groHMM")
        Nrep <- NROW(ans$Start)
        # ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)),
        #           chromStart =ans$Start, chromEnd =ans$End,
        #           name =rep("N", Nrep), score =rep("1", Nrep),
        #           strand =rep("-", Nrep)))
        ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)),
                start=ans$Start, end=ans$End, strand =rep("-", Nrep)))
    }

    #BWem[[4]] <- ANS
    #names(BWem) <- c("EmisParams", "TransParams", "ViterbiStates",
    #                   "Transcripts")

    BWem[[4]] <- GRanges(seqnames = Rle(ANS$chrom),
                    ranges = IRanges(ANS$start, ANS$end-1),
                    strand = Rle(strand(ANS$strand)), type=Rle("tx",NROW(ANS)),
                    ID=paste(ANS$chrom, "_", ANS$start, ANS$strand, sep=""))
        names(BWem) <- c("emisParams", "transParams", "viterbiStates",
                            "transcripts")

    if(debug) {
        print(BWem[[1]])
        print(BWem[[2]])
    }

    return(BWem)
}


S_split =  import(input_f)
hmmResult_AllEM_split <- detectTranscripts_AllEM(S_split)
GRangeTobed(hmmResult_AllEM_split$transcripts, paste(strsplit(input_f, '.bed', T)[[1]][1], "_HMM.bed", sep = ""))

cat("Checking for output .bed file from R script... ")
cat(file.exists(paste(strsplit(input_f, '.bed', T)[[1]][1], "_HMM.bed", sep = "")), "\n")
