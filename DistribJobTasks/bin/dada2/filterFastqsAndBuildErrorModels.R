#! /usr/bin/env R
# This script is designed to be run once per dataset after data files have been provisioned
# Given a folder of amplicon sequencing reads, with one file per sample,
# filters the reads for quality and trim/truncate them as required
# Produces a folder of filtered fastqs
# Also produces error models for the next step of the workflow

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(optparse))

opt = parse_args(OptionParser(option_list=list(
  make_option("--fastqsInDir", type="character", help = "fastqsInDir"),  #Was: fastqsInDir <- args[1] 
  make_option("--fastqsOutDir", type="character", help = "fastqsOutDir"),  #Used to be baked in as fastqsInDir/filtered 
  make_option("--errorsOutDir", type="character", help = "errorsOutDir"),  #Was: errorsOutDir <- args[2] 
  make_option("--samplesInfo", type="character", help = "samplesInfo", default=NULL),  #Was: samplesInfo <- args[9] 
  make_option("--isPaired", type="logical", help = "isPaired"),  #Was: isPaired <- args[10] 
  make_option("--trimLeft", type="integer", help = "trimLeft", default=0),  #Was: trimLeft <- args[11] 
  make_option("--trimLeftR", type="integer", help = "trimLeftR", default=0),  #Was: trimLeftR <- args[12] 
  make_option("--truncLen", type="integer", help = "truncLen - will try to figure it out if not provided", default=NULL),  #Was: truncLen <- args[13] 
  make_option("--truncLenR", type="integer", help = "truncLenR - will try to figure it out if not provided. Danielle said it doesn't work as well for paired reads", default=NULL),  #Was: truncLenR <- args[14] 
  make_option("--readLen", type="integer", help = "readLen - avg for illumina and max allowed for 454"),  #Was: readLen <- as.numeric(args[15]) 
  make_option("--platform", type="character", help = "platform"),  #Was: platform <- args[16] 
  make_option("--mergeTechReps", type="logical", help = "mergeTechReps")  #Was: mergeTechReps <- args[17] 
)))
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!# ARGUMENTS #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

if (is.null(opt$samplesInfo)) {
  samples.info <- NULL
} else {
  samples.info <- fread(opt$samplesInfo)
}

#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!# FUNCTIONS !#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

ilSingleErr <- function(inputFiles, filesOutDir, errFile = NULL, truncLen = NULL, trimLeft = NULL, platform = NULL, readLen = NULL, mergeTechReps) {

             message("ilSingleErr processing ", length(inputFiles), " input files")
             sample.names <- sapply(strsplit(basename(inputFiles), ".fastq"), `[`, 1)

             filts <- file.path(filesOutDir, paste0(sample.names, "_filt.fastq"))

             #unlink(filesOutDir, recursive = TRUE)
             ## TRIM AND FILTER ###
             if (platform == "454") {
               out <- suppressWarnings(filterAndTrim(inputFiles, filts, truncLen=truncLen, trimLeft=trimLeft,
                                       maxEE=2, truncQ=2, rm.phix=TRUE, multithread=1, 
                                       maxLen=readLen, verbose=TRUE))
             } else {
               out <- suppressWarnings(filterAndTrim(inputFiles, filts, truncLen=truncLen, trimLeft=trimLeft,
                                       maxEE=2, truncQ=2, rm.phix=TRUE, multithread=1,
                                       verbose=TRUE))
             }
             filts <- filts[filts %in% list.files(filesOutDir, pattern=".fastq", full.names=TRUE)]
             if(length(filts) == 0) {
               stop("No reads passed the filter (was truncLen longer than the read length?)")
             } else {
               message(length(filts), " read files passed the filter")
             }

           ### LEARN ERROR RATES ###
             # Dereplicate enough samples to get 1 million total unique reads to model errors on
             NREADS <- 0
             drps <- vector("list", length(filts))
             for(i in seq_along(filts)) {
               drps[[i]] <- derepFastq(filts[[i]])
               NREADS <- NREADS + sum(drps[[i]]$uniques)
               if(NREADS > 1000000) { break }
             }
             # Run dada in self-consist mode on those samples
             dds <- vector("list", length(filts))
            
             if (platform == "454") {
               if(i==1) { # breaks list assignment
                 dds[[1]] <- dada(drps[[1]], err=NULL, selfConsist=TRUE, multithread=1,
                                  HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
               } else { # more than one sample, no problem with list assignment
                 dds[1:i] <- dada(drps[1:i], err=NULL, selfConsist=TRUE, multithread=1,
                                  HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
               }
             } else {  
               if(i==1) { # breaks list assignment
                 dds[[1]] <- dada(drps[[1]], err=NULL, selfConsist=TRUE, multithread=1)#,
                           #       VECTORIZED_ALIGNMENT=FALSE)
               } else { # more than one sample, no problem with list assignment
                 dds[1:i] <- dada(drps[1:i], err=NULL, selfConsist=TRUE, multithread=1)#,
                 #VECTORIZED_ALIGNMENT=FALSE)
               }
             }
             #message(dds)
             err <- dds[[1]]$err_out
             message("Saving RDS to ", errFile)
             saveRDS(err,errFile)
}

ilPairedErr <- function(inputFilesF, inputFilesR, filesOutDir, errFile = NULL, truncLen = NULL, trimLeft = NULL,
                     truncLenR = NULL, trimLeftR = NULL, platform = NULL, readLen = NULL, mergeTechReps) {
              if (length(inputFilesF) != length(inputFilesR)){ 
                stop("ilPairedErr found an unequal number of input files: ", length(inputFilesF), " forward and ", length(inputFilesR), " reverse")
              } else {
                message("ilPairedErr processing ", length(inputFilesF), " input files")
              }
              sample.names <- sapply(strsplit(basename(inputFilesF), "_"), `[`, 1)
              filtsF <- file.path(filesOutDir, basename(inputFilesF))
              filtsR <- file.path(filesOutDir, basename(inputFilesR))

              ### TRIM AND FILTER ###
              message("filter and trimming...")
              if (platform == "454") {
                out <- suppressWarnings(filterAndTrim(inputFilesF, filtsF, inputFilesR, filtsR,
                                        truncLen=c(truncLen, truncLenR), trimLeft=c(trimLeft, trimLeftR),
                                        maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                        multithread=1, maxLen = readLen, verbose=TRUE))
              } else {
                out <- suppressWarnings(filterAndTrim(inputFilesF, filtsF, inputFilesR, filtsR,
                                        truncLen=c(truncLen, truncLenR), trimLeft=c(trimLeft, trimLeftR),
                                        maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                        multithread=1, verbose=TRUE))
              }
              filtsF <- Filter(file.exists, filtsF)
              filtsR <- Filter(file.exists, filtsR)
              if(length(filtsF) == 0) { # All reads were filtered out
                stop("No reads passed the filter (were truncLenF/R longer than the read lengths?)")
              }
             message("learning error rates...")
              ### LEARN ERROR RATES ###
              # Dereplicate enough samples to get 1 million total reads
              NREADS <- 0
              drpsF <- vector("list", length(filtsF))
              drpsR <- vector("list", length(filtsR))
              denoisedF <- rep(0, length(filtsF))
              getN <- function(x) sum(getUniques(x))
              for(i in seq_along(filtsF)) {
                drpsF[[i]] <- derepFastq(filtsF[[i]])
                drpsR[[i]] <- derepFastq(filtsR[[i]])
                NREADS <- NREADS + sum(drpsF[[i]]$uniques)
                if(NREADS > 1000000) { break }
              }
              # Run dada in self-consist mode on those samples
              drpsF <- drpsF[1:i]
              drpsR <- drpsR[1:i]
              if (platform == "454") {
                ddsF <- dada(drpsF, err=NULL, selfConsist=TRUE, multithread=1,
                             HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
                ddsR <- dada(drpsR, err=NULL, selfConsist=TRUE, multithread=1,
                             HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
              } else {
                ddsF <- dada(drpsF, err=NULL, selfConsist=TRUE, multithread=1)#,
                             #VECTORIZED_ALIGNMENT=FALSE)
                ddsR <- dada(drpsR, err=NULL, selfConsist=TRUE, multithread=1)#,
                             #VECTORIZED_ALIGNMENT=FALSE)
              }
              if(i==1) {
                errF <- ddsF$err_out
                errR <- ddsR$err_out
              } else {
                errF <- ddsF[[1]]$err_out
                errR <- ddsR[[1]]$err_out
              }
              err <- list(errF,errR)
              message("Saving RDS to ", errFile)
              saveRDS(err,errFile)
}

".getBpParam" <- function(mc.cores=1L){
  return(switch(as.character(mc.cores),
                "1"=SerialParam(),
                MulticoreParam(workers=mc.cores)))
}


findTruncLen = function(forwardReads = NULL, reverseReads = NULL, sample_size = 10000, readLen = NULL, expected_errors = 2, overlap = 20) {
 
  if (!is.null(reverseReads)) {
    lengths <- data.frame("forward"=c(NA),"reverse"=c(NA))
  } else {
    lengths <- data.frame("forward"=c(NA))
  } 

  numSamples <- length(forwardReads)
  if (numSamples < 50) {
    numSamplesToUse <- numSamples
  } else {
    numSamplesToUse <- floor(numSamples/10)
  }

  for (x in sample(x=seq(numSamples),numSamplesToUse)) {
 
    fqFile <- forwardReads[x]
    reads <- readFastq(fqFile)  

    if(sample_size > length(sread(reads))) {
      sample_size <- length(sread(reads))
    }
    
    ## Extract the quality on random subset and calculate the expected error for each sequence
    qual <- as(quality(reads[sample(x = seq(1,length(sread(reads))), size = sample_size, replace = F)]), "matrix")
    qual[is.na(qual)] <- 0 
    error <- 10^(-1*qual/10)
    summed <- apply(error, 1, cumsum)
      
    ## Select reads with less expected errors than the threshold
    trunc_f <- apply(summed, 2, function(x) max(which(x < expected_errors)))

    ## Create the percentile distribution of sequence length that satisfies the threshold
    f_quant <- quantile(trunc_f, probs=c(seq(0,1,0.01)))

    ## The same for reverse if its not NULL
    if (!is.null(reverseReads)) {
      fqFile <- reverseReads[x]
      reads <- readFastq(fqFile)
      qual <- as(quality(reads[sample(x = seq(1,length(sread(reads))), size = sample_size, replace = F)]), "matrix")
      qual[is.na(qual)] <- 0
      error <- 10^(-1*qual/10)
      summed <- apply(error, 1, cumsum)
      trunc_r <- apply(summed, 2, function(x) max(which(x < expected_errors)))
      r_quant <- quantile(trunc_r, probs=c(seq(0,1,0.01)))
 
      i <- 1

      # Goes through the percentile distribution until the forward length and reverse sums more the amplicon length and the minimum required overlap
      while(f_quant[i] + r_quant[i] < readLen + overlap) {
        i <- i+1
      }

      lengths[x,] <- c(f_quant[i], r_quant[i])
      rm(i)
    } else {
      lengths[x,] <- f_quant[1]
    }
  }
  truncLen <- floor(mean(lengths[,1],na.rm=TRUE))
  if (!is.null(reverseReads)) {
    truncLenR <- floor(mean(lengths[,2],na.rm=TRUE))
    return(c(truncLen,truncLenR))
  } else {
    return(truncLen)
  }
}

#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#! WORKFLOW #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
if (is.null(opt$truncLen)){
  message("Determining truncLen... ")
  if (opt$isPaired) {
    fReads <- list.files(opt$fastqsInDir, pattern = "_1.fastq", full.names = TRUE)
    rReads <- list.files(opt$fastqsInDir, pattern = "_2.fastq", full.names = TRUE)
    truncLens <- findTruncLen(forwardReads = fReads, reverseReads = rReads, readLen = opt$readLen)
    opt$truncLen <- truncLens[1]
    opt$truncLenR <- truncLens[2]
  } else {
    reads <- list.files(opt$fastqsInDir, pattern = ".fastq", full.names = TRUE)
    opt$truncLen <- findTruncLen(forwardReads = reads, readLen = opt$readLen)
  }

  #need 50 nts to assign taxonomy so..
  if (opt$truncLen < (opt$trimLeft + 50)) {
    opt$truncLen <- opt$trimLeft + 50
  }
  message("Determined truncLen = " , opt$truncLen)
} 

#run dada2 in selfConsist mode to build error model
if (!is.null(samples.info) && "GROUPS" %in% colnames(samples.info)) {
  groups <- unique(samples.info$GROUPS)
  for (group in groups) {
    myFiles <- samples.info$NAMES[samples.info$GROUPS == group]
    errFileName <- paste0(group, "_err.rds")
    errFile <- file.path(opt$errorsOutDir, errFileName)
    unlink(errFile)
    message("Filtering fastqs and preparing the error file: ", errFile, " for group: ", group)
    if (opt$isPaired) {
      filesF <- file.path(opt$fastqsInDir, paste0(myFiles, "_1.fastq"))
      filesR <- file.path(opt$fastqsInDir, paste0(myFiles, "_2.fastq"))
      ilPairedErr(filesF, filesR, opt$fastqsOutDir, errFile, truncLen = opt$truncLen, trimLeft = opt$trimLeft,
                                         truncLenR = opt$truncLenR, trimLeftR = opt$trimLeftR, platform = opt$platform,
                                         readLen = opt$readLen, mergeTechReps= opt$mergeTechReps)
    } else {
      myFiles <- file.path(opt$fastqsInDir, paste0(myFiles, ".fastq"))
      ilSingleErr(myFiles, opt$fastqsOutDir, errFile, truncLen = opt$truncLen, trimLeft = opt$trimLeft, 
                                         platform = opt$platform, readLen = opt$readLen, mergeTechReps=opt$mergeTechReps)
    }
  }
} else {
  errFileName = "err.rds"
  errFile <- file.path(opt$errorsOutDir, errFileName)
  unlink(errFile)
  message("Filtering fastqs and preparing the error file: ", errFile)
  if (opt$isPaired) {
    filesF <- sort(list.files(opt$fastqsInDir, pattern="_1.fastq", full.names = TRUE))
    filesR <- sort(list.files(opt$fastqsInDir, pattern="_2.fastq", full.names = TRUE))
    ilPairedErr(filesF, filesR, opt$fastqsOutDir, errFile, truncLen = opt$truncLen, trimLeft = opt$trimLeft,
                                       truncLenR = opt$truncLenR, trimLeftR = opt$trimLeftR, platform = opt$platform,
                                       readLen = opt$readLen, mergeTechReps= opt$mergeTechReps)
  } else {
    myFiles <- sort(list.files(opt$fastqsInDir, pattern=".fastq", full.names = TRUE))
    ilSingleErr(myFiles, opt$fastqsOutDir, errFile, truncLen = opt$truncLen, trimLeft = opt$trimLeft, 
                                       platform = opt$platform, readLen = opt$readLen, mergeTechReps=opt$mergeTechReps)
  }
}

