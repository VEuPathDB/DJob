#! /usr/bin/env R
# This script is designed to be run once per dataset after data files have been filtered
# Given a folder of amplicon sequencing reads, with one file per sample,
# produces error models for the next step of the workflow

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(optparse))

opt = parse_args(OptionParser(option_list=list(
  make_option("--fastqsInDir", type="character", help = "fastqsInDir"),
  make_option("--errorsOutDir", type="character", help = "errorsOutDir"),
  make_option("--errorsFileNameSuffix", type="character", help = "file name, or if in groups, it does ${group}_${suffix}"),
  make_option("--samplesInfo", type="character", help = "samplesInfo", default=NULL),
  make_option("--isPaired", type="logical", help = "isPaired"),
  make_option("--platform", type="character", help = "platform"),
  make_option("--verbose", type="integer", help = "0,1,2. Default 1", default=1)
)))
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!# ARGUMENTS #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

if (is.null(opt$samplesInfo)) {
  samples.info <- NULL
} else {
  samples.info <- fread(opt$samplesInfo)
}

#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!# FUNCTIONS !#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

errorModelsSingle <- function(inputFiles, errFile, platform, verbose) {
           ### LEARN ERROR RATES ###
             # Dereplicate enough samples to get 1 million total unique reads to model errors on
             NREADS <- 0
             drps <- vector("list", length(inputFiles))
             for(i in seq_along(inputFiles)) {
               message("Dereplicating fastq: ", inputFiles[[i]])
               drps[[i]] <- derepFastq(inputFiles[[i]])
               NREADS <- NREADS + sum(drps[[i]]$uniques)
               message("Total reads dereplicated so far: ", NREADS)
               if(NREADS > 1000000) { break }
             }
             # Run dada in self-consist mode on those samples
             dds <- vector("list", length(inputFiles))
            
             message("Running dada in selfConsist mode")
             if (platform == "454") {
               if(i==1) { # breaks list assignment
                 dds[[1]] <- dada(drps[[1]], err=NULL, selfConsist=TRUE, MAX_CLUST=1000, verbose=verbose, multithread=1,
                                  HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
               } else { # more than one sample, no problem with list assignment
                 dds[1:i] <- dada(drps[1:i], err=NULL, selfConsist=TRUE, MAX_CLUST=1000, verbose=verbose, multithread=1,
                                  HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
               }
             } else {  
               if(i==1) { # breaks list assignment
                 dds[[1]] <- dada(drps[[1]], err=NULL, selfConsist=TRUE, MAX_CLUST=1000, verbose=verbose, multithread=1)#,
                           #       VECTORIZED_ALIGNMENT=FALSE)
               } else { # more than one sample, no problem with list assignment
                 dds[1:i] <- dada(drps[1:i], err=NULL, selfConsist=TRUE, MAX_CLUST=1000, verbose=verbose, multithread=1)#,
                 #VECTORIZED_ALIGNMENT=FALSE)
               }
             }
             #message(dds)
             err <- dds[[1]]$err_out
             message("Saving RDS to ", errFile)
             saveRDS(err,errFile)
}


errorModelsPaired <- function(inputFilesF, inputFilesR, fastqsInDir, errFile, platform, verbose) {
              if (length(inputFilesF) != length(inputFilesR)){ 
                stop("errorModelsPaired found an unequal number of input files: ", length(inputFilesF), " forward and ", length(inputFilesR), " reverse")
              } 
              ### LEARN ERROR RATES ###
              # Dereplicate enough samples to get 1 million total reads
              NREADS <- 0
              drpsF <- vector("list", length(inputFilesF))
              drpsR <- vector("list", length(inputFilesR))
              denoisedF <- rep(0, length(inputFilesF))
              getN <- function(x) sum(getUniques(x))
              for(i in seq_along(inputFilesF)) {
                message("Dereplicating fastqs: ", inputFilesF[[i]], inputFilesR[[i]])
                drpsF[[i]] <- derepFastq(inputFilesF[[i]])
                drpsR[[i]] <- derepFastq(inputFilesR[[i]])
                NREADS <- NREADS + sum(drpsF[[i]]$uniques)
                message("Total reads dereplicated so far: ", NREADS)
                if(NREADS > 1000000) { break }
              }
              # Run dada in self-consist mode on those samples
              drpsF <- drpsF[1:i]
              drpsR <- drpsR[1:i]
              message("Running dada in selfConsist mode")
              if (platform == "454") {
                ddsF <- dada(drpsF, err=NULL, selfConsist=TRUE, MAX_CLUST=1000, verbose=verbose, multithread=1,
                             HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
                ddsR <- dada(drpsR, err=NULL, selfConsist=TRUE, MAX_CLUST=1000, verbose=verbose, multithread=1,
                             HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
              } else {
                ddsF <- dada(drpsF, err=NULL, selfConsist=TRUE, MAX_CLUST=1000, verbose=verbose, multithread=1)#,
                             #VECTORIZED_ALIGNMENT=FALSE)
                ddsR <- dada(drpsR, err=NULL, selfConsist=TRUE, MAX_CLUST=1000, verbose=verbose, multithread=1)#,
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



#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#! WORKFLOW #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

#run dada2 in selfConsist mode to build error model
if (!is.null(samples.info) && "GROUPS" %in% colnames(samples.info)) {
  groups <- unique(samples.info$GROUPS)
  for (group in groups) {
    myFiles <- samples.info$NAMES[samples.info$GROUPS == group]
    errFileName <- paste0(group, "_", opt$errorsFileNameSuffix)
    errFile <- file.path(opt$errorsOutDir, errFileName)
    unlink(errFile)
    message("Preparing the error file: ", errFile, " for group: ", group)
    if (opt$isPaired) {
      filesF <- file.path(opt$fastqsInDir, paste0(myFiles, "_1.fastq"))
      filesR <- file.path(opt$fastqsInDir, paste0(myFiles, "_2.fastq"))
      errorModelsPaired(filesF, filesR, errFile, platform = opt$platform, verbose = opt$verbose)
    } else {
      myFiles <- file.path(opt$fastqsInDir, paste0(myFiles, ".fastq"))
      errorModelsSingle(myFiles, errFile, platform = opt$platform, verbose = opt$verbose)
    }
  }
} else {
  errFile <- file.path(opt$errorsOutDir, opt$errorsFileNameSuffix)
  unlink(errFile)
  message("Preparing the error file: ", errFile)
  if (opt$isPaired) {
    filesF <- sort(list.files(opt$fastqsInDir, pattern="_1.fastq", full.names = TRUE))
    filesR <- sort(list.files(opt$fastqsInDir, pattern="_2.fastq", full.names = TRUE))
    errorModelsPaired(filesF, filesR, errFile, platform = opt$platform, verbose = opt$verbose)
  } else {
    myFiles <- sort(list.files(opt$fastqsInDir, pattern=".fastq", full.names = TRUE))
    errorModelsSingle(myFiles, errFile, platform = opt$platform, verbose = opt$verbose)
  }
}

