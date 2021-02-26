#! /usr/bin/env R
# This script is designed to be run once per dataset after data files have been provisioned
# Given a folder of amplicon sequencing reads, with one file per sample,
# filters the reads for quality and trim/truncate them as required
# Produces a folder of filtered fastqs

suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(optparse))

opt = parse_args(OptionParser(option_list=list(
  make_option("--fastqsInDir", type="character", help = "fastqsInDir"),
  make_option("--fastqsOutDir", type="character", help = "fastqsOutDir"),
  make_option("--samplesInfo", type="character", help = "samplesInfo", default=NULL),
  make_option("--isPaired", type="logical", help = "isPaired"),
  make_option("--trimLeft", type="integer", help = "trimLeft", default=0),
  make_option("--trimLeftR", type="integer", help = "trimLeftR", default=0),
  make_option("--truncLen", type="integer", help = "truncLen"),
  make_option("--truncLenR", type="integer", help = "truncLenR - required for paired reads", default=NULL),
  make_option("--maxLen", type="integer", help = "max allowed - required for 454", default=NULL),
  make_option("--platform", type="character", help = "platform")
)))
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!# ARGUMENTS #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

if (is.null(opt$samplesInfo)) {
  samples.info <- NULL
} else {
  samples.info <- fread(opt$samplesInfo)
}

#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!# FUNCTIONS !#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

ilSingleErr <- function(inputFiles, filesOutDir, truncLen = NULL, trimLeft = NULL, platform = NULL, maxLen = NULL) {

             message("ilSingleErr processing ", length(inputFiles), " input files")
             sample.names <- sapply(strsplit(basename(inputFiles), ".fastq"), `[`, 1)

             filts <- file.path(filesOutDir, paste0(sample.names, "_filt.fastq"))

             #unlink(filesOutDir, recursive = TRUE)
             ## TRIM AND FILTER ###
             if (platform == "454") {
               out <- suppressWarnings(filterAndTrim(inputFiles, filts, truncLen=truncLen, trimLeft=trimLeft,
                                       maxEE=2, truncQ=2, rm.phix=TRUE, multithread=1, 
                                       maxLen=maxLen, verbose=TRUE))
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
}

ilPairedErr <- function(inputFilesF, inputFilesR, filesOutDir, truncLen = NULL, trimLeft = NULL,
                     truncLenR = NULL, trimLeftR = NULL, platform = NULL, maxLen = NULL) {
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
                                        multithread=1, maxLen = maxLen, verbose=TRUE))
              } else {
                out <- suppressWarnings(filterAndTrim(inputFilesF, filtsF, inputFilesR, filtsR,
                                        truncLen=c(truncLen, truncLenR), trimLeft=c(trimLeft, trimLeftR),
                                        maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                        multithread=1, verbose=TRUE))
              }
}

# Wojtek: Not sure what that is - monkeypatching something?
".getBpParam" <- function(mc.cores=1L){
  return(switch(as.character(mc.cores),
                "1"=SerialParam(),
                MulticoreParam(workers=mc.cores)))
}



#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#! WORKFLOW #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

if (!is.null(samples.info) && "GROUPS" %in% colnames(samples.info)) {
  groups <- unique(samples.info$GROUPS)
  for (group in groups) {
    myFiles <- samples.info$NAMES[samples.info$GROUPS == group]
    message("Filtering fastqs for group: ", group)
    if (opt$isPaired) {
      filesF <- file.path(opt$fastqsInDir, paste0(myFiles, "_1.fastq"))
      filesR <- file.path(opt$fastqsInDir, paste0(myFiles, "_2.fastq"))
      ilPairedErr(filesF, filesR, opt$fastqsOutDir, truncLen = opt$truncLen, trimLeft = opt$trimLeft,
                                         truncLenR = opt$truncLenR, trimLeftR = opt$trimLeftR, platform = opt$platform,
                                         maxLen = opt$maxLen)
    } else {
      myFiles <- file.path(opt$fastqsInDir, paste0(myFiles, ".fastq"))
      ilSingleErr(myFiles, opt$fastqsOutDir, truncLen = opt$truncLen, trimLeft = opt$trimLeft, 
                                         platform = opt$platform, maxLen = opt$maxLen)
    }
  }
} else {
  message("Filtering fastqs")
  if (opt$isPaired) {
    filesF <- sort(list.files(opt$fastqsInDir, pattern="_1.fastq", full.names = TRUE))
    filesR <- sort(list.files(opt$fastqsInDir, pattern="_2.fastq", full.names = TRUE))
    ilPairedErr(filesF, filesR, opt$fastqsOutDir, truncLen = opt$truncLen, trimLeft = opt$trimLeft,
                                       truncLenR = opt$truncLenR, trimLeftR = opt$trimLeftR, platform = opt$platform,
                                       maxLen = opt$maxLen)
  } else {
    myFiles <- sort(list.files(opt$fastqsInDir, pattern=".fastq", full.names = TRUE))
    ilSingleErr(myFiles, opt$fastqsOutDir, truncLen = opt$truncLen, trimLeft = opt$trimLeft, 
                                       platform = opt$platform, maxLen = opt$maxLen)
  }
}

