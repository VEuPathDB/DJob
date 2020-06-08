#!/usr/bin/env R
# This script is designed to run in parallel, once per sample: it takes about a minute for a typical fastq with 5000 reads
# Uses dada2::dada to denoise filtered amplicon reads
# Requires filtered reads and error models
# Produces an RDS file: a data frame of count per ASV

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(optparse))

opt = parse_args(OptionParser(option_list=list(
  make_option("--fastqsInDir", type="character", help="Dir with fastqs having right names"),
  make_option("--errorsRdsPath", type="character", help="Errors R object - for paired end, this is a two-element list"),
  make_option("--outRdsPath", type="character", help="Where to save the feature table to"),
  make_option("--isPaired", type="logical", help=""),
  make_option("--platform", type="character", help="sequencing platform"),
  make_option("--mergeTechReps", type="logical", help="Whether we're merging tech reps - works out file names based on this")
)))

callDada <- function(platform, derep, err){
  if(platform == "454"){
    dada(derep, err=err, multithread=1, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
  } else {
    dada(derep, err=err, multithread=1)  # VECTORIZED_ALIGNMENT=FALSE
  }
}


errModel <- readRDS(opt$errorsRdsPath)

if (opt$isPaired) {
  sampleF <- list.files(opt$fastqsInDir, pattern = "_1.fastq", full.names = TRUE)
  sampleR <- list.files(opt$fastqsInDir, pattern = "_2.fastq", full.names = TRUE)
  message("Input paths, forward reads: ", paste(sampleF, sep=", "), ", reverse reads: ", paste(sampleR, sep=", "))
  if (opt$mergeTechReps) {
    #if true, takes sample name as everything before first point. 
    #this means technical replicates must have identical files names up to first point 
    sample.name <- sapply(strsplit(basename(sampleF), ".", fixed=TRUE), `[`, 1)
  } else { 
    sample.name <- sapply(strsplit(basename(sampleF), "_", fixed = TRUE), `[`, 1)
  }
  drpF <- derepFastq(sampleF)
  ddF <- callDada(opt$platform, drpF, errModel[[1]])
  drpR <- derepFastq(sampleR)
  ddR <- callDada(opt$platform, drpR, errModel[[2]])

  mergers <- mergePairs(ddF, drpF, ddR, drpR)

  seqtab <- makeSequenceTable(mergers)
  rownames(seqtab) <- sample.name
} else {
  sample <- list.files(opt$fastqsInDir, pattern = ".fastq", full.names = TRUE)
  message("Input paths: ", paste(sample, sep=", "))
  if (opt$mergeTechReps) {
    sample.name <- sapply(strsplit(basename(sample), ".", fixed=TRUE), `[`, 1)
  } else {
    sample.name <- sapply(strsplit(basename(sample), "_filt.fastq"), `[`, 1)
  }
  drp <- derepFastq(sample)
  dds <- callDada(opt$platform, drp, errModel)

  seqtab <- makeSequenceTable(dds)
  rownames(seqtab) <- sample.name
}
if (length(seqtab)){
  message("Running removeBimeraDenovo on ", length(seqtab), " ASVs")
  seqtab <- removeBimeraDenovo(seqtab, method="consensus", minFoldParentOverAbundance=1, multithread=1)
}
if (length(seqtab)){
  message("Saving ", length(seqtab), " ASVs in RDS format to ", opt$outRdsPath)
  saveRDS(seqtab, file = opt$outRdsPath)
} else {
  warning("No ASVs. Will not write to ", opt$outRdsPath)
}
