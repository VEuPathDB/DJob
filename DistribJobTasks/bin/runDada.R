library(dada2)

args <- commandArgs(TRUE)

dataDir <- args[1]
isPaired <- args[2]
platform <- args[3]
mergeTechReps <- args[4]

featureFile <- list.files(dataDir, pattern = "featureTable.rds", full.names = TRUE)

#if i am the samples used to build error model do nothing 
if (length(featureFile) == 0) {
  errModel <- list.files(dataDir, pattern="err.rds", full.names = TRUE)
  errModel <- readRDS(errModel)

  if (isPaired) {
    sampleF <- list.files(dataDir, pattern = "_R1_001.fastq", full.names = TRUE)
    sampleR <- list.files(dataDir, pattern = "_R2_001.fastq", full.names = TRUE)
    if (mergeTechReps) {
      #if true, takes sample name as everything before first point. 
      #this means technical replicates must have identical files names up to first point 
      sample.name <- sapply(strsplit(basename(sample), ".", fixed=TRUE), `[`, 1)
    } else { 
      sample.name <- sapply(strsplit(basename(sampleF), "_", fixed = TRUE), `[`, 1)
    }
    drpF <- derepFastq(sampleF)
    if (platform == "454") {
      ddF <- dada(drpF, err=errModel[[1]], multithread=1,
                  HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
      drpR <- derepFastq(sampleR)
      ddR <- dada(drpR, err=errModel[[2]], multithread=1,
                  HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
    } else {
      ddF <- dada(drpF, err=errModel[[1]], multithread=1)#,
                 # VECTORIZED_ALIGNMENT=FALSE)
      drpR <- derepFastq(sampleR)
      ddR <- dada(drpR, err=errModel[[2]], multithread=1)#,
                #VECTORIZED_ALIGNMENT=FALSE)
    }
    mergers <- mergePairs(ddF, drpF, ddR, drpR)
    #denoisedF <- getN(ddF)

    #names(mergers) <- sample.name

    seqtab <- makeSequenceTable(mergers)
    rownames(seqtab) <- sample.name
  } else {
    sample <- list.files(dataDir, pattern = ".fastq", full.names = TRUE)
    if (mergeTechReps) {
      sample.name <- sapply(strsplit(basename(sample), ".", fixed=TRUE), `[`, 1)
    } else {
      sample.name <- sapply(strsplit(basename(sample), "_filt.fastq"), `[`, 1)
    }
    drp <- derepFastq(sample)
    if (platform == "454") {
      dds <- dada(drp, err=errModel, multithread=1,
                  HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
    } else { 
      dds <- dada(drp, err=errModel, multithread=1)#,
                  #VECTORIZED_ALIGNMENT=FALSE)
    }
    seqtab <- makeSequenceTable(dds)
    rownames(seqtab) <- sample.name
  }
  
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                      minFoldParentOverAbundance=1,
                                      multithread=1)

  saveRDS(seqtab.nochim, file = file.path(dataDir, "featureTable.rds"))
}
