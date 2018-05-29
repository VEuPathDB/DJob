library(dada2)

args <- commandArgs(TRUE)

dataDir <- args[1]
isPaired <- args[2]

tabFile <- list.files(dataDir, pattern = ".tab", full.names = TRUE)

#if i am a tab file do nothing 
if (length(tabFile) == 0) {
  errModel <- list.files(dataDir, pattern="err.rds", full.names = TRUE)
  errModel <- readRDS(errModel)

  if (isPaired) {
    sampleF <- list.files(dataDir, pattern = "_R1_001.fastq", full.names = TRUE)
    sampleR <- list.files(dataDir, pattern = "_R2_001.fastq", full.names = TRUE)
    sample.name <- sapply(strsplit(basename(sampleF), "_filt.fastq"), `[`, 1)
    drpF <- derepFastq(sampleF)
    ddF <- dada(drpF, err=errModel[1], multithread=1)#,
               # VECTORIZED_ALIGNMENT=FALSE)
    drpR <- derepFastq(sampleR)
    ddR <- dada(drpR, err=errModel[2], multithread=1)#,
                #VECTORIZED_ALIGNMENT=FALSE)
    mergers <- mergePairs(ddF, drpF, ddR, drpR)
    denoisedF <- getN(ddF)

    names(mergers) <- sample.name

    seqtab <- makeSequenceTable(mergers)
  } else {
    sample <- list.files(dataDir, pattern = ".fastq", full.names = TRUE)
    sample.name <- sapply(strsplit(basename(sample), "_filt.fastq"), `[`, 1)
    drp <- derepFastq(sample)
    dds <- dada(drp, err=errModel, multithread=1)#,
                #VECTORIZED_ALIGNMENT=FALSE)
  
    seqtab <- makeSequenceTable(dds)
    rownames(seqtab) <- sample.name
  }
  
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                      minFoldParentOverAbundance=1,
                                      multithread=1)

  asvTable <- as.data.frame(t(seqtab.nochim))

  write.table(asvTable, file = file.path(dataDir, "featureTable.tab"), quote=FALSE, sep='\t', col.names = NA)
}
