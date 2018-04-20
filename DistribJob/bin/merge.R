library(data.table)
library(dada2)

args <- commandArgs(TRUE)

filesDir <- args[1]
taxonRefFile <- args[2]

myFiles <- list.files(filesDir, pattern = ".tab" , full.names = TRUE)
if (length(myFiles) > 1) {
  myMerge <- function(x,y) {
    #if not a data.table then assume its a string containing a path to a file
    if (!is.data.table(x)) {
      x <- fread(x)
    }
    if (!is.data.table(y)) {
      y <- fread(y)
    }
    out <- merge(x,y, all = TRUE, by = "V1")
  }

  asvTable <- Reduce(myMerge, myFiles)

  asvTable[is.na(asvTable)] <- 0

  file.remove(myFiles)
} else {
  message("only one tab file... renaming to final_featureTable.tab.")
  asvTable <- fread(myFiles[1])
  asvTable[is.na(asvTable)] <- 0
  file.remove(myFiles[1])
}

seqs <- asvTable$V1
rownames(asvTable) <- asvTable$V1
asvTable$V1 <- NULL
taxa <- assignTaxonomy(seqs, taxonRefFile)

#uncomment if we dont want to actually write the seqs themselves, but want unique identifiers
#rownames(taxa) <- unname(sapply(rownames(taxa), digest, algo = "md5", serialize = FALSE))
#rownames(asvTable) <- unname(sapply(rownames(asvTable), digest, algo = "md5", serialize = FALSE))

write.table(taxa, file = file.path(filesDir, "final_taxaMap.tab"), quote=FALSE, sep = '\t', col.names = NA)

write.table(asvTable, file = file.path(filesDir, "final_featureTable.tab"), quote=FALSE, sep = '\t', col.names = NA)
