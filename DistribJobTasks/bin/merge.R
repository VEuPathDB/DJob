library(data.table)
library(dada2)

args <- commandArgs(TRUE)

filesDir <- args[1]
taxonRefFile <- args[2]
mergeTechReps <- args[3]

myFiles <- list.files(filesDir, pattern = ".rds", full.names = TRUE)
if (length(myFiles) > 1) {

  #no longer in use !!!!!!
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

  #asvTable <- Reduce(myMerge, myFiles)

  if (mergeTechReps) {
    seqtab <- mergeSequenceTables(tables = myFiles, repeats = "sum")
  } else {
    seqtab <- mergeSequenceTables(tables = myFiles)
  }

  asvTable <- as.data.table(t(seqtab), keep.rownames=TRUE)
  asvTable[is.na(asvTable)] <- 0

  file.remove(myFiles)
} else {
  message("only one rds file... renaming to final_featureTable.tab.")
  seqtab <- readRDS(myFiles[1])
  asvTable <- as.data.table(t(seqtab), keep.rownames=TRUE)
  asvTable[is.na(asvTable)] <- 0
  file.remove(myFiles[1])
}

seqs <- asvTable$rn
taxa <- assignTaxonomy(seqs, taxonRefFile)
taxa <- as.data.table(taxa, keep.rownames=TRUE)
colnames(taxa)[colnames(taxa) == "rn"] <- "V1"

finalTable <- merge(taxa, asvTable, by = "V1")
finalTable$V1 <- NULL

#remove high level NA values
finalTable <- finalTable[-which(is.na(finalTable$Phylum)),]
#fix NA values
finalTable$Class[is.na(finalTable$Class)] <- "c__"
finalTable$Order[is.na(finalTable$Order)] <- "o__"
finalTable$Family[is.na(finalTable$Family)] <- "f__"
finalTable$Genus[is.na(finalTable$Genus)] <- "g__"
finalTable$Species[is.na(finalTable$Species)] <- "s__"

finalTable$taxaString <- paste0(finalTable$Kingdom, "; ", finalTable$Phylum, "; ", finalTable$Class, "; ", finalTable$Order, "; ", finalTable$Family, "; ", finalTable$Genus, "; ", finalTable$Species)
finalTable[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") := NULL]
#moves taxaString col from last to first
setcolorder(finalTable, c("taxaString", colnames(finalTable)[1:ncol(finalTable)-1]))

#write taxonString to seperate file as well, to be consumed by insertSequenceTaxon plugin
write.table(finalTable$taxaString, file = file.path(filesDir, "final_taxonString.tab"), quote=FALSE, sep = '\t', col.names = FALSE, row.names=FALSE)

write.table(finalTable, file = file.path(filesDir, "final_featureTable.tab"), quote=FALSE, sep = '\t', col.names = TRUE, row.names=FALSE)
