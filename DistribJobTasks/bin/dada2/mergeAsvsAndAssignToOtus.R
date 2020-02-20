library(data.table)
library(dada2)
library(optparse)

opt = parse_args(OptionParser(option_list=list(
  make_option("--asvRdsInDir", type="character", help="Comma separated list of paths to RDS files with ASVs"),
  make_option("--assignTaxonomyRefPath", type="character", help="Fasta for assignTaxonomy"),
  make_option("--addSpeciesRefPath", type="character", help="Fasta for addSpecies"),
  make_option("--outPath", type="character", help="Output file", default = "final_taxonString.tab"),
  make_option("--chunkSize", type="integer", help="Batch size for assignTaxonomy+addSpecies. Limits memory usage, which also depends on ASV length and size of reference", default = 1500)
)))
asvFiles <- list.files(opt$asvRdsInDir, pattern = ".rds", full.names = TRUE)

if(length(asvFiles) > 1) {
# This used to error resultsAggregatedByLineage for repeats when we weren't merging technical replicates, but I got rid of the mergeTechReps flag here so now we always do repeats = "sum".
# We're not merging anything else here, anyway
  message("Merging ", length(asvFiles), " sequence tables")
  seqtab <- mergeSequenceTables(tables = asvFiles, repeats = "sum") 
} else {
  seqtab <- readRDS(asvFiles[1])
}
asvTable <- as.data.table(t(seqtab), keep.rownames=TRUE)
asvTable[is.na(asvTable)] <- 0
seqs <- asvTable$rn

totalChunks=ceiling(length(seqs)/opt$chunkSize)
taxaL <- vector("list", totalChunks)
for(chunkNum in seq(totalChunks)){
  i.lo <-  opt$chunkSize * (chunkNum - 1) + 1
  i.hi <-  min(opt$chunkSize * chunkNum + 1 , length(seqs))
  message("Assigning taxonomy from ", opt$assignTaxonomyRefPath, " for ASVs ", i.lo, ":",i.hi, " out of total ", length(seqs))
  taxa <- assignTaxonomy(seqs[i.lo:i.hi], opt$assignTaxonomyRefPath, tryRC = TRUE)
  message("Adding species from " , opt$addSpeciesRefPath)
  taxa <- addSpecies(taxa, opt$addSpeciesRefPath, tryRC = TRUE)
  taxa <- as.data.table(taxa, keep.rownames=TRUE)
  taxaL[[chunkNum]] <- taxa
}

# Join on the ASV strings
results <- merge(rbindlist(taxaL), asvTable, by = "rn")

message("Writing full results to ", paste0(opt$outPath, ".full"))

if(!dir.exists(dirname(opt$outPath))){
  dir.create(dirname(opt$outPath), recursive=TRUE)
}
write.table(results, file = paste0(opt$outPath, ".full"), quote=FALSE, sep = '\t', col.names = TRUE, row.names=FALSE, na="")

naphylumixs <- which(is.na(results$Phylum))
if (length(naphylumixs)){
  message("Removing ", length(naphylumixs) , " ASVs with no phylum (those were not assigned to anything in the taxonomy)")
  results <- results[-naphylumixs,]
}

# Zip by ;, removing nulls
# dada2 only returns common name for a species - fix that
formatLineage <- function(r){
  species<-r[length(r)]
  genus<-r[length(r)-1]
  lineageToGenus <- head(r, -1)
  s<- Filter(function(i){!is.na(i)},lineageToGenus)

  ifelse(
    is.na(species),
    paste0(s,collapse=";"),
    paste0(c(s, paste0(genus, " ", species)), collapse=";")
  )
}

# Label the rows by lineages, and aggregate by lineage
resultsAggregatedByLineage <-aggregate(
  results[,-c("rn", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")],
  by=list(Lineage = apply(results[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, formatLineage)),
  FUN=sum
)
message("Aggregated results into ", length(rownames(resultsAggregatedByLineage)), " lineages")
rownames(resultsAggregatedByLineage) = resultsAggregatedByLineage$Lineage
resultsAggregatedByLineage$Lineage <- NULL

rs <- rowSums(resultsAggregatedByLineage)
rsixs <- which (rs < 2)
if(length(rsixs)){
  message("Removing ", length(rsixs), " lineages with less than 2 total absolute abundance (across all samples, at most one read mapped to an ASV that was then assigned to this lineage - with only one piece of evidence, can we really claim that this lineage is present in our data ?) ")
  resultsAggregatedByLineage <- resultsAggregatedByLineage[-rsixs,]
}

cs <- colSums(resultsAggregatedByLineage)
csixs <- which(cs < 100)
if(length(csixs)){
  message("Removing  ", length(csixs), " samples with less than 100 total absolute abundance (samples from which less than a 100 reads remained after filtering, denoising, dropping, etc. - we deem them insufficient for meaningful comparisons to other samples)")
  resultsAggregatedByLineage <- resultsAggregatedByLineage[,-csixs]
}

message("Writing results to ", opt$outPath)
write.table(resultsAggregatedByLineage, file = opt$outPath, quote=FALSE, sep = '\t', col.names = NA, row.names=TRUE, na="")

