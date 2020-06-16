
library(dada2)

path<-"../MiSeq_SOP/"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnRs
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

fnFs
strsplit(basename(fnFs),"_")
sapply(strsplit(basename(fnFs), "_"), `[`, 1)


plotQualityProfile(fl =fnRs,aggregate = T)


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))



filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


head(out)


errF <- learnErrors(filtFs, multithread=TRUE,nbases = 3e8)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)



plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


table(nchar(getSequences(seqtab)))
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track2<-track/track[,1]*100

taxa <- assignTaxonomy(seqtab.nochim, "~/projects/databases/silva_nr_v132_train_set.fa.gz", multithread=TRUE)





dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs


load("~/projects/databases/SILVA_SSU_r132_March2018.rdata") # CHANGE TO THE PATH OF YOUR TRAINING SET


ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy


taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; 
rownames(taxid) <- getSequences(seqtab.nochim)


saveRDS(object = taxid,file = "taxid.rds")
saveRDS(object = taxa,file = "taxa.rds")
saveRDS(object = seqtab.nochim,file = "seqtab_nochim.rds")

