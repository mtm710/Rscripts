# Binatur analysis 16S
# R script for DADA2 analysis of DNA metabarcoding data using the 16S primer set

# modified from https://benjjneb.github.io/dada2/tutorial.html

library(dada2); packageVersion("dada2")

# path to each plate (0,1,2,3)

path<- "~/Raw_NGS_Data/binatur/16S_data/Plate0_cut_16S_20221017_mtm/"
path<- "~/Raw_NGS_Data/binatur/16S_data/Plate1_cut_16S_132-143_20220909_mtm/"
path<- "~/Raw_NGS_Data/binatur/16S_data/Plate2_cut_137-158_20221017_hanna_trimmed/"
path <- "~/Raw_NGS_Data/binatur/16S_data/Plate3_cut_16S_158-170_20230209_mtm/" 

#note that these reads were trimmed with cutadapt and filenames are not standard (R1 instead of R001 and fq.gz instead of fastq.gz)

list.files(path) 

fnFs <- sort(list.files(path, pattern="_R1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

plotQualityProfile(fnFs[1:12])

plotQualityProfile(fnRs[1:12])

#filter and trim (trimLeft at 0,0 because cutadapt was used)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(0,0), truncLen=c(250,150), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)

#error rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

# merge paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# sequence table - change names for each plate

seqtab_plate0 <- makeSequenceTable(mergers)
dim(seqtab_plate0)
table(nchar(getSequences(seqtab_plate0)))

seqtab_plate1 <- makeSequenceTable(mergers)
dim(seqtab_plate1)
table(nchar(getSequences(seqtab_plate1)))

seqtab_plate2 <- makeSequenceTable(mergers)
dim(seqtab_plate2)
table(nchar(getSequences(seqtab_plate2)))

seqtab_plate3 <- makeSequenceTable(mergers)
dim(seqtab_plate3)
table(nchar(getSequences(seqtab_plate3)))

# check tables
write.csv(seqtab_plate0, "seqtab_plate0.csv", quote=FALSE)
write.csv(seqtab_plate1, "seqtab_plate1.csv", quote=FALSE)
write.csv(seqtab_plate2, "seqtab_plate2.csv", quote=FALSE)
write.csv(seqtab_plate3, "seqtab_plate3.csv", quote=FALSE)

# merge sequence tables (0,1,2,3)
mergetab0123 <- mergeSequenceTables(seqtab_plate0, seqtab_plate1, seqtab_plate2, seqtab_plate3, repeats = "sum")

# check merged table
write.csv(mergetab0123, "mergetab0123.csv", quote=FALSE)

# remove chimeras from the merged table
seqtab.nochim <- removeBimeraDenovo(mergetab0123, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# collapse no mismatch
# collapsed.nochim <- collapseNoMismatch(seqtab.nochim, identicalOnly = TRUE)

# write ASV table - transposes table
write.table(t(seqtab.nochim), "seqtab_nochim.tsv")

#assign taxonomy locally with rdp trainset
saveRDS(seqtab.nochim, "seqtab.nochim.RDS")
taxa <- assignTaxonomy(seqtab.nochim, "~/scratch/binatur_test/rdp_train_set_18.fa.gz")
saveRDS(taxa, "taxa.RDS")
write.csv(cbind(t(seqtab.nochim), taxa), "seqtab_nochim_taxa_id.csv", quote=FALSE)

#assign species locally
taxaspecies <- assignSpecies(seqtab.nochim, "~/scratch/binatur_test/rdp_species_assignment_18.fa.gz")
saveRDS(taxaspecies, "taxaspecies.RDS")
write.csv(cbind(t(seqtab.nochim), taxaspecies), "seqtab_nochim_taxaspecies_id.csv", quote=FALSE)

# assign taxonomy with HPC using rdp databast
saveRDS(seqtab.nochim, "rdp.seqtab.nochim.RDS") 

#read file from HPC
taxa_hpc<-readRDS("taxa.hpc.20230222_pr2_taxa.RDS")
write.csv(cbind(t(seqtab.nochim), taxa_hpc), "taxa_hpc.rdp.seqtab.nochim.id.csv", quote=FALSE)

# save workspace
save.image("plates0123.RData")
