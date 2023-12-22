# BiNatUr project
# R script for DADA2 analysis of DNA metabarcoding data
# DIV4 primer set

# Script modified from https://benjjneb.github.io/dada2/tutorial.html

# Load dada2 package

library(dada2); packageVersion("dada2")

# Set path to each plate

path <- "~/Raw_NGS_Data/binatur/Div4/Plate1_2021_41_DIV4/" 
path <- "~/Raw_NGS_Data/binatur/Div4/Plate2_2021_41_Div4/" 
path <- "~/Raw_NGS_Data/binatur/Div4/Plate3_2021_41_Div4/" 

# The following steps  (lines 20-68) are repeated for each plate
# Check that the path is correct

list.files(path)

# Read in file names and obtain matched lists of the paired reads

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# Check read quality

plotQualityProfile(fnFs[1:12])

plotQualityProfile(fnRs[1:12])

# Filter and trim
# Trim left removes primers
# TruncLen is based on thorough analysis and optimization of parameters

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(20,20), truncLen=c(260,200), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)

# Error rates.

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

# Merge paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# Create sequence tables - repeat for each plate (1,2,3)

seqtab_plate1 <- makeSequenceTable(mergers)
dim(seqtab_plate1)
table(nchar(getSequences(seqtab_plate1)))

seqtab_plate2 <- makeSequenceTable(mergers)
dim(seqtab_plate2)
table(nchar(getSequences(seqtab_plate2)))

seqtab_plate3 <- makeSequenceTable(mergers)
dim(seqtab_plate3)
table(nchar(getSequences(seqtab_plate3)))

# Merge sequence tables

mergetab123 <- mergeSequenceTables(seqtab_plate1, seqtab_plate2, seqtab_plate3, repeats = "sum")

# Remove chimeras

seqtab.nochim <- removeBimeraDenovo(mergetab123, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Extract ASV table with sequence and abundance per sample
write.table(t(seqtab.nochim), "seqtab_nochim.tsv")

# write to csv
write.csv(t(seqtab.nochim), "seqtab.nochim.csv")

# Assign taxonomy locally with PR2
# PR2 downloaded from https://github.com/pr2database/pr2database/releases/tag/v4.14.0 on 22.02.2023
saveRDS(seqtab.nochim, "seqtab.nochim.RDS")
taxa <- assignTaxonomy(seqtab.nochim, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"), "~/data/R/dada2/pr2_version_4.14.0_SSU_dada2.fasta.gz")
saveRDS(taxa, "taxa.RDS")
write.csv(cbind(t(seqtab.nochim), taxa), "seqtab_nochim_taxa_id.csv", quote=FALSE)

# Assign taxonomy with HPC by saving .RDS for analysis 
saveRDS(collapsed.nochim, "pr2.seqtab.nochim.RDS") 

# Read RDS file from HPC
taxa_hpc <- readRDS("taxa.hpc.20230222_pr2_taxa.RDS")
write.csv(cbind(t(seqtab.nochim), taxa_hpc), "binatur.div4.plates123.taxa.pr2.seqtab.nochim.csv")

# save workspace
save.image("plates123.RData")
