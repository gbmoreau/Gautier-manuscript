#### DADA2 Pipeline for Second Run of Gautier 16S Samples ##########################################
# NAME: G. Brett Moreau
# DATE: November 8, 2021

#### PACKAGES ######################################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
packageVersion("BiocManager") # I'm using version 1.30.16

#BiocManager::install("dada2")
library(dada2)
packageVersion("dada2") # I'm using version 1.22.0

#BiocManager::install("phyloseq")
library(phyloseq)
packageVersion("phyloseq") # I'm using version 1.38.0

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # I'm using version 1.3.1

#install.packages("ggplot2")
library(ggplot2)
packageVersion("ggplot2") # I'm using version 3.3.5


#### INTRODUCTION ##################################################################################
# The Gautier TUMI pilot project is focused on investigating changes in the intestinal microbiome 
# of mice with depressive phenotypes as well as the impact of mucin degradation on this process.
# Samples were previously run, but the resulting fastq files had generally poor quality and resulted
# in significant loss of reads throughout the pipeline. Other samples run around the same time also
# had issues, so we thought this may be a Miseq issue and decided to re-run the samples on the TUMI 
# MiSeq. The sample library was re-diluted and re-run on the TUMI MiSeq, and these .fastq files will
# now be run through the DADA2 pipeline.




#### ORGANIZING READS #####################################################################################
# Before starting the DADA2 pipeline, I checked the .fastq files to see if the primer sequences 
# needed to be trimmed from each file. Sequences at the beginning of each forward .fastq file
# corresponded to the nucleotides immediately downstream of the forward primer sequence, and 
# sequences at the beginning of each reverse .fastq read corresponded to the nucleotides immediately 
# downstream of the reverse primer sequence. In addition, running .fastq files through cutadapt to
# remove primer sequences resulted in no matches. Overall, these results indicate that the primer
# sequences were removed during the sequencing pipeline and that .fastq files are ready for the 
# DADA2 pipeline.

# In order to process the .fastq files using DADA2, I'll first set up the directory paths to find
# the .fastq files.


#setwd("..") # I need to set the working directory up one to the parent directory for the code to 
# work properly. If this has been done previously during the same session then this step should be 
# skipped, or the working directory will be moved out of the parent directory.

path <- "./data/raw_reads" 
list.files(path) # We're in the right directory

# Pull forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)




#### READ QUALITY PROFILES #########################################################################
# The first step of the process is the check the quality read profiles for forward and reverse 
# reads. I'll look at these profiles in aggregate and use these profiles to determine where reads 
# need to be trimmed.



# Set seed for reproducibility
#sample(1:1000, 1) # It selected 624
set.seed(624) 


### FORWARD READ PROFILES ###
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnFs, aggregate = TRUE)

ggsave("./results/figures/Quality_Profile_forward_reads.png", width = 5, height = 3)


# Overall, the forward reads look ok. There is an initial drop in quality in the first 15-20 base
# pairs. After that, quality is good untile around base pair 150, where it starts to drop below a 
# quality score of 30. I'll trim the first 15 base pairs, then trim at base pair 160.

### REVERSE READ PROFILES ### 
plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnRs, aggregate = TRUE)

ggsave("./results/figures/Quality_Profile_reverse_reads.png", width = 5, height = 3)


# The reverse reads are pretty rough, even for reverse reads. The reads average pretty high quality
# for the first 150bp with some dips in quality, but at 150bp they quickly crater. I'll trim the 
# first 10 bp off of these reads, then trim at 150bp.




#### FILTER AND TRIM READS ########################################################################
# Place filtered files in filtered subdirectory
filtFs <- file.path("./data/filtered_reads", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./data/filtered_reads", paste0(sample.names, "_R_filt.fastq.gz"))

getwd()

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(160,150),
                     trimLeft = c(15,10),
                     maxN = 0, 
                     maxEE = c(2,2), 
                     truncQ = 2, 
                     rm.phix = TRUE,
                     compress = TRUE, 
                     multithread = TRUE) 

head(out) # None of these samples look like they lost a lot of reads during the filtering process.




#### DEREPLICATE FILTERED SEQUENCES ################################################################
# The next step of the pipeline is dereplication, which collapses together all reads that have the 
# same sequence. This reduces computational time for downstream processes. Dereplication with DADA2 
# maintains the quality information of the dereplicated sequences.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names




#### LEARN ERROR RATES #######################################################################################
# DADA2 uses a parametric model to estimate error rates for the sequence reads. These error rates are used
# downstream for sample inference.

errF <- learnErrors(derepFs, multithread = TRUE)
errR <- learnErrors(derepRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
ggsave("./results/figures/Error_Rates_forward_reads.png", width = 4, height = 4)

plotErrors(errR, nominalQ = TRUE)
ggsave("./results/figures/Error_Rates_reverse_reads.png", width = 4, height = 4)

# Looking at the error rates, which are outlined in black, they largely follow the dots of the actual
# data, indicating that they are accurate. In addition, they show the same trends as expected for each
# transition, which is outlined by the red lines. This means we can proceed with downstream processes 
# of the pipeline.




#### INFER SAMPLE COMPOSITION ################################################################################
# The DADA2 sample inference algorithm will now be applied to the filtered,trimmed, and dereplicated reads.
# This will count the number of unique sequence variants in each sample.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)



#### MERGE PAIRED READS #####################################################################################
# Now I'll merge forward and reverse reads together to generate full sequences. By default, DADA2 requires
# that forward and reverse reads overlap by at least 12 base pairs, and there must be no mismatches along
# these base pairs. I'll keep these parameters.

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])




#### CONSTRUCT SEQUENCE TABLE ######################################################################
# Finally, I'll construct a table of all Amplicon Sequence Variants (ASVs) present in our samples.
seqtab <- makeSequenceTable(mergers)

dim(seqtab) # There are 29850 unique sequences in 171 samples.

table(nchar(getSequences(seqtab))) 
plot(table(nchar(getSequences(seqtab))), xlab = "read length", ylab="number of reads") # Almost all 
# sequences are 227-229 base pairs in length. The expected product size is ~250bp. 15bp were removed
# from the end of the forward sequence, while 10bp were removed from the reverse sequence, giving a 
# predicted product length of ~220bp. Based on this, I'll keep the most abundant lengths, which are
# those between 227-229bp.

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(227,229)]
table(nchar(getSequences(seqtab2))) 
plot(table(nchar(getSequences(seqtab2))), xlab = "read length", ylab="number of reads") # Now all 
# sequences in the data set are 215 base pairs in length.




#### REMOVE CHIMERAS ###############################################################################
# Now I will remove chimeras (sequencing artifacts that are combinations of multiple unique 
# sequences). This new table will include all sequences after chimeras have been removed.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtab.nochim) # The original data set had 29850 total ASVs in 171 samples. After removing chimeras, 
# the new data set includes only 1795 ASVs in 171 samples, indicating that 93.9% of all ASVs from the 
# original data set were chimeras.

sum(seqtab.nochim)/sum(seqtab2) # While almost all of the original ASVs were chimeras, the chimeras
# made up a smaller (~13%) proportion of the ASVs overall abundance.




#### TRACK PIPELINE READS ##########################################################################
# For quality control I'm going to track the reads over the course of the pipeline. While I'd expect 
# there to be a loss of sequences at the filtering step, there should not be a significantly lower 
# number of reads at any other step of the pipeline. If there is then that step should be examined 
# in more detail.

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) # Looking over the tracked reads, it doesn't look like there was a major drop in reads at 
# any step of the pipeline, which is what we'd expect.

# Now I'll format the data frame to add the percent of input reads that made it through the pipeline.
track <- as.data.frame(track)
track <- plyr::rename(track, replace = c("input" = "Input", 
                                         "filtered" = "Filtered",
                                         "denoisedF" = "DenoisedF",
                                         "denoisedR" = "DenoisedR",
                                         "nonchim" = "Non-Chimera"))

track$Percent.Input <- track$`Non-Chimera`/track$Input * 100

# I'll print the table of tracked reads for future reference.
write.csv(track, file = "./results/tables/tracked read counts.csv")




#### ASSIGN TAXONOMY ##############################################################################
# Now I'll assign taxonomy to the ASV table using an assignment file from the Silva reference 
# database. For this project I'll use version 138.1 for both genus and species assignment.

# NOTES FROM DADA2: The DADA2 package provides a native implementation of the naive Bayesian 
# classifier method for this purpose. The assignTaxonomy function takes as input a set of sequences 
# to be classified and a training set of reference sequences with known taxonomy, and outputs taxonomic 
# assignments with at least minBoot bootstrap confidence.

taxa <- assignTaxonomy(seqtab.nochim, "./data/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa <- addSpecies(taxa, "./data/silva_species_assignment_v138.1.fa.gz") # Instructions for obtaining 
# these reference database files are available on the DADA2 pipeline website.

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

head(taxa.print) 

# There is a significant proportion of ASVs that are missing Genus-level designations, but Family-
# level designations are present for almost all ASVs.




#### HAND OFF THE PHYLOSEQ ########################################################################
# I'll now be transferring the analysis from the DADA2 package to the Phyloseq package. This will 
# allow us to construct a phyloseq object that contains all of the DADA2 pipeline information as 
# well as metadata.

# In order to assemble the phyloseq object, I'll need to import sample metadata and merge it with
# taxanomic and ASV data. I noticed that the dates present in the sample names are slightly 
# different between the first MiSeq run (20210611) and the second MiSeq run (20210603), so I will
# change the sample names in the metadata to ensure they are identical to those in the current
# samples.

metadata <- read.csv("./data/metadata.csv") # This includes the metadata for the samples such as
# experiment name and condition.

# Format to correct sample names
Sample.ID <- metadata$Sample.ID
Sample.ID <- gsub('.{2}$', '', Sample.ID)
metadata$Sample.ID <- Sample.ID
metadata$Sample.ID <- paste0(metadata$Sample.ID, "03")


# Now I'll form the phyloseq object from the metadata, ASV table, and taxa information.
row.names(metadata)<-metadata$Sample.ID
metadata$Sample.ID <- NULL

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxa))


# DADA2 NOTES: It is more convenient to use short names for ASVs (e.g. ASV21) rather than the 
# full DNA sequence when working with some of the tables and visualizations from phyloseq, but 
# we want to keep the full DNA sequences for other purposes like merging with other datasets or 
# indexing into reference databases like the Earth Microbiome Project. For that reason we’ll store 
# the DNA sequences of our ASVs in the “refseq slot” of the phyloseq object, and then rename the 
# taxa to a short string. That way, the short new taxa names will appear in tables and plots, and 
# we can still recover the DNA sequences corresponding to each ASV as needed with refseq(ps).

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# The Phyloseq object is now ready for use. I'll save this R object for subsequent use in the 
# analysis pipeline.

save.image(file = "./results/DADA2 Environment.RData")



