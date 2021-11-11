#### SEPARATION OF GAUTIER MISEQ DATA INTO INDIVIDUAL EXPERIMENTS ##################################
# NAME: G. Brett Moreau
# DATE: November 9, 2021

#### PACKAGES ######################################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
packageVersion("BiocManager") # I'm using version 1.30.16

#BiocManager::install("phyloseq")
library(phyloseq)
packageVersion("phyloseq") # I'm using version 1.38.0


#### INTRODUCTION ##################################################################################
# The Gautier TUMI pilot project is focused on investigating changes in the intestinal microbiome 
# of mice with depressive phenotypes as well as the impact of mucin degradation on this process.
# Samples have been run on the MiSeq and sample reads have been processed using DADA2 and organized
# into a phyloseq object. 

# The MiSeq run consisted of samples from several different experiments, so they will need to be
# separated into these groups before 




#### PROCESSING OF PHYLOSEQ OBJECT #################################################################
#setwd("..") # I need to set the working directory up one to the parent directory for the code to 
# work properly. If this has been done previously during the same session then this step should be 
# skipped, or the working directory will be moved out of the parent directory.


# Before we get started, I'll load the environment from the DADA2 run.
load("./results/DADA2 Environment.RData")


# I'll start by looking at the tracked reads to see if there are any samples with a low number of 
# reads.
View(track)

# Sample #167 had very few reads both at input (537) and at the end of the pipeline (32), so this
# sample will be removed from the data set. Sample #170 was also lower than other samples (about 
# half as many reads as the next highest sample), but the percent recovered through the pipeline
# is in line with other samples. Because of this, I'll go ahead and include this sample in the 
# analysis.

# On the other end of the spectrum, Sample #25 had significantly more reads (2x the next lowest) 
# than other samples. 

# On the other end of the spectrum, Sample # 25 had over 2x more reads than the next highest 
# sample. There are several schools of thought on how to deal with large differences in read count
# between samples, mostly focusing around whether samples should be rarefied/subsampled. This 
# normalizes read count based on the sample with the lowest number of reads, however several 
# publications have come out suggesting that this does more harm than good. Therefore, I'll move
# ahead without subsampling, keeping this sample in the data set with no adjustment.


ps <- prune_samples(sample_names(ps) != "167-20210603", ps)

# This will now be the final phyloseq object which can be used for downstream analyses.




#### SEPARATION OF DIFFERENT EXPERIMENTS FROM THE PHYLOSEQ OBJECT ###############################
# There are six experiments in this MiSeq run. I will separate the phyloseq object into its
# individual experiments based on the experiment metadata.


# Experiment 1: Naive vs Stress mice: changes in microbiota
ps.stress <- subset_samples(ps, Experiment == "Nov Stress Experiment")

# Experiment 2: Mucin Supplementation 
ps.mucin <- subset_samples(ps, Experiment == "Mucin Therapeutic ")


# Save phyloseq objects into a separate .RData file
rm(list = ls()[!ls() %in% c("ps.stress", "ps.mucin", "seqtab.nochim")])

save.image(file = "./results/phyloseq objects by experiment.RData")

