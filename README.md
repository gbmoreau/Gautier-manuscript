# Gautier TUMI Pilot: Microbiota Changes Associated with Stress Treatment

## Project Summary
The following repository contains the data and code used for analysis of the intestinal microbiota in a mouse model of chronic stress. This data is based off of a pilot proposal from the University of Virginia Trans-University Microbiome Initiative and is in preparation for an upcoming manuscript. Mice (n = 24 per group) were either untreated or treated with chronic mild stress for 3 weeks. Stool samples were collected both at Baseline and after stress treatment for analysis of microbiota diversity and community composition.



## Packages Used
The specific R packages used for each analysis are located at the top of each individual code file. A comprehensive list of the packages used for this analysis as well as their release versions are listed below.

* tidyverse (Version 1.3.0)
* dada2 (Version 1.18.0)
* phyloseq (Version 1.34.0)
* microbiome (Version 1.16.0)
* vegan (Version 2.5.7)
* caret (Version 6.0.90)
* randomForest (Version 4.6.14)
* pROC (Version 1.18.0)
* ggplot2 (Version 3.3.3)
* Hmisc (Version 4.6.0)

## Dependencies
This repository includes all of the code files to perform this analysis. Raw FASTQ files for 16S rRNA sequencing are not included in this repositiory, but will be uploaded to a publically availalbe database at the time of manuscript submission. In addition, the reference database files used to assign taxonomy during the DADA2 pipeline are not included in this repository. These files were obtained from the Sivla reference database. Taxa assignment version numbers are included in the DADA2 pipeline code file.

## Setup
For setup of the initial analysis, follow these intial steps:

1. Download the raw FASTQ files for the 16S rRNA seqencing and place them in a "raw reads" directory. This directory should be inside the "data" directory.
2. Download the reference database files for taxonomy assignment and place them in the "data" directory.
3. Run the "DADA2 Pipeline for Second Run of Gautier 16S Samples.R" and "Separation of Gautier MiSeq Data into Individual Experiments.R" code files. These will generate the "DADA2 Environment.RData" and "phyloseq objects by experiment.RData" files, which are required for downstream analysis.

Note that all directory locations are relative for code in this repository, so you have to ensure that you are in the root directory before starting analysis. Each code file has a "setwd()" function which can be used to move to the root directory, and this function should be run once the "src.Rproj" file is loaded. 

Analysis is separated into Stress and Mucin Supplementation experiments. The pipeline of analysis for the Stress experiment consists of "Analysis of Stress Experiment Samples.R" and "Selection of ASVs Associated with Stress Experiment Outcomes.R" code files. The pipeline of analysis for the Mucin experiment consists of "Analysis of Mucin Supplementation Experiment for Gautier Pilot.R" and "Selection of ASVs Associated with Mucin Experiment Outcomes.R".

## Contact Information

G. Brett Moreau - gbm5pn@virginia.edu

