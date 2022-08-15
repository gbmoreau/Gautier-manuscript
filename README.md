# MANUSCRIPT NAME HERE: Baseline vs. Stressed 16S rRNA Sequencing Analysis

## Project Summary
The following repository contains the data and code used for analysis of the intestinal microbiota in a mouse model of chronic stress. This data is based off of a pilot proposal from the University of Virginia Trans-University Microbiome Initiative and has been submitted for publication in ***INSERT JOURNAL NAME HERE***. Mice (n = 24 per group) were treated with chronic mild stress for 3 weeks. Stool samples were collected both at baseline and after stress treatment for analysis of changes in microbiota diversity and community composition. This analysis found significant differences in beta diversity and community composition after stress treatment, and a random forest model clearly distinguished between groups based on these data.


## Abstract
***INSERT MANUSCRIPT ABSTRACT HERE***


## Accessing this Repository (UPDATE THE LINKS IN THIS SECTION ONCE MANUSCRIPT HAS BEEN FINALIZED)
All data and source code from this repository can be downloaded by cloning the git repository:
```
git clone https://github.com/gbmoreau/Gautier-manuscript.git
```
Alternatively, you can download a [zip archive of this repository](https://github.com/gbmoreau/Gautier-manuscript/archive/refs/heads/main.zip
). 


## Software Implementation
All source code for this analysis is stored in the `src` directory. Raw data and sample metadata used for this analysis are stored in the `data` directory. Results generated by this code are stored in the `results` directory.

## Dependencies
Before getting started with this analysis, you'll need to download the raw FASTQ files for 16S rRNA sequencing, which can be found on the Sequence Read Archive at *INSERT IDENTIFYING INFORMATION AND LINK HERE*. You will also need the reference database files used to assign taxonomy during the DADA2 pipeline. These files were obtained from the DADA2 website, and the version numbers for the specific Sivla reference databases are included in the DADA2 pipeline code file.

The specific R packages used for each analysis are located at the top of each individual code file. A comprehensive list of the packages used for this analysis as well as their release versions are listed below.

* BiocManager (Version 1.30.16)
* tidyverse (Version 1.3.1)
* dada2 (Version 1.22.0)
* phyloseq (Version 1.38.0)
* microbiome (Version 1.16.0)
* vegan (Version 2.6.2)
* caret (Version 6.0.90)
* randomForest (Version 4.6.14)
* pROC (Version 1.18.0)
* ggplot2 (Version 3.3.5)


## Setup
For setup of the initial analysis, follow these intial steps:

1. Download the raw FASTQ files for the 16S rRNA seqencing and place them in a `raw_reads` directory. This directory should be inside the `data` directory.
2. Download the reference database files for taxonomy assignment and place them in the `data` directory.
3. Run the "DADA2 Pipeline for Second Run of Gautier 16S Samples.R" code file. This will generate the "DADA2 Environment.RData" and "phyloseq object for analysis.RData" files, which are required for downstream analysis. You are now ready to perform all subsequent analyses with this data.


## Contact Information

G. Brett Moreau - gbm5pn@virginia.edu

