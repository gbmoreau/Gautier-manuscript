#### ANALYSIS OF GAUTIER 16S SAMPLES FOR PILOT PROJECT #############################################
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

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # I'm using version 1.3.1

#install.packages("ggplot2")
library(ggplot2)
packageVersion("ggplot2") # I'm using version 3.3.5

#BiocManager::install("microbiome")
library(microbiome) 
packageVersion("microbiome") # I'm using version 1.16.0

#install.packages("vegan")
library(vegan)
packageVersion("vegan") # I'm using version 2.5.7


#### INTRODUCTION ##################################################################################
# The Gautier TUMI pilot project is focused on investigating changes in the intestinal microbiome 
# of mice with depressive phenotypes as well as the impact of mucin degradation on this process.
# Samples were previously run, but the resulting fastq files had generally poor quality and resulted
# in significant loss of reads throughout the pipeline. Other samples run around the same time also
# had issues, so we thought this may be a Miseq issue and decided to re-run the samples on the TUMI 
# MiSeq. The sample library was re-diluted and re-run on the TUMI MiSeq, and these .fastq files will
# now be run through the DADA2 pipeline.




#### LOAD PHYLOSEQ OBJECT ##########################################################################
# I'll start by loading the initial phyloseq object generated after the DADA2 pipeline. This will be 
# used for analysis.

#setwd("..") # I need to set the working directory up one to the parent directory for the code to 
# work properly. If this has been done previously during the same session then this step should be 
# skipped, or the working directory will be moved out of the parent directory.

load(file = "./results/phyloseq objects by experiment.RData")

# For this analysis, I'll be focusing on the Naive vs Stressed mice experiment (ps.stress), so I'll
# remove other objects from the environment.

rm(ps.mucin)

# We're now ready to start analysis of the data.




#### CHARACTERIZATION OF THE DATA SET ##############################################################
# Before I move on, I'll first characterize the data set as a whole

View(ps.stress@sam_data)

# Overall, there are 24 samples in this experiment. Samples were divided into two groups (Naive and 
# Stressed) and two different timepoints (pre- and post-stress intervention) looked at pre- and post-
# stress. Each treatment group and timepoint is divided into 6 samples each.

dim(ps.stress@otu_table) # From these 24 samples, there were 1795 unique ASVs.

total.reads <- sample_sums(ps.stress)

sum(total.reads) # There are a total of 1,741,609 reads across the 40 samples in the data set.
mean(total.reads) # The average read number is 72,567
median(total.reads) # The median read number is 71,421
range(total.reads) # The range of reads in samples is from 43,957-117,986 total reads.


# Before I perform downstream analyses, I'll first convert the data into relative abundances. I'll 
# also divide by the number of samples because they are identical in each group, resulting in 
# values that represent percent of total for each taxa.

ps.stress@sam_data$Condition <- factor(ps.stress@sam_data$Condition, levels = c("Naïve Pre", "Naïve Post", "Stress Pre", "Stress Post"))
ps.stress.prop <- transform_sample_counts(ps.stress, function(ASV) ASV/sum(ASV)/6*100) 




#### ALPHA DIVERSITY #############################################################################
# I'll first look at alpha diversity. The DADA2 pipeline removes singlets (sequences detected only
# once) and doublets (sequences detected twice) from the data set, as it does not have enough 
# confidence to determine whether these sequences are real or artifacts. Because of this, measures
# of richness (how many ASVs are present in a sample) faile explicitly. I'll focus mainly on 
# measures of evenness but give measures of richness a glance to see if there are any major 
# differences.

### RICHNESS ###
# I'll look at the number of Observed ASVs as my measure of richness.

plot_richness(ps.stress, x = "Condition", measures = c("Observed"), 
              color = "Condition") + 
  geom_point() +
  geom_boxplot() + 
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  scale_color_discrete(name = "Condition", 
                       breaks = c("Naïve Pre", "Naïve Post", "Stress Pre", "Stress Post"), 
                       labels = c("Naïve Pre", "Naïve Post", "Stress Pre", "Stress Post")) 


ggsave("./results/figures/stress/Stress Experiment_richness_observed.png")

# None of the groups appear to have signifcant differences between each other regarding richness. 
# There looks like a small increase in observed ASVs in both Naive and Stressed mice post-
# intervention, but this is likely not significant. To confirm, I'll perform a paired Wilcoxon
# test with a bonferroni correction for multiple comparisons.


# Now I'll look for significant differences using a paired Wilcoxon test.
richness <- estimate_richness(ps.stress)
pairwise.wilcox.test(richness$Observed, p.adjust.method = "bonferroni", sample_data(ps.stress)$Condition)

# There are no significant differences between any of the groups looking at alpha diversity.


# I'll focus more on evenness, which are more robust to noise from rare sequence varients and 
# can be interpretedusing this pipeline. I'll use Pielou's Evenness as my evenness measure.

# First I'll organize the data.
sample.data.stress <- ps.stress@sam_data
sample.data.stress <- cbind("Sample.Names" = rownames(sample.data.stress), sample.data.stress)

evenness<- evenness(ps.stress, index = c("pielou"))
evenness<- cbind("Sample.Names" = rownames(evenness), evenness)

sample.data.stress <- full_join(sample.data.stress, evenness, by = "Sample.Names")


### PIELOU ###
ggplot(sample.data.stress, aes(x = Condition, y = pielou, color = Condition)) +
  geom_point() +
 geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none")

ggsave("./results/figures/stress/Stress Experiment_evenness_pielou.png")

# Overall, there doesn't appear to be a whole lot of separation in evenness between groups. 
# There is a lot of spread between individuals within each group, which makes the data a 
# little difficult to interpret. Now I'll look for statistical differences.


pairwise.wilcox.test(evenness$pielou, p.adjust.method = "bonferroni", sample_data(ps.stress)$Condition)

# No significant differences were observed in evenness between groups, either.




#### BETA DIVERSITY #######################################################################


#### BETA DIVERSITY: BRAY-CURTIS DISSIMILARITY  #####################################################
# I'll now move on to measures of beta diversity, which investigates variation between samples.
# I'll start by looking at Bray-Curtis Dissimilarity.

sample(1:1000, 1) # It selected 763
set.seed(763) # Set seed for reproducibility

# First, let's ordinate Bray-Curtis Dissimilarity using NMDS, then plot this ordination to visualize
# separation between different samples.

ord.nmds.bray <- ordinate(ps.stress.prop, method = "NMDS", distance = "bray")

plot_ordination(ps.stress.prop, ord.nmds.bray, color = "Condition", 
                title = "Bray-Curtis Dissimilarity") +
  theme_bw() +
  scale_color_discrete(name = "Condition", 
                       breaks = c("Naïve Pre", "Naïve Post", "Stress Pre", "Stress Post"), 
                       labels = c("Naïve Pre", "Naïve Post", "Stress Pre", "Stress Post")) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) 

ggsave("./results/figures/stress/Stress Experiment_beta_diversity_bc.png", width = 6, height = 6)


# Bray-Curtis distance is ordinated using Non-metric multidimensional scaling (NMDS). It looks like
# while there is overlap between groups, groups are clustering somewhat distinctly from each other.
# To better explore this, I'll perform a PERMANOVA to see if there are any differences in clustering
# between groups.

adonis(distance(ps.stress.prop, method = "bray") ~sample_data(ps.stress.prop)$Condition)

# The PERMANOVA found significant differences in clustering (p = 0.001), but this doesn't specifiy
# which groups are significantly different from each other. To test this, I'll perform pairwise
# PERMANOVAs with multiple comparisons correction.

permanova.combinations <- combn(x = levels(ps.stress.prop@sam_data$Condition), m = 2)
permanova.p.bray <- vector(length = 6)
for(i in 1:ncol(permanova.combinations)){
  ps.subs.bray <- subset_samples(ps.stress.prop, Condition %in% permanova.combinations[,i])
  metadata.subs.bray <- as(sample_data(ps.subs.bray), "data.frame")
  pairwise.bray <- adonis(distance(ps.subs.bray, method = "bray") ~ Condition, data = metadata.subs.bray)
  permanova.p.bray[i] <- pairwise.bray$aov.tab[1,ncol(pairwise.bray$aov.tab)]
}
combinations <- c("Naïve Pre/Naïve Post", "Naïve Pre/Stress Pre", "Naïve Pre/Stress Post",
                  "Naïve Post/Stress Pre", "Naïve Post/Stress Post", "Stress Pre/Stress Post")

permanova.FDR.bray <- p.adjust(permanova.p.bray, method = "bonferroni")
permanova.bray.table <- cbind(combinations, permanova.p.bray)
permanova.bray.table <- cbind(permanova.bray.table, permanova.FDR.bray)

View(permanova.bray.table)

# Looking at the pairwise comparisons, it looks like clustering is significantly different for
# all comaprisons except for Naive Pre vs. Naive Post and Naive Pre vs Stress Pre. I'll now
# look at Unweighted Unifrac as a measure of beta diversity.




#### BETA DIVERSITY: UNWEIGHTED UNIFRAC #############################################################
# I'll also look at unweighted Unifrac, another measure of beta diversity. While Bray-Curtis 
# dissimilarity looks only at counts between samples, Unifrac looks at both counts and differences in 
# phylogeny between samples, so we'll need to add some taxonomic information.

random_tree <- ape::rtree(ntaxa(ps.stress.prop), rooted = TRUE, tip.label = taxa_names(ps.stress.prop))
samples.out <- rownames(seqtab.nochim)
ps.stress.prop <- merge_phyloseq(ps.stress.prop, samples.out, random_tree)

unifrac_dist <- phyloseq::distance(ps.stress.prop, method = "unifrac", weighted = FALSE)
ordination <- ordinate(ps.stress.prop, method = "PCoA", distance = unifrac_dist)

plot_ordination(ps.stress.prop, ordination, color = "Condition") + 
  theme_bw() +
  scale_color_discrete(name = "Condition", 
                       breaks = c("Naïve Pre", "Naïve Post", "Stress Pre", "Stress Post"), 
                       labels = c("Naïve Pre", "Naïve Post", "Stress Pre", "Stress Post")) +
  labs(title = "Unweighted Unifrac Distance") +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) 

ggsave("./results/figures/stress/Stress Experiment_beta_diversity_unifrac.png", width = 6, height = 6)


# Unweighted Unifrac distance is ordinated using Principal Coordinate Analysis (PCoA). Once again,
# groups appear to overlap but there is somewhat distinct clustering for some groups. I'll once
# again perform a PERMANOVA to see if there are any differences in clustering.

adonis(distance(ps.stress.prop, method = "unifrac", weighted = FALSE) ~sample_data(ps.stress.prop)$Condition)

# The PERMANOVA found significant differences in clustering (p = 0.001), but this doesn't specifiy
# which groups are significantly different from each other. To test this, I'll perform pairwise
# PERMANOVAs with multiple comparisons correction.


permanova.p.unifrac <- vector(length = 6)
for(i in 1:ncol(permanova.combinations)){
  ps.subs.unifrac <- subset_samples(ps.stress.prop, Condition %in% permanova.combinations[,i])
  metadata.subs.unifrac <- as(sample_data(ps.subs.unifrac), "data.frame")
  pairwise.unifrac <- adonis(distance(ps.subs.unifrac, method = "unifrac") ~ Condition, 
                             data = metadata.subs.unifrac)
  permanova.p.unifrac[i] <- pairwise.unifrac$aov.tab[1,ncol(pairwise.unifrac$aov.tab)]
}

permanova.FDR.unifrac <- p.adjust(permanova.p.unifrac, method = "bonferroni")
permanova.unifrac.table <- cbind(combinations, permanova.p.unifrac)
permanova.unifrac.table <- cbind(permanova.unifrac.table, permanova.FDR.unifrac)

View(permanova.unifrac.table)
# Pairwise comparisons of PERMANOVA from Unifrac distance found the same groups were significantly
# different as the groups identified from the Bray-Curtis pairwise comparisons.Overall, beta 
# diversity measures found very similar results, suggesting that there are some differences 
# between groups.

