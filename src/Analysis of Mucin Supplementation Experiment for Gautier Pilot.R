#### ANALYSIS OF MUCIN SUPPLEMENTATION EXPERIMENT FOR GAUTIER PILOT #################################
# NAME: G. Brett Moreau
# DATE: November 11, 2021

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

#install.packages("caret")
library(caret)
packageVersion("caret") # I'm using version 6.0.90

#install.packages("randomForest")
library(randomForest)
packageVersion("randomForest") # I'm using version 4.6.14

#install.packages("pROC")
library(pROC)
packageVersion("pROC") # I'm using version 1.18.0


#### INTRODUCTION ##################################################################################
# The Gautier TUMI pilot project is focused on investigating changes in the intestinal microbiome 
# of mice with depressive phenotypes as well as the impact of mucin degradation on this process.
# Samples were previously run, but the resulting fastq files had generally poor quality and resulted
# in significant loss of reads throughout the pipeline. Other samples run around the same time also
# had issues, so we thought this may be a Miseq issue and decided to re-run the samples on the TUMI 
# MiSeq. The sample library was re-diluted and re-run on the TUMI MiSeq, and these .fastq files will
# now be run through the DADA2 pipeline.


setwd("..")

#### LOAD PHYLOSEQ OBJECT ##########################################################################
# I'll start by loading the initial phyloseq object generated after the DADA2 pipeline. This will be 
# used for analysis.

#setwd("..") # I need to set the working directory up one to the parent directory for the code to 
# work properly. If this has been done previously during the same session then this step should be 
# skipped, or the working directory will be moved out of the parent directory.

load(file = "./results/phyloseq objects by experiment.RData")

# For this analysis, I'll be focusing on the Mucin Supplementation mice experiment (ps.mucin), so I'll
# remove other objects from the environment.

rm(ps.stress)

# We're now ready to start analysis of the data.




#### CHARACTERIZATION OF THE DATA SET ##############################################################
# Before I move on, I'll first characterize the data set as a whole

View(ps.mucin@sam_data)

# Overall, there are 71 samples in this experiment. Samples were divided into two groups (Control and 
# Mucin-supplemented) and three different timepoints (pre-stress, 3 weeks post-stress, and 4 weeks
# post-stress). Importantly, Mucin treatment does not start until Week 4, so Control and Mucin-
# supplemented groups are identical at baseline and 3 weeks post-stress. There are 12 mice for each 
# group/timepoint, but Sample #167 was removed due to low read counts, giving our 71 samples.

dim(ps.mucin@otu_table) # From these 71 samples, there were 1795 unique ASVs.

total.reads <- sample_sums(ps.mucin)

sum(total.reads) # There are a total of 4,784,563 reads across the 71 samples in the data set.
mean(total.reads) # The average read number is 67,388
median(total.reads) # The median read number is 66,966
range(total.reads) # The range of reads in samples is from 39,993-101,514 total reads.


# Before I perform downstream analyses, I'll first convert the data into relative abundances and
# order the sample groups.

ps.mucin.prop <- transform_sample_counts(ps.mucin, function(ASV) ASV/sum(ASV)) 




#### ALPHA DIVERSITY #############################################################################
# I'll first look at alpha diversity. The DADA2 pipeline removes singlets (sequences detected only
# once) and doublets (sequences detected twice) from the data set, as it does not have enough 
# confidence to determine whether these sequences are real or artifacts. Because of this, measures
# of richness (how many ASVs are present in a sample) fail explicitly. I'll focus mainly on 
# measures of evenness but give measures of richness a glance to see if there are any major 
# differences.

### RICHNESS ###
# I'll look at the number of Observed ASVs as my measure of richness.

plot_richness(ps.mucin, x = "Group", measures = c("Observed"), 
              color = "Group") + 
  geom_point() +
  geom_boxplot() + 
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  scale_color_discrete(name = "Group", 
                       breaks = c("Baseline", "3 Wk Stress", "Stress Alone", "Stress + Mucin"), 
                       labels = c("Baseline", "3 Wk Stress", "Stress Alone", "Stress + Mucin")) 


ggsave("./results/figures/mucin/Mucin Experiment_richness_observed.png", width = 5, height = 4)

# There don't appear to be any large changes in richness between groups. I'll perform a paired Wilcoxon test
# to confirm.

richness <- estimate_richness(ps.mucin)
pairwise.wilcox.test(richness$Observed, p.adjust.method = "bonferroni", sample_data(ps.mucin)$Group)

# There are no significant differences between any of the groups looking at Observed ASVs.


# I'll focus more on evenness, which are more robust to noise from rare sequence varients and 
# can be interpreted using this pipeline. I'll use Pielou's Evenness as my evenness measure.

# First I'll organize the data.
sample.data.mucin <- ps.mucin@sam_data
sample.data.mucin <- cbind("Sample.Names" = rownames(sample.data.mucin), sample.data.mucin)

evenness<- evenness(ps.mucin, index = c("pielou"))
evenness<- cbind("Sample.Names" = rownames(evenness), evenness)

sample.data.mucin <- full_join(sample.data.mucin, evenness, by = "Sample.Names")


### PIELOU ###
ggplot(sample.data.mucin, aes(x = Group, y = pielou, color = Group)) +
  geom_point() +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none")

ggsave("./results/figures/mucin/Mucin Experiment_evenness_pielou.png", width = 5, height = 4)

# Overall, there doesn't appear to be a whole lot of separation in evenness between groups, 
# although there is a good amount of spread in the data. It does appear that evenness trends
# upwards around 4 weeks in both the Control and Mucin-supplemented group. Now I'll look 
# for statistical differences.


pairwise.wilcox.test(evenness$pielou, p.adjust.method = "bonferroni", sample_data(ps.mucin)$Group)

# Significant differences are mainly between 4 week samples and all other samples. The 
# Control samples (Stress Alone) at 4 weeks are significantly different from Baseline 
# (p=0.0023) and 3 weeks (p=0.0005). The Mucin samples (Stress + Mucin) at 4 weeks are 
# significantly different from Baseline (p=0.0023) and 3 weeks (p=0.0011). Control and 
# Mucin samples at 4 weeks were not significantly different from each other.

# Overall, this suggests that samples at 4 weeks show more evenness than those pre-stress
# or after 3 weeks stress.




#### BETA DIVERSITY: BRAY-CURTIS DISSIMILARITY  #####################################################
# I'll now move on to measures of beta diversity, which investigates variation between samples.
# I'll start by looking at Bray-Curtis Dissimilarity.

#sample(1:1000, 1) # It selected 966
set.seed(966) # Set seed for reproducibility

# First, let's ordinate Bray-Curtis Dissimilarity using NMDS, then plot this ordination to visualize
# separation between different samples.

ord.nmds.bray <- ordinate(ps.mucin.prop, method = "NMDS", distance = "bray")

plot_ordination(ps.mucin.prop, ord.nmds.bray, color = "Group", 
                title = "Bray-Curtis Dissimilarity") +
  theme_bw() +
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798", "#521b92", "#4d8e00"),
                       breaks = c("Baseline", "3 Wk Stress", "Stress Alone", "Stress + Mucin"), 
                       labels = c("Baseline", "3 Wk Stress", "Stress Alone", "Stress + Mucin")) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) 

ggsave("./results/figures/mucin/Mucin Experiment_beta_diversity_bc.png", width = 4, height = 4)


# Bray-Curtis distance is ordinated using Non-metric multidimensional scaling (NMDS). Based on the
# plot, it echos the results of our evenness measure, which found differences between 4 week samples
# and all other samples. In addition, all 3 timepoints appear to be clustering independently, with
# 4 week experimental conditions (Control vs Mucin-supplemented) separating from each other, as well.

# I'll perform a PERMANOVA to look statistically if there are any differences in clustering between groups.

adonis(distance(ps.mucin.prop, method = "bray") ~sample_data(ps.mucin.prop)$Group)

# The PERMANOVA found significant differences in clustering (p = 0.001), but this doesn't specifiy
# which groups are significantly different from each other. To test this, I'll perform pairwise
# PERMANOVAs with multiple comparisons correction.

permanova.combinations <- combn(x = levels(ps.mucin.prop@sam_data$Group), m = 2)
permanova.p.bray <- vector(length = 6)
for(i in 1:ncol(permanova.combinations)){
  ps.subs.bray <- subset_samples(ps.mucin.prop, Group %in% permanova.combinations[,i])
  metadata.subs.bray <- as(sample_data(ps.subs.bray), "data.frame")
  pairwise.bray <- adonis(distance(ps.subs.bray, method = "bray") ~ Group, data = metadata.subs.bray)
  permanova.p.bray[i] <- pairwise.bray$aov.tab[1,ncol(pairwise.bray$aov.tab)]
}
combinations <- c("Baseline/3 Wk Stress", "Baseline/Stress Alone", "Baseline/Stress + Mucin", 
                  "3 Wk Stress/Stress Alone", "3 Wk Stress/Stress + Mucin",
                  "Stress Alone/Stress + Mucin")

permanova.FDR.bray <- p.adjust(permanova.p.bray, method = "bonferroni")
permanova.bray.table <- cbind(combinations, permanova.p.bray)
permanova.bray.table <- cbind(permanova.bray.table, permanova.FDR.bray)

View(permanova.bray.table)

# According to the pairwise comparisons of the PERMANOVA, all groups are significantly different
# from each other (FDR-adjusted p value = 0.006). This is consistent with initial impressions 
# of the plot when looking by eye.




#### BETA DIVERSITY: UNWEIGHTED UNIFRAC #############################################################
# I'll also look at unweighted Unifrac, another measure of beta diversity. While Bray-Curtis 
# dissimilarity looks only at counts between samples, Unifrac looks at both counts and differences in 
# phylogeny between samples, so we'll need to add some taxonomic information.

random_tree <- ape::rtree(ntaxa(ps.mucin.prop), rooted = TRUE, tip.label = taxa_names(ps.mucin.prop))
samples.out <- rownames(seqtab.nochim)
ps.mucin.prop <- merge_phyloseq(ps.mucin.prop, samples.out, random_tree)

unifrac_dist <- phyloseq::distance(ps.mucin.prop, method = "unifrac", weighted = FALSE)
ordination <- ordinate(ps.mucin.prop, method = "PCoA", distance = unifrac_dist)

plot_ordination(ps.mucin.prop, ordination, color = "Group") + 
  theme_bw() +
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798", "#521b92", "#4d8e00"),
                     breaks = c("Baseline", "3 Wk Stress", "Stress Alone", "Stress + Mucin"), 
                     labels = c("Baseline", "3 Wk Stress", "Stress Alone", "Stress + Mucin")) +
  labs(title = "Unweighted Unifrac Distance") +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) 

ggsave("./results/figures/mucin/Mucin Experiment_beta_diversity_unifrac.png", width = 6, height = 6)


# Unweighted Unifrac distance is ordinated using Principal Coordinate Analysis (PCoA). The trends
# here are similar: all timepoints (Baseline, 3 weeks, and 4 weeks) appear to separate from each
# other, with 4 weeks amples being distinct from the others. In addition, Control and Mucin-
# supplemented samples appear to be separate from each other.




### COMMUNITY COMPOSITION: PHYLUM LEVEL #############################################################
# There appear to be some broad differences in diversity between groups. I'll now see if this is
# reflected by any broad changes in community composition across groups. I'll start by looking at
# the phylum level, looking at all phyla that make up at least 1% of the relative abundance in at
# least one group.

ps.mucin.prop.phylum <- tax_glom(ps.mucin.prop, "Phylum", NArm = FALSE)
phylum.table.mucin <- psmelt(ps.mucin.prop.phylum) # Organize in long format for ggplot.

# There are 24 samples for "Baseline" and "3 Wk Stress" groups, 12
# samples for the "Stress Alone" group, and 11 samples for the "Stress + Mucin" group. I'll 
# separate these groups to generate relative abundances as a percentage by dividing by the total
# number of samples in each group.

phylum.table.mucin.baseline.3wk <- filter(phylum.table.mucin, Group == "Baseline" | Group == "3 Wk Stress")
phylum.table.mucin.baseline.3wk$Relative.Abundance <- NA
phylum.table.mucin.baseline.3wk$Relative.Abundance <- phylum.table.mucin.baseline.3wk$Abundance / 24 *100


phylum.table.stress <- filter(phylum.table.mucin, Group == "Stress Alone")
phylum.table.stressRelative.Abundance <- NA
phylum.table.stress$Relative.Abundance <- phylum.table.stress$Abundance / 12 *100


phylum.table.stress.mucin <- filter(phylum.table.mucin, Group == "Stress + Mucin")
phylum.table.stress.mucinRelative.Abundance <- NA
phylum.table.stress.mucin$Relative.Abundance <- phylum.table.stress.mucin$Abundance / 11 *100


phylum.table.mucin <- rbind(phylum.table.mucin.baseline.3wk, phylum.table.stress)
phylum.table.mucin <- rbind(phylum.table.mucin, phylum.table.stress.mucin)


phylum.table.mucin <- as.data.frame(phylum.table.mucin)


# I'll plot the data to get an overview of differences between groups.
ggplot(phylum.table.mucin, aes(x = Group, y = Relative.Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") 

# This figure shows that the ratio of Bacteroidota to Firmicutes increases over time. However, I 
# want to clean up the figure by removing all Phyla that make up less than 1% of the total in any
# group and placing them in their own category.


# First, I'll summarize relative abundance information per taxonomic group.
summarized.abundance.phylum <- phylum.table.mucin %>%
  group_by(Phylum, Group) %>%
  summarize(Total.Relative.Abundance = (sum(Relative.Abundance)))

summarized.abundance.phylum$Phylum <- as.character(replace_na(summarized.abundance.phylum$Phylum, replace = "Unidentified"))

# Now, I'll collect taxa names for taxa with relative abundance above 1% in any group. These
# names will be preserved in each group.
summarized.abundance.phylum.above1 <- filter(summarized.abundance.phylum, Total.Relative.Abundance >= 1.0)
summarized.abundance.phylum.above1$Phylum <- factor(summarized.abundance.phylum.above1$Phylum, 
                                                    levels = sort(unique(summarized.abundance.phylum.above1$Phylum)))
abundant.families <- as.character(unique(summarized.abundance.phylum.above1$Phylum)) 
table(abundant.families) # There are 3 Phyla that have a Total Relative Abundance above 1%

# Make new data frames either including or excluding abundant families.
summarized.abundance.phylum.no.remainder <- filter(summarized.abundance.phylum, 
                                                   Phylum %in% abundant.families)

# I'll now summarize the remainder samples that compose less than 1% of relative abundance per group.
summarized.abundance.phylum.remainders.only <- anti_join(summarized.abundance.phylum, summarized.abundance.phylum.no.remainder, by = "Phylum")

summarized.abundance.phylum.remainders.only.summary <- summarized.abundance.phylum.remainders.only %>%
  group_by(Group) %>%
  summarize(Total.Relative.Abundance = (sum(Total.Relative.Abundance))) %>%
  mutate(Phylum = "Phyla < 1%", .before = 1)


# Now I'll combine abundant families with the low abundance remainders to account for 100% of total 
# abundance. This data frame will be used to generate the community composition figure.

summarized.abundance.phylum.figure <- rbind(summarized.abundance.phylum.no.remainder, summarized.abundance.phylum.remainders.only.summary)
summarized.abundance.phylum.figure$Phylum <- factor(summarized.abundance.phylum.figure$Phylum, levels = c("Bacteroidota", "Firmicutes", "Verrucomicrobiota", "Phyla < 1%"))
                                                                                                          
ggplot(summarized.abundance.phylum.figure, aes(x = Group, y = Total.Relative.Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_bw()


ggsave("./results/figures/mucin/Mucin Experiment_community_composition_phylum.png", width = 6, height = 4)




### COMMUNITY COMPOSITION: FAMILY LEVEL #############################################################
# There appear to be some differences at the phylum level. I'll now look at the family level to see
# how these differences are reflected here.

ps.mucin.prop.family <- tax_glom(ps.mucin.prop, "Family", NArm = FALSE)
family.table.mucin <- psmelt(ps.mucin.prop.family) # Organize in long format for ggplot.

table(ps.mucin@sam_data$Group) # There are 24 samples for "Baseline" and "3 Wk Stress" groups, 
# 12 samples for the "Stress Alone" group, and 11 samples for the "Stress + Mucin" group. I'll 
# separate these groups to generate relative abundances as a percentage by dividing by the total
# number of samples in each group.

family.table.mucin.baseline.3wk <- filter(family.table.mucin, Group == "Baseline" | Group == "3 Wk Stress")
family.table.mucin.baseline.3wk$Relative.Abundance <- NA
family.table.mucin.baseline.3wk$Relative.Abundance <- family.table.mucin.baseline.3wk$Abundance / 24 *100


family.table.stress <- filter(family.table.mucin, Group == "Stress Alone")
family.table.stressRelative.Abundance <- NA
family.table.stress$Relative.Abundance <- family.table.stress$Abundance / 12 *100


family.table.stress.mucin <- filter(family.table.mucin, Group == "Stress + Mucin")
family.table.stress.mucinRelative.Abundance <- NA
family.table.stress.mucin$Relative.Abundance <- family.table.stress.mucin$Abundance / 11 *100


family.table.mucin <- rbind(family.table.mucin.baseline.3wk, family.table.stress)
family.table.mucin <- rbind(family.table.mucin, family.table.stress.mucin)


family.table.mucin <- as.data.frame(family.table.mucin)

# I'll plot the data to get an overview of differences between groups.
ggplot(family.table.mucin, aes(x = Group, y = Relative.Abundance, fill = Family)) +
  geom_bar(stat = "identity") 

# This figure is pretty cluttered. To clean it up, I'll filter out all Families that make up less than 1
# % of total in any group and place them in their own category.


# First, I'll summarize relative abundance information per taxonomic group.
summarized.abundance.family <- family.table.mucin %>%
  group_by(Family, Group) %>%
  summarize(Total.Relative.Abundance = (sum(Relative.Abundance)))

summarized.abundance.family$Family <- as.character(replace_na(summarized.abundance.family$Family, replace = "Unidentified"))

# Now, I'll collect taxa names for taxa with relative abundance above 1% in any group. These
# names will be preserved in each group.
summarized.abundance.family.above1 <- filter(summarized.abundance.family, Total.Relative.Abundance >= 1.0)
summarized.abundance.family.above1$Family <- factor(summarized.abundance.family.above1$Family, 
                                                    levels = sort(unique(summarized.abundance.family.above1$Family)))
abundant.families <- as.character(unique(summarized.abundance.family.above1$Family)) 
table(abundant.families) # There are 13 Families (and 1 group of unidentified families) that have a Total Relative 
# Abundance above 1%

# Make new data frames either including or excluding abundant families.
summarized.abundance.family.no.remainder <- filter(summarized.abundance.family, 
                                                   Family %in% abundant.families)

# I'll now summarize the remainder samples that compose less than 1% of relative abundance per group.
summarized.abundance.family.remainders.only <- anti_join(summarized.abundance.family, summarized.abundance.family.no.remainder, by = "Family")

summarized.abundance.family.remainders.only.summary <- summarized.abundance.family.remainders.only %>%
  group_by(Group) %>%
  summarize(Total.Relative.Abundance = (sum(Total.Relative.Abundance))) %>%
  mutate(Family = "Families < 1%", .before = 1)


# Now I'll combine abundant families with the low abundance remainders to account for 100% of total 
# abundance. This data frame will be used to generate the community composition figure.

summarized.abundance.family.figure <- rbind(summarized.abundance.family.no.remainder, summarized.abundance.family.remainders.only.summary)
summarized.abundance.family.figure$Family <- factor(summarized.abundance.family.figure$Family, levels = c("[Eubacterium] coprostanoligenes group", "Acholeplasmataceae", "Akkermansiaceae", "Clostridiaceae", "Erysipelotrichaceae", "Lachnospiraceae", 
                                                                                                          "Lactobacillaceae", "Marinifilaceae", "Muribaculaceae", "Oscillospiraceae", "Peptostreptococcaceae", "Rikenellaceae", "Ruminococcaceae", 
                                                                                                          "Unidentified", "Families < 1%"))

ggplot(summarized.abundance.family.figure, aes(x = Group, y = Total.Relative.Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_bw()


ggsave("./results/figures/mucin/Mucin Experiment_community_composition_family.png", width = 6, height = 4)



#### FAMILY LEVEL DIFFERENCES: BASELINE VS 3 WEEK STRESS ############################################
# I now want to plot these differences by Family group. I'll only be plotting the abundant families
# (those with a relative abundance above 1% in at least one Group). I'll start by taking the family
# table I generated previously and filtering it to include only these abundant families.

family.table.mucin.abundant <- filter(family.table.mucin, Family %in% abundant.families)

# I'll look at these comparisons in pairs to keep the figure from getting too complicated. I'll start
# by looking at the effect of stress alone, looking at Baseline vs. 3 week stress samples. I'll
# filter to include only these samples.

family.table.mucin.abundant.stress <- filter(family.table.mucin.abundant, Group == "Baseline" | Group == "3 Wk Stress")
family.table.mucin.abundant.stress$Group <- factor(family.table.mucin.abundant.stress$Group, levels = c("Baseline", "3 Wk Stress"))

# I'll now plot abundance for each sample. Abundances are expressed as a percentage of each sample,
# so I'll format the x axis as a percentage.

ggplot(family.table.mucin.abundant.stress, aes(x = Abundance, y = Family, color = Group)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798"), 
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress")) +
  labs(x = "Relative Abundance", y = NULL) 
 
ggsave("./results/figures/mucin/Mucin Experiment_family abundance_mucin.png", width = 6, height = 5)




#### FAMILY LEVEL DIFFERENCES: STRESS ALONE VS STRESS + MUCIN ########################################
# Now I'll repeat this process, comparing the 4 week samples. I want to see if there are any major
# differences in families between the Control group (Stress Alone) and the Mucin-supplemented group
# (Stress + Mucin). I'll filter to include only these samples.

family.table.mucin.abundant.mucin <- filter(family.table.mucin.abundant, Group == "Stress Alone" | Group == "Stress + Mucin")
family.table.mucin.abundant.mucin$Group <- factor(family.table.mucin.abundant.mucin$Group, levels = c("Stress Alone", "Stress + Mucin"))

# I'll now plot abundance for each sample. Abundances are expressed as a percentage of each sample,
# so I'll format the x axis as a percentage.

ggplot(family.table.mucin.abundant.mucin, aes(x = Abundance, y = Family, color = Group)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(name = "Group", 
                     values = c("#521b92", "#4d8e00"), 
                     breaks = c("Stress Alone", "Stress + Mucin"), 
                     labels = c("Stress Alone", "Stress + Mucin")) +
  labs(x = "Relative Abundance", y = NULL) 

ggsave("./results/figures/mucin/Mucin Experiment_family abundance_stress.png", width = 6, height = 5)




#### COMMUNITY COMPOSITION: ORDER LEVEL #############################################################
# The Gautier lab also wanted to see what relative abundances looked like at the order level. I'll
# now calculate relative abundance values for each order level taxa and organize it into a table.
# For this analysis I won't be grouping Orders that are less than 1% relative abundance together.
# I'll also be using the "Group" column to generate relative abundance.


ps.mucin.prop.order <- tax_glom(ps.mucin.prop, "Order", NArm = FALSE)
order.table <- psmelt(ps.mucin.prop.order) # Organize in long format for ggplot.

table(ps.mucin.prop@sam_data$Group) # There are 24 samples for "Baseline" and "3 Wk Stress" groups, 12
# samples for the "Stress Alone" group, and 11 samples for the "Stress + Mucin" group. I'll 
# separate these groups to generate relative abundances as a percentage by dividing by the total
# number of samples in each group.

order.table.mucin.baseline.3wk <- filter(order.table, Group == "Baseline" | Group == "3 Wk Stress")
order.table.mucin.baseline.3wk$Relative.Abundance <- NA
order.table.mucin.baseline.3wk$Relative.Abundance <- order.table.mucin.baseline.3wk$Abundance / 24 *100


order.table.stress <- filter(order.table, Group == "Stress Alone")
order.table.stressRelative.Abundance <- NA
order.table.stress$Relative.Abundance <- order.table.stress$Abundance / 12 *100


order.table.stress.mucin <- filter(order.table, Group == "Stress + Mucin")
order.table.stress.mucinRelative.Abundance <- NA
order.table.stress.mucin$Relative.Abundance <- order.table.stress.mucin$Abundance / 11 *100


order.table.mucin <- rbind(order.table.mucin.baseline.3wk, order.table.stress)
order.table.mucin <- rbind(order.table.mucin, order.table.stress.mucin)


order.table.mucin <- as.data.frame(order.table.mucin)

# First, I'll summarize relative abundance information per taxonomic group.
summarized.abundance.order <- order.table.mucin %>%
  group_by(Order, Group) %>%
  summarize(Total.Relative.Abundance = (sum(Relative.Abundance)))

summarized.abundance.order$Order <- as.character(replace_na(summarized.abundance.order$Order, replace = "Unidentified"))

# Now I'll organize the data out of a wide format, so that relative abundance values for each 
# condition are their own column.

summarized.abundance.order <- filter(summarized.abundance.order, Total.Relative.Abundance > 0)

summarized.abundance.order.baseline <- filter(summarized.abundance.order, Group == "Baseline")
summarized.abundance.order.baseline <- plyr::rename(summarized.abundance.order.baseline, replace = c("Total.Relative.Abundance" = "Baseline Relative Abundance (%)"))
summarized.abundance.order.baseline <- select(summarized.abundance.order.baseline, -Group)

summarized.abundance.order.3wk <- filter(summarized.abundance.order, Group == "3 Wk Stress")
summarized.abundance.order.3wk <- plyr::rename(summarized.abundance.order.3wk, replace = c("Total.Relative.Abundance" = "3 Wk Stress Relative Abundance (%)"))
summarized.abundance.order.3wk <- select(summarized.abundance.order.3wk, -Group)

summarized.abundance.order.stress <- filter(summarized.abundance.order, Group == "Stress Alone")
summarized.abundance.order.stress <- plyr::rename(summarized.abundance.order.stress, replace = c("Total.Relative.Abundance" = "Stress Alone Relative Abundance (%)"))
summarized.abundance.order.stress <- select(summarized.abundance.order.stress, -Group)

summarized.abundance.order.mucin <- filter(summarized.abundance.order, Group == "Stress + Mucin")
summarized.abundance.order.mucin <- plyr::rename(summarized.abundance.order.mucin, replace = c("Total.Relative.Abundance" = "Stress + Mucin Relative Abundance (%)"))
summarized.abundance.order.mucin <- select(summarized.abundance.order.mucin, -Group)

summarized.abundance.order.wide <- full_join(summarized.abundance.order.baseline, summarized.abundance.order.3wk, by = "Order")
summarized.abundance.order.wide <- full_join(summarized.abundance.order.wide, summarized.abundance.order.stress, by = "Order")
summarized.abundance.order.wide <- full_join(summarized.abundance.order.wide, summarized.abundance.order.mucin, by = "Order")

summarized.abundance.order.wide <- replace_na(summarized.abundance.order.wide, list("Baseline Relative Abundance (%)" = 0, 
                                                                                     "3 Wk Stress Relative Abundance (%)" = 0,
                                                                                     "Stress Alone Relative Abundance (%)" = 0,
                                                                                     "Stress + Mucin Relative Abundance (%)" = 0))


# Get additional taxonomic information
taxa.order <- order.table.mucin %>%
  select(Kingdom, Phylum, Class, Order) %>%
  group_by(Order) %>%
  slice(1L) %>%
  ungroup()

summarized.abundance.order.wide <- left_join(summarized.abundance.order.wide, taxa.order, by = "Order")
summarized.abundance.order.wide <- select(summarized.abundance.order.wide, Kingdom, Phylum, Class, Order, everything())

write.csv(summarized.abundance.order.wide, file = "results/tables/mucin/Order Level Relative Abundance.csv", row.names = FALSE)
