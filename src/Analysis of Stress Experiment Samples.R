#### ANALYSIS OF GAUTIER 16S SAMPLES FOR PILOT PROJECT #############################################
# NAME: G. Brett Moreau
# DATE: November 9, 2021

#### PACKAGES ######################################################################################

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
# Samples have been run through the DADA2 pipeline and are ready for analysis. This analysis will
# focus on the Stress experiment, which examined the role of stress on intestinal microbiota
# diversity and composition.




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
# Stressed) and two different timepoints (pre- and post-stress intervention). Each treatment group 
# and timepoint is divided into 6 samples each.

dim(ps.stress@otu_table) # From these 24 samples, there were 1795 unique ASVs.

total.reads <- sample_sums(ps.stress)

sum(total.reads) # There are a total of 1,741,609 reads across the 24 samples in the data set.
mean(total.reads) # The average read number is 72,567
median(total.reads) # The median read number is 71,421
range(total.reads) # The range of reads in samples is from 43,957-117,986 total reads.


# Before I perform downstream analyses, I'll first convert the data into relative abundances. I'll 
# also divide by the number of samples because they are identical in each group, resulting in 
# values that represent percent of total for each taxa.

ps.stress@sam_data$Condition <- factor(ps.stress@sam_data$Condition, 
                                       levels = c("Naïve Pre", "Stress Pre", "Naïve Post", "Stress Post"))
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
                       breaks = c("Naïve Pre", "Stress Pre", "Naïve Post", "Stress Post"), 
                       labels = c("Naïve Pre", "Stress Pre", "Naïve Post", "Stress Post")) 


#ggsave("./results/figures/stress/Stress Experiment_richness_observed.png", width = 4, height = 4)

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
  labs(title = "Microbial Evenness Across Groups", x = NULL, y = "Pielou Eveness Index") +
  theme(legend.position = "none")

#ggsave("./results/figures/stress/Stress Experiment_evenness_pielou.png", width = 4, height = 4)

# Overall, there doesn't appear to be a whole lot of separation in evenness between groups. 
# There is a lot of spread between individuals within each group, which makes the data a 
# little difficult to interpret. Now I'll look for statistical differences.


pairwise.wilcox.test(evenness$pielou, p.adjust.method = "bonferroni", sample_data(ps.stress)$Condition)

# No significant differences were observed in evenness between groups, either.




#### BETA DIVERSITY: BRAY-CURTIS DISSIMILARITY  #####################################################
# I'll now move on to measures of beta diversity, which investigates variation between samples.
# I'll start by looking at Bray-Curtis Dissimilarity.

#sample(1:1000, 1) # It selected 763
set.seed(763) # Set seed for reproducibility

# First, let's ordinate Bray-Curtis Dissimilarity using NMDS, then plot this ordination to visualize
# separation between different samples.

ord.nmds.bray <- ordinate(ps.stress.prop, method = "NMDS", distance = "bray")

plot_ordination(ps.stress.prop, ord.nmds.bray, color = "Condition", 
                title = "Bray-Curtis Dissimilarity") +
  theme_bw() +
  scale_color_discrete(name = "Condition", 
                       breaks = c("Naïve Pre", "Stress Pre", "Naïve Post", "Stress Post"), 
                       labels = c("Naïve Pre", "Stress Pre", "Naïve Post", "Stress Post")) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) 

#ggsave("./results/figures/stress/Stress Experiment_beta_diversity_bc.png", width = 5, height = 5)


# Bray-Curtis distance is ordinated using Non-metric multidimensional scaling (NMDS). It looks like
# while there is overlap between groups, groups are clustering somewhat distinctly from each other.
# In particular, the Stress Post group seems to be distinct from other groups, while the other 3
# groups are more mixed amongst each other.To better explore this, I'll perform a PERMANOVA to see 
# if there are any differences in clustering between groups.

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
# all comaprisons except for Naive Pre vs. Naive Post and Naive Pre vs Stress Pre. Naive Pre
# and Post should be identical (no differences in treatment between these groups) and Naive Pre
# and Stress Pre should be identical (neither has received stress treatment), so these results
# make sense.




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
                       breaks = c("Naïve Pre", "Stress Pre", "Naïve Post", "Stress Post"), 
                       labels = c("Naïve Pre", "Stress Pre", "Naïve Post", "Stress Post")) +
  labs(title = "Unweighted Unifrac Distance") +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) 

#ggsave("./results/figures/stress/Stress Experiment_beta_diversity_unifrac.png", width = 5, height = 5)


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

# Pairwise comparisons of PERMANOVA from Unifrac distance found similar results to those seen
# with Bray-Curtis dissimilarity. Specifically, the Stress Post group was statistically 
# significant from both pre-intervention groups, as hypothesized. However, it was not significantly
# different from the Naive Post group. Overall, beta diversity measures found very similar results, 
# suggesting that there are some differences between groups.




### COMMUNITY COMPOSITION: PHYLUM LEVEL #############################################################
# There appear to be some broad differences in diversity between groups. I'll now see if this is
# reflected by any broad changes in community composition across groups. I'll start by looking at
# the phylum level, looking at all phyla that make up at least 1% of the relative abundance in at
# least one group.

ps.stress.prop.phylum <- tax_glom(ps.stress.prop, "Phylum", NArm = FALSE)
phylum.table.stress <- psmelt(ps.stress.prop.phylum) # Organize in long format for ggplot.

table(ps.stress.prop@sam_data$Condition)

# There are an equal number of samples (n=6) in each group. Because of this, the 
# Abundance column is already formatted to generate relative abundance as a percent
# of the total condition, so I'll just use this for analysis.

# I'll plot the data to get an overview of differences between groups.
ggplot(phylum.table.stress, aes(x = Condition, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") 

# I want to clean this figure up a little bit by including only phyla that make up
# at least 1% in at least one group.

# First, I'll summarize relative abundance information per taxonomic group.
summarized.abundance.phylum <- phylum.table.stress %>%
  group_by(Phylum, Condition) %>%
  summarize(Total.Relative.Abundance = (sum(Abundance)))

summarized.abundance.phylum$Phylum <- as.character(replace_na(summarized.abundance.phylum$Phylum, replace = "Unidentified"))

# Now, I'll collect taxa names for taxa with relative abundance above 1% in any group. These
# names will be preserved in each group.
summarized.abundance.phylum.above1 <- filter(summarized.abundance.phylum, Total.Relative.Abundance >= 1.0)
summarized.abundance.phylum.above1$Phylum <- factor(summarized.abundance.phylum.above1$Phylum, 
                                                    levels = sort(unique(summarized.abundance.phylum.above1$Phylum)))
abundant.families <- as.character(unique(summarized.abundance.phylum.above1$Phylum)) 
table(abundant.families) # There are 4 Phyla that have a Total Relative Abundance above 1%

# Make new data frames either including or excluding abundant families.
summarized.abundance.phylum.no.remainder <- filter(summarized.abundance.phylum, 
                                                   Phylum %in% abundant.families)

# I'll now summarize the remainder samples that compose less than 1% of relative abundance per group.
summarized.abundance.phylum.remainders.only <- anti_join(summarized.abundance.phylum, summarized.abundance.phylum.no.remainder, by = "Phylum")

summarized.abundance.phylum.remainders.only.summary <- summarized.abundance.phylum.remainders.only %>%
  group_by(Condition) %>%
  summarize(Total.Relative.Abundance = (sum(Total.Relative.Abundance))) %>%
  mutate(Phylum = "Phyla < 1%", .before = 1)


# Now I'll combine abundant families with the low abundance remainders to account for 100% of total 
# abundance. This data frame will be used to generate the community composition figure.

summarized.abundance.phylum.figure <- rbind(summarized.abundance.phylum.no.remainder, summarized.abundance.phylum.remainders.only.summary)
summarized.abundance.phylum.figure$Phylum <- factor(summarized.abundance.phylum.figure$Phylum, levels = c("Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota", "Phyla < 1%"))

ggplot(summarized.abundance.phylum.figure, aes(x = Condition, y = Total.Relative.Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = NULL, y = "Relative Abundance (%)")

#ggsave("./results/figures/stress/stress Experiment_community_composition_phylum.png", width = 5, height = 4)




### COMMUNITY COMPOSITION: FAMILY LEVEL #############################################################
# There appear to be some differences at the phylum level. I'll now look at the family level to see
# how these differences are reflected here.

ps.stress.prop.family <- tax_glom(ps.stress.prop, "Family", NArm = FALSE)
family.table.stress <- psmelt(ps.stress.prop.family) # Organize in long format for ggplot.

table(ps.stress@sam_data$Condition) 

# Since all groups have equal number of samples, I'll just use the Abundance column.

# I'll plot the data to get an overview of differences between groups.
ggplot(family.table.stress, aes(x = Condition, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") 

# This figure is pretty cluttered. To clean it up, I'll filter out all Families that make up less than 1
# % of total in any group and place them in their own category.


# First, I'll summarize relative abundance information per taxonomic group.
summarized.abundance.family <- family.table.stress %>%
  group_by(Family, Condition) %>%
  summarize(Total.Relative.Abundance = (sum(Abundance)))

summarized.abundance.family$Family <- as.character(replace_na(summarized.abundance.family$Family, replace = "Unidentified"))

# Now, I'll collect taxa names for taxa with relative abundance above 1% in any group. These
# names will be preserved in each group.
summarized.abundance.family.above1 <- filter(summarized.abundance.family, Total.Relative.Abundance >= 1.0)
summarized.abundance.family.above1$Family <- factor(summarized.abundance.family.above1$Family, 
                                                    levels = sort(unique(summarized.abundance.family.above1$Family)))
abundant.families <- as.character(unique(summarized.abundance.family.above1$Family)) 
table(abundant.families) # There are 11 Families (and 1 group of unidentified families) that have a Total Relative 
# Abundance above 1%

# Make new data frames either including or excluding abundant families.
summarized.abundance.family.no.remainder <- filter(summarized.abundance.family, 
                                                   Family %in% abundant.families)

# I'll now summarize the remainder samples that compose less than 1% of relative abundance per group.
summarized.abundance.family.remainders.only <- anti_join(summarized.abundance.family, summarized.abundance.family.no.remainder, by = "Family")

summarized.abundance.family.remainders.only.summary <- summarized.abundance.family.remainders.only %>%
  group_by(Condition) %>%
  summarize(Total.Relative.Abundance = (sum(Total.Relative.Abundance))) %>%
  mutate(Family = "Families < 1%", .before = 1)


# Now I'll combine abundant families with the low abundance remainders to account for 100% of total 
# abundance. This data frame will be used to generate the community composition figure.

summarized.abundance.family.figure <- rbind(summarized.abundance.family.no.remainder, summarized.abundance.family.remainders.only.summary)
summarized.abundance.family.figure$Family <- factor(summarized.abundance.family.figure$Family, levels = c("Akkermansiaceae", "Bacteroidaceae","Clostridiaceae", "Erysipelotrichaceae", "Lachnospiraceae", 
                                                                                                          "Lactobacillaceae", "Muribaculaceae", "Oscillospiraceae", "Peptostreptococcaceae", "Ruminococcaceae", 
                                                                                                          "Sutterellaceae", "Unidentified", "Families < 1%"))

ggplot(summarized.abundance.family.figure, aes(x = Condition, y = Total.Relative.Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = NULL, y = "Relative Abundance (%)")


#ggsave("./results/figures/stress/stress Experiment_community_composition_family.png", width = 6, height = 4)



#### FAMILY LEVEL DIFFERENCES ########################################################################
# I now want to plot these differences by Family group. I'll only be plotting the abundant families
# (those with a relative abundance above 1% in at least one Condition).

family.table.stress.abundant <- filter(family.table.stress, Family %in% abundant.families)
family.table.stress.abundant$Abundance <- family.table.stress.abundant$Abundance *10

# I'll now plot abundance for each sample. Abundances are expressed as a percentage of each sample,
# so I'll format the x axis as a percentage.

ggplot(family.table.stress.abundant, aes(x = Abundance, y = Family, color = Condition)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_discrete(limits = rev) +
  labs(x = "Relative Abundance (%)", y = NULL) 

#ggsave("./results/figures/stress/stress Experiment_family abundance_stress.png", width = 6, height = 5)




#### COMMUNITY COMPOSITION: ORDER LEVEL #############################################################
# The Gautier lab also wanted to see what relative abundances looked like at the order level. I'll
# now calculate relative abundance values for each order level taxa and organize it into a table.
# For this analysis I won't be grouping Orders that are less than 1% relative abundance together.
# I'll also be using the "Group" column to generate relative abundance.


ps.stress.prop.order <- tax_glom(ps.stress.prop, "Order", NArm = FALSE)
order.table.stress <- psmelt(ps.stress.prop.order) # Organize in long format for ggplot.

table(ps.stress.prop@sam_data$Condition) 

# Since all groups have equal number of samples, I'll just use the Abundance column.


# First, I'll summarize relative abundance information per taxonomic group.
summarized.abundance.order <- order.table.stress %>%
  group_by(Order, Condition) %>%
  summarize(Total.Relative.Abundance = (sum(Abundance)))

summarized.abundance.order$Order <- as.character(replace_na(summarized.abundance.order$Order, replace = "Unidentified"))

# Now I'll organize the data out of a wide format, so that relative abundance values for each 
# condition are their own column.

summarized.abundance.order <- filter(summarized.abundance.order, Total.Relative.Abundance > 0)

summarized.abundance.order.naive.pre <- filter(summarized.abundance.order, Condition == "Naïve Pre")
summarized.abundance.order.naive.pre <- plyr::rename(summarized.abundance.order.naive.pre, replace = c("Total.Relative.Abundance" = "Naive Pre Relative Abundance (%)"))
summarized.abundance.order.naive.pre <- select(summarized.abundance.order.naive.pre, -Condition)

summarized.abundance.order.stress.pre <- filter(summarized.abundance.order, Condition == "Stress Pre")
summarized.abundance.order.stress.pre <- plyr::rename(summarized.abundance.order.stress.pre, replace = c("Total.Relative.Abundance" = "Stress Pre Relative Abundance (%)"))
summarized.abundance.order.stress.pre <- select(summarized.abundance.order.stress.pre, -Condition)

summarized.abundance.order.naive.post <- filter(summarized.abundance.order, Condition == "Naïve Post")
summarized.abundance.order.naive.post <- plyr::rename(summarized.abundance.order.naive.post, replace = c("Total.Relative.Abundance" = "Naive Post Relative Abundance (%)"))
summarized.abundance.order.naive.post <- select(summarized.abundance.order.naive.post, -Condition)

summarized.abundance.order.stress.post <- filter(summarized.abundance.order, Condition == "Stress Post")
summarized.abundance.order.stress.post <- plyr::rename(summarized.abundance.order.stress.post, replace = c("Total.Relative.Abundance" = "Stress Post Relative Abundance (%)"))
summarized.abundance.order.stress.post <- select(summarized.abundance.order.stress.post, -Condition)

summarized.abundance.order.wide <- full_join(summarized.abundance.order.naive.pre, summarized.abundance.order.stress.pre, by = "Order")
summarized.abundance.order.wide <- full_join(summarized.abundance.order.wide, summarized.abundance.order.naive.post, by = "Order")
summarized.abundance.order.wide <- full_join(summarized.abundance.order.wide, summarized.abundance.order.stress.post, by = "Order")

summarized.abundance.order.wide <- replace_na(summarized.abundance.order.wide, list("Naive Pre Relative Abundance (%)" = 0, 
                                                                                    "Stress Pre Relative Abundance (%)" = 0,
                                                                                    "Naive Post Relative Abundance (%)" = 0,
                                                                                    "Stress Post Relative Abundance (%)" = 0))


# Get additional taxonomic information
taxa.order <- order.table.stress %>%
  select(Kingdom, Phylum, Class, Order) %>%
  group_by(Order) %>%
  slice(1L) %>%
  ungroup()

summarized.abundance.order.wide <- left_join(summarized.abundance.order.wide, taxa.order, by = "Order")
summarized.abundance.order.wide <- select(summarized.abundance.order.wide, Kingdom, Phylum, Class, Order, everything())

#write.csv(summarized.abundance.order.wide, file = "results/tables/stress/Order Level Relative Abundance.csv", row.names = FALSE)

