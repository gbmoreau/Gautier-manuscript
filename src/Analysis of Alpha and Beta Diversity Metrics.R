#### ANALYSIS OF ALPHA AND DIVERSITY METRICS #######################################################
# NAME: G. Brett Moreau
# DATE: February 2, 2022

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
packageVersion("vegan") # I'm using version 2.6.2


#### INTRODUCTION #################################################################################
# This code will outline the analysis of alpha and beta diversity for the 16S samples for the 
# Gautier manuscript.These include richness and evenness as alpha diversity measures and Bray-
# Curtis Dissimilarity for beta diversity.

# Load the data
load(file = "../results/phyloseq object for analysis.RData")




#### ALPHA DIVERSITY METRICS ######################################################################
# PURPOSE: Generate Richness and Evenness plots for this data.

### RICHNESS ###
plot_richness(ps.analysis, x = "Group", measures = c("Observed"), 
              color = "Group") + 
  geom_point() +
  geom_boxplot() + 
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)) +
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress"))


#ggsave("../results/figures/richness-observed.png", width = 5, height = 4)

richness <- estimate_richness(ps.analysis)
pairwise.wilcox.test(richness$Observed, p.adjust.method = "bonferroni", sample_data(ps.analysis)$Group)

# Richness is not significantly different between groups (p = 0.28).


### EVENNESS ###
# First I'll organize the data.
sample.data <- ps.analysis@sam_data
sample.data <- cbind("Sample.Names" = rownames(sample.data), sample.data)

evenness<- evenness(ps.analysis, index = c("pielou"))
evenness<- cbind("Sample.Names" = rownames(evenness), evenness)

sample.data <- full_join(sample.data, evenness, by = "Sample.Names")


### PIELOU ###
ggplot(sample.data, aes(x = Group, y = pielou, color = Group)) +
  geom_point() +
  geom_boxplot() +
  theme_bw() +
  labs(x = NULL, y = "Pielou Eveness Index") +
  theme(legend.position = "none") +
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress"))

#ggsave("../results/figures/evenness-pielou.png", width = 5, height = 4)

pairwise.wilcox.test(evenness$pielou, p.adjust.method = "bonferroni", sample_data(ps.analysis)$Group)

# Evenness is not significantly different between groups (p = 0.27).




#### BETA DIVERSITY ###########################################################################
# PURPOSE: Generate NMDS plot of Bray-Curtis Dissimilarity

set.seed(966) # Set seed for reproducibility

# First, let's ordinate Bray-Curtis Dissimilarity using NMDS, then plot this ordination to 
# visualize separation between different samples.

ord.nmds.bray <- ordinate(ps.analysis, method = "NMDS", distance = "bray")

# NMDS Plot
ordplot <- plot_ordination(ps.analysis, ord.nmds.bray, color = "Group") +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress"))

ordplot

#ggsave("../results/figures/beta-diversity-bc.png", width = 4, height = 3)


# Add ordination ellipses
ordplot +
  stat_ellipse(type = "norm") +
  theme_bw() +
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))  +
  coord_fixed()

#ggsave("../results/figures/beta-diversity-bc-ellipses.png", width = 4, height = 3)

# Stats
adonis2(distance(ps.analysis, method = "bray") ~sample_data(ps.analysis)$Group)

# There is a significant difference between groups by PERMANOVA (p = 0.001)


# Overall, this analysis saw no differences in alpha diversity (richness or evenness) between
# baseline and stressed groups. However, there was a significant difference in Bray-Curtis
# Dissimilarity between groups.