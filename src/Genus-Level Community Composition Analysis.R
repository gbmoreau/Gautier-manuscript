#### GENUS-LEVEL COMMUNITY COMPOSITION ANALYSIS ####################################################
# NAME: G. Brett Moreau
# DATE: April 26, 2023

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


#### INTRODUCTION ##################################################################################
# This code will outline the analysis of community composition for the 16S samples for the Gautier 
# manuscript. This code will calculate relative abundance and statistics at the Genus level.



#### ORGANIZATION OF THE DATA SET ##################################################################
# The first step is to organize the data as relative abundance at the Genus level.

# Load the phyloseq data
load(file = "../results/phyloseq object for analysis.RData")

# Convert data to relative abundances
ps.analysis.prop <-transform_sample_counts(ps.analysis, function(ASV) ASV/sum(ASV))

# Organize data based on Genus level
ps.analysis.prop.genus <- tax_glom(ps.analysis.prop, "Genus", NArm = FALSE)
genus.table <- psmelt(ps.analysis.prop.genus) # Organize in long format for ggplot.
genus.table$Genus[is.na(genus.table$Genus)] <- "Unidentified"

length(unique(genus.table$Genus)) # There are 60 genera (+1 unidentified group) in this data set.




#### SUMMARIZE TAXA <1% IN ANY GROUP ################################################################
# I'll summarize Genera with < 1% abundance in either group

#Collect taxa names for families >1% in any sample
abundant.genera <- genus.table %>%
  group_by(Genus, Group) %>%
  summarize(Abundance = sum(Abundance)) %>%
  filter(Abundance >= .01) %>%
  select(Genus) %>%
  unique()

abundant.genera <- as.character(abundant.genera$Genus) # There are 39 abundant genera

# Add Genera <1% designation as a separate column.
genus.table$Genus.per.group <- genus.table$Genus
genus.table$Genus.per.group[! genus.table$Genus %in% abundant.genera] <- "Genera <1%"




#### GENERATE PER-SAMPLE ABUNDANCE TABLES ###########################################################

# Summarize Genus abundance per sample
genus.table.sample.export <- genus.table %>%
  group_by(Sample, Genus) %>%
  summarize(Abundance = sum(Abundance) * 100)

# Reorganize into a long format, with Samples as rows and Genera as columns.
genus.table.sample.export <- pivot_wider(genus.table.sample.export, names_from = Genus, values_from = Abundance)

# Add back sample name information
sample.info <- select(genus.table, Sample, Group, Animal.Number)
sample.info <- unique(sample.info)
genus.table.sample.export <- left_join(genus.table.sample.export, sample.info, by = "Sample")
genus.table.sample.export <- relocate(genus.table.sample.export, c(Group, Animal.Number), .after = Sample)

#write.csv(genus.table.sample.export, "../results/tables/per-sample relative abundance table_genus.csv", row.names = FALSE)




#### VISUALIZE DIFFERENCES IN GENUS ABUNDANCE BETWEEN GROUPS ########################################

### BOXPLOT OF GENUS-LEVEL COMMUNITY COMPOSITION ###
ggplot(genus.table, aes(x = Abundance, y = Genus, color = Group)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(limits = rev) +
  labs(x = "Relative Abundance (%)", y = NULL, color = "Group") +
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress"))

#ggsave("../results/figures/per-group boxplot_genus.png", width = 6, height = 6, dpi = 300)



# Perform Wilcoxon rank sum test on all genera in the data set.
pvalue <- lapply(genus.table.sample.export[, 4:64], function(x) wilcox.test(x ~ genus.table.sample.export$Group)$p.value)

# Extract names
names <- names(pvalue)
names <- names[1:61]

# Extract values
values <- as.numeric(unname(unlist(pvalue)))
values <- values[! is.na(values)]

# Make combined table and add correction for multiple comparisons.
Wilcox.all <- cbind(names, values)
Wilcox.all <- as.data.frame(Wilcox.all)
Wilcox.all$Adjusted_p <- p.adjust(Wilcox.all$values, method = "bonferroni", n = 61)
Wilcox.all <- plyr::rename(Wilcox.all, replace = c("names" = "Genus",
                                                   "values" = "p value",
                                                   "Adjusted_p" = "Bonferroni-adjusted p value"))

#write.csv(Wilcox.all, file = "../results/tables/wilcoxon rank comparison-genus composition.csv", row.names = FALSE)
