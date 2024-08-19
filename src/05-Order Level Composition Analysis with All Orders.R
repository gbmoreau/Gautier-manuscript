#### ORDER-LEVEL COMPOSITION ANALYSIS WITH ALL ORDERS ##############################################
# NAME: G. Brett Moreau
# DATE: March 23, 2022

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


#### INTRODUCTION #################################################################################
# In addition to the order-level community composition I performed previously, I also want to look
# at this with all orders, not just those that make up at least 1% of the total relative abundance.


# Load the data
load(file = "../results/phyloseq object for analysis.RData")




#### ORDER LEVEL DISTRIBUTION #####################################################################

# Organize phyloseq object based on Order-level classifications
ps.analysis.prop <-transform_sample_counts(ps.analysis, function(ASV) ASV/sum(ASV))

ps.analysis.prop.order <- tax_glom(ps.analysis.prop, "Order", NArm = FALSE)
order.table <- psmelt(ps.analysis.prop.order) # Organize in long format for ggplot.

# First, I'll compare the number of Orders in the data set with and without excluding Orders < 1%
summarized.abundance.order <- order.table %>%
  group_by(Order, Group) %>%
  summarize(Relative.Abundance = (sum(Abundance) / 24))

summarized.abundance.order$Order <- as.character(replace_na(summarized.abundance.order$Order, replace = "Unidentified"))


summarized.abundance.order.above1 <- filter(summarized.abundance.order, Relative.Abundance >= 0.01)
summarized.abundance.order.above1$Order <- factor(summarized.abundance.order.above1$Order, 
                                                  levels = sort(unique(summarized.abundance.order.above1$Order)))
abundant.orders <- as.character(unique(summarized.abundance.order.above1$Order)) 
length(abundant.orders) # There are 8 Orders that have a Total Relative Abundance above 1%

abundant.orders.all <- as.character(unique(summarized.abundance.order$Order))
length(abundant.orders.all) # There are 23 Orders in the data set overall.


# Make new data frame with sample-level data.
summarized.abundance.order.sample <- order.table %>%
  group_by(Order, Sample) %>%
  summarize(Relative.Abundance = (sum(Abundance)))

summarized.abundance.order.sample$Order <- as.character(replace_na(summarized.abundance.order.sample$Order, 
                                                                   replace = "Unidentified"))

# Add back sample and group metadata
Sample.Group <- select(order.table, Sample, Group, Animal.Number)
Sample.Group <- Sample.Group %>%
  group_by(Sample) %>%
  slice(1L)

summarized.abundance.order.sample.table <- full_join(summarized.abundance.order.sample, Sample.Group, by = "Sample")
summarized.abundance.order.sample.table <- ungroup(summarized.abundance.order.sample.table)




#### STATISTICAL DIFFERENCES BETWEEN ORDERS IN CONTROL AND STRESSED SAMPLES  #######################

# Here is the list of the 23 Orders in the data set.
table(abundant.orders.all)

# I'll now use a Wilcoxon rank test and Bonferroni correction for multiple comparisons to identify 
#significant differences between each of these order groups. 


# Calculate Wilcoxon Ranked Sum p values for each Order.

# ACHOLEPLASMATALES
wilcox.ach <- filter(summarized.abundance.order.sample.table, Order == "Acholeplasmatales")
wilcox.ach.baseine <- filter(wilcox.ach, Group == "Baseline")
wilcox.ach.stress <- filter(wilcox.ach, Group == "3 Wk Stress")
wilcox.test.ach <- wilcox.test(wilcox.ach.baseine$Relative.Abundance, 
                               wilcox.ach.stress$Relative.Abundance)$p.value


# BACILLALES
wilcox.bac <- filter(summarized.abundance.order.sample.table, Order == "Bacillales")
wilcox.bac.baseine <- filter(wilcox.bac, Group == "Baseline")
wilcox.bac.stress <- filter(wilcox.bac, Group == "3 Wk Stress")
wilcox.test.bac <- wilcox.test(wilcox.bac.baseine$Relative.Abundance, 
                               wilcox.bac.stress$Relative.Abundance)$p.value


# BACTEROIDALES
wilcox.bact <- filter(summarized.abundance.order.sample.table, Order == "Bacteroidales")
wilcox.bact.baseine <- filter(wilcox.bact, Group == "Baseline")
wilcox.bact.stress <- filter(wilcox.bact, Group == "3 Wk Stress")
wilcox.test.bact <- wilcox.test(wilcox.bact.baseine$Relative.Abundance, 
                               wilcox.bact.stress$Relative.Abundance)$p.value


# BURKHOLDERIALES
wilcox.burk <- filter(summarized.abundance.order.sample.table, Order == "Burkholderiales")
wilcox.burk.baseine <- filter(wilcox.burk, Group == "Baseline")
wilcox.burk.stress <- filter(wilcox.burk, Group == "3 Wk Stress")
wilcox.test.burk <- wilcox.test(wilcox.burk.baseine$Relative.Abundance, 
                               wilcox.burk.stress$Relative.Abundance)$p.value


# CHLAMYDIALES
wilcox.chla <- filter(summarized.abundance.order.sample.table, Order == "Chlamydiales")
wilcox.chla.baseine <- filter(wilcox.chla, Group == "Baseline")
wilcox.chla.stress <- filter(wilcox.chla, Group == "3 Wk Stress")
wilcox.test.chla <- wilcox.test(wilcox.chla.baseine$Relative.Abundance, 
                                wilcox.chla.stress$Relative.Abundance)$p.value


# CHRISTENSENELLALES
wilcox.chris <- filter(summarized.abundance.order.sample.table, Order == "Christensenellales")
wilcox.chris.baseine <- filter(wilcox.chris, Group == "Baseline")
wilcox.chris.stress <- filter(wilcox.chris, Group == "3 Wk Stress")
wilcox.test.chris <- wilcox.test(wilcox.chris.baseine$Relative.Abundance, 
                                wilcox.chris.stress$Relative.Abundance)$p.value


# CLOSTRIDIALES
wilcox.clos <- filter(summarized.abundance.order.sample.table, Order == "Clostridiales")
wilcox.clos.baseine <- filter(wilcox.clos, Group == "Baseline")
wilcox.clos.stress <- filter(wilcox.clos, Group == "3 Wk Stress")
wilcox.test.clos <- wilcox.test(wilcox.clos.baseine$Relative.Abundance, 
                                wilcox.clos.stress$Relative.Abundance)$p.value


# CORIOBACTERIALES
wilcox.cori <- filter(summarized.abundance.order.sample.table, Order == "Coriobacteriales")
wilcox.cori.baseine <- filter(wilcox.cori, Group == "Baseline")
wilcox.cori.stress <- filter(wilcox.cori, Group == "3 Wk Stress")
wilcox.test.cori <- wilcox.test(wilcox.cori.baseine$Relative.Abundance, 
                                wilcox.cori.stress$Relative.Abundance)$p.value


# CORYNEBACTERIALES
wilcox.cory <- filter(summarized.abundance.order.sample.table, Order == "Corynebacteriales")
wilcox.cory.baseine <- filter(wilcox.cory, Group == "Baseline")
wilcox.cory.stress <- filter(wilcox.cory, Group == "3 Wk Stress")
wilcox.test.cory <- wilcox.test(wilcox.cory.baseine$Relative.Abundance, 
                                wilcox.cory.stress$Relative.Abundance)$p.value


# DESULFOVIBRIONALES
wilcox.desulf <- filter(summarized.abundance.order.sample.table, Order == "Desulfovibrionales")
wilcox.desulf.baseine <- filter(wilcox.desulf, Group == "Baseline")
wilcox.desulf.stress <- filter(wilcox.desulf, Group == "3 Wk Stress")
wilcox.test.desulf <- wilcox.test(wilcox.desulf.baseine$Relative.Abundance, 
                                wilcox.desulf.stress$Relative.Abundance)$p.value


# ERYSIPELOTRICHALES
wilcox.erys <- filter(summarized.abundance.order.sample.table, Order == "Erysipelotrichales")
wilcox.erys.baseine <- filter(wilcox.erys, Group == "Baseline")
wilcox.erys.stress <- filter(wilcox.erys, Group == "3 Wk Stress")
wilcox.test.erys <- wilcox.test(wilcox.erys.baseine$Relative.Abundance, 
                                wilcox.erys.stress$Relative.Abundance)$p.value


# EXIGUOBACTERALES
wilcox.exig <- filter(summarized.abundance.order.sample.table, Order == "Exiguobacterales")
wilcox.exig.baseine <- filter(wilcox.exig, Group == "Baseline")
wilcox.exig.stress <- filter(wilcox.exig, Group == "3 Wk Stress")
wilcox.test.exig <- wilcox.test(wilcox.exig.baseine$Relative.Abundance, 
                                wilcox.exig.stress$Relative.Abundance)$p.value


# LACHNOSPIRALES
wilcox.lachno <- filter(summarized.abundance.order.sample.table, Order == "Lachnospirales")
wilcox.lachno.baseine <- filter(wilcox.lachno, Group == "Baseline")
wilcox.lachno.stress <- filter(wilcox.lachno, Group == "3 Wk Stress")
wilcox.test.lachno <- wilcox.test(wilcox.lachno.baseine$Relative.Abundance, 
                                  wilcox.lachno.stress$Relative.Abundance)$p.value


# LACTOBACILLALES
wilcox.lacto <- filter(summarized.abundance.order.sample.table, Order == "Lactobacillales")
wilcox.lacto.baseine <- filter(wilcox.lacto, Group == "Baseline")
wilcox.lacto.stress <- filter(wilcox.lacto, Group == "3 Wk Stress")
wilcox.test.lacto <- wilcox.test(wilcox.lacto.baseine$Relative.Abundance, 
                                 wilcox.lacto.stress$Relative.Abundance)$p.value

# MICROCOCCALES
wilcox.micro <- filter(summarized.abundance.order.sample.table, Order == "Micrococcales")
wilcox.micro.baseine <- filter(wilcox.micro, Group == "Baseline")
wilcox.micro.stress <- filter(wilcox.micro, Group == "3 Wk Stress")
wilcox.test.micro <- wilcox.test(wilcox.micro.baseine$Relative.Abundance, 
                                wilcox.micro.stress$Relative.Abundance)$p.value


# MONOGLOBALES
wilcox.mono <- filter(summarized.abundance.order.sample.table, Order == "Monoglobales")
wilcox.mono.baseine <- filter(wilcox.mono, Group == "Baseline")
wilcox.mono.stress <- filter(wilcox.mono, Group == "3 Wk Stress")
wilcox.test.mono <- wilcox.test(wilcox.mono.baseine$Relative.Abundance, 
                                 wilcox.mono.stress$Relative.Abundance)$p.value


# OSCILLOSPIRALES
wilcox.osci <- filter(summarized.abundance.order.sample.table, Order == "Oscillospirales")
wilcox.osci.baseine <- filter(wilcox.osci, Group == "Baseline")
wilcox.osci.stress <- filter(wilcox.osci, Group == "3 Wk Stress")
wilcox.test.osci <- wilcox.test(wilcox.osci.baseine$Relative.Abundance, 
                                wilcox.osci.stress$Relative.Abundance)$p.value


# PEPTOCOCCALES
wilcox.pepto <- filter(summarized.abundance.order.sample.table, Order == "Peptococcales")
wilcox.pepto.baseine <- filter(wilcox.pepto, Group == "Baseline")
wilcox.pepto.stress <- filter(wilcox.pepto, Group == "3 Wk Stress")
wilcox.test.pepto <- wilcox.test(wilcox.pepto.baseine$Relative.Abundance, 
                                 wilcox.pepto.stress$Relative.Abundance)$p.value


# PEPTOSTREPTOCOCCALES-TISSIERELLALES
wilcox.PT <- filter(summarized.abundance.order.sample.table, Order == "Peptostreptococcales-Tissierellales")
wilcox.PT.baseine <- filter(wilcox.PT, Group == "Baseline")
wilcox.PT.stress <- filter(wilcox.PT, Group == "3 Wk Stress")
wilcox.test.PT <- wilcox.test(wilcox.PT.baseine$Relative.Abundance, 
                                 wilcox.PT.stress$Relative.Abundance)$p.value


# STAPHYLOCOCCALES
wilcox.staph <- filter(summarized.abundance.order.sample.table, Order == "Staphylococcales")
wilcox.staph.baseine <- filter(wilcox.staph, Group == "Baseline")
wilcox.staph.stress <- filter(wilcox.staph, Group == "3 Wk Stress")
wilcox.test.staph <- wilcox.test(wilcox.staph.baseine$Relative.Abundance, 
                                 wilcox.staph.stress$Relative.Abundance)$p.value


# STREPTOSPORANGIALES
wilcox.strep <- filter(summarized.abundance.order.sample.table, Order == "Streptosporangiales")
wilcox.strep.baseine <- filter(wilcox.strep, Group == "Baseline")
wilcox.strep.stress <- filter(wilcox.strep, Group == "3 Wk Stress")
wilcox.test.strep <- wilcox.test(wilcox.strep.baseine$Relative.Abundance, 
                                 wilcox.strep.stress$Relative.Abundance)$p.value


# THERMOMICROBIALES
wilcox.thermo <- filter(summarized.abundance.order.sample.table, Order == "Thermomicrobiales")
wilcox.thermo.baseine <- filter(wilcox.thermo, Group == "Baseline")
wilcox.thermo.stress <- filter(wilcox.thermo, Group == "3 Wk Stress")
wilcox.test.thermo <- wilcox.test(wilcox.thermo.baseine$Relative.Abundance, 
                                 wilcox.thermo.stress$Relative.Abundance)$p.value


# VERRUMICROBIALES
wilcox.verru <- filter(summarized.abundance.order.sample.table, Order == "Verrucomicrobiales")
wilcox.verru.baseine <- filter(wilcox.verru, Group == "Baseline")
wilcox.verru.stress <- filter(wilcox.verru, Group == "3 Wk Stress")
wilcox.test.verru <- wilcox.test(wilcox.verru.baseine$Relative.Abundance, 
                                 wilcox.verru.stress$Relative.Abundance)$p.value


Wilcox.Orders <- c("Acholeplasmatales", "Bacillales", "Bacteroidales", "Burkholderiales", "Chlamydiales", 
                   "Christensenellales", "Clostridiales", "Coriobacteriales", "Corynebacteriales", 
                   "Desulfovibrionales", "Erysipelotrichales", "Exiguobacterales", "Lachnospirales", 
                   "Lactobacillales", "Micrococcales", "Monoglobales", "Oscillospirales", "Peptococcales", 
                   "Peptostreptococcales-Tissierellales", "Staphylococcales", "Streptosporangiales", 
                   "Thermomicrobiales", "Verrucomicrobiales")

Wilcox.pvalues <- c(wilcox.test.ach, wilcox.test.bac, wilcox.test.bact, wilcox.test.burk, wilcox.test.chla, 
                    wilcox.test.chris, wilcox.test.clos, wilcox.test.cori, wilcox.test.cory, wilcox.test.desulf, 
                    wilcox.test.erys, wilcox.test.exig, wilcox.test.lachno, wilcox.test.lacto, wilcox.test.micro,
                    wilcox.test.mono, wilcox.test.osci, wilcox.test.pepto, wilcox.test.PT, wilcox.test.staph, 
                    wilcox.test.strep, wilcox.test.thermo, wilcox.test.verru)

Wilcox.pvalues.Bonferonni <- p.adjust(Wilcox.pvalues, method = "bonferroni", n = 23)

Wilcox.all <- cbind(Wilcox.Orders, Wilcox.pvalues)
Wilcox.all <- cbind(Wilcox.all, Wilcox.pvalues.Bonferonni)
Wilcox.all <- as.data.frame(Wilcox.all)
Wilcox.all <- plyr::rename(Wilcox.all, replace = c("Wilcox.Orders" = "Order", 
                                                   "Wilcox.pvalues" = "Raw p value", 
                                                   "Wilcox.pvalues.Bonferonni" = "Bonferroni-adjusted p value"))

#write.csv(Wilcox.all, file = "../results/tables/wilcoxon rank comparison-order composition-all orders.csv", row.names = FALSE)

