#### FAMILY-LEVEL COMMUNITY COMPOSITION ANALYSIS ###################################################
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

#install.packages("caret")
library(caret)
packageVersion("caret") # I'm using version 6.0.90

#install.packages("randomForest")
library(randomForest)
packageVersion("randomForest") # I'm using version 4.6.14

#install.packages("pROC")
library(pROC)
packageVersion("pROC") # I'm using version 1.18.0


#### INTRODUCTION #################################################################################
# This code will outline the analysis of community composition for the 16S samples for the Gautier 
# manuscript. I'll focus on family-level composition differences in this analysis, along with a 
# random forest model with individual ASVs.

# Load the data
load(file = "../results/phyloseq object for analysis.RData")




#### FAMILY LEVEL DISTRIBUTION #################################################################
# PURPOSE: Look at Family level distribution on a per-sample level.
ps.analysis.prop <-transform_sample_counts(ps.analysis, function(ASV) ASV/sum(ASV))

ps.analysis.prop.family <- tax_glom(ps.analysis.prop, "Family", NArm = FALSE)
family.table <- psmelt(ps.analysis.prop.family) # Organize in long format for ggplot.

# Summarize Familiess for each Group
summarized.abundance.family <- family.table %>%
  group_by(Family, Group) %>%
  summarize(Relative.Abundance = (sum(Abundance) / 24))


# Replace NAs with "Unidentified"
summarized.abundance.family$Family <- as.character(replace_na(summarized.abundance.family$Family, replace = "Unidentified"))


#Collect taxa names for Families > 1% in any sample
summarized.abundance.family.above1 <- filter(summarized.abundance.family, Relative.Abundance >= 0.01)
summarized.abundance.family.above1$Family <- factor(summarized.abundance.family.above1$Family, 
                                                    levels = sort(unique(summarized.abundance.family.above1$Family)))
abundant.families <- as.character(unique(summarized.abundance.family.above1$Family)) 
table(abundant.families) # There are 10 Families that have a Total Relative Abundance above 1%


# Make new data frames either including or excluding abundant families with sample-level data.
summarized.abundance.family.sample <- family.table %>%
  group_by(Family, Sample) %>%
  summarize(Relative.Abundance = (sum(Abundance)))

summarized.abundance.family.sample$Family <- as.character(replace_na(summarized.abundance.family.sample$Family, 
                                                                     replace = "Unidentified"))

summarized.abundance.family.sample.abundant <- filter(summarized.abundance.family.sample, 
                                                      Family %in% abundant.families)

summarized.abundance.family.sample.not.abundant <- filter(summarized.abundance.family.sample, 
                                                          ! Family %in% abundant.families)

summarized.abundance.family.sample.not.abundant <- summarized.abundance.family.sample.not.abundant %>%
  group_by(Sample) %>%
  summarize(Relative.Abundance = (sum(Relative.Abundance))) %>%
  mutate(Family = "Family < 1%", .before = 1)


# Bring data back together and add Group metadata
summarized.abundance.family.sample.table <- rbind(summarized.abundance.family.sample.abundant, 
                                                  summarized.abundance.family.sample.not.abundant)


Sample.Group <- select(family.table, Sample, Group, Animal.Number)
Sample.Group <- Sample.Group %>%
  group_by(Sample) %>%
  slice(1L)

summarized.abundance.family.sample.table <- full_join(summarized.abundance.family.sample.table, Sample.Group, by = "Sample")
summarized.abundance.family.sample.table <- ungroup(summarized.abundance.family.sample.table)

# Organize data into a wide format

# EUBACTERIA
per.sample.abundance.eub <- filter(summarized.abundance.family.sample.table, Family == "[Eubacterium] coprostanoligenes group")
per.sample.abundance.eub$"Eubacterium Relative Abundance (%)" <- per.sample.abundance.eub$Relative.Abundance *100
per.sample.abundance.eub <- select(per.sample.abundance.eub, Sample, Animal.Number, "Eubacterium Relative Abundance (%)")

# ACHOLEPLASMATACEAE
per.sample.abundance.ach <- filter(summarized.abundance.family.sample.table, Family == "Acholeplasmataceae")
per.sample.abundance.ach$"Acholeplasmataceae Relative Abundance (%)" <- per.sample.abundance.ach$Relative.Abundance *100
per.sample.abundance.ach <- select(per.sample.abundance.ach, Sample, "Acholeplasmataceae Relative Abundance (%)")

# CLOSTRIDIACEAE
per.sample.abundance.clos <- filter(summarized.abundance.family.sample.table, Family == "Clostridiaceae")
per.sample.abundance.clos$"Clostridiaceae Relative Abundance (%)" <- per.sample.abundance.clos$Relative.Abundance *100
per.sample.abundance.clos <- select(per.sample.abundance.clos, Sample, "Clostridiaceae Relative Abundance (%)")

# ERYSIPELOTRICHACEAE
per.sample.abundance.erys <- filter(summarized.abundance.family.sample.table, Family == "Erysipelotrichaceae")
per.sample.abundance.erys$"Erysipelotrichaceae Relative Abundance (%)" <- per.sample.abundance.erys$Relative.Abundance *100
per.sample.abundance.erys <- select(per.sample.abundance.erys, Sample, "Erysipelotrichaceae Relative Abundance (%)")

# LACHNOSPIRACEAE
per.sample.abundance.lachno <- filter(summarized.abundance.family.sample.table, Family == "Lachnospiraceae")
per.sample.abundance.lachno$"Lachnospiraceae Relative Abundance (%)" <- per.sample.abundance.lachno$Relative.Abundance *100
per.sample.abundance.lachno <- select(per.sample.abundance.lachno, Sample, "Lachnospiraceae Relative Abundance (%)")

# LACTOBACILLACEAE
per.sample.abundance.lacto <- filter(summarized.abundance.family.sample.table, Family == "Lactobacillaceae")
per.sample.abundance.lacto$"Lactobacillaceae Relative Abundance (%)" <- per.sample.abundance.lacto$Relative.Abundance *100
per.sample.abundance.lacto <- select(per.sample.abundance.lacto, Sample, "Lactobacillaceae Relative Abundance (%)")

# MURIBACULACEAE
per.sample.abundance.muri <- filter(summarized.abundance.family.sample.table, Family == "Muribaculaceae")
per.sample.abundance.muri$"Muribaculaceae Relative Abundance (%)" <- per.sample.abundance.muri$Relative.Abundance *100
per.sample.abundance.muri <- select(per.sample.abundance.muri, Sample, "Muribaculaceae Relative Abundance (%)")

# OSCILLOSPIRACEAE
per.sample.abundance.osci <- filter(summarized.abundance.family.sample.table, Family == "Oscillospiraceae")
per.sample.abundance.osci$"Oscillospiraceae Relative Abundance (%)" <- per.sample.abundance.osci$Relative.Abundance *100
per.sample.abundance.osci <- select(per.sample.abundance.osci, Sample, "Oscillospiraceae Relative Abundance (%)")

# PEPTOSTREPTOCOCCACEAE
per.sample.abundance.pepto <- filter(summarized.abundance.family.sample.table, Family == "Peptostreptococcaceae")
per.sample.abundance.pepto$"Peptostreptococcaceae Relative Abundance (%)" <- per.sample.abundance.pepto$Relative.Abundance *100
per.sample.abundance.pepto <- select(per.sample.abundance.pepto, Sample, "Peptostreptococcaceae Relative Abundance (%)")

# RUMINOCOCCACEAE
per.sample.abundance.rumin <- filter(summarized.abundance.family.sample.table, Family == "Ruminococcaceae")
per.sample.abundance.rumin$"Ruminococcaceae Relative Abundance (%)" <- per.sample.abundance.rumin$Relative.Abundance *100
per.sample.abundance.rumin <- select(per.sample.abundance.rumin, Sample, "Ruminococcaceae Relative Abundance (%)")

#FAMILY <1%
per.sample.abundance.1p <- filter(summarized.abundance.family.sample.table, Family == "Family < 1%")
per.sample.abundance.1p$"Family < 1% Relative Abundance (%)" <- per.sample.abundance.1p$Relative.Abundance *100
per.sample.abundance.1p <- select(per.sample.abundance.1p, Sample, "Family < 1% Relative Abundance (%)")

per.sample.abundance <- full_join(per.sample.abundance.eub, per.sample.abundance.ach, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.clos, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.erys, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.lachno, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.lacto, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.muri, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.osci, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.pepto, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.rumin, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.1p, by = "Sample")

#write.csv(per.sample.abundance, file = "../results/tables/composition per sample-family.csv", row.names = FALSE)

summarized.abundance.family.sample.table.figure <- summarized.abundance.family.sample.table
summarized.abundance.family.sample.table.figure$Relative.Abundance.Percentage <- summarized.abundance.family.sample.table.figure$Relative.Abundance * 100


# Plot community composition per sample
ggplot(summarized.abundance.family.sample.table.figure, aes(x = Animal.Number, y = Relative.Abundance.Percentage, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_x_discrete(limits = rev) +
  theme(axis.text = element_text(size = 5), 
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5)) +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 90))

#ggsave("../results/figures/per sample community composition-family.png", width = 6, height = 4)

# Plot with no Animal number labels.
ggplot(summarized.abundance.family.sample.table.figure, aes(x = Animal.Number, y = Relative.Abundance.Percentage, fill = Family)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_x_discrete(limits = rev) +
  theme(axis.text = element_text(size = 5), 
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5)) +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())

#ggsave("../results/figures/per sample community composition-family-no names.png", width = 6, height = 4)


# I'll now plot abundance for each sample. Abundances are expressed as a percentage of each sample, so I'll format the 
# x axis as a percentage.

summarized.abundance.family.sample.table$Family <- factor(summarized.abundance.family.sample.table$Family, levels = c("[Eubacterium] coprostanoligenes group", "Acholeplasmataceae", "Clostridiaceae", "Erysipelotrichaceae",
                                                                                                                      "Lachnospiraceae", "Lactobacillaceae", "Muribaculaceae", "Oscillospiraceae", "Peptostreptococcaceae", 
                                                                                                                      "Ruminococcaceae", "Family < 1%"))

summarized.abundance.family.sample.table$Group <- factor(summarized.abundance.family.sample.table$Group, 
                                                         levels = c("3 Wk Stress", "Baseline"), ordered = TRUE)

ggplot(summarized.abundance.family.sample.table, aes(x = Relative.Abundance, y = Family, color = Group)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(limits = rev) +
  labs( x = "Relative Abundance", y = NULL, color = "Group") + 
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress"))

#ggsave("../results/figures/family abundance-boxplot.png", width = 6, height = 6)

summarized.abundance.family.sample.table$Group <- factor(summarized.abundance.family.sample.table$Group, 
                                                         levels = c("Baseline", "3 Wk Stress"))


### SIGNIFICANT DIFFERENCES: HIGHER VS LOWER GA AT FIRST TIMEPOINT ###
# I'll now use a Wilcoxon rank test to look for significant differences between each of these Family groups. First, I'll 
# filter the per-sample table for only abundant Families. I'll then perform a Wilcoxon rank test and Bonferroni correction
# for multiple comparisons on each Family to identify significant differences.

wilcox.family <- filter(summarized.abundance.family.sample.table, Family != "Family < 1%" )

# Now I'll calculate Wilcoxon Ranked Sum p values for each Family.

# EUBACTERIUM
wilcox.eub <- filter(wilcox.family, Family == "[Eubacterium] coprostanoligenes group")
wilcox.eub.baseline <- filter(wilcox.eub, Group == "Baseline")
wilcox.eub.stress <- filter(wilcox.eub, Group == "3 Wk Stress")
wilcox.test.eub <- wilcox.test(wilcox.eub.baseline$Relative.Abundance, 
                               wilcox.eub.stress$Relative.Abundance)$p.value

# ACHOLEPLASMATACEAE
wilcox.acho <- filter(wilcox.family, Family == "Acholeplasmataceae")
wilcox.acho.baseline <- filter(wilcox.acho, Group == "Baseline")
wilcox.acho.stress <- filter(wilcox.acho, Group == "3 Wk Stress")
wilcox.test.acho <- wilcox.test(wilcox.acho.baseline$Relative.Abundance, 
                                wilcox.acho.stress$Relative.Abundance)$p.value

# CLOSTRIDIACEAE
wilcox.clos <- filter(wilcox.family, Family == "Clostridiaceae")
wilcox.clos.baseline <- filter(wilcox.clos, Group == "Baseline")
wilcox.clos.stress <- filter(wilcox.clos, Group == "3 Wk Stress")
wilcox.test.clos <- wilcox.test(wilcox.clos.baseline$Relative.Abundance, 
                                wilcox.clos.stress$Relative.Abundance)$p.value

# ERYSIPELOTRICHACEAE
wilcox.erys <- filter(wilcox.family, Family == "Erysipelotrichaceae")
wilcox.erys.baseline <- filter(wilcox.erys, Group == "Baseline")
wilcox.erys.stress <- filter(wilcox.erys, Group == "3 Wk Stress")
wilcox.test.erys <- wilcox.test(wilcox.erys.baseline$Relative.Abundance, 
                                wilcox.erys.stress$Relative.Abundance)$p.value

# LACHNOSPIRACEAE
wilcox.lachno <- filter(wilcox.family, Family == "Lachnospiraceae")
wilcox.lachno.baseline <- filter(wilcox.lachno, Group == "Baseline")
wilcox.lachno.stress <- filter(wilcox.lachno, Group == "3 Wk Stress")
wilcox.test.lachno <- wilcox.test(wilcox.lachno.baseline$Relative.Abundance, 
                                  wilcox.lachno.stress$Relative.Abundance)$p.value

# LACTOBACILLACEAE
wilcox.lacto <- filter(wilcox.family, Family == "Lactobacillaceae")
wilcox.lacto.baseline <- filter(wilcox.lacto, Group == "Baseline")
wilcox.lacto.stress <- filter(wilcox.lacto, Group == "3 Wk Stress")
wilcox.test.lacto <- wilcox.test(wilcox.lacto.baseline$Relative.Abundance, 
                                 wilcox.lacto.stress$Relative.Abundance)$p.value

# MURIBACULACEAE
wilcox.muri <- filter(wilcox.family, Family == "Muribaculaceae")
wilcox.muri.baseline <- filter(wilcox.muri, Group == "Baseline")
wilcox.muri.stress <- filter(wilcox.muri, Group == "3 Wk Stress")
wilcox.test.muri <- wilcox.test(wilcox.muri.baseline$Relative.Abundance, 
                                wilcox.muri.stress$Relative.Abundance)$p.value

# OSCILLOSPIRACEAE
wilcox.osci <- filter(wilcox.family, Family == "Oscillospiraceae")
wilcox.osci.baseline <- filter(wilcox.osci, Group == "Baseline")
wilcox.osci.stress <- filter(wilcox.osci, Group == "3 Wk Stress")
wilcox.test.osci <- wilcox.test(wilcox.osci.baseline$Relative.Abundance, 
                                wilcox.osci.stress$Relative.Abundance)$p.value

# PEPTOSTREPTOCOCCACEAE
wilcox.pepto <- filter(wilcox.family, Family == "Peptostreptococcaceae")
wilcox.pepto.baseline <- filter(wilcox.pepto, Group == "Baseline")
wilcox.pepto.stress <- filter(wilcox.pepto, Group == "3 Wk Stress")
wilcox.test.pepto <- wilcox.test(wilcox.pepto.baseline$Relative.Abundance, 
                                 wilcox.pepto.stress$Relative.Abundance)$p.value

# RUMINOCOCCACEAE
wilcox.rumin <- filter(wilcox.family, Family == "Ruminococcaceae")
wilcox.rumin.baseline <- filter(wilcox.rumin, Group == "Baseline")
wilcox.rumin.stress <- filter(wilcox.rumin, Group == "3 Wk Stress")
wilcox.test.rumin <- wilcox.test(wilcox.rumin.baseline$Relative.Abundance, 
                                 wilcox.rumin.stress$Relative.Abundance)$p.value


Wilcox.Families <- c("[Eubacterium] coprostanoligenes group", "Acholeplasmataceae", "Clostridiaceae", "Erysipelotrichaceae",
                     "Lachnospiraceae", "Lactobacillaceae", "Muribaculaceae", "Oscillospiraceae", "Peptostreptococcaceae", 
                     "Ruminococcaceae")



Wilcox.pvalues <- c(wilcox.test.eub, wilcox.test.acho, wilcox.test.clos, wilcox.test.erys, wilcox.test.lachno, wilcox.test.lacto,
                    wilcox.test.muri, wilcox.test.osci, wilcox.test.pepto, wilcox.test.rumin)

Wilcox.pvalues.Bonferonni <- p.adjust(Wilcox.pvalues, method = "bonferroni", n = 10)

Wilcox.all <- cbind(Wilcox.Families, Wilcox.pvalues)
Wilcox.all <- cbind(Wilcox.all, Wilcox.pvalues.Bonferonni)
Wilcox.all <- as.data.frame(Wilcox.all)
Wilcox.all <- plyr::rename(Wilcox.all, replace = c("Wilcox.Families" = "Family", 
                                                   "Wilcox.pvalues" = "Raw p value", 
                                                   "Wilcox.pvalues.Bonferonni" = "Bonferroni-adjusted p value"))

#write.csv(Wilcox.all, file = "../results/tables/wilcoxon rank comparison-family composition.csv", row.names = FALSE)




#### RANDOM FOREST PLOT OF ASVs ##################################################################
# PURPOSE: Make figure of Random Forest Importance

ASV.table.stress <- psmelt(ps.analysis)

ASV.table.summary.stress <- ASV.table.stress %>%
  group_by(OTU, Group) %>%
  summarize(Abundance.per.Group = (sum(Abundance))) %>%
  filter(any(Abundance.per.Group > 0)) # Removes any ASVs with an abundance of 0 in both groups.

length(unique(ASV.table.summary.stress$OTU)) # There are a total of 430 unique ASVs in these samples.


# Now I'll collect the ASV names for these ASVs. These will be used to merge only the ASV Abundance
# values for ASVs present in at least one group.
unique.ASVs.stress <- as.data.frame(unique(ASV.table.summary.stress$OTU))
unique.ASVs.stress <- plyr::rename(unique.ASVs.stress, replace = c("unique(ASV.table.summary.stress$OTU)" = "OTU"))

ASV.table.samples.stress <- left_join(unique.ASVs.stress, ASV.table.stress, by = "OTU")

length(unique(ASV.table.samples.stress$OTU)) # 430 unique ASVs as expected.


# I'll reshape the data into a wide format for analysis.
ASV.table.wide.stress <- select(ASV.table.samples.stress, OTU, Sample, Abundance, Group, Animal.Number)
ASV.table.wide.stress <- spread(ASV.table.wide.stress, key = OTU, value = Abundance)

# Now I'll split the data into Baseline and 3 Week Stress groups
ASV.table.wide.Baseline <- filter(ASV.table.wide.stress, Group == "Baseline")
ASV.table.wide.Baseline$Group <- factor(ASV.table.wide.Baseline$Group, levels = c("Baseline"))

ASV.table.wide.3Wk.Stress <- filter(ASV.table.wide.stress, Group == "3 Wk Stress")
ASV.table.wide.3Wk.Stress$Group <- factor(ASV.table.wide.3Wk.Stress$Group, levels = c("3 Wk Stress"))


### TRAINING AND TEST SET GENERATION: STRESS TREATMENT ###############################################
# There are 24 mice per group, which is a pretty good sample size for mouse experiments. Because of 
# this, I'll separate the data into training and test sets to better optimize the random forest 
# model and check model performance in an outside population.


# Set seed for reproducibility
#sample(1:1000, 1) # It selected 572
set.seed(572) 

# Now I'll separate each group into training and test sets.

### BASELINE ###
testIndex.Baseline <- createDataPartition(ASV.table.wide.Baseline$Group,
                                          p = 0.30,
                                          list = FALSE,
                                          times = 1)

featureTest.Baseline <- ASV.table.wide.Baseline[testIndex.Baseline,]
featureTrain.Baseline  <- ASV.table.wide.Baseline[-testIndex.Baseline,]


### 3 WEEK STRESS ###
testIndex.3Wk.Stress <- createDataPartition(ASV.table.wide.3Wk.Stress$Group,
                                            p = 0.30,
                                            list = FALSE,
                                            times = 1)

featureTest.3Wk.Stress <- ASV.table.wide.3Wk.Stress[testIndex.3Wk.Stress,]
featureTrain.3Wk.Stress  <- ASV.table.wide.3Wk.Stress[-testIndex.3Wk.Stress,]



# Finally, I'll combine each training and test partition together to make a complete training and 
# test set.

training.set.stress <- rbind(featureTrain.Baseline, featureTrain.3Wk.Stress) # Combine training sets together
test.set.stress <- rbind(featureTest.Baseline, featureTest.3Wk.Stress) # Combine test sets together


# Predictors
predictors.stress.Train <- select(training.set.stress, -Sample, -Group, -Animal.Number) # Remove metadata
predictors.stress.Test <- select(test.set.stress, -Sample, -Group, -Animal.Number) # Remove metadata

# Outcomes
outcome.stress.Train <- as.factor(training.set.stress$Group)
outcome.stress.Test <- as.factor(test.set.stress$Group)




#### CONSTRUCTING THE RANDOM FOREST MODEL ###########################################################
# The parameter I'm going to tune for the random forest model is mtry, the number of features sampled 
# at each node of the decision tree. For this random forest, I'll use ntree (the number of trees in 
# the forest) = 1000, which is a pretty large number but not too computationally taxing. Then I'll 
# look at a variety of different mtry values to see which minimizes the Out-of-Bag error rate.

# By default mtry is set to the square root of the number of features. There are 516 features in 
# this data set, meaning mtry be default is set to 23. I'll test mtry values around 23 to see if 
# this improves model error.


### TUNING MTRY ###
model.rf.mtry19 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 19)
model.rf.mtry19 # The OOB error rate is 3.12%

model.rf.mtry21 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 21)
model.rf.mtry21 # The OOB error rate is 0%

model.rf.mtry23 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 23)
model.rf.mtry23 # The OOB error rate is 0%

model.rf.mtry25 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 25)
model.rf.mtry25 # The OOB error rate is 3.12%

model.rf.mtry27 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 27)
model.rf.mtry27 # The OOB error rate is 3.12%

# OOB error rates were low for all models. Models with mtry = 21 and 23 perfectly classified samples.
# I'll select mtry = 23 because it is the standard recommended mtry value.

model.rf.stress.final <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 23)
model.rf.stress.final




#### TESTING THE RANDOM FOREST MODEL ################################################################
# I now want to test the model on the test set data to see how well it can predict a data set it 
# hasn't seen before. 

pred <- predict(model.rf.stress.final, newdata = predictors.stress.Test)

prediction.stress.test <- select(test.set.stress, Sample, Group)
prediction.stress.test <- cbind(prediction.stress.test, pred)
prediction.stress.test <- plyr::rename(prediction.stress.test, replace = c("Group" = "Actual", 
                                                                           "pred" = "Predicted"))

View(prediction.stress.test)

# Overall, 8/8 (100%) of Baseline samples were correctly classified, while 7/8 (87.5%)  of 3 week
# stress samples were correctly classified. This suggests that the model is not really overfit
# and performs well on data it has not been trained on.

# I'll now generate ROC curves for both training and test sets to summarize model performance.

### TRAINING SET ###
roc(outcome.stress.Train, model.rf.stress.final$votes[,1], plot = TRUE, legacy.axes = TRUE, percent = TRUE, main = "Training Set ROC Curve",
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# The RF model performs perfectly on the training set, with a 100% AUC. This is expected from
# the performance of the model on the training set data.


### TEST SET ###
pred.prob.stress <- predict(model.rf.stress.final, newdata = predictors.stress.Test, type = "prob") # Re-run prediction 
# looking at probabilities instead of classification.

pred.stress.data.frame <- as.data.frame(pred.prob.stress)
roc(outcome.stress.Test, pred.stress.data.frame$`3 Wk Stress`, plot = TRUE, legacy.axes = TRUE, percent = TRUE, main = "Test Set ROC Curve",
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# The RF model for the test set performs perfectly, as well, with a 100% AUC value.




#### IMPORTANCE VALUES FROM RANDOM FOREST MODEL ####################################################
# Now that I've generated the model, I'll look at which features best discriminate between groups as
# selected by the model. I'll do this by checking the importance values from the Random Forest model, 
# using Gini index as the measure of importance.

importance.rf.stress <- importance(model.rf.stress.final, type = 2)

# Format the importance data as a data frame
importance.rf.stress <- cbind(OTU = rownames(importance.rf.stress), importance.rf.stress)
rownames(importance.rf.stress) <- NULL
importance.rf.stress <- as.data.frame(importance.rf.stress)
importance.rf.stress$MeanDecreaseGini <- as.numeric(importance.rf.stress$MeanDecreaseGini)
str(importance.rf.stress$MeanDecreaseGini)

# I now want to combine the RF ranking data with the taxonomic information for identified ASVs. I'll
# pull this information from the phyloseq object metadata.

ASV.taxonomy <- as.data.frame(ps.analysis@tax_table)
ASV.taxonomy <- cbind(OTU = rownames(ASV.taxonomy), ASV.taxonomy)
rownames(ASV.taxonomy) <- NULL

importance.rf.stress.taxonomy <- left_join(importance.rf.stress, ASV.taxonomy, by = "OTU")

# This is the table of ASVs listed by mean decrease in Gini Index. I'll rank the features according
# to Gini index to make the data easier to interpret.

RF.ranking.stress <- importance.rf.stress.taxonomy %>%
  arrange(desc(MeanDecreaseGini))

RF.ranking.stress <- plyr::rename(RF.ranking.stress, replace = c("MeanDecreaseGini" = "Stress.Mean.Decrease.Gini"))
RF.ranking.stress$Stress.RF.Rank <- rank(-RF.ranking.stress$Stress.Mean.Decrease.Gini)

# I also want to add whether ASVs are enriched at baseline or after 3 weeks of stress treatment.
ASV.table.enrichment.stress <- ASV.table.stress %>%
  group_by(OTU, Group) %>%
  summarize(Average.Abundance.per.Group = (mean(Abundance))) %>%
  filter(any(Average.Abundance.per.Group > 0)) # Removes any ASVs with an abundance of 0 in both groups.

ASV.enrichment.baseline <- filter(ASV.table.enrichment.stress, Group == "Baseline")
ASV.enrichment.baseline <- plyr::rename(ASV.enrichment.baseline, replace = c("Average.Abundance.per.Group" = "Baseline.Average.Abundance.per.Group"))
ASV.enrichment.baseline <- select(ASV.enrichment.baseline, OTU, Baseline.Average.Abundance.per.Group)

ASV.enrichment.3Wk.stress <- filter(ASV.table.enrichment.stress, Group == "3 Wk Stress")
ASV.enrichment.3Wk.stress <- plyr::rename(ASV.enrichment.3Wk.stress, replace = c("Average.Abundance.per.Group" = "3Wk.stress.Average.Abundance.per.Group"))
ASV.enrichment.3Wk.stress <- select(ASV.enrichment.3Wk.stress, OTU, '3Wk.stress.Average.Abundance.per.Group')

ASV.enrichment.stress <- full_join(ASV.enrichment.baseline, ASV.enrichment.3Wk.stress, by = "OTU")
ASV.enrichment.stress$Enrichment.Stress <- NA
ASV.enrichment.stress$Enrichment.Stress[ASV.enrichment.stress$Baseline.Average.Abundance.per.Group > ASV.enrichment.stress$`3Wk.stress.Average.Abundance.per.Group`] <- "Baseline"
ASV.enrichment.stress$Enrichment.Stress[ASV.enrichment.stress$Baseline.Average.Abundance.per.Group < ASV.enrichment.stress$`3Wk.stress.Average.Abundance.per.Group`] <- "3 Wk Stress"

ASV.enrichment <- select(ASV.enrichment.stress, OTU, Enrichment.Stress)

RF.ranking.stress <- full_join(RF.ranking.stress, ASV.enrichment, by = "OTU")

RF.ranking.stress$OTU <- factor(RF.ranking.stress$OTU, levels = RF.ranking.stress$OTU)

# I'll make new names for each figure, which will concatenate the ASV number and the Family
# designation for each ASV.
RF.ranking.stress$ASV.Family <- NA
RF.ranking.stress$ASV.Family <- paste(RF.ranking.stress$OTU, " (", RF.ranking.stress$Family, ")")
RF.ranking.stress$ASV.Family <- factor(RF.ranking.stress$ASV.Family, levels = RF.ranking.stress$ASV.Family)


ggplot(RF.ranking.stress, aes(x = Stress.Mean.Decrease.Gini, y = ASV.Family, color = Enrichment.Stress)) +
  geom_point(size = 3) +
  labs(title = "Top Ranked Random Forest Predictors", x = "Mean Decrease in Node Impurity (Gini)", y = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# To reduce the number of ASVs graphed on the figure, I'm going to set an arbitrarily cutoff 
# for the Gini Index to include only the most important ASVs.

RF.ranking.stress.figure <- filter(RF.ranking.stress, Stress.Mean.Decrease.Gini > 0.201)

ggplot(RF.ranking.stress.figure, aes(x = Stress.Mean.Decrease.Gini, y = ASV.Family, color = Enrichment.Stress)) +
  geom_point(size = 1) +
  labs(x = "Mean Decrease in Node Impurity (Gini)", y = NULL, color = "Enrichment") +
  scale_color_manual(name = "Enrichment", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  theme(axis.text = element_text(size = 5), 
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(face = "bold", size = 5),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))

#ggsave("../results/figures/RF Rankings-Baseline vs Stress-ASV.png", width = 3, height = 4)




#### RANDOM FOREST PLOT OF FAMILIES ##############################################################
# PURPOSE: Make figure of Random Forest Importance, this time at the Family level rather than the
# ASV level.

ps.analysis.prop <-transform_sample_counts(ps.analysis, function(ASV) ASV/sum(ASV))

ps.analysis.prop.family <- tax_glom(ps.analysis.prop, "Family", NArm = FALSE)
family.table <- psmelt(ps.analysis.prop.family) # Organize in long format for ggplot.

# I'll reshape the data into a wide format for analysis.
family.table.wide <- select(family.table, Family, Sample, Abundance, Group, Animal.Number)
family.table.wide <- spread(family.table.wide, key = Family, value = Abundance)


# Now I'll split the data into Baseline and 3 Week Stress groups
family.table.wide.Baseline <- filter(family.table.wide, Group == "Baseline")
family.table.wide.Baseline$Group <- factor(ASV.table.wide.Baseline$Group, levels = c("Baseline"))

family.table.wide.3Wk.Stress <- filter(family.table.wide, Group == "3 Wk Stress")
family.table.wide.3Wk.Stress$Group <- factor(family.table.wide.3Wk.Stress$Group, levels = c("3 Wk Stress"))




### TRAINING AND TEST SET GENERATION: STRESS TREATMENT ###############################################
# There are 24 mice per group, which is a pretty good sample size for mouse experiments. Because of 
# this, I'll separate the data into training and test sets to better optimize the random forest 
# model and check model performance in an outside population.


# Set seed for reproducibility
#sample(1:1000, 1) # It selected 763
set.seed(763) 

# Now I'll separate each group into training and test sets.

### BASELINE ###
testIndex.Baseline <- createDataPartition(family.table.wide.Baseline$Group,
                                          p = 0.30,
                                          list = FALSE,
                                          times = 1)

featureTest.Baseline <- family.table.wide.Baseline[testIndex.Baseline,]
featureTrain.Baseline  <- family.table.wide.Baseline[-testIndex.Baseline,]


### 3 WEEK STRESS ###
testIndex.3Wk.Stress <- createDataPartition(family.table.wide.3Wk.Stress$Group,
                                            p = 0.30,
                                            list = FALSE,
                                            times = 1)

featureTest.3Wk.Stress <- family.table.wide.3Wk.Stress[testIndex.3Wk.Stress,]
featureTrain.3Wk.Stress  <- family.table.wide.3Wk.Stress[-testIndex.3Wk.Stress,]



# Finally, I'll combine each training and test partition together to make a complete training and 
# test set.

training.set.stress <- rbind(featureTrain.Baseline, featureTrain.3Wk.Stress) # Combine training sets together
test.set.stress <- rbind(featureTest.Baseline, featureTest.3Wk.Stress) # Combine test sets together


# Predictors
predictors.stress.Train <- select(training.set.stress, -Sample, -Group, -Animal.Number) # Remove metadata
predictors.stress.Test <- select(test.set.stress, -Sample, -Group, -Animal.Number) # Remove metadata

# Outcomes
outcome.stress.Train <- as.factor(training.set.stress$Group)
outcome.stress.Test <- as.factor(test.set.stress$Group)




#### CONSTRUCTING THE RANDOM FOREST MODEL ###########################################################
# The parameter I'm going to tune for the random forest model is mtry, the number of features sampled 
# at each node of the decision tree. For this random forest, I'll use ntree (the number of trees in 
# the forest) = 1000, which is a pretty large number but not too computationally taxing. Then I'll 
# look at a variety of different mtry values to see which minimizes the Out-of-Bag error rate.

# By default mtry is set to the square root of the number of features. There are 33 features in 
# this data set, meaning mtry be default is set to 5.7. I'll test mtry values around 6 to see if 
# this improves model error.


### TUNING MTRY ###
model.rf.mtry4 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 4)
model.rf.mtry4 # The OOB error rate is 3.12%

model.rf.mtry5 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 5)
model.rf.mtry5 # The OOB error rate is 3.12%

model.rf.mtry6 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 6)
model.rf.mtry6 # The OOB error rate is 3.12%

model.rf.mtry7 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 7)
model.rf.mtry7 # The OOB error rate is 3.12%

model.rf.mtry8 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 8)
model.rf.mtry8 # The OOB error rate is 3.12

# All random forest models performed very well at correctly classifying Baseline and 3 week stress 
# samples, regardless of mtry value. Because of this, I'll just go with the standard mtry value
# of 6.

model.rf.stress.final <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 6)
model.rf.stress.final




#### TESTING THE RANDOM FOREST MODEL ################################################################
# I now want to test the model on the test set data to see how well it can predict a data set it 
# hasn't seen before. 

pred <- predict(model.rf.stress.final, newdata = predictors.stress.Test)

prediction.stress.test <- select(test.set.stress, Sample, Group)
prediction.stress.test <- cbind(prediction.stress.test, pred)
prediction.stress.test <- plyr::rename(prediction.stress.test, replace = c("Group" = "Actual", 
                                                                           "pred" = "Predicted"))

View(prediction.stress.test)

# Overall, 8/8 (100%) of Baseline samples were correctly classified, while 6/8 (75%)  of 3 week
# stress samples were correctly classified. This suggests that the model is not really overfit
# and performs well on data it has not been trained on.

# I'll now generate ROC curves for both training and test sets to summarize model performance.

### TRAINING SET ###
roc(outcome.stress.Train, model.rf.stress.final$votes[,1], plot = TRUE, legacy.axes = TRUE, percent = TRUE, main = "Training Set ROC Curve",
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# The RF model performs very well on the training set, with a 99.6% AUC. This is expected from
# the performance of the model on the training set data.


### TEST SET ###
pred.prob.stress <- predict(model.rf.stress.final, newdata = predictors.stress.Test, type = "prob") # Re-run prediction 
# looking at probabilities instead of classification.

pred.stress.data.frame <- as.data.frame(pred.prob.stress)
roc(outcome.stress.Test, pred.stress.data.frame$`3 Wk Stress`, plot = TRUE, legacy.axes = TRUE, percent = TRUE, main = "Test Set ROC Curve",
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# The RF model for the test set performs perfectly with a 100% AUC value.




#### IMPORTANCE VALUES FROM RANDOM FOREST MODEL ####################################################
# Now that I've generated the model, I'll look at which features best discriminate between groups as
# selected by the model. I'll do this by checking the importance values from the Random Forest model, 
# using Gini index as the measure of importance.

importance.rf.stress <- importance(model.rf.stress.final, type = 2)

# Format the importance data as a data frame
importance.rf.stress <- cbind(Family = rownames(importance.rf.stress), importance.rf.stress)
rownames(importance.rf.stress) <- NULL
importance.rf.stress <- as.data.frame(importance.rf.stress)
importance.rf.stress$MeanDecreaseGini <- as.numeric(importance.rf.stress$MeanDecreaseGini)
str(importance.rf.stress$MeanDecreaseGini)

# This is the table of familiess listed by mean decrease in Gini Index. I'll rank the features 
# according to Gini index to make the data easier to interpret.

RF.ranking.stress <- importance.rf.stress %>%
  arrange(desc(MeanDecreaseGini))

RF.ranking.stress <- plyr::rename(RF.ranking.stress, replace = c("MeanDecreaseGini" = "Stress.Mean.Decrease.Gini"))
RF.ranking.stress$Stress.RF.Rank <- rank(-RF.ranking.stress$Stress.Mean.Decrease.Gini)

# I also want to add whether families are enriched at baseline or after 3 weeks of stress treatment.
family.table.enrichment.stress <- family.table %>%
  group_by(Family, Group) %>%
  summarize(Average.Abundance.per.Group = (mean(Abundance))) %>%
  filter(any(Average.Abundance.per.Group > 0)) # Removes any familys with an abundance of 0 in both groups.

family.enrichment.baseline <- filter(family.table.enrichment.stress, Group == "Baseline")
family.enrichment.baseline <- plyr::rename(family.enrichment.baseline, replace = c("Average.Abundance.per.Group" = "Baseline.Average.Abundance.per.Group"))
family.enrichment.baseline <- select(family.enrichment.baseline, Family, Baseline.Average.Abundance.per.Group)

family.enrichment.3Wk.stress <- filter(family.table.enrichment.stress, Group == "3 Wk Stress")
family.enrichment.3Wk.stress <- plyr::rename(family.enrichment.3Wk.stress, replace = c("Average.Abundance.per.Group" = "3Wk.stress.Average.Abundance.per.Group"))
family.enrichment.3Wk.stress <- select(family.enrichment.3Wk.stress, Family, '3Wk.stress.Average.Abundance.per.Group')

family.enrichment.stress <- full_join(family.enrichment.baseline, family.enrichment.3Wk.stress, by = "Family")
family.enrichment.stress$Enrichment.Stress <- NA
family.enrichment.stress$Enrichment.Stress[family.enrichment.stress$Baseline.Average.Abundance.per.Group > family.enrichment.stress$`3Wk.stress.Average.Abundance.per.Group`] <- "Baseline"
family.enrichment.stress$Enrichment.Stress[family.enrichment.stress$Baseline.Average.Abundance.per.Group < family.enrichment.stress$`3Wk.stress.Average.Abundance.per.Group`] <- "3 Wk Stress"

family.enrichment <- select(family.enrichment.stress, Family, Enrichment.Stress)

RF.ranking.stress <- full_join(RF.ranking.stress, family.enrichment, by = "Family")

RF.ranking.stress$Family <- factor(RF.ranking.stress$Family, levels = RF.ranking.stress$Family)


ggplot(RF.ranking.stress, aes(x = Stress.Mean.Decrease.Gini, y = Family, color = Enrichment.Stress)) +
  geom_point(size = 1) +
  labs(x = "Mean Decrease in Node Impurity (Gini)", y = NULL, color = "Enrichment") +
  scale_color_manual(name = "Enrichment", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  theme(axis.text = element_text(size = 5), 
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(face = "bold", size = 5),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))

#ggsave("../results/figures/RF Rankings-Baseline vs Stress-family.png", width = 3, height = 4)
