#### SELECTION OF ASVs ASSOCIATED WITH MUCIN EXPERIMENT OUTCOMES ###################################
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


#### INTRODUCTION ###################################################################################
# I've broadly looked at differences in diversity and community composition with the Mucin 
# Supplementation experiment for the Gautier TUMI Pilot. # I'll now move on to identify ASVs that best 
# discriminate between groups. There are two comparisons we are interested in for this experiment:
## 1) Baseline (n=24) vs. 3 Wk Stress samples (n=24)
## 2) Stress only (n=12) vs. Stress + Mucin (n=11)


# I'll start with the Baseline vs 3 Wk Stress comparison, then move on to the Control vs. Mucin 
# comparison. I'll be using a Random Forest model to select ASVs that best discriminate between these 
# conditions.




#### ORGANIZATION OF ASVs FROM BASELINE AND 3 WEEK STRESS SAMPLES ####################################
# I'll be comparing only the Baseline and 3 week stress groups for this analysis, so I'll start by 
# filtering out all the 4 week samples. 

#setwd("..") # I need to set the working directory up one to the parent directory for the code to 
# work properly. If this has been done previously during the same session then this step should be 
# skipped, or the working directory will be moved out of the parent directory.

load(file = "./results/phyloseq objects by experiment.RData")
rm(ps.stress)

ps.mucin.stress.ASV <- subset_samples(ps.mucin, Group == "Baseline" | Group == "3 Wk Stress")
ASV.table.stress <- psmelt(ps.mucin.stress.ASV)

length(unique(ASV.table.stress$OTU)) #1795 unique ASVs across all samples.

# I want to limit the data set to only ASVs that are present in at least one group. To do this, I'll
# first summarize the Abundance of each ASV in each group, then select all ASVs that are >0 in any
# group.

ASV.table.summary.stress <- ASV.table.stress %>%
  group_by(OTU, Group) %>%
  summarize(Abundance.per.Group = (sum(Abundance))) %>%
  filter(any(Abundance.per.Group > 0)) # Removes any ASVs with an abundance of 0 in both groups.

length(unique(ASV.table.summary.stress$OTU)) # There are a total of 516 unique ASVs in these samples.


# Now I'll collect the ASV names for these ASVs. These will be used to merge only the ASV Abundance
# values for ASVs present in at least one group.
unique.ASVs.stress <- as.data.frame(unique(ASV.table.summary.stress$OTU))
unique.ASVs.stress <- plyr::rename(unique.ASVs.stress, replace = c("unique(ASV.table.summary.stress$OTU)" = "OTU"))

ASV.table.samples.stress <- left_join(unique.ASVs.stress, ASV.table.stress, by = "OTU")

length(unique(ASV.table.samples.stress$OTU)) # 516 unique ASVs as expected.


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
model.rf.mtry19 # The OOB error rate is 0%

model.rf.mtry21 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 21)
model.rf.mtry21 # The OOB error rate is 0%

model.rf.mtry23 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 23)
model.rf.mtry23 # The OOB error rate is 0%

model.rf.mtry25 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 25)
model.rf.mtry25 # The OOB error rate is 0%

model.rf.mtry27 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 27)
model.rf.mtry27 # The OOB error rate is 0%

# All random forest models performed perfectly at correctly classifying Baseline and 3 week stress 
# samples, regardless of mtry value. Because of this, I'll just go with the standard mtry value
# of 23.

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

ASV.taxonomy <- as.data.frame(ps.mucin@tax_table)
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

RF.ranking.stress.figure <- filter(RF.ranking.stress, Stress.Mean.Decrease.Gini > 0.12)

ggplot(RF.ranking.stress.figure, aes(x = Stress.Mean.Decrease.Gini, y = ASV.Family, color = Enrichment.Stress)) +
  geom_point(size = 3) +
  labs(title = "Top Ranked Random Forest Predictors", x = "Mean Decrease in Node Impurity (Gini)", y = NULL, color = "Enrichment") +
  scale_color_manual(name = "Enrichment", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("./results/figures/mucin/RF Rankings_Baseline_vs_Stress.png", width = 7, height = 5)


# Print Random Forest importance rankings
write.csv(RF.ranking.stress, file = "./results/tables/mucin/RF Rankings_Baseline_vs_Stress.csv", row.names = FALSE)




#### RANDOM FOREST MODEL WITH TOP 5 ASVs ONLY #######################################################
# The top 5 ASVs selected by the random forest model look clearly separated from other ASVs. Because
# the model descriminates so well between Baseline and 3 Week samples, I want to see whether these 
# ASVs alone are sufficient to separate groups. To test this, I'll filter only these ASVs for each
# sample and use this data to generate a new random forest model.

# Set seed for reproducibility
#sample(1:1000, 1) # It selected 182
set.seed(182) 

predictors.stress.Train.subset <- select(predictors.stress.Train, ASV195, ASV159, ASV50, ASV73, ASV104)
predictors.stress.Test.subset <- select(predictors.stress.Test, ASV195, ASV159, ASV50, ASV73, ASV104)


model.rf.stress.subset <- randomForest(x = predictors.stress.Train.subset, y = outcome.stress.Train, ntree = 1000)

model.rf.stress.subset # The model has 3.12% error rate, misclassifying 1 sample.

pred.subset <- predict(model.rf.stress.subset, newdata = predictors.stress.Test.subset)

prediction.stress.test <- select(test.set.stress, Sample, Group)
prediction.stress.test.subset <- cbind(prediction.stress.test, pred.subset)
prediction.stress.test.subset <- plyr::rename(prediction.stress.test.subset, 
                                              replace = c("Group" = "Actual", "pred.subset" = "Predicted"))

View(prediction.stress.test.subset) # 15/16 Test set samples were correctly predicted.

# These data indicate that the Top 5 ranked ASVs are sufficient to clearly separate Baseline and
# 3 week Stress samples.




#### ORGANIZATION OF ASVs FROM CONTROL AND MUCIN-SUPPLEMENTED SAMPLES ################################
# I'll now move on to the second half of the experiment. Here I'll be classifying 4 week samples only
# as either Control or Mucin-supplemented groups.

ps.mucin.ASV <- subset_samples(ps.mucin, Group == "Stress Alone" | Group == "Stress + Mucin")
ASV.table.mucin <- psmelt(ps.mucin.ASV)

length(unique(ASV.table.stress$OTU)) #1795 unique ASVs across all samples.

# I want to limit the data set to only ASVs that are present in at least one group. To do this, I'll
# first summarize the Abundance of each ASV in each group, then select all ASVs that are >0 in any
# group.

ASV.table.summary.mucin <- ASV.table.mucin %>%
  group_by(OTU, Group) %>%
  summarize(Abundance.per.Group = (sum(Abundance))) %>%
  filter(any(Abundance.per.Group > 0)) # Removes any ASVs with an abundance of 0 in both groups.

length(unique(ASV.table.summary.mucin$OTU)) # There are a total of 498 unique ASVs in these samples.


# Now I'll collect the ASV names for these ASVs. These will be used to merge only the ASV Abundance
# values for ASVs present in at least one group.
unique.ASVs.mucin <- as.data.frame(unique(ASV.table.summary.mucin$OTU))
unique.ASVs.mucin <- plyr::rename(unique.ASVs.mucin, replace = c("unique(ASV.table.summary.mucin$OTU)" = "OTU"))

ASV.table.samples.mucin <- left_join(unique.ASVs.mucin, ASV.table.mucin, by = "OTU")

length(unique(ASV.table.samples.mucin$OTU)) # 498 unique ASVs as expected.


# I'll reshape the data into a wide format for analysis.
ASV.table.wide.mucin <- select(ASV.table.samples.mucin, OTU, Sample, Abundance, Group, Animal.Number)
ASV.table.wide.mucin <- spread(ASV.table.wide.mucin, key = OTU, value = Abundance)

# Now I'll split the data into Control and Mucin groups
ASV.table.wide.Control <- filter(ASV.table.wide.mucin, Group == "Stress Alone")
ASV.table.wide.Control$Group <- factor(ASV.table.wide.Control$Group, levels = c("Stress Alone"))

ASV.table.wide.Mucin <- filter(ASV.table.wide.mucin, Group == "Stress + Mucin")
ASV.table.wide.Mucin$Group <- factor(ASV.table.wide.Mucin$Group, levels = c("Stress + Mucin"))


### TRAINING AND TEST SET GENERATION: STRESS TREATMENT ###############################################
# There are 12 mice in the Control group and 11 mice in the Mucin-supplemented group. This is smaller
# than the sample size for the Baseline vs Stress comparison, but should still be ok for model
# building. I'll now separate the data into training and test sets to better optimize the random forest 
# model and check model performance in an outside population.


# Set seed for reproducibility
#sample(1:1000, 1) # It selected 389
set.seed(389) 

# Now I'll separate each group into training and test sets.

### CONTROL ###
testIndex.Control <- createDataPartition(ASV.table.wide.Control$Group,
                                          p = 0.30,
                                          list = FALSE,
                                          times = 1)

featureTest.Control <- ASV.table.wide.Control[testIndex.Control,]
featureTrain.Control  <- ASV.table.wide.Control[-testIndex.Control,]


### MUCIN ###
testIndex.Mucin <- createDataPartition(ASV.table.wide.Mucin$Group,
                                            p = 0.30,
                                            list = FALSE,
                                            times = 1)

featureTest.Mucin <- ASV.table.wide.Mucin[testIndex.Mucin,]
featureTrain.Mucin  <- ASV.table.wide.Mucin[-testIndex.Mucin,]



# Finally, I'll combine each training and test partition together to make a complete training and 
# test set.

training.set.mucin <- rbind(featureTrain.Control, featureTrain.Mucin) # Combine training sets together
test.set.mucin <- rbind(featureTest.Control, featureTest.Mucin) # Combine test sets together


# Predictors
predictors.mucin.Train <- select(training.set.mucin, -Sample, -Group, -Animal.Number) # Remove metadata
predictors.mucin.Test <- select(test.set.mucin, -Sample, -Group, -Animal.Number) # Remove metadata

# Outcomes
outcome.mucin.Train <- as.factor(training.set.mucin$Group)
outcome.mucin.Test <- as.factor(test.set.mucin$Group)




#### CONSTRUCTING THE RANDOM FOREST MODEL ###########################################################
# The parameter I'm going to tune for the random forest model is mtry, the number of features sampled 
# at each node of the decision tree. For this random forest, I'll use ntree (the number of trees in 
# the forest) = 1000, which is a pretty large number but not too computationally taxing. Then I'll 
# look at a variety of different mtry values to see which minimizes the Out-of-Bag error rate.

# By default mtry is set to the square root of the number of features. There are 498 features in 
# this data set, meaning mtry be default is set to 22. I'll test mtry values around 22 to see if 
# this improves model error.


### TUNING MTRY ###
model.rf.mtry18 <- randomForest(x = predictors.mucin.Train, y = outcome.mucin.Train, ntree = 1000, mtry = 18)
model.rf.mtry18 # The OOB error rate is 13.33%

model.rf.mtry20 <- randomForest(x = predictors.mucin.Train, y = outcome.mucin.Train, ntree = 1000, mtry = 20)
model.rf.mtry20 # The OOB error rate is 13.33%

model.rf.mtry22 <- randomForest(x = predictors.mucin.Train, y = outcome.mucin.Train, ntree = 1000, mtry = 22)
model.rf.mtry22 # The OOB error rate is 20%

model.rf.mtry24 <- randomForest(x = predictors.mucin.Train, y = outcome.mucin.Train, ntree = 1000, mtry = 24)
model.rf.mtry24 # The OOB error rate is 13.33%

model.rf.mtry26 <- randomForest(x = predictors.mucin.Train, y = outcome.mucin.Train, ntree = 1000, mtry = 26)
model.rf.mtry26 # The OOB error rate is 6.67%

# These models performed more poorly than the earlier models, but still performed well, with OOB 
# error rates between 7-20%. Repeating these models several times, it looks like there 
# heterogeneity in the error rate for each model, so they're all roughly similar. I'll go with 
# mtry == 26, because it was slightly more consistent with a lower OOB error rate.


model.rf.mucin.final <- randomForest(x = predictors.mucin.Train, y = outcome.mucin.Train, ntree = 1000, mtry = 26)
model.rf.mucin.final




#### TESTING THE RANDOM FOREST MODEL ################################################################
# I now want to test the model on the test set data to see how well it can predict a data set it 
# hasn't seen before. 

pred <- predict(model.rf.mucin.final, newdata = predictors.mucin.Test)

prediction.mucin.test <- select(test.set.mucin, Sample, Group)
prediction.mucin.test <- cbind(prediction.mucin.test, pred)
prediction.mucin.test <- plyr::rename(prediction.mucin.test, replace = c("Group" = "Actual", 
                                                                           "pred" = "Predicted"))

View(prediction.mucin.test)

# Overall, 4/4 (100%) of Control (Stress Alone) samples were correctly classified, while 4/4 (100%)  
# of Mucin-supplemented (Stress + Mucin) samples were correctly classified. This suggests that the 
# model is not really overfit and performs well on data it has not been trained on.

# I'll now generate ROC curves for both training and test sets to summarize model performance.

### TRAINING SET ###
roc(outcome.mucin.Train, model.rf.mucin.final$votes[,1], plot = TRUE, legacy.axes = TRUE, percent = TRUE, main = "Training Set ROC Curve",
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# The RF model performs very well on the training set, with a 98.2% AUC. This is expected from
# the performance of the model on the training set data.


### TEST SET ###
pred.prob.mucin <- predict(model.rf.mucin.final, newdata = predictors.mucin.Test, type = "prob") # Re-run prediction 
# looking at probabilities instead of classification.

pred.mucin.data.frame <- as.data.frame(pred.prob.mucin)
roc(outcome.mucin.Test, pred.mucin.data.frame$`Stress + Mucin`, plot = TRUE, legacy.axes = TRUE, percent = TRUE, main = "Test Set ROC Curve",
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# The RF model for the test set performs perfectly, with a 100% AUC value.




#### IMPORTANCE VALUES FROM RANDOM FOREST MODEL ####################################################
# Now that I've generated the model, I'll look at which features best discriminate between groups as
# selected by the model. I'll do this by checking the importance values from the Random Forest model, 
# using Gini index as the measure of importance.

importance.rf.mucin <- importance(model.rf.mucin.final, type = 2)

# Format the importance data as a data frame
importance.rf.mucin <- cbind(OTU = rownames(importance.rf.mucin), importance.rf.mucin)
rownames(importance.rf.mucin) <- NULL
importance.rf.mucin <- as.data.frame(importance.rf.mucin)
importance.rf.mucin$MeanDecreaseGini <- as.numeric(importance.rf.mucin$MeanDecreaseGini)
str(importance.rf.mucin$MeanDecreaseGini)

# I now want to combine the RF ranking data with the taxonomic information for identified ASVs. I'll
# pull this information from the phyloseq object metadata.

importance.rf.mucin.taxonomy <- left_join(importance.rf.mucin, ASV.taxonomy, by = "OTU")

# This is the table of ASVs listed by mean decrease in Gini Index. I'll rank the features according
# to Gini index to make the data easier to interpret.

RF.ranking.mucin <- importance.rf.mucin.taxonomy %>%
  arrange(desc(MeanDecreaseGini))

RF.ranking.mucin <- plyr::rename(RF.ranking.mucin, replace = c("MeanDecreaseGini" = "Mucin.Mean.Decrease.Gini"))
RF.ranking.mucin$Mucin.RF.Rank <- rank(-RF.ranking.mucin$Mucin.Mean.Decrease.Gini)


# I also want to add whether ASVs are enriched in the Control group or in the Mucin-supplemented group.
ASV.table.enrichment.mucin <- ASV.table.mucin %>%
  group_by(OTU, Group) %>%
  summarize(Average.Abundance.per.Group = (mean(Abundance))) %>%
  filter(any(Average.Abundance.per.Group > 0)) # Removes any ASVs with an abundance of 0 in both groups.

ASV.enrichment.Control <- filter(ASV.table.enrichment.mucin, Group == "Stress Alone")
ASV.enrichment.Control <- plyr::rename(ASV.enrichment.Control, replace = c("Average.Abundance.per.Group" = "Control.Average.Abundance.per.Group"))
ASV.enrichment.Control <- select(ASV.enrichment.Control, OTU, Control.Average.Abundance.per.Group)

ASV.enrichment.Mucin <- filter(ASV.table.enrichment.mucin, Group == "Stress + Mucin")
ASV.enrichment.Mucin <- plyr::rename(ASV.enrichment.Mucin, replace = c("Average.Abundance.per.Group" = "Mucin.Average.Abundance.per.Group"))
ASV.enrichment.Mucin <- select(ASV.enrichment.Mucin, OTU, 'Mucin.Average.Abundance.per.Group')

ASV.enrichment.mucin <- full_join(ASV.enrichment.Control, ASV.enrichment.Mucin, by = "OTU")
ASV.enrichment.mucin$Enrichment.Mucin <- NA
ASV.enrichment.mucin$Enrichment.Mucin[ASV.enrichment.mucin$Control.Average.Abundance.per.Group > ASV.enrichment.mucin$Mucin.Average.Abundance.per.Group] <- "Stress Alone"
ASV.enrichment.mucin$Enrichment.Mucin[ASV.enrichment.mucin$Control.Average.Abundance.per.Group < ASV.enrichment.mucin$Mucin.Average.Abundance.per.Group] <- "Stress + Mucin"

ASV.enrichment <- select(ASV.enrichment.mucin, OTU, Enrichment.Mucin)

RF.ranking.mucin <- full_join(RF.ranking.mucin, ASV.enrichment, by = "OTU")

RF.ranking.mucin$OTU <- factor(RF.ranking.mucin$OTU, levels = RF.ranking.mucin$OTU)

# I'll make new names for each figure, which will concatenate the ASV number and the Family
# designation for each ASV.
RF.ranking.mucin$ASV.Family <- NA
RF.ranking.mucin$ASV.Family <- paste(RF.ranking.mucin$OTU, " (", RF.ranking.mucin$Family, ")")
RF.ranking.mucin$ASV.Family <- factor(RF.ranking.mucin$ASV.Family, levels = RF.ranking.mucin$ASV.Family)


ggplot(RF.ranking.mucin, aes(x = Mucin.Mean.Decrease.Gini, y = ASV.Family, color = Enrichment.Mucin)) +
  geom_point(size = 3) +
  labs(title = "Top Ranked Random Forest Predictors", x = "Mean Decrease in Node Impurity (Gini)", y = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# To reduce the number of ASVs graphed on the figure, I'm going to set an arbitrarily cutoff 
# for the Gini Index to include only the most important ASVs.

RF.ranking.mucin.figure <- filter(RF.ranking.mucin, Mucin.Mean.Decrease.Gini > 0.07)

RF.ranking.mucin.figure$Enrichment.Mucin <- factor(RF.ranking.mucin.figure$Enrichment.Mucin, levels = c("Stress Alone", "Stress + Mucin"))

ggplot(RF.ranking.mucin.figure, aes(x = Mucin.Mean.Decrease.Gini, y = ASV.Family, color = Enrichment.Mucin)) +
  geom_point(size = 3) +
  labs(title = "Top Ranked Random Forest Predictors", x = "Mean Decrease in Node Impurity (Gini)", y = NULL, color = "Enrichment") +
  scale_color_manual(name = "Enrichment", 
                     values = c("#521b92", "#4d8e00"),
                     breaks = c("Stress Alone", "Stress + Mucin"), 
                     labels = c("Stress Alone", "Stress + Mucin")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("./results/figures/mucin/RF Rankings_Control_vs_Mucin.png", width = 7, height = 5)


# Print Random Forest importance rankings
write.csv(RF.ranking.mucin, file = "./results/tables/mucin/RF Rankings_Control_vs_Mucin.csv", row.names = FALSE)




#### RANDOM FOREST MODEL WITH TOP 5 ASVs ONLY #######################################################
# The top 5 ASVs selected by the random forest model look clearly separated from other ASVs. Because
# the model descriminates so well between Control and Mucin-supplemented samples, I want to see whether 
# these ASVs alone are sufficient to separate groups. To test this, I'll filter only these ASVs for 
# each sample and use this data to generate a new random forest model.

# Set seed for reproducibility
#sample(1:1000, 1) # It selected 972
set.seed(972) 

predictors.mucin.Train.subset <- select(predictors.mucin.Train, ASV85, ASV63, ASV176, ASV186, ASV49)
predictors.mucin.Test.subset <- select(predictors.mucin.Test, ASV85, ASV63, ASV176, ASV186, ASV49)


model.rf.mucin.subset <- randomForest(x = predictors.mucin.Train.subset, y = outcome.mucin.Train, ntree = 1000)

model.rf.mucin.subset # The model has 0% OOB error rate.

pred.subset <- predict(model.rf.mucin.subset, newdata = predictors.mucin.Test.subset)

prediction.mucin.test <- select(test.set.mucin, Sample, Group)
prediction.mucin.test.subset <- cbind(prediction.mucin.test, pred.subset)
prediction.mucin.test.subset <- plyr::rename(prediction.mucin.test.subset, 
                                              replace = c("Group" = "Actual", "pred.subset" = "Predicted"))

View(prediction.mucin.test.subset) # 3/4 (75%) Control samples were correctly classified, while 3/4 (75%)
# Mucin-supplemented samples were correctly classified in the Test set. samples were correctly predicted.
# These data indicate that the Top 5 ranked ASVs still do a pretty good job at separating Control and 
# Mucin supplemented samples.




#### OVERLAP BETWEEN STRESS AND MUCIN PREDICTIVE ASVs #####################################################
# Finally, I want to check if any of the ASVs identified as discriminating between Baseline and Stress
# conditions also discreminate between Stress Alone and Stress + Mucin conditions. To do this, I'll combine
# ranking data sets and compare ASVs.

RF.ranking.combined <- select(RF.ranking.stress, OTU, Stress.Mean.Decrease.Gini, Stress.RF.Rank, Enrichment.Stress)

RF.ranking.combined <- full_join(RF.ranking.combined, RF.ranking.mucin, by = "OTU")
RF.ranking.combined <- select(RF.ranking.combined, OTU, Stress.Mean.Decrease.Gini, Mucin.Mean.Decrease.Gini,
                              Stress.RF.Rank, Mucin.RF.Rank, Enrichment.Stress, Enrichment.Mucin, everything() )                     

View(RF.ranking.combined)
# It doesn't look like there's a whole lot of overlap between data sets. To look for sure, I'll arbitrarily
# limit the data set to ASVs that were in the top 50 for each outcome by mean decrease in Gini Index.

RF.ranking.combined.important <- filter(RF.ranking.combined, Stress.RF.Rank < 50 & Mucin.RF.Rank < 50)                              

View(RF.ranking.combined.important)

# Only 7 ASVs were ranked at least 50 in both outcomes. All 7 were from the Lachnospiraceae Family.
# Enrichment for 4/7 was as hypothesized (I.e., enriched at both Baseline and in Stress + Mucin, or 
# enriched both after 3 Wk Stress and in Stress Alone).

# I'll now print this data set for subsequent analysis.
write.csv(RF.ranking.combined.important, file = "./results/tables/mucin/RF Rankings_Combined.csv", row.names = FALSE)

