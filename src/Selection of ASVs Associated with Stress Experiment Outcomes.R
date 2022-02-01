#### SELECTION OF ASVs ASSOCIATED WITH STRESS EXPERIMENT OUTCOMES ##################################
# NAME: G. Brett Moreau
# DATE: November 29, 2021

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

#install.packages("randomForest")
library(randomForest)
packageVersion("randomForest") # I'm using version 4.6.14

#install.packages("pROC")
library(pROC)
packageVersion("pROC") # I'm using version 1.18.0


#### INTRODUCTION ###################################################################################
# I've broadly looked at differences in diversity and community composition with the Stress 
# experiment for the Gautier TUMI Pilot. I'll now move on to identify ASVs that best discriminate 
# between groups. This experiment looks at mice before and after Stress treatment, so the major 
# comparisons will be Naive vs Stressed mice in the Post timepoint and Pre vs. Post timepoints for 
# the Stressed mice. There are 6 mice/group.




#### ORGANIZATION OF ASVs FROM NAIVE AND STRESSED POST-STRESS SAMPLES ###############################

#setwd("..") # I need to set the working directory up one to the parent directory for the code to 
# work properly. If this has been done previously during the same session then this step should be 
# skipped, or the working directory will be moved out of the parent directory.

load(file = "./results/phyloseq objects by experiment.RData")
rm(ps.mucin, seqtab.nochim)

ps.stress.ASV <- subset_samples(ps.stress, Condition == "Naïve Post" | Condition == "Stress Post")
ASV.table.stress <- psmelt(ps.stress.ASV)

length(unique(ASV.table.stress$OTU)) #1795 unique ASVs across all samples.

# I want to limit the data set to only ASVs that are present in at least one group. To do this, I'll
# first summarize the Abundance of each ASV in each group, then select all ASVs that are >0 in any
# group.

ASV.table.summary.stress <- ASV.table.stress %>%
  group_by(OTU, Condition) %>%
  summarize(Abundance.per.Condition = (sum(Abundance))) %>%
  filter(any(Abundance.per.Condition > 0)) # Removes any ASVs with an abundance of 0 in both groups.

length(unique(ASV.table.summary.stress$OTU)) # There are a total of 392 unique ASVs in these samples.


# Now I'll collect the ASV names for these ASVs. These will be used to merge only the ASV Abundance
# values for ASVs present in at least one group.
unique.ASVs.stress <- as.data.frame(unique(ASV.table.summary.stress$OTU))
unique.ASVs.stress <- plyr::rename(unique.ASVs.stress, replace = c("unique(ASV.table.summary.stress$OTU)" = "OTU"))

ASV.table.samples.stress <- left_join(unique.ASVs.stress, ASV.table.stress, by = "OTU")

length(unique(ASV.table.samples.stress$OTU)) # 392 unique ASVs as expected.


# I'll reshape the data into a wide format for analysis.
ASV.table.wide.stress <- select(ASV.table.samples.stress, OTU, Sample, Abundance, Condition, Animal.Number)
ASV.table.wide.stress <- spread(ASV.table.wide.stress, key = OTU, value = Abundance)
ASV.table.wide.stress$Condition <- factor(ASV.table.wide.stress$Condition, levels = c("Naïve Post", "Stress Post"))




### RANDOM FOREST MODEL GENERATION ##################################################################
# There are 6 mice per group, which is a pretty small sample size. Because of this, I'll forego 
# using a training and test set and instead generate the random forest model on the entire data set.

# First, I'll separate the data into predictors and outcomes.
predictors <- select(ASV.table.wide.stress, -Sample, -Condition, -Animal.Number) # Remove metadata
outcome <- as.factor(ASV.table.wide.stress$Condition)


# The parameter I'll  tune for the random forest model is mtry, the number of features sampled 
# at each node of the decision tree. For this random forest, I'll use ntree (the number of trees in 
# the forest) = 1000, which is a pretty large number but not too computationally taxing. Then I'll 
# look at a variety of different mtry values to see which minimizes the Out-of-Bag error rate.

# By default mtry is set to the square root of the number of features. There are 392 features in 
# this data set, meaning mtry be default is set to 20. I'll test mtry values around 20 to see if 
# this improves model error.

# Set seed for reproducibility
#sample(1:1000, 1) # It selected 763
set.seed(763) 

### TUNING MTRY ###
model.rf.mtry16 <- randomForest(x = predictors, y = outcome, ntree = 1000, mtry = 16)
model.rf.mtry16 # The OOB error rate is 8.33%

model.rf.mtry18 <- randomForest(x = predictors, y = outcome, ntree = 1000, mtry = 18)
model.rf.mtry18 # The OOB error rate is 8.33%

model.rf.mtry20 <- randomForest(x = predictors, y = outcome, ntree = 1000, mtry = 20)
model.rf.mtry20 # The OOB error rate is 8.33%

model.rf.mtry22 <- randomForest(x = predictors, y = outcome, ntree = 1000, mtry = 22)
model.rf.mtry22 # The OOB error rate is 8.33%

model.rf.mtry24 <- randomForest(x = predictors, y = outcome, ntree = 1000, mtry = 24)
model.rf.mtry24 # The OOB error rate is 8.33%

# All random forest models performed the same, with an OOB error rate of 8.3%. This result was
# very consistent for all models, so I'll just go with the default mtry value of mtry == 20 for
# the final model.

model.rf.stress.final <- randomForest(x = predictors, y = outcome, ntree = 1000, mtry = 20)
model.rf.stress.final


roc(outcome, model.rf.stress.final$votes[,1], plot = TRUE, legacy.axes = TRUE, percent = TRUE, main = "Training Set ROC Curve",
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# The RF model performs well on training data, with a 88.9% AUC. 




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

# I now want to combine the RF ranking data with the taxonomic information for identified ASVs. I'll
# pull this information from the phyloseq object metadata.

ASV.taxonomy <- as.data.frame(ps.stress@tax_table)
ASV.taxonomy <- cbind(OTU = rownames(ASV.taxonomy), ASV.taxonomy)
rownames(ASV.taxonomy) <- NULL

importance.rf.stress.taxonomy <- left_join(importance.rf.stress, ASV.taxonomy, by = "OTU")

# This is the table of ASVs listed by mean decrease in Gini Index. I'll rank the features according
# to Gini index to make the data easier to interpret.

RF.ranking.stress <- importance.rf.stress.taxonomy %>%
  arrange(desc(MeanDecreaseGini))

RF.ranking.stress <- plyr::rename(RF.ranking.stress, replace = c("MeanDecreaseGini" = "Stress.Mean.Decrease.Gini"))
RF.ranking.stress$Stress.RF.Rank <- rank(-RF.ranking.stress$Stress.Mean.Decrease.Gini)

# I also want to add whether ASVs are enriched in the Naive or Stressed group.
ASV.table.enrichment.stress <- ASV.table.stress %>%
  group_by(OTU, Condition) %>%
  summarize(Average.Abundance.per.Group = (mean(Abundance))) %>%
  filter(any(Average.Abundance.per.Group > 0)) # Removes any ASVs with an abundance of 0 in both groups.


ASV.enrichment.Naive <- filter(ASV.table.enrichment.stress, Condition == "Naïve Post")
ASV.enrichment.Naive <- plyr::rename(ASV.enrichment.Naive, replace = c("Average.Abundance.per.Group" = "Naive.Average.Abundance.per.Group"))
ASV.enrichment.Naive <- select(ASV.enrichment.Naive, OTU, Naive.Average.Abundance.per.Group)

ASV.enrichment.Stress <- filter(ASV.table.enrichment.stress, Condition == "Stress Post")
ASV.enrichment.Stress <- plyr::rename(ASV.enrichment.Stress, replace = c("Average.Abundance.per.Group" = "Stress.Average.Abundance.per.Group"))
ASV.enrichment.Stress <- select(ASV.enrichment.Stress, OTU, 'Stress.Average.Abundance.per.Group')

ASV.enrichment.stress <- full_join(ASV.enrichment.Naive, ASV.enrichment.Stress, by = "OTU")
ASV.enrichment.stress$Enrichment.Stress <- NA
ASV.enrichment.stress$Enrichment.Stress[ASV.enrichment.stress$Naive.Average.Abundance.per.Group > ASV.enrichment.stress$Stress.Average.Abundance.per.Group] <- "Naive Post"
ASV.enrichment.stress$Enrichment.Stress[ASV.enrichment.stress$Naive.Average.Abundance.per.Group < ASV.enrichment.stress$Stress.Average.Abundance.per.Group] <- "Stress Post"

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

RF.ranking.stress.figure <- filter(RF.ranking.stress, Stress.Mean.Decrease.Gini > 0.025)

ggplot(RF.ranking.stress.figure, aes(x = Stress.Mean.Decrease.Gini, y = ASV.Family, color = Enrichment.Stress)) +
  geom_point(size = 3) +
  labs(title = "Top Ranked Random Forest Predictors", x = "Mean Decrease in Node Impurity (Gini)", y = NULL) +
  scale_color_manual(name = "Enrichment", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Naive Post", "Stress Post"), 
                     labels = c("Naive Post", "Stress Post")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#ggsave("./results/figures/stress/RF Rankings_Naive_Post_vs_Stress_Post.png", width = 7, height = 5)


# Print Random Forest importance rankings
write.csv(RF.ranking.stress, file = "./results/tables/stress/RF Rankings_Naive_Post_vs_Stress_Post.csv", row.names = FALSE)




#### ORGANIZATION OF ASVs FROM STRESSED PRE VS POST SAMPLES #####################################
# I'll now look at differences in the Stressed mice Pre-stress vs Post-stress. 

ps.stress.ASV.time <- subset_samples(ps.stress, Condition == "Stress Pre" | Condition == "Stress Post")
ASV.table.stress.time <- psmelt(ps.stress.ASV.time)
View(ASV.table.stress.time)
length(unique(ASV.table.stress.time$OTU)) #1795 unique ASVs across all samples.

# I want to limit the data set to only ASVs that are present in at least one group. To do this, I'll
# first summarize the Abundance of each ASV in each group, then select all ASVs that are >0 in any
# group.

ASV.table.summary.stress.time <- ASV.table.stress.time %>%
  group_by(OTU, Condition) %>%
  summarize(Abundance.per.Condition = (sum(Abundance))) %>%
  filter(any(Abundance.per.Condition > 0)) # Removes any ASVs with an abundance of 0 in both groups.

length(unique(ASV.table.summary.stress.time$OTU)) # There are a total of 370 unique ASVs in these samples.


# Now I'll collect the ASV names for these ASVs. These will be used to merge only the ASV Abundance
# values for ASVs present in at least one group.
unique.ASVs.stress.time <- as.data.frame(unique(ASV.table.summary.stress.time$OTU))
unique.ASVs.stress.time <- plyr::rename(unique.ASVs.stress.time, replace = c("unique(ASV.table.summary.stress.time$OTU)" = "OTU"))

ASV.table.samples.stress.time <- left_join(unique.ASVs.stress.time, ASV.table.stress.time, by = "OTU")

length(unique(ASV.table.samples.stress.time$OTU)) # 370 unique ASVs as expected.


# I'll reshape the data into a wide format for analysis.
ASV.table.wide.stress.time <- select(ASV.table.samples.stress.time, OTU, Sample, Abundance, Condition, Animal.Number)
ASV.table.wide.stress.time <- spread(ASV.table.wide.stress.time, key = OTU, value = Abundance)
ASV.table.wide.stress.time$Condition <- factor(ASV.table.wide.stress.time$Condition, levels = c("Stress Pre", "Stress Post"))

View(ASV.table.wide.stress.time)


### RANDOM FOREST MODEL GENERATION ##################################################################
# There are 6 mice per group, which is a pretty small sample size. Because of this, I'll forego 
# using a training and test set and instead generate the random forest model on the entire data set.

# First, I'll separate the data into predictors and outcomes.
predictors.time <- select(ASV.table.wide.stress.time, -Sample, -Condition, -Animal.Number) # Remove metadata
outcome.time <- as.factor(ASV.table.wide.stress.time$Condition)

# The parameter I'll  tune for the random forest model is mtry, the number of features sampled 
# at each node of the decision tree. For this random forest, I'll use ntree (the number of trees in 
# the forest) = 1000, which is a pretty large number but not too computationally taxing. Then I'll 
# look at a variety of different mtry values to see which minimizes the Out-of-Bag error rate.

# By default mtry is set to the square root of the number of features. There are 370 features in 
# this data set, meaning mtry be default is set to 19. I'll test mtry values around 19 to see if 
# this improves model error.

# Set seed for reproducibility
#sample(1:1000, 1) # It selected 745
set.seed(745) 

### TUNING MTRY ###
model.rf.mtry15 <- randomForest(x = predictors.time, y = outcome.time, ntree = 1000, mtry = 15)
model.rf.mtry15 # The OOB error rate is 8.33%

model.rf.mtry17 <- randomForest(x = predictors.time, y = outcome.time, ntree = 1000, mtry = 17)
model.rf.mtry17 # The OOB error rate is 8.33%

model.rf.mtry19 <- randomForest(x = predictors.time, y = outcome.time, ntree = 1000, mtry = 19)
model.rf.mtry19 # The OOB error rate is 8.33%

model.rf.mtry21 <- randomForest(x = predictors.time, y = outcome.time, ntree = 1000, mtry = 21)
model.rf.mtry21 # The OOB error rate is 8.33%

model.rf.mtry23 <- randomForest(x = predictors.time, y = outcome.time, ntree = 1000, mtry = 23)
model.rf.mtry23 # The OOB error rate is 8.33%

# All random forest models performed the same, with an OOB error rate of 8.3%. This result was
# pretty consistent for all models, so I'll just go with the default mtry value of mtry == 19 
# for the final model.

model.rf.stress.time.final <- randomForest(x = predictors.time, y = outcome.time, ntree = 1000, mtry = 19)
model.rf.stress.time.final


roc(outcome.time, model.rf.stress.time.final$votes[,1], plot = TRUE, legacy.axes = TRUE, percent = TRUE, main = "ROC Curve",
    xlab = "False Positive Percentage", ylab = "True Positive Percentage", col = "#377eb8", lwd = 4, print.auc = TRUE)

# The RF model performs pretty well on training data, but not as well as previous models. It has 
# an AUC of 83.3% 




#### IMPORTANCE VALUES FROM RANDOM FOREST MODEL ####################################################
# Now that I've generated the model, I'll look at which features best discriminate between groups as
# selected by the model. I'll do this by checking the importance values from the Random Forest model, 
# using Gini index as the measure of importance.

importance.rf.stress.time <- importance(model.rf.stress.time.final, type = 2)

# Format the importance data as a data frame
importance.rf.stress.time <- cbind(OTU = rownames(importance.rf.stress.time), importance.rf.stress.time)
rownames(importance.rf.stress.time) <- NULL
importance.rf.stress.time <- as.data.frame(importance.rf.stress.time)
importance.rf.stress.time$MeanDecreaseGini <- as.numeric(importance.rf.stress.time$MeanDecreaseGini)

importance.rf.stress.time.taxonomy <- left_join(importance.rf.stress.time, ASV.taxonomy, by = "OTU")

# This is the table of ASVs listed by mean decrease in Gini Index. I'll rank the features according
# to Gini index to make the data easier to interpret.

RF.ranking.stress.time <- importance.rf.stress.time.taxonomy %>%
  arrange(desc(MeanDecreaseGini))

RF.ranking.stress.time <- plyr::rename(RF.ranking.stress.time, replace = c("MeanDecreaseGini" = "Pre.vs.Post.Mean.Decrease.Gini"))
RF.ranking.stress.time$Pre.Post.RF.Rank <- rank(-RF.ranking.stress.time$Pre.vs.Post.Mean.Decrease.Gini)

# I also want to add whether ASVs are enriched in the Naive or Stressed group.
ASV.table.enrichment.stress.time <- ASV.table.stress.time %>%
  group_by(OTU, Condition) %>%
  summarize(Average.Abundance.per.Group = (mean(Abundance))) %>%
  filter(any(Average.Abundance.per.Group > 0)) # Removes any ASVs with an abundance of 0 in both groups.

ASV.enrichment.Pre <- filter(ASV.table.enrichment.stress.time, Condition == "Stress Pre")
ASV.enrichment.Pre <- plyr::rename(ASV.enrichment.Pre, replace = c("Average.Abundance.per.Group" = "Pre.Average.Abundance.per.Group"))
ASV.enrichment.Pre <- select(ASV.enrichment.Pre, OTU, Pre.Average.Abundance.per.Group)

ASV.enrichment.Post <- filter(ASV.table.enrichment.stress.time, Condition == "Stress Post")
ASV.enrichment.Post <- plyr::rename(ASV.enrichment.Post, replace = c("Average.Abundance.per.Group" = "Post.Average.Abundance.per.Group"))
ASV.enrichment.Post <- select(ASV.enrichment.Post, OTU, 'Post.Average.Abundance.per.Group')

ASV.enrichment.stress.time <- full_join(ASV.enrichment.Pre, ASV.enrichment.Post, by = "OTU")
ASV.enrichment.stress.time$Enrichment.Pre.Post <- NA
ASV.enrichment.stress.time$Enrichment.Pre.Post[ASV.enrichment.stress.time$Pre.Average.Abundance.per.Group > ASV.enrichment.stress.time$Post.Average.Abundance.per.Group] <- "Stress Pre"
ASV.enrichment.stress.time$Enrichment.Pre.Post[ASV.enrichment.stress.time$Pre.Average.Abundance.per.Group < ASV.enrichment.stress.time$Post.Average.Abundance.per.Group] <- "Stress Post"

ASV.enrichment.time <- select(ASV.enrichment.stress.time, OTU, Enrichment.Pre.Post)

RF.ranking.stress.time <- full_join(RF.ranking.stress.time, ASV.enrichment.time, by = "OTU")

RF.ranking.stress.time$OTU <- factor(RF.ranking.stress.time$OTU, levels = RF.ranking.stress.time$OTU)

# I'll make new names for each figure, which will concatenate the ASV number and the Family
# designation for each ASV.
RF.ranking.stress.time$ASV.Family <- NA
RF.ranking.stress.time$ASV.Family <- paste(RF.ranking.stress.time$OTU, " (", RF.ranking.stress.time$Family, ")")
RF.ranking.stress.time$ASV.Family <- factor(RF.ranking.stress.time$ASV.Family, levels = RF.ranking.stress.time$ASV.Family)


ggplot(RF.ranking.stress.time, aes(x = Pre.vs.Post.Mean.Decrease.Gini, y = ASV.Family, color = Enrichment.Pre.Post)) +
  geom_point(size = 3) +
  labs(title = "Top Ranked Random Forest Predictors", x = "Mean Decrease in Node Impurity (Gini)", y = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# To reduce the number of ASVs graphed on the figure, I'm going to set an arbitrarily cutoff 
# for the Gini Index to include only the most important ASVs.

RF.ranking.stress.time.figure <- filter(RF.ranking.stress.time, Pre.vs.Post.Mean.Decrease.Gini > 0.025)

ggplot(RF.ranking.stress.time.figure, aes(x = Pre.vs.Post.Mean.Decrease.Gini, y = ASV.Family, color = Enrichment.Pre.Post)) +
  geom_point(size = 3) +
  labs(title = "Top Ranked Random Forest Predictors", x = "Mean Decrease in Node Impurity (Gini)", y = NULL) +
  scale_color_manual(name = "Enrichment", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Stress Pre", "Stress Post"), 
                     labels = c("Stress Pre", "Stress Post")) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#ggsave("./results/figures/stress/RF Rankings_Stress_Pre_vs_Stress_Post.png", width = 7, height = 5)


# Print Random Forest importance rankings
#write.csv(RF.ranking.stress.time, file = "./results/tables/stress/RF Rankings_Stress_Pre_vs_Stress_Post.csv", row.names = FALSE)




#### OVERLAP BETWEEN ASVs FROM STRESS AND MUCIN EXPERIMENTS #########################################
# While the Mucin supplementation experiment looks at the effect of mucin supplementation on stress,
# the "Baseline" and "3 Weeks" samples are identical to the "Stress Pre" and "Stress Post" groups in
# this experiment. Therefore, I'll look for overlap between the ASVs that discriminate in each 
# experiment.


# Read in data from mucin experiemnt.
RF.ranking.baseline.vs.stress <- read.csv("./results/tables/mucin/RF Rankings_Baseline_vs_Stress.csv")


RF.ranking.stress.all <- select(RF.ranking.stress.time, OTU, Pre.vs.Post.Mean.Decrease.Gini, Pre.Post.RF.Rank, Enrichment.Pre.Post)

RF.ranking.stress.all <- full_join(RF.ranking.stress.all, RF.ranking.baseline.vs.stress, by = "OTU")
RF.ranking.stress.all <- select(RF.ranking.stress.all, OTU, Stress.Mean.Decrease.Gini, Pre.vs.Post.Mean.Decrease.Gini,
                                Stress.RF.Rank, Pre.Post.RF.Rank, Enrichment.Stress, Enrichment.Pre.Post, everything())

View(RF.ranking.stress.all)

# It doesn't look like there's a whole lot of overlap between data sets. To look for sure, I'll arbitrarily
# limit the data set to ASVs that were in the top 50 for each outcome by mean decrease in Gini Index.

RF.ranking.stress.all.important <- filter(RF.ranking.stress.all, Stress.RF.Rank < 50 & Pre.Post.RF.Rank < 50)                              

View(RF.ranking.stress.all.important)

# There are 10 ASVs that are ranked in the top 50 for the stress arms of both the Mucin and Stress 
# experiments. All 10 ASVs are enriched consistently (i.e., enriched both at Baseline (Mucin) 
# and Pre-Stress (Stress) or enriched both after 3 weeks (Mucin) and Post-Stress (Stress)). 

# 6/10 are from Family Lachnospiraceae. Of these, 5/6 were enriched post-stress treatment. In addition,
# one of these ASVs is from Family Lactobacillus, which has been implicated previously in response to
# stress. In both experiments, the Lactobacillus ASV was enriched pre-stress, which is consistent with
# previous reports.

# I'll now print these importance rankings.
#write.csv(RF.ranking.stress.all.important, file = "./results/tables/stress/RF Rankings_Stress and Mucin Experiment ASVs.csv")
