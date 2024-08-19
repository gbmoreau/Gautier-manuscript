#### ORDER-LEVEL COMMUNITY COMPOSITION ANALYSIS ###################################################
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
# manuscript. I'll focus on order-level composition differences in this analysis.

# Load the data
load(file = "../results/phyloseq object for analysis.RData")




#### ORDER LEVEL DISTRIBUTION #####################################################################
# PURPOSE: Look at Order level distribution on a per-sample level.
ps.analysis.prop <-transform_sample_counts(ps.analysis, function(ASV) ASV/sum(ASV))

ps.analysis.prop.order <- tax_glom(ps.analysis.prop, "Order", NArm = FALSE)
order.table <- psmelt(ps.analysis.prop.order) # Organize in long format for ggplot.

# Summarize Orders for each Group
summarized.abundance.order <- order.table %>%
  group_by(Order, Group) %>%
  summarize(Relative.Abundance = (sum(Abundance) / 24))


# Replace NAs with "Unidentified"
summarized.abundance.order$Order <- as.character(replace_na(summarized.abundance.order$Order, replace = "Unidentified"))


#Collect taxa names for orders > 1% in any sample
summarized.abundance.order.above1 <- filter(summarized.abundance.order, Relative.Abundance >= 0.01)
summarized.abundance.order.above1$Order <- factor(summarized.abundance.order.above1$Order, 
                                                    levels = sort(unique(summarized.abundance.order.above1$Order)))
abundant.orders <- as.character(unique(summarized.abundance.order.above1$Order)) 
table(abundant.orders) # There are 8 Orders that have a Total Relative Abundance above 1%


# Make new data frames either including or excluding abundant orders with sample-level data.
summarized.abundance.order.sample <- order.table %>%
  group_by(Order, Sample) %>%
  summarize(Relative.Abundance = (sum(Abundance)))

summarized.abundance.order.sample$Order <- as.character(replace_na(summarized.abundance.order.sample$Order, 
                                                                     replace = "Unidentified"))

summarized.abundance.order.sample.abundant <- filter(summarized.abundance.order.sample, 
                                                      Order %in% abundant.orders)

summarized.abundance.order.sample.not.abundant <- filter(summarized.abundance.order.sample, 
                                                          ! Order %in% abundant.orders)

summarized.abundance.order.sample.not.abundant <- summarized.abundance.order.sample.not.abundant %>%
  group_by(Sample) %>%
  summarize(Relative.Abundance = (sum(Relative.Abundance))) %>%
  mutate(Order = "Order < 1%", .before = 1)


# Bring data back together and add Group metadata
summarized.abundance.order.sample.table <- rbind(summarized.abundance.order.sample.abundant, 
                                                  summarized.abundance.order.sample.not.abundant)


Sample.Group <- select(order.table, Sample, Group, Animal.Number)
Sample.Group <- Sample.Group %>%
  group_by(Sample) %>%
  slice(1L)

summarized.abundance.order.sample.table <- full_join(summarized.abundance.order.sample.table, Sample.Group, by = "Sample")
summarized.abundance.order.sample.table <- ungroup(summarized.abundance.order.sample.table)


# Organize data into a wide format

# ACHOLEPLASMATALES
per.sample.abundance.ach <- filter(summarized.abundance.order.sample.table, Order == "Acholeplasmatales")
per.sample.abundance.ach$"Acholeplasmatales Relative Abundance (%)" <- per.sample.abundance.ach$Relative.Abundance *100
per.sample.abundance.ach <- select(per.sample.abundance.ach, Sample, Animal.Number, "Acholeplasmatales Relative Abundance (%)")


# BACTEROIDALES
per.sample.abundance.bac <- filter(summarized.abundance.order.sample.table, Order == "Bacteroidales")
per.sample.abundance.bac$"Bacteroidales Relative Abundance (%)" <- per.sample.abundance.bac$Relative.Abundance *100
per.sample.abundance.bac <- select(per.sample.abundance.bac, Sample, "Bacteroidales Relative Abundance (%)")


# CLOSTRIDIALES
per.sample.abundance.clos <- filter(summarized.abundance.order.sample.table, Order == "Clostridiales")
per.sample.abundance.clos$"Clostridiales Relative Abundance (%)" <- per.sample.abundance.clos$Relative.Abundance *100
per.sample.abundance.clos <- select(per.sample.abundance.clos, Sample, "Clostridiales Relative Abundance (%)")


# ERYSIPELOTRIACHALES
per.sample.abundance.erys <- filter(summarized.abundance.order.sample.table, Order == "Erysipelotrichales")
per.sample.abundance.erys$"Erysipelotrichales Relative Abundance (%)" <- per.sample.abundance.erys$Relative.Abundance *100
per.sample.abundance.erys <- select(per.sample.abundance.erys, Sample, "Erysipelotrichales Relative Abundance (%)")


# LACHNOSPIRALES
per.sample.abundance.lachno <- filter(summarized.abundance.order.sample.table, Order == "Lachnospirales")
per.sample.abundance.lachno$"Lachnospirales Relative Abundance (%)" <- per.sample.abundance.lachno$Relative.Abundance *100
per.sample.abundance.lachno <- select(per.sample.abundance.lachno, Sample, "Lachnospirales Relative Abundance (%)")


# LACTOBACILLALES
per.sample.abundance.lacto <- filter(summarized.abundance.order.sample.table, Order == "Lactobacillales")
per.sample.abundance.lacto$"Lactobacillales Relative Abundance (%)" <- per.sample.abundance.lacto$Relative.Abundance *100
per.sample.abundance.lacto <- select(per.sample.abundance.lacto, Sample, "Lactobacillales Relative Abundance (%)")


# OSCILLOSPIRALES
per.sample.abundance.osci <- filter(summarized.abundance.order.sample.table, Order == "Oscillospirales")
per.sample.abundance.osci$"Oscillospirales Relative Abundance (%)" <- per.sample.abundance.osci$Relative.Abundance *100
per.sample.abundance.osci <- select(per.sample.abundance.osci, Sample, "Oscillospirales Relative Abundance (%)")


# PEPTOSTREPTOCOCCALES-TISSIERELLALES
per.sample.abundance.pepto <- filter(summarized.abundance.order.sample.table, Order == "Peptostreptococcales-Tissierellales")
per.sample.abundance.pepto$"Peptostreptococcales-Tissierellales Relative Abundance (%)" <- per.sample.abundance.pepto$Relative.Abundance *100
per.sample.abundance.pepto <- select(per.sample.abundance.pepto, Sample, "Peptostreptococcales-Tissierellales Relative Abundance (%)")


#ORDER <1%
per.sample.abundance.1p <- filter(summarized.abundance.order.sample.table, Order == "Order < 1%")
per.sample.abundance.1p$"Order < 1% Relative Abundance (%)" <- per.sample.abundance.1p$Relative.Abundance *100
per.sample.abundance.1p <- select(per.sample.abundance.1p, Sample, "Order < 1% Relative Abundance (%)")


per.sample.abundance <- full_join(per.sample.abundance.ach, per.sample.abundance.bac, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.clos, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.erys, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.lachno, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.lacto, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.osci, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.pepto, by = "Sample")
per.sample.abundance <- full_join(per.sample.abundance, per.sample.abundance.1p, by = "Sample")

#write.csv(per.sample.abundance, file = "../results/tables/composition per sample-order.csv", row.names = FALSE)


summarized.abundance.order.sample.table.figure <- summarized.abundance.order.sample.table
summarized.abundance.order.sample.table.figure$Relative.Abundance.Percentage <- summarized.abundance.order.sample.table.figure$Relative.Abundance * 100


# Plot community composition per sample
ggplot(summarized.abundance.order.sample.table.figure, aes(x = Animal.Number, y = Relative.Abundance.Percentage, fill = Order)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_x_discrete(limits = rev) +
  theme(axis.text = element_text(size = 5), 
        axis.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5)) +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 90))

#ggsave("../results/figures/per sample community composition-order.png", width = 6, height = 4)

# Plot with no Animal number labels.
ggplot(summarized.abundance.order.sample.table.figure, aes(x = Animal.Number, y = Relative.Abundance.Percentage, fill = Order)) +
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

#ggsave("../results/figures/per sample community composition-order-no names.png", width = 6, height = 4)


# I'll now plot abundance for each sample. Abundances are expressed as a percentage of each sample, so I'll format the 
# x axis as a percentage.

summarized.abundance.order.sample.table$Order <- factor(summarized.abundance.order.sample.table$Order, levels = c("Acholeplasmatales", "Bacteroidales", "Clostridiales", "Erysipelotrichales", "Lachnospirales", 
                                                                                                                  "Lactobacillales", "Oscillospirales", "Peptostreptococcales-Tissierellales","Order < 1%"))

summarized.abundance.order.sample.table$Group <- factor(summarized.abundance.order.sample.table$Group, 
                                                         levels = c("3 Wk Stress", "Baseline"), ordered = TRUE)

ggplot(summarized.abundance.order.sample.table, aes(x = Relative.Abundance, y = Order, color = Group)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(limits = rev) +
  labs( x = "Relative Abundance", y = NULL, color = "Group") + 
  scale_color_manual(name = "Group", 
                     values = c("#941650", "#1c7798"),
                     breaks = c("Baseline", "3 Wk Stress"), 
                     labels = c("Baseline", "3 Wk Stress"))

#ggsave("../results/figures/order abundance-boxplot.png", width = 6, height = 6)

summarized.abundance.order.sample.table$Group <- factor(summarized.abundance.order.sample.table$Group, 
                                                         levels = c("Baseline", "3 Wk Stress"))




#### SIGNIFICANT DIFFERENCES: BASELINE VS STRESSED SAMPLES ##########################################
# I'll now use a Wilcoxon rank test to look for significant differences between each of these order 
# groups. First, I'll filter the per-sample table for only abundant orders. I'll then perform a 
# Wilcoxon rank test and Bonferroni correction for multiple comparisons on each Order to identify 
# significant differences.

wilcox.order <- filter(summarized.abundance.order.sample.table, Order != "Order < 1%" )

# Now I'll calculate Wilcoxon Ranked Sum p values for each Order.

# ACHOLEPLASMATALES
wilcox.ach <- filter(wilcox.order, Order == "Acholeplasmatales")
wilcox.ach.baseine <- filter(wilcox.ach, Group == "Baseline")
wilcox.ach.stress <- filter(wilcox.ach, Group == "3 Wk Stress")
wilcox.test.ach <- wilcox.test(wilcox.ach.baseine$Relative.Abundance, 
                               wilcox.ach.stress$Relative.Abundance)$p.value


# BACTEROIDALES
wilcox.bac <- filter(wilcox.order, Order == "Bacteroidales")
wilcox.bac.baseine <- filter(wilcox.bac, Group == "Baseline")
wilcox.bac.stress <- filter(wilcox.bac, Group == "3 Wk Stress")
wilcox.test.bac <- wilcox.test(wilcox.bac.baseine$Relative.Abundance, 
                               wilcox.bac.stress$Relative.Abundance)$p.value


# CLOSTRIDIALES
wilcox.clos <- filter(wilcox.order, Order == "Clostridiales")
wilcox.clos.baseine <- filter(wilcox.clos, Group == "Baseline")
wilcox.clos.stress <- filter(wilcox.clos, Group == "3 Wk Stress")
wilcox.test.clos <- wilcox.test(wilcox.clos.baseine$Relative.Abundance, 
                               wilcox.clos.stress$Relative.Abundance)$p.value


# ERYSIPELOTRICHALES
wilcox.erys <- filter(wilcox.order, Order == "Erysipelotrichales")
wilcox.erys.baseine <- filter(wilcox.erys, Group == "Baseline")
wilcox.erys.stress <- filter(wilcox.erys, Group == "3 Wk Stress")
wilcox.test.erys <- wilcox.test(wilcox.erys.baseine$Relative.Abundance, 
                                wilcox.erys.stress$Relative.Abundance)$p.value


# LACHNOSPIRALES
wilcox.lachno <- filter(wilcox.order, Order == "Lachnospirales")
wilcox.lachno.baseine <- filter(wilcox.lachno, Group == "Baseline")
wilcox.lachno.stress <- filter(wilcox.lachno, Group == "3 Wk Stress")
wilcox.test.lachno <- wilcox.test(wilcox.lachno.baseine$Relative.Abundance, 
                                wilcox.lachno.stress$Relative.Abundance)$p.value


# LACTOBACILLALES
wilcox.lacto <- filter(wilcox.order, Order == "Lactobacillales")
wilcox.lacto.baseine <- filter(wilcox.lacto, Group == "Baseline")
wilcox.lacto.stress <- filter(wilcox.lacto, Group == "3 Wk Stress")
wilcox.test.lacto <- wilcox.test(wilcox.lacto.baseine$Relative.Abundance, 
                                wilcox.lacto.stress$Relative.Abundance)$p.value


# OSCILLOSPIRALES
wilcox.osci <- filter(wilcox.order, Order == "Oscillospirales")
wilcox.osci.baseine <- filter(wilcox.osci, Group == "Baseline")
wilcox.osci.stress <- filter(wilcox.osci, Group == "3 Wk Stress")
wilcox.test.osci <- wilcox.test(wilcox.osci.baseine$Relative.Abundance, 
                                wilcox.osci.stress$Relative.Abundance)$p.value


# PEPTOSTREPTOCOCCALES-TISSIERELLALES
wilcox.pepto <- filter(wilcox.order, Order == "Peptostreptococcales-Tissierellales")
wilcox.pepto.baseine <- filter(wilcox.pepto, Group == "Baseline")
wilcox.pepto.stress <- filter(wilcox.pepto, Group == "3 Wk Stress")
wilcox.test.pepto <- wilcox.test(wilcox.pepto.baseine$Relative.Abundance, 
                                wilcox.pepto.stress$Relative.Abundance)$p.value

Wilcox.Orders <- c("Acholeplasmatales", "Bacteroidales", "Clostridiales", "Erysipelotrichales", "Lachnospirales", 
                     "Lactobacillales", "Oscillospirales", "Peptostreptococcales-Tissierellales")




Wilcox.pvalues <- c(wilcox.test.ach, wilcox.test.bac, wilcox.test.clos, wilcox.test.erys, wilcox.test.lachno,
                    wilcox.test.lacto, wilcox.test.osci, wilcox.test.pepto)

Wilcox.pvalues.Bonferonni <- p.adjust(Wilcox.pvalues, method = "bonferroni", n = 8)

Wilcox.all <- cbind(Wilcox.Orders, Wilcox.pvalues)
Wilcox.all <- cbind(Wilcox.all, Wilcox.pvalues.Bonferonni)
Wilcox.all <- as.data.frame(Wilcox.all)
Wilcox.all <- plyr::rename(Wilcox.all, replace = c("Wilcox.Orders" = "Order", 
                                                   "Wilcox.pvalues" = "Raw p value", 
                                                   "Wilcox.pvalues.Bonferonni" = "Bonferroni-adjusted p value"))

#write.csv(Wilcox.all, file = "../results/tables/wilcoxon rank comparison-order composition.csv", row.names = FALSE)





#### RANDOM FOREST PLOT OF FAMILIES ##############################################################
# PURPOSE: Make figure of Random Forest Importance looking at the Order level

ps.analysis.prop <-transform_sample_counts(ps.analysis, function(ASV) ASV/sum(ASV))

ps.analysis.prop.order<- tax_glom(ps.analysis.prop, "Order", NArm = FALSE)
order.table <- psmelt(ps.analysis.prop.order) # Organize in long format for ggplot.

# I'll reshape the data into a wide format for analysis.
order.table.wide <- select(order.table, Order, Sample, Abundance, Group, Animal.Number)
order.table.wide <- spread(order.table.wide, key = Order, value = Abundance)


# Now I'll split the data into Baseline and 3 Week Stress groups
order.table.wide.Baseline <- filter(order.table.wide, Group == "Baseline")
order.table.wide.Baseline$Group <- factor(order.table.wide.Baseline$Group, levels = c("Baseline"))

order.table.wide.3Wk.Stress <- filter(order.table.wide, Group == "3 Wk Stress")
order.table.wide.3Wk.Stress$Group <- factor(order.table.wide.3Wk.Stress$Group, levels = c("3 Wk Stress"))




### TRAINING AND TEST SET GENERATION: STRESS TREATMENT ###############################################
# There are 24 mice per group, which is a pretty good sample size for mouse experiments. Because of 
# this, I'll separate the data into training and test sets to better optimize the random forest 
# model and check model performance in an outside population.


# Set seed for reproducibility
#sample(1:1000, 1) # It selected 794
set.seed(794) 

# Now I'll separate each group into training and test sets.

### BASELINE ###
testIndex.Baseline <- createDataPartition(order.table.wide.Baseline$Group,
                                          p = 0.30,
                                          list = FALSE,
                                          times = 1)

featureTest.Baseline <- order.table.wide.Baseline[testIndex.Baseline,]
featureTrain.Baseline  <- order.table.wide.Baseline[-testIndex.Baseline,]


### 3 WEEK STRESS ###
testIndex.3Wk.Stress <- createDataPartition(order.table.wide.3Wk.Stress$Group,
                                            p = 0.30,
                                            list = FALSE,
                                            times = 1)

featureTest.3Wk.Stress <- order.table.wide.3Wk.Stress[testIndex.3Wk.Stress,]
featureTrain.3Wk.Stress  <- order.table.wide.3Wk.Stress[-testIndex.3Wk.Stress,]



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

# By default mtry is set to the square root of the number of features. There are 23 features in 
# this data set, meaning mtry be default is set to 4.8. I'll test mtry values around 5 to see if 
# this improves model error.


### TUNING MTRY ###
model.rf.mtry3 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 3)
model.rf.mtry3 # The OOB error rate is 3.12%

model.rf.mtry4 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 4)
model.rf.mtry4 # The OOB error rate is 3.12%

model.rf.mtry5 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 5)
model.rf.mtry5 # The OOB error rate is 3.12%

model.rf.mtry6 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 6)
model.rf.mtry6 # The OOB error rate is 3.12%

model.rf.mtry7 <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 7)
model.rf.mtry7 # The OOB error rate is 6.25%

# All random forest models performed very well at correctly classifying Baseline and 3 week stress 
# samples, regardless of mtry value. I'll chose the standard mtry value  of 5.

model.rf.stress.final <- randomForest(x = predictors.stress.Train, y = outcome.stress.Train, ntree = 1000, mtry = 5)
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

# Overall, 6/8 (75%) of Baseline samples were correctly classified, while 8/8 (100%)  of 3 week
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

# The RF model performs well on test set data, with a 93.8% AUC value.




#### IMPORTANCE VALUES FROM RANDOM FOREST MODEL ####################################################
# Now that I've generated the model, I'll look at which features best discriminate between groups as
# selected by the model. I'll do this by checking the importance values from the Random Forest model, 
# using Gini index as the measure of importance.

importance.rf.stress <- importance(model.rf.stress.final, type = 2)

# Format the importance data as a data frame
importance.rf.stress <- cbind(Order = rownames(importance.rf.stress), importance.rf.stress)
rownames(importance.rf.stress) <- NULL
importance.rf.stress <- as.data.frame(importance.rf.stress)
importance.rf.stress$MeanDecreaseGini <- as.numeric(importance.rf.stress$MeanDecreaseGini)
str(importance.rf.stress$MeanDecreaseGini)

# This is the table of orders listed by mean decrease in Gini Index. I'll rank the features according
# to Gini index to make the data easier to interpret.

RF.ranking.stress <- importance.rf.stress %>%
  arrange(desc(MeanDecreaseGini))

RF.ranking.stress <- plyr::rename(RF.ranking.stress, replace = c("MeanDecreaseGini" = "Stress.Mean.Decrease.Gini"))
RF.ranking.stress$Stress.RF.Rank <- rank(-RF.ranking.stress$Stress.Mean.Decrease.Gini)

# I also want to add whether orders are enriched at baseline or after 3 weeks of stress treatment.
order.table.enrichment.stress <- order.table %>%
  group_by(Order, Group) %>%
  summarize(Average.Abundance.per.Group = (mean(Abundance))) %>%
  filter(any(Average.Abundance.per.Group > 0)) # Removes any orders with an abundance of 0 in both groups.

order.enrichment.baseline <- filter(order.table.enrichment.stress, Group == "Baseline")
order.enrichment.baseline <- plyr::rename(order.enrichment.baseline, replace = c("Average.Abundance.per.Group" = "Baseline.Average.Abundance.per.Group"))
order.enrichment.baseline <- select(order.enrichment.baseline, Order, Baseline.Average.Abundance.per.Group)

order.enrichment.3Wk.stress <- filter(order.table.enrichment.stress, Group == "3 Wk Stress")
order.enrichment.3Wk.stress <- plyr::rename(order.enrichment.3Wk.stress, replace = c("Average.Abundance.per.Group" = "3Wk.stress.Average.Abundance.per.Group"))
order.enrichment.3Wk.stress <- select(order.enrichment.3Wk.stress, Order, '3Wk.stress.Average.Abundance.per.Group')

order.enrichment.stress <- full_join(order.enrichment.baseline, order.enrichment.3Wk.stress, by = "Order")
order.enrichment.stress$Enrichment.Stress <- NA
order.enrichment.stress$Enrichment.Stress[order.enrichment.stress$Baseline.Average.Abundance.per.Group > order.enrichment.stress$`3Wk.stress.Average.Abundance.per.Group`] <- "Baseline"
order.enrichment.stress$Enrichment.Stress[order.enrichment.stress$Baseline.Average.Abundance.per.Group < order.enrichment.stress$`3Wk.stress.Average.Abundance.per.Group`] <- "3 Wk Stress"

order.enrichment <- select(order.enrichment.stress, Order, Enrichment.Stress)

RF.ranking.stress <- full_join(RF.ranking.stress, order.enrichment, by = "Order")

RF.ranking.stress$Order <- factor(RF.ranking.stress$Order, levels = RF.ranking.stress$Order)


ggplot(RF.ranking.stress, aes(x = Stress.Mean.Decrease.Gini, y = Order, color = Enrichment.Stress)) +
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

#ggsave("../results/figures/RF Rankings-Baseline vs Stress-order.png", width = 3, height = 4)

