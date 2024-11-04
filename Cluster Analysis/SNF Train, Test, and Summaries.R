

# """
#
# Filename: SNF Train, Test, and Summaries.R
# Version: 2024.09.0+375 RStudio 
# Author: Sanaz Nazari
# Date Created: Nov 3, 2024
# Last Modified: Nov 3, 2024
# Description: This script performs SNF on train and test data and creates summaries.
#
# """

#-----------------------
# Load required library
library(dplyr)
library(SNFtool)

###################################################################################################

# SNF

# Hyperparameters
k <- 20       # number of neighbors, usually (10~30)
alpha <- 0.5  # hyperparameter, usually (0.3~0.8)
t <- 20       # Number of Iterations, usually (10~20)
c <- 3        # number of clusters

# Function to standardize and set dimnames for data layers
prepare_data_layer <- function(data, columns, dim_names) {
  layer <- data[, columns, drop = FALSE]
  dimnames(layer) <- list(1:nrow(layer), dim_names)
  layer <- standardNormalization(layer)
  return(layer)
}

# Train data
data.orig <- read.csv("train.labels.asd.td2.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Set layers
data1 <- prepare_data_layer(data.orig, c("coso", "soc"), c("coso", "soc"))
data2 <- prepare_data_layer(data.orig, c("mtr", "fm"), c("mtr", "fm"))
data3 <- prepare_data_layer(data.orig, c("com", "rl", "el"), c("com", "rl", "el"))

data.tr.list <- list(data1, data2, data3)

# Pairwise squared Euclidean distances
dist_list <- lapply(data.tr.list, function(layer) (dist2(as.matrix(layer), as.matrix(layer)))^(1/2))

# Construct similarity graphs
w_list <- lapply(dist_list, function(dist_matrix) affinityMatrix(dist_matrix, k, alpha))

# Fusing all graphs
overall_matrix <- SNF(w_list, k, t)

# Get cluster labels for each data point by spectral clustering
labels.orig <- spectralClustering(overall_matrix, c)

# Merge labels to data
data.orig <- data.frame(data.orig, labels.orig)

######################################################################################

# Predict test labels from train labels

test <- read.csv("data.test.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Set layers for test data
data1_test <- prepare_data_layer(test, c(7, 14), c("coso", "soc"))
data2_test <- prepare_data_layer(test, c(13, 17), c("mtr", "fm"))
data3_test <- prepare_data_layer(test, c(11, 19:20), c("com", "rl", "el"))

data.ts.list <- list(data1_test, data2_test, data3_test)

# Groups from train labels
groups <- data.orig$labels.orig

# Predict clusters for unseen data
pred.labels <- groupPredict(data.tr.list, data.ts.list, groups, k, alpha, t, method = 1)

# Merge predicted labels to test data
data.ts <- data.frame(test, pred.labels[(nrow(data.orig) + 1):length(pred.labels)])
colnames(data.ts)[ncol(data.ts)] <- "pred.labels"

# Save results to CSV for future use
#write.csv(data.ts, "predicted_test_labels.csv", row.names = FALSE)

#####################################################################################

#Summary

#------------------------------
# Create summary for train data

# Read in data
data.orig <- read.csv("train.labels.asd.td2.csv", header = TRUE, na.strings = c("NULL", "NA", ""))
data.orig$Age <- rowMeans(data.orig[, c("a.age", "v.age", "m.age")], na.rm = TRUE)

# Define a function to generate summaries
generate_summary <- function(data, cluster_label) {
  data %>%
    group_by(diag, gender) %>%
    summarise(Age = round(mean(Age), 1),
              SA = round(mean(coso), 1), SOC = round(mean(soc), 1),
              MTR = round(mean(mtr, na.rm = TRUE), 1), FM = round(mean(fm), 1),
              COM = round(mean(com), 1), RL = round(mean(rl, na.rm = TRUE), 1), EL = round(mean(el, na.rm = TRUE), 1),
              n = n(), .groups = "drop") %>%
    mutate(Cluster = cluster_label)
}

# Separate ASD and TD data
data.orig.asd <- subset(data.orig, diag == "ASD")
data.orig.td <- subset(data.orig, diag == "TD")

# Create summaries for ASD
summary.asd <- generate_summary(data.orig.asd, "All data")
summary.asd.c1 <- generate_summary(subset(data.orig.asd, labels.orig == 1), "High (C1)")
summary.asd.c2 <- generate_summary(subset(data.orig.asd, labels.orig == 2), "Medium (C2)")
summary.asd.c3 <- generate_summary(subset(data.orig.asd, labels.orig == 3), "Low (C3)")

# Combine ASD summaries
summary.asd <- rbind(summary.asd, summary.asd.c1, summary.asd.c2, summary.asd.c3) %>%
  select(Cluster, everything())

# Create summaries for TD
summary.td <- generate_summary(data.orig.td, "All data")
summary.td.c1 <- generate_summary(subset(data.orig.td, labels.orig == 1), "High (C1)")
summary.td.c2 <- generate_summary(subset(data.orig.td, labels.orig == 2), "Medium (C2)")

# Combine TD summaries
summary.td <- rbind(summary.td, summary.td.c1, summary.td.c2) %>%
  select(Cluster, everything())

# Print ASD and TD tables
print(summary.asd)
print(summary.td)

#-----------------------------
#create summary for test data

# Read in test data
data <- read.csv("test.labels.asd.td.csv", header = TRUE, na.strings = c("NULL", "NA", ""))
data$Age <- rowMeans(data[, c("a.age", "v.age", "m.age")], na.rm = TRUE)

# Define a function to generate summaries
generate_summary <- function(data, cluster_label) {
  data %>%
    group_by(diag, gender) %>%
    summarise(Age = round(mean(Age), 1),
              SA = round(mean(coso), 1), SOC = round(mean(soc), 1),
              MTR = round(mean(mtr, na.rm = TRUE), 1), FM = round(mean(fm), 1),
              COM = round(mean(com), 1), RL = round(mean(rl, na.rm = TRUE), 1), EL = round(mean(el, na.rm = TRUE), 1),
              n = n(), .groups = "drop") %>%
    mutate(Cluster = cluster_label)
}

# Separate ASD and TD data
data.asd <- subset(data, diag == "ASD")
data.td <- subset(data, diag == "TD")

# Create summaries for ASD
summary.asd <- generate_summary(data.asd, "All data")
summary.asd.c1 <- generate_summary(subset(data.asd, pred.labels == 1), "High (C1)")
summary.asd.c2 <- generate_summary(subset(data.asd, pred.labels == 2), "Medium (C2)")
summary.asd.c3 <- generate_summary(subset(data.asd, pred.labels == 3), "Low (C3)")

# Combine ASD summaries
summary.asd <- rbind(summary.asd, summary.asd.c1, summary.asd.c2, summary.asd.c3) %>%
  select(Cluster, everything())

# Create summaries for TD
summary.td <- generate_summary(data.td, "All data")
summary.td.c1 <- generate_summary(subset(data.td, pred.labels == 1), "High (C1)")
summary.td.c2 <- generate_summary(subset(data.td, pred.labels == 2), "Medium (C2)")

# Combine TD summaries
summary.td <- rbind(summary.td, summary.td.c1, summary.td.c2) %>%
  select(Cluster, everything())

# Print ASD and TD tables separately
print("ASD Summary:")
print(summary.asd)

print("TD Summary:")
print(summary.td)





