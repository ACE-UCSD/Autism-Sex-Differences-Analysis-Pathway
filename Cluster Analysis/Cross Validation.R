
# """
#
# Filename: Cross Validation.R
# Version: 2024.09.0+375 RStudio 
# Author: Sanaz Nazari and Javad Zahiri
# Date Created: Nov 3, 2024
# Last Modified: Nov 3, 2024
# Description: This script performs 5-fold cross validation.
#
# """

#-------------------------
# Load required libraries
library(caret)
library(SNFtool)

# Load data
SNF_cluster <- read.csv("train.labels.asd.td2.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Train set preparation
# Repeated 5-fold cross-validation analysis
lay1 <- c("coso", "soc")
lay2 <- c("mtr", "fm")
lay3 <- c("com", "rl", "el")

cluster.index.vctr <- SNF_cluster$labels.orig
layers.list.outcome <- list(SNF_cluster[, lay1], SNF_cluster[, lay2], SNF_cluster[, lay3])

#-----------------------------------------------------------
# Function to run repeated k-fold cross-validation analysis
SNF.repeated.cross.validation <- function(no.of.repeat, fold.No, layers.list, cluster.index.vctr) {
  # INPUT:
  # no.of.repeat = number of iterations
  # fold.No = number of folds, i.e., k
  # layers.list = list of data frames for SNF
  # cluster.index.vctr = cluster labels
  # OUTPUT:
  # accuracy of prediction
  
  accuracy.vctr.repeated.k.fold.CV <- numeric(no.of.repeat * fold.No)
  
  for (i in seq_len(no.of.repeat)) {
    cat("================================================\n")
    cat(paste0("Iteration: ", i, "\n"))
    
    accuracy.vctr.k.fold.CV <- SNF.cross.validation(fold.No, layers.list, cluster.index.vctr)
    
    for (j in seq_along(accuracy.vctr.k.fold.CV)) {
      cat(paste0("Iteration ", i, " - Fold ", j, " Accuracy: ", accuracy.vctr.k.fold.CV[j], "\n"))
    }
    
    accuracy.vctr.repeated.k.fold.CV[((i - 1) * fold.No + 1):(i * fold.No)] <- accuracy.vctr.k.fold.CV
  }
  
  return(accuracy.vctr.repeated.k.fold.CV)
}

# Function to run k-fold cross-validation
SNF.cross.validation <- function(fold.No, layers.list, cluster.index.vctr) {
  flds <- createFolds(seq_along(cluster.index.vctr), k = fold.No, list = TRUE, returnTrain = FALSE)
  accuracy.vctr <- numeric(fold.No)
  
  for (i in seq_len(fold.No)) {
    cat(paste0("Fold ", i, "\n"))
    
    # Prediction of the subtype of new subjects
    testSampleIndexVctr <- flds[[i]]
    
    train <- lapply(layers.list, function(x) x[-testSampleIndexVctr, ]) # Training set
    test <- lapply(layers.list, function(x) x[testSampleIndexVctr, ]) # Test set
    groups <- cluster.index.vctr[-testSampleIndexVctr] # Labels of the training data
    
    # Set parameters for SNF
    K <- 20
    alpha <- 0.5
    t <- 20
    method <- 1 # 1 means using label propagation
    
    # Apply prediction function
    newcluster.index.vctr.from.SNF <- groupPredict(train, test, groups, K, alpha, t, method = method)
    
    # Compute prediction accuracy
    accuracy.for.this.fold <- sum(cluster.index.vctr[testSampleIndexVctr] == 
                                    newcluster.index.vctr.from.SNF[-seq_len(length(groups))]) /
      length(testSampleIndexVctr)
    
    accuracy.vctr[i] <- accuracy.for.this.fold
    cat(paste0("Accuracy for Fold ", i, ": ", accuracy.for.this.fold, "\n"))
  }
  
  return(accuracy.vctr)
}

# Run repeated cross-validation and measure time
start <- proc.time()
accuracy.vctr.repeated.k.fold.CV <- SNF.repeated.cross.validation(
  no.of.repeat = 10,
  fold.No = 5,
  layers.list = layers.list.outcome,
  cluster.index.vctr = cluster.index.vctr
)
end <- proc.time()

time_taken <- end - start
cat("Time taken:\n")
print(time_taken)

# Output mean accuracy and save results
mean_accuracy <- mean(accuracy.vctr.repeated.k.fold.CV)
cat("Mean Accuracy:", mean_accuracy, "\n")

#write.table(accuracy.vctr.repeated.k.fold.CV, "accuracy-CV.csv", col.names = FALSE, row.names = FALSE, quote = FALSE)
