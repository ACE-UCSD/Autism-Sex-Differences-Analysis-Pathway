
# """
#
# Filename: NMI.R
# Version: 2024.09.0+375 RStudio 
# Author: Sanaz Nazari
# Date Created: Nov 3, 2024
# Last Modified: Nov 3, 2024
# Description: This script performs Mormalized Mutual Information on random removals.
#
# """

#--------------------------
# Load necessary libraries
library(qpcR)
library(SNFtool)
library(dplyr)
library(purrr)

# Hyperparameters
k <- 20       # number of neighbors, usually (10~30)
alpha <- 0.5  # hyperparameter, usually (0.3~0.8)
t <- 20       # Number of Iterations, usually (10~20)
c <- 3        # number of clusters

# Set random seed for reproducibility
set.seed(123)

# Load data
data.orig <- read.csv("train.labels.asd.td2.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Function to standardize and set dimnames for data layers
prepare_data_layer <- function(data, columns, dim_names) {
  layer <- data[, columns, drop = FALSE]
  dimnames(layer) <- list(1:nrow(layer), dim_names)
  layer <- standardNormalization(layer)
  return(layer)
}

# Start timing execution
start <- proc.time()

# Initialize variables for storing results
nmi_results <- list()
random_percents <- c(0.95, 0.90, 0.80, 0.70, 0.60, 0.50)

# Main loop for different percentages
nmi_results <- map(random_percents, function(r) {
  map_dbl(1:100, function(i) {
    # Randomly sample a percentage of the original dataset
    n <- floor(r * nrow(data.orig))
    ind <- sample(seq_len(nrow(data.orig)), size = n)
    data.rr <- data.orig[ind, ]
    rownames(data.rr) <- NULL
    
    # Set layers using the prepare_data_layer function
    data1 <- prepare_data_layer(data.rr, c("coso", "soc"), c("coso", "soc"))
    data2 <- prepare_data_layer(data.rr, c("mtr", "fm"), c("mtr", "fm"))
    data3 <- prepare_data_layer(data.rr, c("com", "rl", "el"), c("com", "rl", "el"))
    
    # Combine layers into a list
    data.tr.list <- list(data1, data2, data3)
    
    # Calculate the pair-wise distance and similarity graphs
    dist_list <- lapply(data.tr.list, function(layer) (dist2(as.matrix(layer), as.matrix(layer)))^(1/2))
    
    w_list <- lapply(dist_list, function(dist_matrix) affinityMatrix(dist_matrix, k, alpha))
    
    # Fuse all the graphs using similarity network fusion (SNF)
    w <- SNF(w_list, k, t)
    
    # Get cluster labels by spectral clustering
    labels.rr <- spectralClustering(w, c)
    
    # Add labels to the dataset
    data.rr <- data.frame(data.rr, labels.rr)
    
    # Calculate NMI
    nmi_value <- calNMI(data.rr$labels.orig, data.rr$labels.rr)
    
    # Print progress
    print(paste("i=", i, ", r=", r, ", nmi=", nmi_value))
    
    return(nmi_value)
  })
})

# Combine NMI results into a data frame
nmi <- do.call(qpcR:::cbind.na, nmi_results)

# End timing execution
end <- proc.time()
print(end - start)

# Print NMI results
nmi


