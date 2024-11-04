
# """
#
# Filename: Subtype Separation in Train, Test, and External Variables.R
# Version: 2024.09.0+375 RStudio 
# Author: Sanaz Nazari
# Date Created: Nov 3, 2024
# Last Modified: Nov 3, 2024
# Description: This script performs ANOVA and pairwise comparisons on train data, test data and 
#              external variables to examine subtype separation.
#
# """

#------------------------
# Load required libraries
library(ez)
library(dplyr)

##############################################################################################

#Train data

# Read in the dataset
data <- read.csv("train.labels.asd.td2.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Prepare dataset for ANOVA
data.anova <- data.frame(data$coso, data$soc, data$mtr, data$fm, data$com, data$rl, data$el, data$labels.orig, data$gender)
data.anova$id <- seq(1, nrow(data.anova), by = 1)
colnames(data.anova)[8:9] <- c("labels", "gender")

# Rename variables to match descriptions
colnames(data.anova)[1:7] <- c("Social Affect", "Socialization", "Motor Skills", "Fine Motor", 
                               "Communication", "Receptive Language", "Expressive Language")

# Convert factors
data.anova$gender <- as.factor(data.anova$gender)
data.anova$id <- as.factor(data.anova$id)
data.anova$labels <- as.factor(data.anova$labels)

# Initialize result container for ANOVA
res.anova <- data.frame(Domain = character(), Subscale = character(), Effect = character(),
                        F_ratio = numeric(), DFn = numeric(), DFd = numeric(),
                        GES = numeric(), P_value = character(), stringsAsFactors = FALSE)

# Run ANOVA for each variable
for (i in 1:7) {
  data_copy <- data.anova
  colnames(data_copy)[i] <- "var"
  
  anv <- ezANOVA(data = data_copy, dv = var, wid = id, type = 3, between = .(labels), detailed = TRUE, return_aov = TRUE)
  
  # Extract necessary information for the "labels" effect (i.e., "Clusters")
  labels_row <- anv$ANOVA[anv$ANOVA$Effect == "labels", ]
  
  # Extract necessary information for the table
  res_row <- data.frame(
    Domain = ifelse(i <= 2, "Social", ifelse(i <= 4, "Motor", "Language")),  # Adjust domain names based on the variable order
    Subscale = colnames(data.anova)[i],
    Effect = "Clusters",
    F_ratio = round(labels_row$F, 2),
    DFn = labels_row$DFn,
    DFd = labels_row$DFd,
    GES = round(labels_row$ges, 2),
    P_value = ifelse(labels_row$p < 0.001, "< .001*", ifelse(labels_row$p < 0.05, "< .05*", as.character(round(labels_row$p, 3))))
  )
  
  res.anova <- rbind(res.anova, res_row)
}

# View result
res.anova

#---------------------------------------
#pairwise comparisons of three clusters

# Function to perform pairwise comparisons for each variable
perform_pairwise_comparisons <- function(var_name, data) {
  # Initialize result container for each variable
  res_var <- data.frame(Variable = character(), Cluster_1 = numeric(), Cluster_2 = numeric(), P_value = character(), stringsAsFactors = FALSE)
  
  # Perform Shapiro-Wilk test for normality for each cluster
  p_values_normality <- sapply(1:3, function(cluster) shapiro.test(data[[var_name]][data$labels == cluster])$p.value)
  
  # Determine which test to use based on normality
  use_kruskal <- any(p_values_normality < 0.05)
  
  for (c1 in 1:2) {
    for (c2 in (c1 + 1):3) {
      if (use_kruskal) {
        # Use Kruskal-Wallis test if normality is violated in any group
        test <- kruskal.test(list(data[[var_name]][data$labels == c1], data[[var_name]][data$labels == c2]))
      } else {
        # Use t-test if normality is not violated in any group
        test <- t.test(data[[var_name]][data$labels == c1], data[[var_name]][data$labels == c2])
      }
      
      # Store the result
      res_row <- data.frame(Variable = var_name, Cluster_1 = c1, Cluster_2 = c2, P_value = round(test$p.value, 3))
      res_var <- rbind(res_var, res_row)
    }
  }
  
  # Apply FDR correction
  res_var$Adjusted_P_value <- p.adjust(res_var$P_value, method = "fdr")
  
  # Format p-values for significance notation
  res_var$P_value <- ifelse(res_var$Adjusted_P_value < 0.001, "< .001*", ifelse(res_var$Adjusted_P_value < 0.05, "< .05*", as.character(res_var$Adjusted_P_value)))
  
  return(res_var)
}

# Initialize overall result container
res_pairwise <- data.frame(Variable = character(), Cluster_1 = numeric(), Cluster_2 = numeric(), P_value = character(), stringsAsFactors = FALSE)

# Run pairwise comparisons for each variable in the dataset
variables <- c("Social Affect", "Socialization", "Motor Skills", "Fine Motor", "Communication", "Receptive Language", "Expressive Language")

for (var in variables) {
  res_pairwise <- rbind(res_pairwise, perform_pairwise_comparisons(var, data.anova))
}


# View result
res_pairwise

#############################################################################################

# Test data

# Read in the dataset
data <- read.csv("test.labels.asd.td.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Prepare dataset for ANOVA
data.anova <- data.frame(data$coso, data$soc, data$mtr, data$fm, data$com, data$rl, data$el, data$pred.labels, data$gender)
data.anova$id <- seq(1, nrow(data.anova), by = 1)
colnames(data.anova)[8:9] <- c("labels", "gender")

# Rename variables to match descriptions
colnames(data.anova)[1:7] <- c("Social Affect", "Socialization", "Motor Skills", "Fine Motor", 
                               "Communication", "Receptive Language", "Expressive Language")

# Convert factors
data.anova$gender <- as.factor(data.anova$gender)
data.anova$id <- as.factor(data.anova$id)
data.anova$labels <- as.factor(data.anova$labels)

# Initialize result container for ANOVA
res.anova <- data.frame(Domain = character(), Subscale = character(), Effect = character(),
                        F_ratio = numeric(), DFn = numeric(), DFd = numeric(),
                        GES = numeric(), P_value = character(), stringsAsFactors = FALSE)

# Run ANOVA for each variable
for (i in 1:7) {
  data_copy <- data.anova
  colnames(data_copy)[i] <- "var"
  
  anv <- ezANOVA(data = data_copy, dv = var, wid = id, type = 3, between = .(labels), detailed = TRUE, return_aov = TRUE)
  
  # Extract necessary information for the "labels" effect (i.e., "Clusters")
  labels_row <- anv$ANOVA[anv$ANOVA$Effect == "labels", ]
  
  # Extract necessary information for the table
  res_row <- data.frame(
    Domain = ifelse(i <= 2, "Social", ifelse(i <= 4, "Motor", "Language")),  # Adjust domain names based on the variable order
    Subscale = colnames(data.anova)[i],
    Effect = "Clusters",
    F_ratio = round(labels_row$F, 2),
    DFn = labels_row$DFn,
    DFd = labels_row$DFd,
    GES = round(labels_row$ges, 2),
    P_value = ifelse(labels_row$p < 0.001, "< .001*", ifelse(labels_row$p < 0.05, "< .05*", as.character(round(labels_row$p, 3))))
  )
  
  res.anova <- rbind(res.anova, res_row)
}

# View result
res.anova

#---------------------------------------
# Pairwise comparisons of three clusters

# Function to perform pairwise comparisons for each variable
perform_pairwise_comparisons <- function(var_name, data) {
  # Initialize result container for each variable
  res_var <- data.frame(Variable = character(), Cluster_1 = numeric(), Cluster_2 = numeric(), P_value = character(), stringsAsFactors = FALSE)
  
  # Perform Shapiro-Wilk test for normality for each cluster
  p_values_normality <- sapply(1:3, function(cluster) shapiro.test(data[[var_name]][data$labels == cluster])$p.value)
  
  # Determine which test to use based on normality
  use_kruskal <- any(p_values_normality < 0.05)
  
  for (c1 in 1:2) {
    for (c2 in (c1 + 1):3) {
      if (use_kruskal) {
        # Use Kruskal-Wallis test if normality is violated in any group
        test <- kruskal.test(list(data[[var_name]][data$labels == c1], data[[var_name]][data$labels == c2]))
      } else {
        # Use t-test if normality is not violated in any group
        test <- t.test(data[[var_name]][data$labels == c1], data[[var_name]][data$labels == c2])
      }
      
      # Store the result
      res_row <- data.frame(Variable = var_name, Cluster_1 = c1, Cluster_2 = c2, P_value = round(test$p.value, 3))
      res_var <- rbind(res_var, res_row)
    }
  }
  
  # Apply FDR correction
  res_var$Adjusted_P_value <- p.adjust(res_var$P_value, method = "fdr")
  
  # Format p-values for significance notation
  res_var$P_value <- ifelse(res_var$Adjusted_P_value < 0.001, "< .001*", ifelse(res_var$Adjusted_P_value < 0.05, "< .05*", as.character(res_var$Adjusted_P_value)))
  
  return(res_var)
}

# Initialize overall result container
res_pairwise <- data.frame(Variable = character(), Cluster_1 = numeric(), Cluster_2 = numeric(), P_value = character(), stringsAsFactors = FALSE)

# Run pairwise comparisons for each variable in the dataset
variables <- c("Social Affect", "Socialization", "Motor Skills", "Fine Motor", "Communication", "Receptive Language", "Expressive Language")

for (var in variables) {
  res_pairwise <- rbind(res_pairwise, perform_pairwise_comparisons(var, data.anova))
}

# View result
res_pairwise


#####################################################################################################

# External Variables

# Read in the dataset
data <- read.csv("train.labels.asd.td2-new-ext.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Prepare dataset for ANOVA
data.anova <- data.frame(data$w.prod, data$w.und, data$early.ges, data$later.ges, data$tot.ges, data$w.prod.ws, 
                         data$percent.fixation.social, data$rr,  data$labels.orig, data$gender)
data.anova$id <- seq(1, nrow(data.anova), by = 1)
colnames(data.anova)[9:10] <- c("labels", "gender")

# Rename variables to match descriptions
colnames(data.anova)[1:8] <- c("WG-Words Produced", "WG-Words Understood", "WG-Early Gestures", "WG-Later Gestures", 
                               "WG-Total Gestures", "WS-Words Produced", "Geopref-Social % Fixation",
                               "ADOS-Restricted and Repetitive Behavior")

# Convert factors
data.anova$gender <- as.factor(data.anova$gender)
data.anova$id <- as.factor(data.anova$id)
data.anova$labels <- as.factor(data.anova$labels)

# Initialize result container for ANOVA
res.anova <- data.frame(Subscale = character(), Effect = character(),
                        F_ratio = numeric(), DFn = numeric(), DFd = numeric(),
                        GES = numeric(), P_value = character(), stringsAsFactors = FALSE)

# Run ANOVA for each variable
for (i in 1:8) {
  data_copy <- data.anova[!is.na(data.anova[[i]]), ]  # Remove rows with NA for the current variable
  colnames(data_copy)[i] <- "var"
  
  anv <- ezANOVA(data = data_copy, dv = var, wid = id, type = 3, between = .(labels), detailed = TRUE, return_aov = TRUE)
  
  # Extract necessary information for the "labels" effect (i.e., "Clusters")
  labels_row <- anv$ANOVA[anv$ANOVA$Effect == "labels", ]
  
  # Extract necessary information for the table
  res_row <- data.frame(
    Subscale = colnames(data.anova)[i],
    Effect = "Clusters",
    F_ratio = round(labels_row$F, 2),
    DFn = labels_row$DFn,
    DFd = labels_row$DFd,
    GES = round(labels_row$ges, 2),
    P_value = ifelse(labels_row$p < 0.001, "< .001*", ifelse(labels_row$p < 0.05, "< .05*", as.character(round(labels_row$p, 3))))
  )
  
  res.anova <- rbind(res.anova, res_row)
}

# View result
res.anova

#---------------------------------------
# Pairwise comparisons of three clusters

# Function to perform pairwise comparisons for each variable
perform_pairwise_comparisons <- function(var_name, data) {
  # Remove rows with NA for the current variable
  data <- data[!is.na(data[[var_name]]), ]
  
  # Initialize result container for each variable
  res_var <- data.frame(Variable = character(), Cluster_1 = numeric(), Cluster_2 = numeric(), P_value = character(), Test_Type = character(), stringsAsFactors = FALSE)
  
  # Perform Shapiro-Wilk test for normality for each cluster
  p_values_normality <- sapply(1:3, function(cluster) shapiro.test(data[[var_name]][data$labels == cluster])$p.value)
  
  for (c1 in 1:2) {
    for (c2 in (c1 + 1):3) {
      if (p_values_normality[c1] < 0.05 || p_values_normality[c2] < 0.05) {
        # Use Kruskal-Wallis test if normality is violated in any group
        test <- kruskal.test(list(data[[var_name]][data$labels == c1], data[[var_name]][data$labels == c2]))
        test_type <- "Kruskal-Wallis"
      } else {
        # Use t-test if normality is not violated in both groups
        test <- t.test(data[[var_name]][data$labels == c1], data[[var_name]][data$labels == c2])
        test_type <- "t-test"
      }
      
      # Store the result
      res_row <- data.frame(Variable = var_name, Cluster_1 = c1, Cluster_2 = c2, P_value = round(test$p.value, 3), Test_Type = test_type)
      res_var <- rbind(res_var, res_row)
    }
  }
  
  # Apply FDR correction
  res_var$Adjusted_P_value <- round(p.adjust(res_var$P_value, method = "fdr"),3)
  
  # Format p-values for significance notation
  res_var$P_value <- ifelse(res_var$Adjusted_P_value < 0.001, "< .001*", ifelse(res_var$Adjusted_P_value < 0.05, "< .05*", as.character(res_var$Adjusted_P_value)))
  
  return(res_var)
}

# Initialize overall result container
res_pairwise <- data.frame(Variable = character(), Cluster_1 = numeric(), Cluster_2 = numeric(), P_value = character(), Test_Type = character(), stringsAsFactors = FALSE)

# Run pairwise comparisons for each variable in the dataset
variables <- c("WG-Words Produced", "WG-Words Understood", "WG-Early Gestures", "WG-Later Gestures", 
               "WG-Total Gestures", "WS-Words Produced", "Geopref-Social % Fixation","ADOS-Restricted and Repetitive Behavior")

for (var in variables) {
  res_pairwise <- rbind(res_pairwise, perform_pairwise_comparisons(var, data.anova))
}

# View result
res_pairwise


