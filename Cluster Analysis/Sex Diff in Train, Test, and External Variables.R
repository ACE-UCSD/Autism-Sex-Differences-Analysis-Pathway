
# """
#
# Filename: Sex Diff in Train, Test, and External Variables.R
# Version: 2024.09.0+375 RStudio 
# Author: Sanaz Nazari
# Date Created: Nov 3, 2024
# Last Modified: Nov 3, 2024
# Description: This script examines sex differences in train data, test data and external variables
#              in each cluster.
#
# """

###########################################################################################################

# Train data

# Load dataset
data <- read.csv("train.labels.asd.td2.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Define variables of interest
variables <- c("coso", "soc", "mtr", "fm", "com", "rl", "el")
variable_names <- c("Social Affect", "Socialization", "Motor Skills", "Fine Motor", "Communication", 
                    "Receptive Language", "Expressive Language")

# Initialize an empty data frame to store the results
results <- data.frame(Group = character(), Subscale = character(), Cluster = integer(), 
                      Statistic = numeric(), P_value = numeric(), Effect_size = character(), stringsAsFactors = FALSE)

# Function to calculate effect size (Cohen's d for t-test or eta squared for Kruskal-Wallis)
calculate_effect_size <- function(test_result, test_type, data_subset, var) {
  if (test_type == "t") {
    group_means <- tapply(data_subset[[var]], data_subset$gender, mean)
    group_sds <- tapply(data_subset[[var]], data_subset$gender, sd)
    n_f <- sum(data_subset$gender == "F")
    n_m <- sum(data_subset$gender == "M")
    pooled_sd <- sqrt(((n_f - 1) * group_sds["F"]^2 + (n_m - 1) * group_sds["M"]^2) / (n_f + n_m - 2))
    effect_size <- (group_means["F"] - group_means["M"]) / pooled_sd
    return(round(effect_size, 2))
  } else if (test_type == "kruskal") {
    n <- length(data_subset[[var]])
    h_stat <- test_result$statistic
    effect_size <- h_stat / (n - 1)
    return(round(effect_size, 2))
  }
  return("-")
}

# Function to test normality and perform the appropriate test
test_gender_effect <- function(data, cluster_num, diag_group) {
  data_subset <- subset(data, labels.orig == cluster_num & diag == diag_group)
  
  for (i in seq_along(variables)) {
    var <- variables[i]
    subscale <- variable_names[i]
    
    # Normality test
    shapiro_f <- shapiro.test(data_subset[[var]][data_subset$gender == "F"])
    shapiro_m <- shapiro.test(data_subset[[var]][data_subset$gender == "M"])
    
    # Choose the appropriate test based on normality results
    if (shapiro_f$p.value > 0.05 && shapiro_m$p.value > 0.05) {
      t_result <- t.test(data_subset[[var]] ~ data_subset$gender)
      statistic <- t_result$statistic
      p_value <- t_result$p.value
      effect_size <- ifelse(p_value < 0.05, calculate_effect_size(t_result, "t", data_subset, var), "-")
    } else {
      k_result <- kruskal.test(data_subset[[var]] ~ data_subset$gender)
      statistic <- k_result$statistic
      p_value <- k_result$p.value
      effect_size <- ifelse(p_value < 0.05, calculate_effect_size(k_result, "kruskal", data_subset, var), "-")
    }
    
    # Append the results to the summary data frame
    results <<- rbind(results, data.frame(Group = diag_group, Subscale = subscale, Cluster = cluster_num, 
                                          Statistic = round(statistic, 2), P_value = round(p_value, 3), 
                                          Effect_size = effect_size, stringsAsFactors = FALSE))
  }
}

# Analyze sex effects for each cluster and diagnosis group
for (cluster in 1:3) {
  for (group in c("ASD", "TD")) {
    test_gender_effect(data, cluster, group)
  }
}

# Display the summary table
print(results)


################################################################################################


# Test Data

# Load dataset (Test)
data_test <- read.csv("test.labels.asd.td.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Define variables of interest
variables <- c("coso", "soc", "mtr", "fm", "com", "rl", "el")
variable_names <- c("Social Affect", "Socialization", "Motor Skills", "Fine Motor", "Communication", 
                    "Receptive Language", "Expressive Language")

# Initialize an empty data frame to store the results for test data
results_test <- data.frame(Group = character(), Subscale = character(), Cluster = integer(), 
                           Statistic = numeric(), P_value = numeric(), Effect_size = character(), stringsAsFactors = FALSE)

# Function to calculate effect size (Cohen's d for t-test or eta squared for Kruskal-Wallis)
calculate_effect_size <- function(test_result, test_type, data_subset, var) {
  if (test_type == "t") {
    group_means <- tapply(data_subset[[var]], data_subset$gender, mean)
    group_sds <- tapply(data_subset[[var]], data_subset$gender, sd)
    n_f <- sum(data_subset$gender == "F")
    n_m <- sum(data_subset$gender == "M")
    pooled_sd <- sqrt(((n_f - 1) * group_sds["F"]^2 + (n_m - 1) * group_sds["M"]^2) / (n_f + n_m - 2))
    effect_size <- (group_means["F"] - group_means["M"]) / pooled_sd
    return(round(effect_size, 2))
  } else if (test_type == "kruskal") {
    n <- length(data_subset[[var]])
    h_stat <- test_result$statistic
    effect_size <- h_stat / (n - 1)
    return(round(effect_size, 2))
  }
  return("-")
}

# Function to test normality and perform the appropriate test
test_gender_effect_test <- function(data, cluster_num, diag_group) {
  data_subset <- subset(data, pred.labels == cluster_num & diag == diag_group)
  
  # Check if there are any data points in the subset
  if (nrow(data_subset) == 0) {
    cat("\nNo data available for Cluster", cluster_num, "and Group", diag_group, "\n")
    return()
  }
  
  for (i in seq_along(variables)) {
    var <- variables[i]
    subscale <- variable_names[i]
    
    # Normality test
    shapiro_f <- tryCatch(shapiro.test(data_subset[[var]][data_subset$gender == "F"]), error = function(e) NULL)
    shapiro_m <- tryCatch(shapiro.test(data_subset[[var]][data_subset$gender == "M"]), error = function(e) NULL)
    
    # Choose the appropriate test based on normality results
    if (!is.null(shapiro_f) && !is.null(shapiro_m) && shapiro_f$p.value > 0.05 && shapiro_m$p.value > 0.05) {
      t_result <- t.test(data_subset[[var]] ~ data_subset$gender)
      statistic <- t_result$statistic
      p_value <- t_result$p.value
      effect_size <- ifelse(p_value < 0.05, calculate_effect_size(t_result, "t", data_subset, var), "-")
    } else {
      k_result <- kruskal.test(data_subset[[var]] ~ data_subset$gender)
      statistic <- k_result$statistic
      p_value <- k_result$p.value
      effect_size <- ifelse(p_value < 0.05, calculate_effect_size(k_result, "kruskal", data_subset, var), "-")
    }
    
    # Append the results to the summary data frame
    results_test <<- rbind(results_test, data.frame(Group = diag_group, Subscale = subscale, Cluster = cluster_num, 
                                                    Statistic = round(statistic, 2), P_value = round(p_value, 3), 
                                                    Effect_size = effect_size, stringsAsFactors = FALSE))
  }
}

# Analyze sex effects for each cluster and diagnosis group
for (cluster in 1:3) {
  for (group in c("ASD", "TD")) {
    test_gender_effect_test(data_test, cluster, group)
  }
}

# Display the summary table for test data
cat("\nTest Data Results:\n")
print(results_test)


#####################################################################################################

# External Variables

# Load dataset
data_ext <- read.csv("train.labels.asd.td2-new-ext.csv", header = TRUE, na.strings = c("NULL", "NA", ""))

# Define external variables of interest
variables_ext <- c("rr","w.prod", "w.und", "early.ges", "later.ges", "tot.ges", "w.prod.ws", "percent.fixation.social")
variable_names_ext <- c("ADOS-Restricted and Repetitive Behavior", "WG-Words Produced", "WG-Words Understood", 
                        "WG-Early Gestures", "WG-Later Gestures", "WG-Total Gestures", "WS-Words Produced", "Geopref-Social % Fixation")

# Initialize an empty data frame to store the results for external variables
results_ext <- data.frame(Group = character(), Subscale = character(), Cluster = integer(), 
                          Statistic = numeric(), P_value = numeric(), Effect_size = character(), stringsAsFactors = FALSE)

# Function to calculate effect size (Cohen's d for t-test or eta squared for Kruskal-Wallis)
calculate_effect_size_ext <- function(test_result, test_type, data_subset, var) {
  if (test_type == "t") {
    group_means <- tapply(data_subset[[var]], data_subset$gender, mean, na.rm = TRUE)
    group_sds <- tapply(data_subset[[var]], data_subset$gender, sd, na.rm = TRUE)
    n_f <- sum(data_subset$gender == "F", na.rm = TRUE)
    n_m <- sum(data_subset$gender == "M", na.rm = TRUE)
    if (is.na(group_sds["F"]) || is.na(group_sds["M"]) || n_f < 2 || n_m < 2) {
      return("-")
    }
    pooled_sd <- sqrt(((n_f - 1) * group_sds["F"]^2 + (n_m - 1) * group_sds["M"]^2) / (n_f + n_m - 2))
    effect_size <- (group_means["F"] - group_means["M"]) / pooled_sd
    return(round(effect_size, 2))
  } else if (test_type == "kruskal") {
    n <- length(data_subset[[var]])
    h_stat <- test_result$statistic
    effect_size <- h_stat / (n - 1)
    return(round(effect_size, 2))
  }
  return("-")
}

# Function to test normality and perform the appropriate test for external variables
test_gender_effect_ext <- function(data, cluster_num, diag_group) {
  data_subset <- subset(data, labels.orig == cluster_num & diag == diag_group)
  
  # Check if there are any data points in the subset
  if (nrow(data_subset) == 0) {
    cat("\nNo data available for Cluster", cluster_num, "and Group", diag_group, "\n")
    return()
  }
  
  for (i in seq_along(variables_ext)) {
    var <- variables_ext[i]
    subscale <- variable_names_ext[i]
    
    # Normality test
    shapiro_f <- tryCatch(shapiro.test(data_subset[[var]][data_subset$gender == "F"]), error = function(e) NULL)
    shapiro_m <- tryCatch(shapiro.test(data_subset[[var]][data_subset$gender == "M"]), error = function(e) NULL)
    
    # Choose the appropriate test based on normality results
    if (!is.null(shapiro_f) && !is.null(shapiro_m) && shapiro_f$p.value > 0.05 && shapiro_m$p.value > 0.05) {
      t_result <- t.test(data_subset[[var]] ~ data_subset$gender)
      statistic <- t_result$statistic
      p_value <- t_result$p.value
      effect_size <- ifelse(p_value < 0.05, calculate_effect_size_ext(t_result, "t", data_subset, var), "-")
    } else {
      k_result <- kruskal.test(data_subset[[var]] ~ data_subset$gender)
      statistic <- k_result$statistic
      p_value <- k_result$p.value
      effect_size <- ifelse(p_value < 0.05, calculate_effect_size_ext(k_result, "kruskal", data_subset, var), "-")
    }
    
    # Append the results to the summary data frame
    results_ext <<- rbind(results_ext, data.frame(Group = diag_group, Subscale = subscale, 
                                                  Cluster = cluster_num, Statistic = round(statistic, 2), 
                                                  P_value = round(p_value, 3), Effect_size = effect_size, stringsAsFactors = FALSE))
  }
}

# Analyze sex effects for each cluster and diagnosis group
for (cluster in 1:3) {
  for (group in c("ASD", "TD")) {
    test_gender_effect_ext(data_ext, cluster, group)
  }
}

# Display the summary table for external variables
cat("\nTrain Data with External Variables Results:\n")
print(results_ext)









