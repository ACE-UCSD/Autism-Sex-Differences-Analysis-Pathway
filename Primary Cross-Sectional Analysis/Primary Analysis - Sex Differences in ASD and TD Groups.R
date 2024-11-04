

# """
#
# Filename: Primary Analysis - Sex Differences in ASD and TD Groups.R
# Version: 2024.09.0+375 RStudio 
# Author: Sanaz Nazari
# Date Created: Nov 3, 2024
# Last Modified: Nov 3, 2024
# Description: This script examines sex differences in ASD and TD groups across several tests and subscales.
#
# """

# ------------------------
# Load required libraries
library(dplyr)
library(tidyverse)
library(rstatix)
library(lsr)

######################################################################################################

# ADOS
# -----------------------
# Step 1: Load ADOS Data

data.long <- read.csv("ADOS.long.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Check the number of unique participants
length(unique(data.long$id)) #2138

# Calculate mean age for ASD and TD groups
mean(data.long$age[data.long$diag == "ASD"]) # 28.9
mean(data.long$age[data.long$diag == "TD"])  # 25.7

# ---------------------------
# Step 2: Matching Procedure

# Explanation: The matching process was initially conducted using the MatchIt package.
# However, due to fundamental updates in the package and lack of backward compatibility, 
# the matching code is no longer functional. Below is the code previously used for matching.

# Uncomment if using an older version of MatchIt (not recommended):
# library(MatchIt)
# library(Rglpk)
# data.unq <- data.long %>% group_by(sid) %>% slice_sample(n=1)
# t.test(data.unq$age[data.unq$diag == "ASD"], data.unq$age[data.unq$diag == "TD"])
# data.unq[,10] <- data.unq$diag == "TD"
# colnames(data.unq)[10] <- "group"
# match.it <- matchit(group ~ age, data = data.unq, method = "cardinality", ratio = 1)
# df.match <- match.data(match.it)[1:ncol(data.unq)]

# Current Solution: Instead, the pre-saved matched dataset will be loaded for analysis.
df.match <- read.csv("df.match.ados2.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Mean age post-matching
mean(df.match$age[df.match$diag == "ASD"]) # 26.12
mean(df.match$age[df.match$diag == "TD"])  # 25.52

# T-test on matched data
t.test(df.match$age[df.match$diag == "ASD"], df.match$age[df.match$diag == "TD"])

# Separate ASD and TD datasets for further analysis
d.asd <- subset(df.match, diag == "ASD")
d.td <- subset(df.match, diag == "TD")

# ---------------------------------------------------------
# Step 3: Statistical Tests Based on Normality Assumptions

# Function to check normality using Shapiro-Wilk test
check_normality <- function(data, variable, group_col) {
  groups <- unique(data[[group_col]])
  normality_results <- map(groups, ~ shapiro.test(data[[variable]][data[[group_col]] == .x])$p.value)
  all(map_lgl(normality_results, ~ .x > 0.05))  # Returns TRUE if normality is met in both groups
}

# Function to run the appropriate statistical test (t-test or Kruskal-Wallis) based on normality
run_stat_test <- function(data_asd, data_td, variable, group) {
  # Check normality for both groups
  normal_asd <- check_normality(data_asd, variable, group)
  normal_td <- check_normality(data_td, variable, group)
  
  # If normality is met in both groups, perform t-test
  if (normal_asd && normal_td) {
    asd_test <- t.test(as.formula(paste(variable, "~", group)), data = data_asd)
    td_test <- t.test(as.formula(paste(variable, "~", group)), data = data_td)
  } else {
    # Perform Kruskal-Wallis test if normality is not met in either group
    asd_test <- kruskal.test(as.formula(paste(variable, "~", group)), data = data_asd)
    td_test <- kruskal.test(as.formula(paste(variable, "~", group)), data = data_td)
  }
  
  list(asd = asd_test, td = td_test, test_type = ifelse(normal_asd && normal_td, "t-test", "Kruskal-Wallis"))
}

# Function to calculate Cohen's d
calculate_cohens_d <- function(data, variable, group) {
  cohensD(as.formula(paste(variable, "~", group)), data = data)
}

# Variables of interest
var_list <- c("coso", "rr", "tot")

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(d.asd, d.td, .x, "gender"))

# Collect p-values for each variable (ASD and TD together)
pvalue_list <- map(test_results, ~ c(.x$asd$p.value, .x$td$p.value))

# Apply FDR adjustment within each variable
p.adj_list <- map(pvalue_list, ~ p.adjust(.x, method = "fdr"))

# Extract and combine the test statistics and p-values for each variable
stat_list <- map(test_results, ~ c(.x$asd$statistic, .x$td$statistic))

# Combine statistics, p-values, and adjusted p-values into a table
results <- data.frame(
  Subscale = c("ASD-Social Affect", "TD-Social Affect", "ASD-Restricted & Repetitive Behavior", 
               "TD-Restricted & Repetitive Behavior", "ASD-Overall Total", "TD-Overall Total"),
  Statistic = round(unlist(stat_list), 2),
  P_value = round(unlist(pvalue_list), 3),
  Adj_P_value = round(unlist(p.adj_list), 3),
  Test_Type = rep(map_chr(test_results, ~ .x$test_type), each = 2)
)

# --------------------------------------------------------------------------------
# Step 4: Effect Size Calculation (Cohen's d for t-test, or Kruskal's effect size)

# Add a column for Effect Size, calculate Cohen's d for t-test, or Kruskal-Wallis effect size where applicable
results$Effect_size <- NA

# Calculate effect sizes where applicable
for (i in 1:nrow(results)) {
  if (results$Test_Type[i] == "t-test" && results$Adj_P_value[i] < 0.05) {
    # Calculate Cohen's d if the adjusted p-value is < 0.05 for t-tests
    if (i %% 2 == 0) {  # TD group
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[(i + 1) / 2], "~ gender")), d.td), 2)
    } else {  # ASD group
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[(i + 1) / 2], "~ gender")), d.asd), 2)
    }
  } else if (results$Test_Type[i] == "Kruskal-Wallis" && results$Adj_P_value[i] < 0.05) {
    # Calculate effect size if the adjusted p-value is < 0.05 for Kruskal-Wallis
    if (i %% 2 == 0) {  # TD group
      effect_size <- d.td %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    } else {  # ASD group
      effect_size <- d.asd %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    }
    results$Effect_size[i] <- round(effect_size$effsize, 2)
  }
}

# --------------------
# Display final table 
print(results)

################################################################################################

# Vineland
# ---------------------------
# Step 1: Load Vineland Data

data.long <- read.csv("Vineland.long.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Check the number of unique participants
length(unique(data.long$sid)) #2140

# Calculate mean age for ASD and TD groups
mean(data.long$age[data.long$diag == "ASD"]) # 28.63
mean(data.long$age[data.long$diag == "TD"])  # 25.75

# ---------------------------
# Step 2: Matching Procedure

# Uncomment if using an older version of MatchIt (not recommended):
# library(MatchIt)
# library(Rglpk)
# data.unq <- data.long %>% group_by(sid) %>% slice_sample(n=1)
# t.test(data.unq$age[data.unq$diag == "ASD"], data.unq$age[data.unq$diag == "TD"])
# data.unq[,12] <- data.unq$diag == "TD"
# colnames(data.unq)[12] <- "group"
# match.it <- matchit(group ~ age, data = data.unq, method = "cardinality", ratio = 2, time = 180)
# df.match <- match.data(match.it)[1:ncol(data.unq)]

# Current Solution: Use the pre-saved matched dataset for analysis
df.match <- read.csv("df.match.vine2.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Mean age post-matching
mean(df.match$age[df.match$diag == "ASD"]) # 25.9
mean(df.match$age[df.match$diag == "TD"])  # 25.2

# T-test on matched data to ensure no significant age difference
t.test(df.match$age[df.match$diag == "ASD"], df.match$age[df.match$diag == "TD"])

# Separate ASD and TD datasets for further analysis
d.asd <- subset(df.match, diag == "ASD")
d.td <- subset(df.match, diag == "TD")

# --------------------------------------------------------
# Step 3: Statistical Tests Based on Normality Assumptions

# Variables of interest
var_list <- c("com", "dly", "mtr", "soc", "adap")

# Function to check normality using Shapiro-Wilk test
check_normality <- function(data, variable, group_col) {
  groups <- unique(data[[group_col]])
  normality_results <- map(groups, ~ shapiro.test(data[[variable]][data[[group_col]] == .x])$p.value)
  all(map_lgl(normality_results, ~ .x > 0.05))  # Returns TRUE if normality is met in both groups
}

# Function to run the appropriate statistical test (t-test or Kruskal-Wallis) based on normality
run_stat_test <- function(data_asd, data_td, variable, group) {
  # Check normality for both groups
  normal_asd <- check_normality(data_asd, variable, group)
  normal_td <- check_normality(data_td, variable, group)
  
  # If normality is met in both groups, perform t-test
  if (normal_asd && normal_td) {
    asd_test <- t.test(as.formula(paste(variable, "~", group)), data = data_asd)
    td_test <- t.test(as.formula(paste(variable, "~", group)), data = data_td)
  } else {
    # Perform Kruskal-Wallis test if normality is not met in either group
    asd_test <- kruskal.test(as.formula(paste(variable, "~", group)), data = data_asd)
    td_test <- kruskal.test(as.formula(paste(variable, "~", group)), data = data_td)
  }
  
  list(asd = asd_test, td = td_test, test_type = ifelse(normal_asd && normal_td, "t-test", "Kruskal-Wallis"))
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(d.asd, d.td, .x, "gender"))

# Collect p-values for each variable (ASD and TD together)
pvalue_list <- map(test_results, ~ c(.x$asd$p.value, .x$td$p.value))

# Apply FDR adjustment within each variable (for both ASD and TD groups)
p.adj_list <- map(pvalue_list, ~ p.adjust(.x, method = "fdr"))

# Extract and combine the test statistics and p-values for each variable
stat_list <- map(test_results, ~ c(.x$asd$statistic, .x$td$statistic))

# Combine statistics, p-values, and adjusted p-values into a table
results <- data.frame(
  Subscale = c("ASD-Com", "TD-Com", "ASD-Dly", "TD-Dly", "ASD-Mtr", "TD-Mtr", "ASD-Soc", "TD-Soc", "ASD-Adap", "TD-Adap"),
  Statistic = round(unlist(stat_list), 2),
  P_value = round(unlist(pvalue_list), 3),
  Adj_P_value = round(unlist(p.adj_list), 3),
  Test_Type = rep(map_chr(test_results, ~ .x$test_type), each = 2)
)

# --------------------------------------------------------------------------------
# Step 4: Effect Size Calculation (Cohen's d for t-test, or Kruskal's effect size)

# Add a column for Effect Size, calculate Cohen's d for t-test, or Kruskal-Wallis effect size where applicable
results$Effect_size <- NA

# Calculate effect sizes where applicable
for (i in 1:nrow(results)) {
  if (results$Test_Type[i] == "t-test" && results$Adj_P_value[i] < 0.05) {
    # Calculate Cohen's d if the adjusted p-value is < 0.05 for t-tests
    if (i %% 2 == 0) {  # TD group
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[(i + 1) / 2], "~ gender")), d.td), 2)
    } else {  # ASD group
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[(i + 1) / 2], "~ gender")), d.asd), 2)
    }
  } else if (results$Test_Type[i] == "Kruskal-Wallis" && results$Adj_P_value[i] < 0.05) {
    # Calculate effect size if the adjusted p-value is < 0.05 for Kruskal-Wallis
    if (i %% 2 == 0) {  # TD group
      effect_size <- d.td %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    } else {  # ASD group
      effect_size <- d.asd %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    }
    results$Effect_size[i] <- round(effect_size$effsize, 3)
  }
}

# --------------------
# Display final table
print(results)

#########################################################################################

# Mullen
# ------------------------
# Step 1: Load Mullen Data

data.long <- read.csv("Mullen.long.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Check the number of unique participants
length(unique(data.long$sid)) # 2108 unique participants

mean(data.long$age[data.long$diag=="ASD"]) #28.7
mean(data.long$age[data.long$diag=="TD"]) #23.9

# ---------------------------
# Step 2: Matching Procedure

# Uncomment if using an older version of MatchIt (not recommended):
# library(MatchIt)
# library(Rglpk)
#data.unq <- data.long %>% group_by(sid) %>% slice_sample(n=1)
#t.test(data.unq$age[data.unq$diag=="ASD"], data.unq$age[data.unq$diag=="TD"])
#data.unq[,11] <- data.unq$diag=="TD"
#colnames(data.unq)[11]<- "group"
#match.it <- matchit(group ~ age, data = data.unq, method="cardinality", ratio=1.5, time=300)
#df.match <- match.data(match.it)[1:ncol(data.unq)]

# Current Solution: Use the pre-saved matched dataset for analysis
df.match <- read.csv("df.match.mul2.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Mean age post-matching
mean(df.match$age[df.match$diag == "ASD"]) # ASD mean age
mean(df.match$age[df.match$diag == "TD"])  # TD mean age

# T-test to ensure no significant age difference between ASD and TD
t.test(df.match$age[df.match$diag == "ASD"], df.match$age[df.match$diag == "TD"])

# Separate ASD and TD datasets for analysis
d.asd <- subset(df.match, diag == "ASD")
d.td <- subset(df.match, diag == "TD")

# --------------------------------------------------------
# Step 3: Statistical Tests Based on Normality Assumptions

# Variables of interest
var_list <- c("fm", "vr", "rl", "el")

# Function to check normality using Shapiro-Wilk test for each gender group
check_normality_per_gender <- function(data, variable, group_col) {
  genders <- unique(data[[group_col]])
  
  normality_results <- map_lgl(genders, function(gender) {
    p_value <- shapiro.test(data[[variable]][data[[group_col]] == gender])$p.value
    p_value > 0.05
  })
  
  all(normality_results)  # Returns TRUE if normality is met in both gender groups
}

# Function to run the appropriate statistical test (t-test or Kruskal-Wallis) based on normality
run_stat_test <- function(data_asd, data_td, variable, group) {
  # Check normality for ASD and TD groups separately for each gender
  normal_asd <- check_normality_per_gender(data_asd, variable, group)
  normal_td <- check_normality_per_gender(data_td, variable, group)
  
  # Perform appropriate test for ASD group
  asd_test <- if (normal_asd) {
    t.test(as.formula(paste(variable, "~", group)), data = data_asd)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_asd)
  }
  
  # Perform appropriate test for TD group
  td_test <- if (normal_td) {
    t.test(as.formula(paste(variable, "~", group)), data = data_td)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_td)
  }
  
  list(
    asd = asd_test, 
    td = td_test, 
    asd_test_type = ifelse(normal_asd, "t-test", "Kruskal-Wallis"),
    td_test_type = ifelse(normal_td, "t-test", "Kruskal-Wallis")
  )
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(d.asd, d.td, .x, "gender"))

# Collect p-values for each variable (ASD and TD together)
pvalue_list <- map(test_results, ~ c(.x$asd$p.value, .x$td$p.value))

# Apply FDR adjustment within each variable (for both ASD and TD groups)
p.adj_list <- map(pvalue_list, ~ p.adjust(.x, method = "fdr"))

# Extract and combine the test statistics and p-values for each variable
stat_list <- map(test_results, ~ c(.x$asd$statistic, .x$td$statistic))

# Define Subscale labels
subscale_labels <- c("ASD-FM", "TD-FM", "ASD-VR", "TD-VR", "ASD-RL", "TD-RL", "ASD-EL", "TD-EL")

# Extract ASD and TD test types separately
asd_test_types <- map_chr(test_results, ~ .x$asd_test_type)
td_test_types <- map_chr(test_results, ~ .x$td_test_type)

# Create vectors for ASD and TD test types for the final table
asd_test_type_vec <- rep(NA, length(subscale_labels))
td_test_type_vec <- rep(NA, length(subscale_labels))

# Assign ASD test types to ASD rows and TD test types to TD rows
asd_test_type_vec[grepl("^ASD", subscale_labels)] <- asd_test_types
td_test_type_vec[grepl("^TD", subscale_labels)] <- td_test_types

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(unlist(stat_list), 2),
  P_value = round(unlist(pvalue_list), 3),
  Adj_P_value = round(unlist(p.adj_list), 3),
  ASD_Test_Type = asd_test_type_vec,
  TD_Test_Type = td_test_type_vec
)

# --------------------------------------------------------------------------------------
# Step 4: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size, calculate Cohen's d for t-test, or Kruskal-Wallis effect size where applicable
results$Effect_size <- NA

# Calculate effect sizes where applicable
for (i in 1:nrow(results)) {
  if (results$TD_Test_Type[i] == "t-test" && results$Adj_P_value[i] < 0.05) {
    # Calculate Cohen's d if the adjusted p-value is < 0.05 for t-tests
    if (i %% 2 == 0) {  # TD group
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[(i + 1) / 2], "~ gender")), d.td), 2)
    } else {  # ASD group
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[(i + 1) / 2], "~ gender")), d.asd), 2)
    }
  } else if (results$TD_Test_Type[i] == "Kruskal-Wallis" && results$Adj_P_value[i] < 0.05) {
    # Calculate effect size if the adjusted p-value is < 0.05 for Kruskal-Wallis
    if (i %% 2 == 0) {  # TD group
      effect_size <- d.td %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    } else {  # ASD group
      effect_size <- d.asd %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    }
    results$Effect_size[i] <- round(effect_size$effsize, 2)
  }
}

# --------------------
# Display final table
print(results)


#######################################################################################

# MacArthur-Bates Communicative Development Inventories (CDI)
# words and gestures (wg)
# ------------------------
# Step 1: Load WG Data

data.long <- read.csv("CDI-WG.long.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Check the number of unique participants
length(unique(data.long$sid)) # 554 unique participants

# Mean age per diagnostic group
mean(data.long$age[data.long$diag == "ASD"]) #15.1
mean(data.long$age[data.long$diag == "TD"])  #14.2

# ---------------------------
# Step 2: Matching Procedure

# Uncomment if using an older version of MatchIt (not recommended):
# library(MatchIt)
# library(Rglpk)
# data.unq <- data.long %>% group_by(sid) %>% slice_sample(n=1)
# t.test(data.unq$age[data.unq$diag=="ASD"], data.unq$age[data.unq$diag=="TD"])
#data.unq[,12] <- data.unq$diag=="ASD"
#colnames(data.unq)[12] <- "group"
# match.it <- matchit(group ~ age, data = data.unq, method="cardinality", ratio= 1,time=120)
# df.match <- match.data(match.it)[1:ncol(data.unq)]

# Current Solution: Use the pre-saved matched dataset for analysis
df.match <- read.csv("df.match.wg.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Mean age post-matching
mean(df.match$age[df.match$diag == "ASD"]) # ASD mean age
mean(df.match$age[df.match$diag == "TD"])  # TD mean age

# T-test to ensure no significant age difference post-matching
t.test(df.match$age[df.match$diag == "ASD"], df.match$age[df.match$diag == "TD"])

# Separate ASD and TD datasets for analysis
d.asd <- subset(df.match, diag == "ASD")
d.td <- subset(df.match, diag == "TD")

# --------------------------------------------------------
# Step 3: Statistical Tests Based on Normality Assumptions

# Variables of interest
var_list <- c("w.prod", "w.und", "early.ges", "later.ges", "tot.ges")

# Function to check normality using Shapiro-Wilk test for each gender group
check_normality_per_gender <- function(data, variable, group_col) {
  genders <- unique(data[[group_col]])
  
  normality_results <- map_lgl(genders, function(gender) {
    p_value <- shapiro.test(data[[variable]][data[[group_col]] == gender])$p.value
    p_value > 0.05
  })
  
  all(normality_results)  # Returns TRUE if normality is met in both gender groups
}

# Function to run the appropriate statistical test (t-test or Kruskal-Wallis) based on normality
run_stat_test <- function(data_asd, data_td, variable, group) {
  # Check normality for ASD and TD groups separately for each gender
  normal_asd <- check_normality_per_gender(data_asd, variable, group)
  normal_td <- check_normality_per_gender(data_td, variable, group)
  
  # Perform appropriate test for ASD group
  asd_test <- if (normal_asd) {
    t.test(as.formula(paste(variable, "~", group)), data = data_asd)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_asd)
  }
  
  # Perform appropriate test for TD group
  td_test <- if (normal_td) {
    t.test(as.formula(paste(variable, "~", group)), data = data_td)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_td)
  }
  
  list(
    asd = asd_test, 
    td = td_test, 
    asd_test_type = ifelse(normal_asd, "t-test", "Kruskal-Wallis"),
    td_test_type = ifelse(normal_td, "t-test", "Kruskal-Wallis")
  )
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(d.asd, d.td, .x, "gender"))

# Collect p-values for each variable (ASD and TD together)
pvalue_list <- map(test_results, ~ c(.x$asd$p.value, .x$td$p.value))

# Apply FDR adjustment within each variable (for both ASD and TD groups)
p.adj_list <- map(pvalue_list, ~ p.adjust(.x, method = "fdr"))

# Extract and combine the test statistics and p-values for each variable
stat_list <- map(test_results, ~ c(.x$asd$statistic, .x$td$statistic))

# Define Subscale labels
subscale_labels <- c("ASD-w.prod", "TD-w.prod", "ASD-w.und", "TD-w.und", "ASD-early.ges", "TD-early.ges", "ASD-later.ges", "TD-later.ges", "ASD-tot.ges", "TD-tot.ges")

# Extract ASD and TD test types separately
asd_test_types <- map_chr(test_results, ~ .x$asd_test_type)
td_test_types <- map_chr(test_results, ~ .x$td_test_type)

# Create vectors for ASD and TD test types for the final table
asd_test_type_vec <- rep(NA, length(subscale_labels))
td_test_type_vec <- rep(NA, length(subscale_labels))

# Assign ASD test types to ASD rows and TD test types to TD rows
asd_test_type_vec[grepl("^ASD", subscale_labels)] <- asd_test_types
td_test_type_vec[grepl("^TD", subscale_labels)] <- td_test_types

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(unlist(stat_list), 2),
  P_value = round(unlist(pvalue_list), 3),
  Adj_P_value = round(unlist(p.adj_list), 3),
  ASD_Test_Type = asd_test_type_vec,
  TD_Test_Type = td_test_type_vec
)

# --------------------------------------------------------------------------------------
# Step 4: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size, calculate Cohen's d for t-test, or Kruskal-Wallis effect size where applicable
results$Effect_size <- NA

# Calculate effect sizes where applicable
for (i in 1:nrow(results)) {
  if (results$TD_Test_Type[i] == "t-test" && results$Adj_P_value[i] < 0.05) {
    # Calculate Cohen's d if the adjusted p-value is < 0.05 for t-tests
    if (i %% 2 == 0) {  # TD group
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[(i + 1) / 2], "~ gender")), d.td), 2)
    } else {  # ASD group
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[(i + 1) / 2], "~ gender")), d.asd), 2)
    }
  } else if (results$TD_Test_Type[i] == "Kruskal-Wallis" && results$Adj_P_value[i] < 0.05) {
    # Calculate effect size if the adjusted p-value is < 0.05 for Kruskal-Wallis
    if (i %% 2 == 0) {  # TD group
      effect_size <- d.td %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    } else {  # ASD group
      effect_size <- d.asd %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    }
    results$Effect_size[i] <- round(effect_size$effsize, 2)
  }
}

# --------------------
# Display final table
print(results)

################################################################

# MacArthur-Bates Communicative Development Inventories (CDI)
# words and sentences (ws)
# ------------------------
# Step 1: Load WS Data

data.long <- read.csv("CDI-WS.long.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Check the number of unique participants
length(unique(data.long$sid)) # 892

# Mean age per diagnostic group
mean(data.long$age[data.long$diag == "ASD"]) # 24.6
mean(data.long$age[data.long$diag == "TD"])  # 24.8

# --------------------------------
# Step 2: Checking Age Difference

# Unique participant data
data.unq <- data.long %>% group_by(sid) %>% slice_sample(n = 1)

# T-test to check age difference between ASD and TD
t.test(data.unq$age[data.unq$diag == "ASD"], data.unq$age[data.unq$diag == "TD"])

# No matching required

# Load saved unique data
data.unq <- read.csv("unq.ws2.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# ---------------------------------------------------------
# Step 3: Statistical Tests Based on Normality Assumptions

# Separate ASD and TD datasets for analysis
d.asd <- subset(data.unq, diag == "ASD")
d.td <- subset(data.unq, diag == "TD")

# Variables of interest
var_list <- c("w.prod")

# Function to check normality using Shapiro-Wilk test for each gender group
check_normality_per_gender <- function(data, variable, group_col) {
  genders <- unique(data[[group_col]])
  
  normality_results <- map_lgl(genders, function(gender) {
    p_value <- shapiro.test(data[[variable]][data[[group_col]] == gender])$p.value
    p_value > 0.05
  })
  
  all(normality_results)  # Returns TRUE if normality is met in both gender groups
}

# Function to run the appropriate statistical test (t-test or Kruskal-Wallis) based on normality
run_stat_test <- function(data_asd, data_td, variable, group) {
  # Check normality for ASD and TD groups separately for each gender
  normal_asd <- check_normality_per_gender(data_asd, variable, group)
  normal_td <- check_normality_per_gender(data_td, variable, group)
  
  # Perform appropriate test for ASD group
  asd_test <- if (normal_asd) {
    t.test(as.formula(paste(variable, "~", group)), data = data_asd)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_asd)
  }
  
  # Perform appropriate test for TD group
  td_test <- if (normal_td) {
    t.test(as.formula(paste(variable, "~", group)), data = data_td)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_td)
  }
  
  list(
    asd = asd_test, 
    td = td_test, 
    asd_test_type = ifelse(normal_asd, "t-test", "Kruskal-Wallis"),
    td_test_type = ifelse(normal_td, "t-test", "Kruskal-Wallis")
  )
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(d.asd, d.td, .x, "gender"))

# Collect p-values for each variable (ASD and TD together)
pvalue_list <- map(test_results, ~ c(.x$asd$p.value, .x$td$p.value))

# Apply FDR adjustment within each variable (for both ASD and TD groups)
p.adj_list <- map(pvalue_list, ~ p.adjust(.x, method = "fdr"))

# Extract and combine the test statistics and p-values for each variable
stat_list <- map(test_results, ~ c(.x$asd$statistic, .x$td$statistic))

# Define Subscale labels
subscale_labels <- c("ASD-w.prod", "TD-w.prod")

# Extract ASD and TD test types separately
asd_test_types <- map_chr(test_results, ~ .x$asd_test_type)
td_test_types <- map_chr(test_results, ~ .x$td_test_type)

# Create vectors for ASD and TD test types for the final table
asd_test_type_vec <- rep(NA, length(subscale_labels))
td_test_type_vec <- rep(NA, length(subscale_labels))

# Assign ASD test types to ASD rows and TD test types to TD rows
asd_test_type_vec[grepl("^ASD", subscale_labels)] <- asd_test_types
td_test_type_vec[grepl("^TD", subscale_labels)] <- td_test_types

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(unlist(stat_list), 2),
  P_value = round(unlist(pvalue_list), 3),
  Adj_P_value = round(unlist(p.adj_list), 3),
  ASD_Test_Type = asd_test_type_vec,
  TD_Test_Type = td_test_type_vec
)


# ------------------------------------------------------------
# Step 4: Effect Size Calculation (Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes where applicable
for (i in 1:nrow(results)) {
  if (results$TD_Test_Type[i] == "Kruskal-Wallis" && results$Adj_P_value[i] < 0.05) {
    # Calculate effect size if the adjusted p-value is < 0.05 for Kruskal-Wallis
    if (i %% 2 == 0) {  # TD group
      effect_size <- d.td %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    } else {  # ASD group
      effect_size <- d.asd %>% kruskal_effsize(as.formula(paste(var_list[(i + 1) / 2], "~ gender")))
    }
    results$Effect_size[i] <- round(effect_size$effsize, 3)
  }
}

# --------------------
# Display final table
print(results)

########################################################################################

#CSBS
# ------------------------
# Step 1: Load CSBS Data

data.long <- read.csv("CSBS.long.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Check the number of unique participants
length(unique(data.long$sid)) # 1363 unique participants

# Mean age per diagnostic group
mean(data.long$age[data.long$diag == "ASD"]) # 17.6
mean(data.long$age[data.long$diag == "TD"])  # 16.3

# ---------------------------
# Step 2: Matching Procedure

# Uncomment if using an older version of MatchIt (not recommended):
# library(MatchIt)
# library(Rglpk)
# data.unq <- data.long %>% group_by(sid) %>% slice_sample(n=1)
# t.test(data.unq$age[data.unq$diag=="ASD"], data.unq$age[data.unq$diag=="TD"])
# data.unq[,23] <- data.unq$diag=="TD"
# colnames(data.unq)[23] <- "group"
# match.it <- matchit(group ~ age, data = data.unq, method="cardinality", ratio= 1, time=200)
# df.match <- match.data(match.it)[1:ncol(data.unq)]

# Current Solution: Use the pre-saved matched dataset for analysis
df.match <- read.csv("df.match.csbs2.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Mean age post-matching
mean(df.match$age[df.match$diag == "ASD"]) # ASD mean age
mean(df.match$age[df.match$diag == "TD"])  # TD mean age

# T-test to ensure no significant age difference between ASD and TD
t.test(df.match$age[df.match$diag == "ASD"], df.match$age[df.match$diag == "TD"])

# Separate ASD and TD datasets for analysis
d.asd <- subset(df.match, diag == "ASD")
d.td <- subset(df.match, diag == "TD")

# --------------------------------------------------------
# Step 3: Statistical Tests Based on Normality Assumptions

# Variables of interest
var_list <- c("soc.pc", "exp.pc", "sym.pc", "tot.pc")

# Function to check normality using Shapiro-Wilk test for each gender group
check_normality_per_gender <- function(data, variable, group_col) {
  genders <- unique(data[[group_col]])
  
  normality_results <- map_lgl(genders, function(gender) {
    p_value <- shapiro.test(data[[variable]][data[[group_col]] == gender])$p.value
    p_value > 0.05
  })
  
  all(normality_results)  # Returns TRUE if normality is met in both gender groups
}

# Function to run the appropriate statistical test (t-test or Kruskal-Wallis) based on normality
run_stat_test <- function(data_asd, data_td, variable, group) {
  # Check normality for ASD and TD groups separately for each gender
  normal_asd <- check_normality_per_gender(data_asd, variable, group)
  normal_td <- check_normality_per_gender(data_td, variable, group)
  
  # Perform appropriate test for ASD group
  asd_test <- if (normal_asd) {
    t.test(as.formula(paste(variable, "~", group)), data = data_asd)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_asd)
  }
  
  # Perform appropriate test for TD group
  td_test <- if (normal_td) {
    t.test(as.formula(paste(variable, "~", group)), data = data_td)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_td)
  }
  
  list(
    asd = asd_test, 
    td = td_test, 
    asd_test_type = ifelse(normal_asd, "t-test", "Kruskal-Wallis"),
    td_test_type = ifelse(normal_td, "t-test", "Kruskal-Wallis")
  )
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(d.asd, d.td, .x, "gender"))

# Collect p-values for each variable (ASD and TD together)
pvalue_list <- map(test_results, ~ c(.x$asd$p.value, .x$td$p.value))

# Apply FDR adjustment within each variable (for both ASD and TD groups)
p.adj_list <- map(pvalue_list, ~ p.adjust(.x, method = "fdr"))

# Extract and combine the test statistics and p-values for each variable
stat_list <- map(test_results, ~ c(.x$asd$statistic, .x$td$statistic))

# Define Subscale labels
subscale_labels <- c("ASD-soc.pc", "TD-soc.pc", "ASD-exp.pc", "TD-exp.pc", "ASD-sym.pc", "TD-sym.pc", "ASD-tot.pc", "TD-tot.pc")

# Extract ASD and TD test types separately
asd_test_types <- map_chr(test_results, ~ .x$asd_test_type)
td_test_types <- map_chr(test_results, ~ .x$td_test_type)

# Create vectors for ASD and TD test types for the final table
asd_test_type_vec <- rep(NA, length(subscale_labels))
td_test_type_vec <- rep(NA, length(subscale_labels))

# Assign ASD test types to ASD rows and TD test types to TD rows
asd_test_type_vec[grepl("^ASD", subscale_labels)] <- asd_test_types
td_test_type_vec[grepl("^TD", subscale_labels)] <- td_test_types

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(unlist(stat_list), 2),
  P_value = round(unlist(pvalue_list), 3),
  Adj_P_value = round(unlist(p.adj_list), 3),
  ASD_Test_Type = asd_test_type_vec,
  TD_Test_Type = td_test_type_vec
)

# -------------------------------------------------------------------------------------
# Step 4: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes where applicable
for (i in 1:nrow(results)) {
  if (!is.na(results$Adj_P_value[i]) && results$Adj_P_value[i] < 0.05) {
    if (results$ASD_Test_Type[i] == "Kruskal-Wallis" || results$TD_Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      if (grepl("^TD", results$Subscale[i])) {
        effect_size <- d.td %>% kruskal_effsize(as.formula(paste(var_list[ceiling(i / 2)], "~ gender")))
      } else if (grepl("^ASD", results$Subscale[i])) {
        effect_size <- d.asd %>% kruskal_effsize(as.formula(paste(var_list[ceiling(i / 2)], "~ gender")))
      }
      results$Effect_size[i] <- round(effect_size$effsize, 2)
    } else if (results$ASD_Test_Type[i] == "t-test" || results$TD_Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      if (grepl("^TD", results$Subscale[i])) {
        results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[ceiling(i / 2)], "~ gender")), d.td), 3)
      } else if (grepl("^ASD", results$Subscale[i])) {
        results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[ceiling(i / 2)], "~ gender")), d.asd), 3)
      }
    }
  }
}

# --------------------
# Display final table
print(results)

###############################################################################

# GeoPref
# --------------------------
# Step 1: Load GeoPref Data

data.g <- read.csv("GeoPref.long.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Check the number of unique participants
length(unique(data.g$sid)) # 1482 unique participants

# Mean age per diagnostic group
mean(data.g$et.age[data.g$diag == "ASD"]) # ASD mean age
mean(data.g$et.age[data.g$diag == "TD"])  # TD mean age

# ---------------------------
# Step 2: Matching Procedure

# Uncomment if using an older version of MatchIt (not recommended):
# library(MatchIt)
# library(Rglpk)
# data.unq <- data.g %>% group_by(sid) %>% slice_sample(n=1)
# t.test(data.unq$et.age[data.unq$diag=="ASD"], data.unq$et.age[data.unq$diag=="TD"])
# data.unq[,9] <- data.unq$diag=="TD"
# colnames(data.unq)[9]<- "group"
# match.it <- matchit(group ~  et.age, data = data.unq, method="cardinality", ratio=1.2, time=150)
# df.match <- match.data(match.it)[1:ncol(data.unq)]

# Current Solution: Use the pre-saved matched dataset for analysis
df.match <- read.csv("df.match.geopref.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# Mean age post-matching
mean(df.match$et.age[df.match$diag=="ASD"]) #27
mean(df.match$et.age[df.match$diag=="TD"]) #26

# T-test to ensure no significant age difference between ASD and TD
t.test(df.match$et.age[df.match$diag=="ASD"], df.match$et.age[df.match$diag=="TD"])

# Separate ASD and TD datasets for analysis
d.asd <- subset(df.match, diag == "ASD")
d.td <- subset(df.match, diag == "TD")

# --------------------------------------------------------
# Step 3: Statistical Tests Based on Normality Assumptions

# Variable of interest
var_list <- c("perc.soc")

# Function to check normality using Shapiro-Wilk test for each gender group
check_normality_per_gender <- function(data, variable, group_col) {
  genders <- unique(data[[group_col]])
  
  normality_results <- map_lgl(genders, function(gender) {
    p_value <- shapiro.test(data[[variable]][data[[group_col]] == gender])$p.value
    p_value > 0.05
  })
  
  all(normality_results)  # Returns TRUE if normality is met in both gender groups
}

# Function to run the appropriate statistical test (t-test or Kruskal-Wallis) based on normality
run_stat_test <- function(data_asd, data_td, variable, group) {
  # Check normality for ASD and TD groups separately for each gender
  normal_asd <- check_normality_per_gender(data_asd, variable, group)
  normal_td <- check_normality_per_gender(data_td, variable, group)
  
  # Perform appropriate test for ASD group
  asd_test <- if (normal_asd) {
    t.test(as.formula(paste(variable, "~", group)), data = data_asd)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_asd)
  }
  
  # Perform appropriate test for TD group
  td_test <- if (normal_td) {
    t.test(as.formula(paste(variable, "~", group)), data = data_td)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data_td)
  }
  
  list(
    asd = asd_test, 
    td = td_test, 
    asd_test_type = ifelse(normal_asd, "t-test", "Kruskal-Wallis"),
    td_test_type = ifelse(normal_td, "t-test", "Kruskal-Wallis")
  )
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(d.asd, d.td, .x, "gender"))

# Collect p-values for each variable (ASD and TD together)
pvalue_list <- map(test_results, ~ c(.x$asd$p.value, .x$td$p.value))

# Apply FDR adjustment within each variable (for both ASD and TD groups)
p.adj_list <- map(pvalue_list, ~ p.adjust(.x, method = "fdr"))

# Extract and combine the test statistics and p-values for each variable
stat_list <- map(test_results, ~ c(.x$asd$statistic, .x$td$statistic))

# Define Subscale labels
subscale_labels <- c("ASD-perc.soc", "TD-perc.soc")

# Extract ASD and TD test types separately
asd_test_types <- map_chr(test_results, ~ .x$asd_test_type)
td_test_types <- map_chr(test_results, ~ .x$td_test_type)

# Create vectors for ASD and TD test types for the final table
asd_test_type_vec <- rep(NA, length(subscale_labels))
td_test_type_vec <- rep(NA, length(subscale_labels))

# Assign ASD test types to ASD rows and TD test types to TD rows
asd_test_type_vec[grepl("^ASD", subscale_labels)] <- asd_test_types
td_test_type_vec[grepl("^TD", subscale_labels)] <- td_test_types

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(unlist(stat_list), 2),
  P_value = round(unlist(pvalue_list), 3),
  Adj_P_value = round(unlist(p.adj_list), 3),
  ASD_Test_Type = asd_test_type_vec,
  TD_Test_Type = td_test_type_vec
)

# --------------------------------------------------------------------------------------
# Step 4: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes where applicable
for (i in 1:nrow(results)) {
  if (!is.na(results$Adj_P_value[i]) && results$Adj_P_value[i] < 0.05) {
    if (results$ASD_Test_Type[i] == "Kruskal-Wallis" || results$TD_Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      if (grepl("^TD", results$Subscale[i])) {
        effect_size <- d.td %>% kruskal_effsize(as.formula(paste(var_list[ceiling(i / 2)], "~ gender")))
      } else if (grepl("^ASD", results$Subscale[i])) {
        effect_size <- d.asd %>% kruskal_effsize(as.formula(paste(var_list[ceiling(i / 2)], "~ gender")))
      }
      results$Effect_size[i] <- round(effect_size$effsize, 3)
    } else if (results$ASD_Test_Type[i] == "t-test" || results$TD_Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      if (grepl("^TD", results$Subscale[i])) {
        results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[ceiling(i / 2)], "~ gender")), d.td), 3)
      } else if (grepl("^ASD", results$Subscale[i])) {
        results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[ceiling(i / 2)], "~ gender")), d.asd), 3)
      }
    }
  }
}

# --------------------
# Display final table
print(results)



