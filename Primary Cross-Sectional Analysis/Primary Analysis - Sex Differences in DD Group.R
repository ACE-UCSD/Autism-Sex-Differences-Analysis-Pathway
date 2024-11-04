

# """
#
# Filename: Primary Analysis - Sex Differences in DD Group.R
# Version: 2024.09.0+375 RStudio 
# Author: Sanaz Nazari
# Date Created: Nov 3, 2024
# Last Modified: Nov 3, 2024
# Description: This script examines sex differences in DD group across several tests and subscales.
#
# """

# ------------------------
# Load required libraries
library(tidyverse)
library(rstatix)

###############################################################################################################

# ADOS
# ------------------------
# Step 1: Load ADoS Data
ados <- read.csv("ados.dd.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# --------------------------------------------------------
# Step 2: Statistical Tests Based on Normality Assumptions

# Variables of interest
var_list <- c("coso", "rr", "tot")

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
run_stat_test <- function(data, variable, group) {
  # Check normality for each gender
  normal <- check_normality_per_gender(data, variable, group)
  
  # Perform appropriate test
  if (normal) {
    t.test(as.formula(paste(variable, "~", group)), data = data)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  }
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(ados, .x, "gender"))

# Collect p-values for each variable
pvalue_list <- map_dbl(test_results, ~ .x$p.value)

# Extract and combine the test statistics and p-values for each variable
stat_list <- map_dbl(test_results, ~ .x$statistic)

# Define Subscale labels
subscale_labels <- c("coso", "rr", "tot")

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(stat_list, 2),
  P_value = round(pvalue_list, 3),
  Test_Type = map_chr(test_results, ~ ifelse(class(.x) == "htest" && .x$method == "Welch Two Sample t-test", "t-test", "Kruskal-Wallis"))
)

# -------------------------------------------------------------------------------------
# Step 3: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes for variables where p < 0.05
for (i in 1:nrow(results)) {
  if (results$P_value[i] < 0.05) {
    if (results$Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      effect_size <- ados %>% kruskal_effsize(as.formula(paste(var_list[i], "~ gender")))
      results$Effect_size[i] <- round(effect_size$effsize, 3)
    } else if (results$Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[i], "~ gender")), ados), 3)
    }
  }
}

# --------------------
# Display final table
print(results)

####################################################################################################

# Vineland
# ---------------------------
# Step 1: Load Vineland Data

vine <- read.csv("vineland.dd.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# --------------------------------------------------------
# Step 2: Statistical Tests Based on Normality Assumptions

# Variables of interest
var_list <- c("com", "dly", "mtr", "soc", "adap")

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
run_stat_test <- function(data, variable, group) {
  # Check normality for each gender
  normal <- check_normality_per_gender(data, variable, group)
  
  # Perform appropriate test
  if (normal) {
    t.test(as.formula(paste(variable, "~", group)), data = data)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  }
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(vine, .x, "gender"))

# Collect p-values for each variable
pvalue_list <- map_dbl(test_results, ~ .x$p.value)

# Extract and combine the test statistics and p-values for each variable
stat_list <- map_dbl(test_results, ~ .x$statistic)

# Define Subscale labels
subscale_labels <- c("com", "dly", "mtr", "soc", "adap")

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(stat_list, 2),
  P_value = round(pvalue_list, 3),
  Test_Type = map_chr(test_results, ~ ifelse(class(.x) == "htest" && .x$method == "Welch Two Sample t-test", "t-test", "Kruskal-Wallis"))
)

# --------------------------------------------------------------------------------------
# Step 3: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes for variables where p < 0.05
for (i in 1:nrow(results)) {
  if (results$P_value[i] < 0.05) {
    if (results$Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      effect_size <- vine %>% kruskal_effsize(as.formula(paste(var_list[i], "~ gender")))
      results$Effect_size[i] <- round(effect_size$effsize, 3)
    } else if (results$Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[i], "~ gender")), vine), 3)
    }
  }
}

# --------------------
# Display final table
print(results)

#############################################################################################

# Mullen
# ------------------------
# Step 1: Load Mullen Data

mul <- read.csv("mullen.dd.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# --------------------------------------------------------
# Step 2: Statistical Tests Based on Normality Assumptions

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
run_stat_test <- function(data, variable, group) {
  # Check normality for each gender
  normal <- check_normality_per_gender(data, variable, group)
  
  # Perform appropriate test
  if (normal) {
    t.test(as.formula(paste(variable, "~", group)), data = data)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  }
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(mul, .x, "gender"))

# Collect p-values for each variable
pvalue_list <- map_dbl(test_results, ~ .x$p.value)

# Extract and combine the test statistics and p-values for each variable
stat_list <- map_dbl(test_results, ~ .x$statistic)

# Define Subscale labels
subscale_labels <- c("fm", "vr", "rl", "el")

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(stat_list, 2),
  P_value = round(pvalue_list, 3),
  Test_Type = map_chr(test_results, ~ ifelse(class(.x) == "htest" && .x$method == "Welch Two Sample t-test", "t-test", "Kruskal-Wallis"))
)

# -------------------------------------------------------------------------------------
# Step 3: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes for variables where p < 0.05
for (i in 1:nrow(results)) {
  if (results$P_value[i] < 0.05) {
    if (results$Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      effect_size <- mul %>% kruskal_effsize(as.formula(paste(var_list[i], "~ gender")))
      results$Effect_size[i] <- round(effect_size$effsize, 3)
    } else if (results$Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[i], "~ gender")), mul), 3)
    }
  }
}

# --------------------
# Display final table
print(results)

###################################################################################

# MacArthur-Bates Communicative Development Inventories (CDI)
# words and gestures (wg)
# ------------------------
# Step 1: Load WG Data

wg <- read.csv("wg.dd.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# --------------------------------------------------------
# Step 2: Statistical Tests Based on Normality Assumptions

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
run_stat_test <- function(data, variable, group) {
  # Check normality for each gender
  normal <- check_normality_per_gender(data, variable, group)
  
  # Perform appropriate test
  if (normal) {
    t.test(as.formula(paste(variable, "~", group)), data = data)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  }
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(wg, .x, "gender"))

# Collect p-values for each variable
pvalue_list <- map_dbl(test_results, ~ .x$p.value)

# Extract and combine the test statistics and p-values for each variable
stat_list <- map_dbl(test_results, function(x) {
  if (inherits(x, "htest") && x$method == "Welch Two Sample t-test") {
    x$statistic  # Extract the t-value
  } else {
    x$statistic
  }
})

# Define Subscale labels
subscale_labels <- c("w.prod", "w.und", "early.ges", "later.ges", "tot.ges")

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(stat_list, 2),
  P_value = round(pvalue_list, 3),
  Test_Type = map_chr(test_results, ~ ifelse(inherits(.x, "htest") && .x$method == "Welch Two Sample t-test", "t-test", "Kruskal-Wallis"))
)

# --------------------------------------------------------------------------------------
# Step 3: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes for variables where p < 0.05
for (i in 1:nrow(results)) {
  if (results$P_value[i] < 0.05) {
    if (results$Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      effect_size <- wg %>% kruskal_effsize(as.formula(paste(var_list[i], "~ gender")))
      results$Effect_size[i] <- round(effect_size$effsize, 3)
    } else if (results$Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[i], "~ gender")), wg), 3)
    }
  }
}

# --------------------
# Display final table
print(results)

###################################################################################

# MacArthur-Bates Communicative Development Inventories (CDI)
# words and sentences (ws)
# ------------------------
# Step 1: Load WS Data

ws <- read.csv("ws.dd.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# --------------------------------------------------------
# Step 2: Statistical Tests Based on Normality Assumptions

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
run_stat_test <- function(data, variable, group) {
  # Check normality for each gender
  normal <- check_normality_per_gender(data, variable, group)
  
  # Perform appropriate test
  if (normal) {
    t.test(as.formula(paste(variable, "~", group)), data = data)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  }
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(ws, .x, "gender"))

# Collect p-values for each variable
pvalue_list <- map_dbl(test_results, ~ .x$p.value)

# Extract and combine the test statistics and p-values for each variable
stat_list <- map_dbl(test_results, function(x) {
  if (inherits(x, "htest") && x$method == "Welch Two Sample t-test") {
    x$statistic  # Extract the t-value
  } else {
    x$statistic
  }
})

# Define Subscale labels
subscale_labels <- c("w.prod")

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(stat_list, 2),
  P_value = round(pvalue_list, 3),
  Test_Type = map_chr(test_results, ~ ifelse(inherits(.x, "htest") && .x$method == "Welch Two Sample t-test", "t-test", "Kruskal-Wallis"))
)

# --------------------------------------------------------------------------------------
# Step 3: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes for variables where p < 0.05
for (i in 1:nrow(results)) {
  if (results$P_value[i] < 0.05) {
    if (results$Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      effect_size <- ws %>% kruskal_effsize(as.formula(paste(var_list[i], "~ gender")))
      results$Effect_size[i] <- round(effect_size$effsize, 3)
    } else if (results$Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[i], "~ gender")), ws), 3)
    }
  }
}

# --------------------
# Display final table
print(results)

##############################################################################################

#CSBS
# ------------------------
# Step 1: Load CSBS Data

csbs3 <- read.csv("csbs.dd.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# --------------------------------------------------------
# Step 2: Statistical Tests Based on Normality Assumptions

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
run_stat_test <- function(data, variable, group) {
  # Check normality for each gender
  normal <- check_normality_per_gender(data, variable, group)
  
  # Perform appropriate test
  if (normal) {
    t.test(as.formula(paste(variable, "~", group)), data = data)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  }
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(csbs3, .x, "gender"))

# Collect p-values for each variable
pvalue_list <- map_dbl(test_results, ~ .x$p.value)

# Extract and combine the test statistics and p-values for each variable
stat_list <- map_dbl(test_results, function(x) {
  if (inherits(x, "htest") && x$method == "Welch Two Sample t-test") {
    x$statistic  # Extract the t-value
  } else {
    x$statistic
  }
})

# Define Subscale labels
subscale_labels <- c("soc.pc", "exp.pc", "sym.pc", "tot.pc")

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(stat_list, 2),
  P_value = round(pvalue_list, 3),
  Test_Type = map_chr(test_results, ~ ifelse(inherits(.x, "htest") && .x$method == "Welch Two Sample t-test", "t-test", "Kruskal-Wallis"))
)

# ---------------------------------------------------------------------------------------
# Step 3: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes for variables where p < 0.05
for (i in 1:nrow(results)) {
  if (results$P_value[i] < 0.05) {
    if (results$Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      effect_size <- csbs3 %>% kruskal_effsize(as.formula(paste(var_list[i], "~ gender")))
      results$Effect_size[i] <- round(effect_size$effsize, 3)
    } else if (results$Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[i], "~ gender")), csbs3), 3)
    }
  }
}

# --------------------
# Display final table
print(results)

##############################################################################################################

# GeoPref
# ------------------------
# Step 1: Load GeoPref Data

geo3 <- read.csv("geopref.dd.csv", header = TRUE, na.strings = c("NULL", "NA",""))

# --------------------------------------------------------
# Step 2: Statistical Tests Based on Normality Assumptions

# Variables of interest
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
run_stat_test <- function(data, variable, group) {
  # Check normality for each gender
  normal <- check_normality_per_gender(data, variable, group)
  
  # Perform appropriate test
  if (normal) {
    t.test(as.formula(paste(variable, "~", group)), data = data)
  } else {
    kruskal.test(as.formula(paste(variable, "~", group)), data = data)
  }
}

# Run tests for each variable
test_results <- map(var_list, ~ run_stat_test(geo3, .x, "gender"))

# Collect p-values for each variable
pvalue_list <- map_dbl(test_results, ~ .x$p.value)

# Extract and combine the test statistics and p-values for each variable
stat_list <- map_dbl(test_results, function(x) {
  if (inherits(x, "htest") && x$method == "Welch Two Sample t-test") {
    x$statistic  # Extract the t-value
  } else {
    x$statistic
  }
})

# Define Subscale labels
subscale_labels <- c("perc.soc")

# Create the results dataframe
results <- data.frame(
  Subscale = subscale_labels,
  Statistic = round(stat_list, 2),
  P_value = round(pvalue_list, 3),
  Test_Type = map_chr(test_results, ~ ifelse(inherits(.x, "htest") && .x$method == "Welch Two Sample t-test", "t-test", "Kruskal-Wallis"))
)

# --------------------------------------------------------------------------------------
# Step 3: Effect Size Calculation (Cohen's d for t-test, or Kruskal-Wallis effect size)

# Add a column for Effect Size
results$Effect_size <- NA

# Calculate effect sizes for variables where p < 0.05
for (i in 1:nrow(results)) {
  if (results$P_value[i] < 0.05) {
    if (results$Test_Type[i] == "Kruskal-Wallis") {
      # Calculate Kruskal-Wallis effect size (eta squared)
      effect_size <- geo3 %>% kruskal_effsize(as.formula(paste(var_list[i], "~ gender")))
      results$Effect_size[i] <- round(effect_size$effsize, 3)
    } else if (results$Test_Type[i] == "t-test") {
      # Calculate Cohen's d
      results$Effect_size[i] <- round(cohensD(as.formula(paste(var_list[i], "~ gender")), geo3), 3)
    }
  }
}

# --------------------
# Display final table
print(results)


