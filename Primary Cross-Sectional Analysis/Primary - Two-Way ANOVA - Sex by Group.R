

# """
#
# Filename: Primary - Two-Way ANOVA - Sex by Group.R
# Version: 2024.09.0+375 RStudio 
# Author: Sanaz Nazari
# Date Created: Nov 3, 2024
# Last Modified: Nov 3, 2024
# Description: # Two-Way ANOVA with Sex and Group factors for all datasets.
#
# """

#-----------------------
# Load required library
library(ez)


# Define datasets and their respective configurations
all_datasets <- list(
  list(filename = "df.match.ados2.csv", id_seq = 1805, cols = 7:9, output_csv = "anova.ados.csv"),
  list(filename = "df.match.vine2.csv", id_seq = 1815, cols = 7:11, output_csv = "anova.vine.csv"),
  list(filename = "df.match.mul2.csv", id_seq = 1451, cols = 7:10, output_csv = "anova.mul.csv"),
  list(filename = "df.match.wg.csv", id_seq = 465, cols = 7:11, output_csv = "anova.wg.csv"),
  list(filename = "unq.ws2.csv", id_seq = 892, cols = 7, output_csv = "anova.ws.csv"),
  list(filename = "df.match.csbs2.csv", id_seq = 914, cols = 7:10, output_csv = "anova.csbs.csv", preprocess = TRUE),
  list(filename = "df.match.geopref.csv", id_seq = 968, cols = 8, output_csv = "anova.geopref.csv")
)


# Initialize consolidated results
consolidated_results <- NULL

# Function to perform ANOVA for a dataset
perform_anova <- function(dataset_info) {
  # Load dataset
  df.match <- read.csv(dataset_info$filename, header = TRUE, na.strings = c("NULL", "NA"))
  
  # Assign id sequence and remove unnecessary columns if needed
  id <- seq(1, dataset_info$id_seq, by = 1)
  df.match$id <- id
  df.match$group <- NULL
  
  # Specific preprocessing for CSBS dataset
  if (!is.null(dataset_info$preprocess) && dataset_info$preprocess) {
    df.match <- df.match[, 1:10]  # Select specific columns for CSBS
  }
  
  # Remove NA values
  df.match <- na.omit(df.match)
  
  # T-test for age comparison
  if ("age" %in% colnames(df.match)) {
    if (length(df.match$age[df.match$diag == "ASD"]) > 1 && length(df.match$age[df.match$diag == "TD"]) > 1) {
      t_test_result <- t.test(df.match$age[df.match$diag == "ASD"], df.match$age[df.match$diag == "TD"])
      print(t_test_result)
    } else {
      message("Not enough observations for age comparison.")
    }
  }
  
  # Convert columns to factors
  df.match$gender <- as.factor(df.match$gender)
  df.match$id <- as.factor(df.match$id)
  df.match$diag <- as.factor(df.match$diag)
  
  # Initialize results container
  res.anova <- NULL
  name_list <- colnames(df.match)
  data_copy <- df.match
  
  # Loop over specified columns for ANOVA
  for (i in dataset_info$cols) {
    # Rename column to "var" for ezANOVA
    colnames(data_copy)[i] <- "var"
    
    # Perform ANOVA with error handling
    tryCatch({
      anv <- ez::ezANOVA(data = data_copy, dv = var, wid = id, type = 3, between = .(gender, diag), detailed = TRUE, return_aov = TRUE)
      colnames(data_copy)[i] <- colnames(df.match)[i]
      temp <- cbind(dataset_info$filename, colnames(df.match)[i], round(anv$ANOVA$ges, 3), anv$ANOVA$`p<.05`, anv$ANOVA$Effect, round(anv$ANOVA$F, 2),
                    anv$ANOVA$DFn, anv$ANOVA$DFd, round(anv$ANOVA$p, 3))
      res.anova <- rbind(res.anova, temp)
    }, error = function(e) {
      message(paste("Error in ANOVA for variable", name_list[i], "in dataset", dataset_info$filename, ":", e$message))
    })
    
    # Restore original column name
    colnames(data_copy)[i] <- name_list[i]
  }
  
  # Write individual results to CSV
  if (!is.null(res.anova)) {
    write.csv(res.anova, dataset_info$output_csv)
    message(paste("ANOVA results written to", dataset_info$output_csv))
  } else {
    message(paste("No ANOVA results were generated for the dataset", dataset_info$filename))
  }
  
  return(res.anova)
}

# Loop through each dataset and perform ANOVA, consolidating the results
for (dataset_info in all_datasets) {
  res.anova <- perform_anova(dataset_info)
  
  # Append results to consolidated table
  if (!is.null(res.anova)) {
    consolidated_results <- rbind(consolidated_results, res.anova)
  } else {
    message(paste("No ANOVA results were generated for the dataset", dataset_info$filename))
  }
}

# Check and write the final consolidated results
if (!is.null(consolidated_results)) {
  colnames(consolidated_results) <- c("Dataset", "Variable", "GES", "P<.05", "Effect", "F_value", "DFn", "DFd", "P_value")
  write.csv(consolidated_results, "consolidated_anova_results.csv")
  print(consolidated_results)
} else {
  message("No consolidated ANOVA results were generated.")
}







