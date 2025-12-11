#!/usr/bin/env Rscript

# Install packages if needed (uncomment if required)
# install.packages(c("tidyverse", "ggplot2", "corrplot", "pheatmap", "viridis"))

library(tidyverse)
library(ggplot2)
library(corrplot)
library(pheatmap)
library(viridis)

# -----------------------------------------------------------------------------
# Configuration 
# -----------------------------------------------------------------------------

# Path to folder containing LDSC output .log files
RESULTS_FOLDER <- "path/to/ldsc/results"

# Output directory for plots and results
OUTPUT_DIR <- "output"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

extract_ldsc_results <- function(file_path) {
  
  # Read the file
  lines <- readLines(file_path)
  

  # Extract trait names from the --rg command
  rg_line <- grep("^--rg", lines, value = TRUE)
  traits <- str_extract(rg_line, "(?<=--rg ).*") %>%
    str_split(",") %>%
    unlist() %>%
    str_remove_all("MUNGE_|\\.sumstats\\.gz")
  
  # Extract genetic correlation statistics
  rg_val <- as.numeric(str_extract(
    grep("^Genetic Correlation:", lines, value = TRUE), 
    "-?\\d+\\.\\d+"
  ))
  
  se_val <- as.numeric(str_extract(
    grep("^Genetic Correlation:", lines, value = TRUE), 
    "(?<=\\()\\d+\\.\\d+"
  ))
  
  z_val <- as.numeric(str_extract(
    grep("^Z-score:", lines, value = TRUE), 
    "-?\\d+\\.\\d+"
  ))
  
  p_val <- as.numeric(str_extract(
    grep("^P:", lines, value = TRUE), 
    "\\d+\\.\\d+"
  ))
  
  # Extract heritability values
  h2_lines <- grep("^Total Observed scale h2:", lines, value = TRUE)
  h2_trait1 <- as.numeric(str_extract(h2_lines[1], "\\d+\\.\\d+"))
  h2_trait2 <- as.numeric(str_extract(h2_lines[2], "\\d+\\.\\d+"))
  
  # Return results as a data frame
  data.frame(
    trait1 = traits[1],
    trait2 = traits[2],
    rg = rg_val,
    se = se_val,
    z = z_val,
    p = p_val,
    h2_trait1 = h2_trait1,
    h2_trait2 = h2_trait2,
    stringsAsFactors = FALSE
  )
}

ldsc_files <- list.files(RESULTS_FOLDER, full.names = TRUE, pattern = "\\.log$")

if (length(ldsc_files) == 0) {
  stop("No LDSC output files found in: ", RESULTS_FOLDER)
}


# Extract results from all files
all_results <- map_df(ldsc_files, extract_ldsc_results)

# Add significance categories and confidence intervals
all_results <- all_results %>%
  mutate(
    significant = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    rg_ci_lower = rg - 1.96 * se,
    rg_ci_upper = rg + 1.96 * se
  )




write.csv(all_results, file.path(OUTPUT_DIR, "ldsc_combined_results.csv"), row.names = FALSE)
