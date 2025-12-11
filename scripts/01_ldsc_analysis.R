#!/usr/bin/env Rscript
# ==============================================================================
# LDSC Genetic Correlation Analysis
# ==============================================================================
# Description: Processes LDSC output files to extract genetic correlations 
#              between cleft lip/palate and neurodevelopmental traits.
#              Generates forest plots, heatmaps, and summary statistics.
#
# Input:  LDSC .log output files from ldsc.py --rg analyses
# Output: Forest plot, heatmap, summary CSV
#
# Required packages: tidyverse, ggplot2, corrplot, pheatmap, viridis
# ==============================================================================

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

# Install packages if needed (uncomment if required)
# install.packages(c("tidyverse", "ggplot2", "corrplot", "pheatmap", "viridis"))

library(tidyverse)
library(ggplot2)
library(corrplot)
library(pheatmap)
library(viridis)

# -----------------------------------------------------------------------------
# Configuration - UPDATE THESE PATHS
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

#' Extract results from LDSC output file
#' 
#' @param file_path Path to LDSC .log output file
#' @return Data frame with genetic correlation statistics
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

# -----------------------------------------------------------------------------
# Main Analysis
# -----------------------------------------------------------------------------

# Read all LDSC output files
ldsc_files <- list.files(RESULTS_FOLDER, full.names = TRUE, pattern = "\\.log$")

if (length(ldsc_files) == 0) {
  stop("No LDSC output files found in: ", RESULTS_FOLDER)
}

cat("Processing", length(ldsc_files), "LDSC output files...\n")

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

# Print summary
cat("\n=== Summary of LDSC Results ===\n\n")
print(all_results %>% select(trait1, trait2, rg, se, p, significant))

# -----------------------------------------------------------------------------
# Visualisation 1: Forest Plot
# -----------------------------------------------------------------------------

p_forest <- ggplot(all_results, aes(x = rg, y = paste(trait1, "-", trait2))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = rg_ci_lower, xmax = rg_ci_upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = significant, x = rg_ci_upper + 0.05), size = 5) +
  theme_minimal() +
  labs(
    x = "Genetic Correlation (rg)",
    y = "Trait Pairs",
    title = "Genetic Correlations from LD Score Regression",
    subtitle = "Error bars represent 95% confidence intervals"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 10)
  )

ggsave(
  file.path(OUTPUT_DIR, "genetic_correlations_forest.png"), 
  p_forest, 
  width = 10, 
  height = 8,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Visualisation 2: Correlation Heatmap
# -----------------------------------------------------------------------------

# Get unique traits
unique_traits <- unique(c(all_results$trait1, all_results$trait2))
n_traits <- length(unique_traits)

# Initialize correlation matrix
cor_matrix <- matrix(
  NA, n_traits, n_traits,
  dimnames = list(unique_traits, unique_traits)
)
diag(cor_matrix) <- 1

# Fill the matrix (symmetric)
for (i in 1:nrow(all_results)) {
  row_idx <- which(unique_traits == all_results$trait1[i])
  col_idx <- which(unique_traits == all_results$trait2[i])
  cor_matrix[row_idx, col_idx] <- all_results$rg[i]
  cor_matrix[col_idx, row_idx] <- all_results$rg[i]
}

# Generate heatmap
png(file.path(OUTPUT_DIR, "genetic_correlations_heatmap.png"), 
    width = 10, height = 8, units = "in", res = 300)

pheatmap(
  cor_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Genetic Correlation Heatmap",
  fontsize = 10,
  fontsize_number = 8
)

dev.off()

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------

summary_stats <- all_results %>%
  summarise(
    n_comparisons = n(),
    n_significant_p05 = sum(p < 0.05),
    n_significant_p01 = sum(p < 0.01),
    mean_abs_rg = mean(abs(rg)),
    median_abs_rg = median(abs(rg)),
    range_rg = paste(round(min(rg), 3), "to", round(max(rg), 3))
  )

cat("\n=== Summary Statistics ===\n")
print(summary_stats)

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

write.csv(all_results, file.path(OUTPUT_DIR, "ldsc_combined_results.csv"), row.names = FALSE)
cat("\nResults saved to:", file.path(OUTPUT_DIR, "ldsc_combined_results.csv"), "\n")

cat("\n=== Analysis Complete ===\n")
