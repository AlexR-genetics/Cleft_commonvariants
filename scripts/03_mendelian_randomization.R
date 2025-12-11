#!/usr/bin/env Rscript
# Mendelian Randomization Analysis


library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(data.table)
library(tidyverse)
library(MRPRESSO)

# -----------------------------------------------------------------------------
# Configuration - UPDATE THESE PATHS
# -----------------------------------------------------------------------------

# Working directory
WORKING_DIR <- "path/to/working/directory"
setwd(WORKING_DIR)

# Cleft GWAS summary statistics file
CLEFT_GWAS_FILE <- "cleft_gwas_sumstats.txt"

# Directory containing outcome GWAS files
OUTCOME_GWAS_DIR <- "outcome_gwas"

# Output directory
OUTPUT_DIR <- "MR_Results"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Sample sizes for cleft GWAS (update with your values)
N_CASES <- 2268
N_CONTROLS <- 7913



OUTCOME_FILES <- c(
  "ADHD_gwas.txt",
  "bipolar_gwas.txt",
  "educational_attainment_gwas.txt",
  "anxiety_gwas.txt",
  "autism_gwas.txt",
  "depression_gwas.txt",
  "intelligence_gwas.txt",
  "schizophrenia_gwas.txt"
)

# Clean names for plotting
OUTCOME_LABELS <- c(
  "ADHD",
  "Bipolar",
  "Education",
  "Anxiety",
  "Autism",
  "Depression",
  "Intelligence",
  "Schizophrenia"
)

# ==============================================================================
# PART 1: PREPARE EXPOSURE DATA (CLEFT)
# ==============================================================================


# Load cleft GWAS
cleft_gwas <- fread(CLEFT_GWAS_FILE)

# Calculate beta and SE from OR if needed
# Adjust column names based on your file structure
if (!"beta" %in% tolower(names(cleft_gwas))) {
  cleft_gwas$beta <- log(cleft_gwas$OR)
}

if (!"se" %in% tolower(names(cleft_gwas))) {
  cleft_gwas$z <- qnorm(cleft_gwas$P / 2, lower.tail = FALSE)
  cleft_gwas$se <- abs(cleft_gwas$beta) / cleft_gwas$z
}

# Calculate effect allele frequency if needed
if (!"eaf" %in% tolower(names(cleft_gwas))) {
  cleft_gwas$eaf <- (cleft_gwas$MAF_cases * N_CASES + 
                     cleft_gwas$MAF_controls * N_CONTROLS) / (N_CASES + N_CONTROLS)
}

cleft_gwas$samplesize <- N_CASES + N_CONTROLS

# Extract genome-wide significant SNPs (p < 5e-8)
cleft_gw <- subset(cleft_gwas, P < 5e-8)
cleft_gw <- as.data.frame(cleft_gw)

cat("Found", nrow(cleft_gw), "genome-wide significant SNPs\n")

# Format as exposure for TwoSampleMR
# Adjust column names to match your data
exposure_dat <- format_data(
  cleft_gw,
  type = "exposure",
  snp_col = "rsid",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  beta_col = "beta",
  se_col = "se",
  pval_col = "P",
  eaf_col = "eaf",
  ncase_col = "N_CASES",
  ncontrol_col = "N_CONTROLS"
)

exposure_dat_clumped <- clump_data(
  exposure_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  pop = "EUR"
)


# Calculate F-statistic for instrument strength
exposure_dat_clumped$F_stat <- (exposure_dat_clumped$beta.exposure^2) / 
                               (exposure_dat_clumped$se.exposure^2)
mean_F <- mean(exposure_dat_clumped$F_stat)

if (mean_F < 10) {
  warning("Mean F-statistic < 10, instruments may be weak")
}

# ==============================================================================
# PART 2: PREPARE OUTCOME DATA
# ==============================================================================


# Get list of instrumental SNPs
instrument_snps <- exposure_dat_clumped$SNP

# Process each outcome GWAS
harmonized_results <- list()

for (i in seq_along(OUTCOME_FILES)) {
  
  file <- OUTCOME_FILES[i]
  outcome_name <- OUTCOME_LABELS[i]
  
  cat("Processing:", outcome_name, "\n")
  
  # Load outcome GWAS
  gwas_data <- fread(file.path(OUTCOME_GWAS_DIR, file))
  
  # Filter to instrumental SNPs
  gwas_filtered <- gwas_data %>%
    filter(SNP %in% instrument_snps | rsid %in% instrument_snps)
  
  if (nrow(gwas_filtered) == 0) {
    cat("  Warning: No matching SNPs found for", outcome_name, "\n")
    next
  }
  
  # Calculate beta from OR if needed
  if (!"BETA" %in% names(gwas_filtered) & "OR" %in% names(gwas_filtered)) {
    gwas_filtered$BETA <- log(gwas_filtered$OR)
  }
  
  # Calculate SE if missing
  if (!"SE" %in% names(gwas_filtered) & "P" %in% names(gwas_filtered)) {
    gwas_filtered$SE <- abs(gwas_filtered$BETA / qnorm(gwas_filtered$P / 2))
  }
  
  gwas_filtered <- as.data.frame(gwas_filtered)
  
  # Format for TwoSampleMR
  outcome_dat <- format_data(
    gwas_filtered,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P"
  )
  
  outcome_dat$outcome <- outcome_name
  
  # Harmonize with exposure
  dat <- harmonise_data(
    exposure_dat = exposure_dat_clumped,
    outcome_dat = outcome_dat
  )
  
  harmonized_results[[outcome_name]] <- dat
  
  cat("  Harmonized:", nrow(dat), "SNPs\n")
}

# ==============================================================================
# PART 3: RUN MR ANALYSES
# ==============================================================================


# Run MR for each outcome
mr_results <- lapply(harmonized_results, function(dat) {
  if (nrow(dat) > 0) {
    mr(dat, method_list = c("mr_ivw", "mr_egger_regression", 
                            "mr_weighted_median", "mr_weighted_mode"))
  } else {
    NULL
  }
}) %>% bind_rows(.id = "outcome")

# Run sensitivity analyses
pleiotropy_results <- lapply(harmonized_results, mr_pleiotropy_test) %>%
  bind_rows(.id = "outcome")

heterogeneity_results <- lapply(harmonized_results, mr_heterogeneity) %>%
  bind_rows(.id = "outcome")

# ==============================================================================
# PART 4: GENERATE PLOTS
# ==============================================================================


# Individual SNP results
snp_results <- lapply(names(harmonized_results), function(outcome_name) {
  dat <- harmonized_results[[outcome_name]]
  if (nrow(dat) > 0) {
    res <- mr_singlesnp(dat)
    res$outcome <- outcome_name
    fwrite(res, file.path(OUTPUT_DIR, paste0(outcome_name, "_SNPresults.csv")))
    return(res)
  }
}) %>% bind_rows()

# Generate diagnostic plots for each outcome
for (outcome_name in names(harmonized_results)) {
  
  dat <- harmonized_results[[outcome_name]]
  if (nrow(dat) < 1) next
  
  # Scatter plot
  p_scatter <- mr_scatter_plot(mr(dat), dat)[[1]] +
    ggtitle(paste("MR Analysis:", outcome_name))
  ggsave(file.path(OUTPUT_DIR, paste0(outcome_name, "_scatter.png")),
         p_scatter, width = 8, height = 6, dpi = 300)
  
  # Forest plot
  p_forest <- mr_forest_plot(mr_singlesnp(dat))[[1]] +
    ggtitle(paste("SNP Effects:", outcome_name))
  ggsave(file.path(OUTPUT_DIR, paste0(outcome_name, "_forest.png")),
         p_forest, width = 10, height = 8, dpi = 300)
  
  # Funnel plot
  p_funnel <- mr_funnel_plot(mr_singlesnp(dat))[[1]] +
    ggtitle(paste("Funnel Plot:", outcome_name))
  ggsave(file.path(OUTPUT_DIR, paste0(outcome_name, "_funnel.png")),
         p_funnel, width = 8, height = 6, dpi = 300)
  
  # Leave-one-out analysis
  loo <- mr_leaveoneout(dat)
  p_loo <- mr_leaveoneout_plot(loo)[[1]] +
    ggtitle(paste("Leave-One-Out Analysis:", outcome_name))
  ggsave(file.path(OUTPUT_DIR, paste0(outcome_name, "_leaveoneout.png")),
         p_loo, width = 10, height = 8, dpi = 300)
}

# ==============================================================================
# PART 5: SUMMARY VOLCANO PLOT
# ==============================================================================


# Calculate Bonferroni threshold
n_tests <- nrow(mr_results)
bonf_threshold <- 0.05 / n_tests

# Prepare data for volcano plot
method_abbrev <- c(
  "Inverse variance weighted" = "IVW",
  "MR Egger" = "Egger",
  "Weighted median" = "WM",
  "Weighted mode" = "Mode"
)

volcano_data <- mr_results %>%
  mutate(
    log_p = -log10(pval),
    sig = ifelse(pval < bonf_threshold, "sig", "ns"),
    method_abbrev = method_abbrev[method],
    custom_label = paste(outcome, method_abbrev)
  )

bonf_line <- -log10(bonf_threshold)

p_volcano <- ggplot(volcano_data, aes(x = b, y = log_p, color = sig)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_hline(yintercept = bonf_line, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = custom_label), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("ns" = "grey60", "sig" = "firebrick")) +
  labs(
    x = "Effect Size (beta)",
    y = "-log10(p-value)",
    title = "Mendelian Randomization Results",
    subtitle = paste0("Bonferroni threshold (p = ", 
                      format(bonf_threshold, scientific = TRUE, digits = 2),
                      ") shown in red")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "MR_volcano_plot.png"),
       p_volcano, width = 10, height = 8, dpi = 300)

# ==============================================================================
# PART 6: SAVE RESULTS
# ==============================================================================

cat("\n=== Saving Results ===\n")

# Save main results
fwrite(mr_results, file.path(OUTPUT_DIR, "full_mr_results.csv"))
fwrite(pleiotropy_results, file.path(OUTPUT_DIR, "pleiotropy_results.csv"))
fwrite(heterogeneity_results, file.path(OUTPUT_DIR, "heterogeneity_results.csv"))
fwrite(snp_results, file.path(OUTPUT_DIR, "individual_snp_results.csv"))

# Create summary table
summary_table <- mr_results %>%
  filter(method == "Inverse variance weighted") %>%
  left_join(
    pleiotropy_results %>% select(outcome, egger_intercept, pval),
    by = "outcome",
    suffix = c("", "_pleiotropy")
  ) %>%
  left_join(
    heterogeneity_results %>% 
      filter(method == "Inverse variance weighted") %>%
      select(outcome, Q, Q_pval),
    by = "outcome"
  ) %>%
  select(outcome, nsnp, b, se, pval, Q, Q_pval, egger_intercept, pval_pleiotropy)

fwrite(summary_table, file.path(OUTPUT_DIR, "MR_summary_table.csv"))

