#!/usr/bin/env Rscript
# Polygenic Risk Score Analysis

library(data.table)
library(tidyverse)
library(fmsb)
library(haven)
library(officer)
library(flextable)

# -----------------------------------------------------------------------------
# Configuration - UPDATE THESE PATHS
# -----------------------------------------------------------------------------

# Working directory containing PRSice output and phenotype files
WORKING_DIR <- "path/to/working/directory"
setwd(WORKING_DIR)

# Output directory for results
OUTPUT_DIR <- "results"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# File paths - UPDATE THESE
# -----------------------------------------------------------------------------

# PCA and phenotype files
PCA_FILE <- "CC_PCs.eigenvec"
EIGENVAL_FILE <- "CC_PCs.eigenval"
FAM_FILE <- "CC_QC1.fam"

# PRSice output files (p < 0.05 threshold scores)
PRS_FILES <- list(
  ADHD = "ADHD_PRS_CC.all_score",
  SCZ = "sz_PRS_CC.all_score",
  BP = "bipolar_PRS_eur_CC.all_score",
  ASD = "autism_PRS_CC.all_score",
  ANX = "anxiety_PRS_CC.all_score",
  EA = "educ_PRS_CC.all_score",
  DEP = "depression_PRS_CC.all_score",
  INT = "int_PRS_CC.all_score",
  CLEFT = "Cleft_PRS_meta.all_score",
)

# Phenotype data files
PHENOTYPE_CLASSIFIED_FILE <- "phenotypes_classified.csv"
STUDY_CHILDREN_FILE <- "studychildren_clean.csv"
MAIN_RELEASE_FILE <- "release_file.dta"

# -----------------------------------------------------------------------------
# Load and Prepare Data
# -----------------------------------------------------------------------------


# Load PCA data
pc1 <- read.table(PCA_FILE, header = FALSE)
names(pc1)[3:ncol(pc1)] <- paste0("PC", 1:(ncol(pc1) - 2))
names(pc1)[1:2] <- c("FID", "IID")

# Load eigenvalues
eigenval <- scan(EIGENVAL_FILE)

# Load phenotype from FAM file
pheno <- read.table(FAM_FILE, header = FALSE)
pheno <- pheno[, c(1, 2, 6)]
names(pheno) <- c("FID", "IID", "Phenotype")

# Merge PCA with phenotype
pca <- merge(pc1, pheno, by = c("FID", "IID"))
pca <- as_tibble(data.frame(pca))

# Code phenotype as factor
pca$Phenotype <- factor(
  pca$Phenotype,
  levels = c(1, 2),
  labels = c("control", "case")
)

# -----------------------------------------------------------------------------
# Load PRS Data
# -----------------------------------------------------------------------------


# Function to load and rename PRS file
load_prs <- function(filepath, score_name) {
  dat <- fread(filepath)
  # Assuming column 3 contains the p < 0.05 threshold score
  colnames(dat)[3] <- score_name
  return(dat[, c(1, 2, 3)])
}

# Load all PRS files
prs_adhd <- load_prs(PRS_FILES$ADHD, "PRS_0.05_ADHD")
prs_scz <- load_prs(PRS_FILES$SCZ, "PRS_0.05_SCZ")
prs_bp <- load_prs(PRS_FILES$BP, "PRS_0.05_BP")
prs_asd <- load_prs(PRS_FILES$ASD, "PRS_0.05_ASD")
prs_anx <- load_prs(PRS_FILES$ANX, "PRS_0.05_anx")
prs_ea <- load_prs(PRS_FILES$EA, "PRS_0.05_EA")
prs_dep <- load_prs(PRS_FILES$DEP, "PRS_0.05_DEP")
prs_int <- load_prs(PRS_FILES$INT, "PRS_0.05_INT")

# Cleft PRS has multiple thresholds - load columns 3 and 6
prs_cleft <- fread(PRS_FILES$CLEFT)
colnames(prs_cleft)[3] <- "PRS_Genomewide_CL"
colnames(prs_cleft)[6] <- "PRS_0.05_CL"
prs_cleft <- prs_cleft[, c(1, 2, 3, 6)]


# -----------------------------------------------------------------------------
# Merge All PRS Data
# -----------------------------------------------------------------------------


# Sequential merging of all PRS with PCA/phenotype data
scores_allsample <- pca
for (prs_df in list(prs_int, prs_scz, prs_bp, prs_asd, prs_ea, 
                    prs_anx, prs_adhd, prs_dep, prs_cleft, prs_dd)) {
  scores_allsample <- merge(scores_allsample, prs_df, 
                            by.x = c("FID", "IID"), 
                            by.y = c("FID", "IID"))
}

# -----------------------------------------------------------------------------
# Standardise PRS
# -----------------------------------------------------------------------------

prs_vars <- c('PRS_0.05_ADHD', 'PRS_0.05_DEP', 'PRS_0.05_anx', 'PRS_0.05_EA',
              'PRS_0.05_ASD', 'PRS_0.05_BP', 'PRS_0.05_INT', 'PRS_0.05_SCZ',
              'PRS_0.05_CL', 'PRS_Genomewide_CL', 'PRS_0.05_DD', 'PRS_Genomewide_DD')

scores_allsample <- scores_allsample %>%
  mutate(across(all_of(prs_vars), ~ scale(.) %>% as.vector))

# Note: ASD PRS direction may need inverting depending on GWAS coding
# scores_allsample$PRS_0.05_ASD <- -(scores_allsample$PRS_0.05_ASD)

# -----------------------------------------------------------------------------
# Load Additional Phenotype Data (for subtype and outcome analyses)
# -----------------------------------------------------------------------------

#Adjust on the basis of phenotype file structure
#end with scores_all_sample2 = merge (score_all_sample, phenotype_file by.x = (IID, FID),
# by y= (IID,FID))

# -----------------------------------------------------------------------------
# Create Cleft Subtype Variables
# -----------------------------------------------------------------------------


scores_all_sample2 <- scores_all_sample2 %>%
  mutate(
    # Cleft lip only vs control
    cleft_lip_only_vs_control = case_when(
      cleft_lip_only == TRUE & Phenotype == "case" ~ 1,
      Phenotype == "control" ~ 0,
      TRUE ~ NA_real_
    ),
    # Cleft lip only vs cleft palate only
    cleft_lip_only_vs_cleft_palate_only = case_when(
      cleft_lip_only == TRUE & Phenotype == "case" ~ 1,
      CPO == TRUE & Phenotype == "case" ~ 0,
      TRUE ~ NA_real_
    ),
    # Cleft lip only vs cleft lip and palate
    cleft_lip_only_vs_cleft_lip_and_palate = case_when(
      cleft_lip_only == TRUE & Phenotype == "case" ~ 1,
      cleft_lip_with_without_palate == TRUE & Phenotype == "case" & 
        cleft_lip_only == FALSE ~ 0,
      TRUE ~ NA_real_
    ),
    # Cleft palate only vs control
    cleft_palate_only_vs_control = case_when(
      CPO == TRUE & Phenotype == "case" ~ 1,
      Phenotype == "control" ~ 0,
      TRUE ~ NA_real_
    ),
    # Cleft lip and palate vs control
    cleft_lip_and_palate_vs_control = case_when(
      cleft_lip_with_without_palate == TRUE & cleft_lip_only == FALSE & 
        Phenotype == "case" ~ 1,
      Phenotype == "control" ~ 0,
      TRUE ~ NA_real_
    ),
    # Cleft lip and palate vs cleft palate only
    cleft_lip_and_palate_vs_cleft_palate_only = case_when(
      cleft_lip_with_without_palate == TRUE & cleft_lip_only == FALSE ~ 1,
      CPO == TRUE & Phenotype == "case" ~ 0,
      TRUE ~ NA_real_
    ),
    # Any cleft lip vs cleft palate only
    any_cleft_lip_vs_cleft_palate_only = case_when(
      cleft_lip_with_without_palate == TRUE ~ 1,
      CPO == TRUE & Phenotype == "case" ~ 0,
      TRUE ~ NA_real_
    )
  )

# ==============================================================================
# ANALYSIS 1: CASE-CONTROL COMPARISONS
# ==============================================================================


variables <- c('PRS_0.05_ADHD', 'PRS_0.05_DEP', 'PRS_0.05_anx', 'PRS_0.05_EA',
               'PRS_0.05_ASD', 'PRS_0.05_BP', 'PRS_0.05_INT', 'PRS_0.05_SCZ',
               'PRS_0.05_CL', 'PRS_Genomewide_CL')

results_case_control <- data.frame(
  Variable = character(),
  Lower_CI = numeric(),
  OR = numeric(),
  Upper_CI = numeric(),
  P_Value = numeric(),
  R2_Change = numeric(),
  stringsAsFactors = FALSE
)

for (var in variables) {
  
  # Full model with PRS
  model_full <- glm(
    as.formula(paste0("Phenotype ~ ", var, 
                      " + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),
    data = scores_allsample,
    family = "binomial"
  )
  
  # Null model (PCs only)
  model_null <- glm(
    Phenotype ~  gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = scores_allsample,
    family = "binomial"
  )
  
  # Extract statistics
  lower_ci <- exp(confint.default(model_full)[2, 1])
  or_value <- exp(model_full$coefficients[2])
  upper_ci <- exp(confint.default(model_full)[2, 2])
  p_value <- coef(summary(model_full))[2, 'Pr(>|z|)']
  r2_change <- (NagelkerkeR2(model_full)$R2 - NagelkerkeR2(model_null)$R2) * 100
  
  results_case_control <- rbind(
    results_case_control,
    data.frame(
      Variable = var,
      Lower_CI = lower_ci,
      OR = or_value,
      Upper_CI = upper_ci,
      P_Value = p_value,
      R2_Change = r2_change
    )
  )
}

# Save results
doc <- read_docx()
doc <- body_add_flextable(doc, flextable(results_case_control))
print(doc, target = file.path(OUTPUT_DIR, "results_CASE_CONTROL.docx"))


# ==============================================================================
# ANALYSIS 2: CLEFT SUBTYPE COMPARISONS
# ==============================================================================


phenotypes <- c(
  'cleft_lip_only_vs_control',
  'cleft_lip_only_vs_cleft_palate_only',
  'cleft_lip_only_vs_cleft_lip_and_palate',
  'cleft_palate_only_vs_control',
  'cleft_lip_and_palate_vs_control',
  'cleft_lip_and_palate_vs_cleft_palate_only'
)

results_subtypes <- data.frame()

for (phenotype in phenotypes) {
  for (var in variables) {
    
    # Full model
    model_full <- glm(
      as.formula(paste0(phenotype, " ~ ", var,
                        " +  gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),
      data = scores_all_sample2,
      family = "binomial"
    )
    
    # Null model
    model_null <- glm(
      as.formula(paste0(phenotype, 
                        " ~   gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),
      data = scores_all_sample2,
      family = "binomial"
    )
    
    # Extract statistics
    lower_ci <- exp(confint.default(model_full)[2, 1])
    or_value <- exp(model_full$coefficients[2])
    upper_ci <- exp(confint.default(model_full)[2, 2])
    p_value <- coef(summary(model_full))[2, 'Pr(>|z|)']
    r2_change <- (NagelkerkeR2(model_full)$R2 - NagelkerkeR2(model_null)$R2) * 100
    
    results_subtypes <- rbind(
      results_subtypes,
      data.frame(
        Phenotype = phenotype,
        Variable = var,
        Lower_CI = lower_ci,
        OR = or_value,
        Upper_CI = upper_ci,
        P_Value = p_value,
        R2_Change = r2_change
      )
    )
  }
}

doc <- read_docx()
doc <- body_add_flextable(doc, flextable(results_subtypes))
print(doc, target = file.path(OUTPUT_DIR, "results_CleftSubtypes.docx"))


# ==============================================================================
# ANALYSIS 3: BEHAVIOURAL/DEVELOPMENTAL OUTCOMES
# ==============================================================================

cat("\n=== Running Behavioural Outcome Analyses ===\n")

# Function to combine maternal and paternal ratings
combine_ratings <- function(maternal, paternal) {
  ifelse(!is.na(maternal) & !is.na(paternal),
         (maternal + paternal) / 2,
         ifelse(!is.na(maternal), maternal, paternal))
}

scores_all_sample2 <- scores_all_sample2 %>%
  mutate(
    cc5y_10y_sdq_total_combined = combine_ratings(cc5ym_10y_sdq_total, cc5yf_10y_sdq_total),
    cc5y_10y_mfqtotal_combined = combine_ratings(cc5ym_10y_mfqtotal, cc5yf_10y_mfqtotal),
    cc5y_10y_scaredadjtot_combined = combine_ratings(cc5ym_10y_scaredadjtot, cc5yf_10y_scaredadjtot),
    cc5y_8yr_sdq_total_combined = combine_ratings(cc5ym_8yr_sdq_total, cc5yf_8yr_sdq_total),
    ccpn_5yr_sdq_total_combined = combine_ratings(ccpnm_5yr_sdq_total, ccpnf_5yr_sdq_total),
    ccpn_5yr_asqse_adjscr_combined = combine_ratings(ccpnm_5yr_asqse_adjscr, ccpnf_5yr_asqse_adjscr),
    ccpn_3yr_asqse_adjscr_combined = combine_ratings(ccpnm_3yr_asqse_adjscr, ccpnf_3yr_asqse_adjscr),
    ccpn_18m_asqse_adjscr_combined = combine_ratings(ccpnm_18m_asqse_adjscr, ccpnf_18m_asqse_adjscr),
    ccpn_18m_asqse_probsolv_combined = combine_ratings(ccpnm_18m_problemsolv_adjscore, ccpnf_18m_problemsolv_adjscore),
    ccpn_18m_asqse_persoc_combined = combine_ratings(ccpnm_18m_persoc_adjscore, ccpnf_18m_persoc_adjscore),
    ccpn_18m_asqse_finemotor_combined = combine_ratings(ccpnm_18m_finemotor_adjscore, ccpnf_18m_finemotor_adjscore),
    ccpn_18m_asqse_grossmotor_combined = combine_ratings(ccpnm_18m_grossmotor_adjscore, ccpnf_18m_grossmotor_adjscore),
    ccpn_18m_asqse_communication_combined = combine_ratings(ccpnm_18m_communication_adjscore, ccpnf_18m_communication_adjscore),
    ccpn_3yr_asqse_probsolv_combined = combine_ratings(ccpnm_3yr_problemsolv_adjscr, ccpnf_3yr_problemsolv_adjscr),
    ccpn_3yr_asqse_persoc_combined = combine_ratings(ccpnm_3yr_persoc_adjscore, ccpnf_3yr_persoc_adjscore),
    ccpn_3yr_asqse_finemotor_combined = combine_ratings(ccpnm_3yr_finemotor_adjscore, ccpnf_3yr_finemotor_adjscore),
    ccpn_3yr_asqse_grossmotor_combined = combine_ratings(ccpnm_3yr_grossmotor_adjscore, ccpnf_3yr_grossmotor_adjscore),
    ccpn_3yr_asqse_communication_combined = combine_ratings(ccpnm_3yr_communication_adjscore, ccpnf_3yr_communication_adjscore),
    ccpn_5yr_asqse_probsolv_combined = combine_ratings(ccpnm_5yr_problemsolv_adjscr, ccpnf_5yr_problemsolv_adjscr),
    ccpn_5yr_asqse_persoc_combined = combine_ratings(ccpnm_5yr_persoc_adjscore, ccpnf_5yr_persoc_adjscore),
    ccpn_5yr_asqse_finemotor_combined = combine_ratings(ccpnm_5yr_finemotor_adjscore, ccpnf_5yr_finemotor_adjscore),
    ccpn_5yr_asqse_grossmotor_combined = combine_ratings(ccpnm_5yr_grossmotor_adjscore, ccpnf_5yr_grossmotor_adjscore),
    ccpn_5yr_asqse_communication_combined = combine_ratings(ccpnm_5yr_communication_adjscore, ccpnf_5yr_communication_adjscore),
    ccpn_5yr_sdq_emotional_combined = combine_ratings(ccpnm_5yr_sdq_emotional, ccpnf_5yr_sdq_emotional),
    ccpn_5yr_sdq_conduct_combined = combine_ratings(ccpnm_5yr_sdq_conduct, ccpnf_5yr_sdq_conduct),
    ccpn_5yr_sdq_hyperactivity_combined = combine_ratings(ccpnm_5yr_sdq_hyperactivity, ccpnf_5yr_sdq_hyperactivity),
    ccpn_5yr_sdq_peer_combined = combine_ratings(ccpnm_5yr_sdq_peer, ccpnf_5yr_sdq_peer),
    ccpn_5yr_sdq_prosocial_combined = combine_ratings(ccpnm_5yr_sdq_prosocial, ccpnf_5yr_sdq_prosocial),
    cc5y_8yr_sdq_emotional_combined = combine_ratings(cc5ym_8yr_sdq_emotional, cc5yf_8yr_sdq_emotional),
    cc5y_8yr_sdq_conduct_combined = combine_ratings(cc5ym_8yr_sdq_conduct, cc5yf_8yr_sdq_conduct),
    cc5y_8yr_sdq_hyperactivity_combined = combine_ratings(cc5ym_8yr_sdq_hyperactivity, cc5yf_8yr_sdq_hyperactivity),
    cc5y_8yr_sdq_peer_combined = combine_ratings(cc5ym_8yr_sdq_peer, cc5yf_8yr_sdq_peer),
    cc5y_8yr_sdq_prosocial_combined = combine_ratings(cc5ym_8yr_sdq_prosocial, cc5yf_8yr_sdq_prosocial),
    cc5y_10y_sdq_emotional_combined = combine_ratings(cc5ym_10y_sdq_emotional, cc5yf_10y_sdq_emotional),
    cc5y_10y_sdq_conduct_combined = combine_ratings(cc5ym_10y_sdq_conduct, cc5yf_10y_sdq_conduct),
    cc5y_10y_sdq_hyperactivity_combined = combine_ratings(cc5ym_10y_sdq_hyperactivity, cc5yf_10y_sdq_hyperactivity),
    cc5y_10y_sdq_peer_combined = combine_ratings(cc5ym_10y_sdq_peer, cc5yf_10y_sdq_peer),
    cc5y_10y_sdq_prosocial_combined = combine_ratings(cc5ym_10y_sdq_prosocial, cc5yf_10y_sdq_prosocial)
    
    
     )

# Define the new continuous outcome variables
outcomes <- c(
  "cc5y_10y_sdq_total_combined",
  "cc5y_10y_mfqtotal_combined",
  "cc5y_10y_scaredadjtot_combined",
  "cc5y_8yr_sdq_total_combined",
  "ccpn_5yr_sdq_total_combined",
  "ccpn_5yr_asqse_adjscr_combined",
  "ccpn_3yr_asqse_adjscr_combined",
  "ccpn_18m_asqse_adjscr_combined",
  "ccpn_18m_asqse_probsolv_combined",
  "ccpn_18m_asqse_persoc_combined",
  "ccpn_18m_asqse_finemotor_combined",
  "ccpn_18m_asqse_grossmotor_combined",
  "ccpn_18m_asqse_communication_combined",
  "ccpn_3yr_asqse_probsolv_combined",
  "ccpn_3yr_asqse_persoc_combined",
  "ccpn_3yr_asqse_finemotor_combined",
  "ccpn_3yr_asqse_grossmotor_combined",
  "ccpn_3yr_asqse_communication_combined",
  "ccpn_5yr_asqse_probsolv_combined",
  "ccpn_5yr_asqse_persoc_combined",
  "ccpn_5yr_asqse_finemotor_combined",
  "ccpn_5yr_asqse_grossmotor_combined",
  "ccpn_5yr_asqse_communication_combined",
  "ccpn_5yr_sdq_emotional_combined",
  "ccpn_5yr_sdq_conduct_combined",
  "ccpn_5yr_sdq_hyperactivity_combined",
  "ccpn_5yr_sdq_peer_combined",
  "ccpn_5yr_sdq_prosocial_combined",
  "cc5y_8yr_sdq_emotional_combined",
  "cc5y_8yr_sdq_conduct_combined",
  "cc5y_8yr_sdq_hyperactivity_combined",
  "cc5y_8yr_sdq_peer_combined",
  "cc5y_8yr_sdq_prosocial_combined",
  "cc5y_10y_sdq_emotional_combined",
  "cc5y_10y_sdq_conduct_combined",
  "cc5y_10y_sdq_hyperactivity_combined",
  "cc5y_10y_sdq_peer_combined",
  "cc5y_10y_sdq_prosocial_combined"
)




results_continuous <- data.frame()

for (outcome in outcomes) {
  for (var in variables) {
    
    # Linear regression with PRS
    model_full <- lm(
      as.formula(paste0(outcome, " ~ ", var,
                        " + gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),
      data = scores_all_sample2
    )
    
    # Null model
    model_null <- lm(
      as.formula(paste0(outcome,
                        " ~ gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),
      data = scores_all_sample2
    )
    
    # Extract statistics
    effect_size <- coef(model_full)[2]
    lower_ci <- confint(model_full)[2, 1]
    upper_ci <- confint(model_full)[2, 2]
    p_value <- summary(model_full)$coefficients[2, 'Pr(>|t|)']
    r2_diff <- (summary(model_full)$r.squared - summary(model_null)$r.squared) * 100
    n_obs <- nobs(model_full)
    
    results_continuous <- rbind(
      results_continuous,
      data.frame(
        Outcome = outcome,
        Variable = var,
        Effect_Size = effect_size,
        Lower_CI = lower_ci,
        Upper_CI = upper_ci,
        P_Value = p_value,
        R2_Value = r2_diff,
        N = n_obs
      )
    )
  }
}

doc <- read_docx()
doc <- body_add_flextable(doc, flextable(results_continuous))
print(doc, target = file.path(OUTPUT_DIR, "results_SCALES.docx"))


# ==============================================================================
# ANALYSIS 4: ND-CNV CARRIER ANALYSIS
# ==============================================================================

cat("\n=== Running ND-CNV Carrier Analyses ===\n")

# Analysis comparing PRS between CNV carriers and non-carriers within cases
results_cnv <- data.frame()

for (var in variables) {
  
  model_full <- glm(
    as.formula(paste0("NDD_CNV ~ ", var,
                      " + gender + + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")),
    data = scores_all_sample2,
    family = "binomial"
  )
  
  model_null <- glm(
    NDD_CNV ~ gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
    data = scores_all_sample2,
    family = "binomial"
  )
  
  lower_ci <- exp(confint.default(model_full)[2, 1])
  or_value <- exp(model_full$coefficients[2])
  upper_ci <- exp(confint.default(model_full)[2, 2])
  p_value <- coef(summary(model_full))[2, 'Pr(>|z|)']
  r2_change <- (NagelkerkeR2(model_full)$R2 - NagelkerkeR2(model_null)$R2) * 100
  
  results_cnv <- rbind(
    results_cnv,
    data.frame(
      Variable = var,
      Lower_CI = lower_ci,
      OR = or_value,
      Upper_CI = upper_ci,
      P_Value = p_value,
      R2_Change = r2_change
    )
  )
}

doc <- read_docx()
doc <- body_add_flextable(doc, flextable(results_cnv))
print(doc, target = file.path(OUTPUT_DIR, "results_CNV_carriers.docx"))

# Additional analysis controlling for cleft type
results_cnv_adj <- data.frame()

for (var in variables) {
  
  model_full <- glm(
    as.formula(paste0("NDD_CNV ~ ", var,
                      " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10",
                      " + any_cleft_lip_vs_cleft_palate_only")),
    data = scores_all_sample2,
    family = "binomial"
  )
  
  model_null <- glm(
    as.formula(paste0("NDD_CNV ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10",
                      " + any_cleft_lip_vs_cleft_palate_only")),
    data = scores_all_sample2,
    family = "binomial"
  )
  
  lower_ci <- exp(confint.default(model_full)[2, 1])
  or_value <- exp(model_full$coefficients[2])
  upper_ci <- exp(confint.default(model_full)[2, 2])
  p_value <- coef(summary(model_full))[2, 'Pr(>|z|)']
  r2_change <- (NagelkerkeR2(model_full)$R2 - NagelkerkeR2(model_null)$R2) * 100
  
  results_cnv_adj <- rbind(
    results_cnv_adj,
    data.frame(
      Variable = var,
      Lower_CI = lower_ci,
      OR = or_value,
      Upper_CI = upper_ci,
      P_Value = p_value,
      R2_Change = r2_change
    )
  )
}

doc <- read_docx()
doc <- body_add_flextable(doc, flextable(results_cnv_adj))
print(doc, target = file.path(OUTPUT_DIR, "results_CNV_carriers_adjusted.docx"))
