#!/bin/bash
# ==============================================================================
# PRSice-2 Polygenic Risk Score Calculation
# ==============================================================================
# Description: Calculates polygenic risk scores for neurodevelopmental traits
#              and cleft using PRSice-2.
#
# Prerequisites:
#   - PRSice-2 installed (https://www.prsice.info/)
#   - Quality-controlled target genotype data (PLINK format)
#   - GWAS summary statistics for each trait
#   - Phenotype file
#
# Usage: bash 00b_prsice_commands.sh
# ==============================================================================

# -----------------------------------------------------------------------------
# Configuration - UPDATE THESE PATHS
# -----------------------------------------------------------------------------

# Path to PRSice files
PRSICE_R="PRSice.R"
PRSICE_BINARY="PRSice_linux"

# Target data (QC'd genotypes in PLINK format)
TARGET="CC_QC1"

# Phenotype file
PHENO_FILE="CC_pheno1.txt"
PHENO_COL="PHENOTYPE"

# Output directory
OUT_DIR="prs_results"

mkdir -p ${OUT_DIR}

# MHC region to exclude (GRCh37 coordinates)
MHC_EXCLUDE="chr6:28510120-33480577"

# -----------------------------------------------------------------------------
# GWAS Summary Statistics Files - UPDATE THESE
# -----------------------------------------------------------------------------

# Discovery GWAS files (ensure correct column headers: SNP, A1, A2, CHR, BP, P, BETA/OR)
SCZ_GWAS="scz_pgc3_gwas.txt"
ADHD_GWAS="ADHD_ipsych_decode_PGC.txt"
ASD_GWAS="iPSYCH_autism.txt"
BIP_GWAS="BIP_2024_euronly.txt"
DEP_GWAS="iPSYCH_depression.txt"
ANX_GWAS="iPSYCH_anxiety.txt"
EA_GWAS="EA_gwas.txt"
INT_GWAS="savage_intelligence.txt"
CLEFT_GWAS="Cleft_meta.txt"

# ==============================================================================
# PRSice PARAMETERS
# ==============================================================================

# Standard parameters used across all analyses
# --bar-levels: p-value thresholds for score calculation
# --clump-kb: clumping window (250kb)
# --clump-r2: LD threshold for clumping (0.2)
# --fastscore: only calculate scores at specified thresholds
# --all-score: output scores at all thresholds
# --x-range: exclude MHC region

# ==============================================================================
# PRS CALCULATIONS
# ==============================================================================

echo "=== Calculating Polygenic Risk Scores ==="

# -----------------------------------------------------------------------------
# Schizophrenia PRS
# -----------------------------------------------------------------------------
echo "Calculating Schizophrenia PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${SCZ_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat BETA \
    --pvalue P \
    --bar-levels 0.001,0.01,0.05,0.25 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/sz_PRS_CC

# -----------------------------------------------------------------------------
# ADHD PRS
# -----------------------------------------------------------------------------
echo "Calculating ADHD PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${ADHD_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat OR \
    --pvalue P \
    --bar-levels 0.001,0.01,0.05,0.25 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/ADHD_PRS_CC

# -----------------------------------------------------------------------------
# Autism PRS
# -----------------------------------------------------------------------------
echo "Calculating Autism PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${ASD_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat OR \
    --pvalue P \
    --bar-levels 0.001,0.01,0.05,0.25 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/autism_PRS_CC

# -----------------------------------------------------------------------------
# Bipolar Disorder PRS
# -----------------------------------------------------------------------------
echo "Calculating Bipolar Disorder PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${BIP_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat OR \
    --pvalue P \
    --bar-levels 0.001,0.01,0.05,0.25 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/bipolar_PRS_eur_CC

# -----------------------------------------------------------------------------
# Depression PRS
# -----------------------------------------------------------------------------
echo "Calculating Depression PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${DEP_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat OR \
    --pvalue P \
    --bar-levels 0.001,0.01,0.05,0.25 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/depression_PRS_CC

# -----------------------------------------------------------------------------
# Anxiety PRS
# -----------------------------------------------------------------------------
echo "Calculating Anxiety PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${ANX_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat OR \
    --pvalue P \
    --bar-levels 0.001,0.01,0.05,0.25 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/anxiety_PRS_CC

# -----------------------------------------------------------------------------
# Educational Attainment PRS
# -----------------------------------------------------------------------------
echo "Calculating Educational Attainment PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${EA_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat BETA \
    --pvalue P \
    --bar-levels 0.001,0.01,0.05,0.25 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/educ_PRS_CC

# -----------------------------------------------------------------------------
# Intelligence PRS
# -----------------------------------------------------------------------------
echo "Calculating Intelligence PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${INT_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat BETA \
    --pvalue P \
    --bar-levels 0.001,0.01,0.05,0.25 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/int_PRS_CC

# -----------------------------------------------------------------------------
# Cleft PRS (using independent GWAS to avoid overfitting)
# Note: Different thresholds including genome-wide significant
# -----------------------------------------------------------------------------
echo "Calculating Cleft PRS..."

Rscript ${PRSICE_R} \
    --prsice ${PRSICE_BINARY} \
    --base ${CLEFT_GWAS} \
    --target ${TARGET} \
    --pheno ${PHENO_FILE} \
    --pheno-col ${PHENO_COL} \
    --binary-target T \
    --A1 A1 \
    --A2 A2 \
    --chr CHR \
    --snp SNP \
    --stat BETA \
    --pvalue P \
    --bar-levels 5e-08,0.0001,0.001,0.05,0.25,1 \
    --clump-kb 250 \
    --clump-r2 0.2 \
    --clump-p 1 \
    --fastscore \
    --all-score \
    --model add \
    --x-range ${MHC_EXCLUDE} \
    --thread 1 \
    --out ${OUT_DIR}/Cleft_PRS_meta

echo "=== PRS Calculation Complete ==="
echo "Results saved to: ${OUT_DIR}"
echo "Process .all_score files with 02_prs_analysis.R"
