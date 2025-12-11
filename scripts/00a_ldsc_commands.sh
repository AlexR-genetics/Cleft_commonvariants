#!/bin/bash
# ==============================================================================
# LDSC Analysis Pipeline
# ==============================================================================
# Description: Runs LD Score Regression to estimate genetic correlations
#              between cleft lip/palate and neurodevelopmental traits.
#
# Prerequisites:
#   - LDSC installed (https://github.com/bulik/ldsc)
#   - Munge'd summary statistics for each trait
#   - LD score reference panel (eur_w_ld_chr/)
#   - HapMap3 SNP list (w_hm3.snplist)
#
# Usage: bash 00a_ldsc_commands.sh
# ==============================================================================

# -----------------------------------------------------------------------------
# Configuration - UPDATE THESE PATHS
# -----------------------------------------------------------------------------

# Path to LDSC installation
LDSC_DIR="/path/to/ldsc"

# Path to LD score reference files
LD_REF="eur_w_ld_chr/"

# Path to HapMap3 SNP list for munging
HM3_SNPS="w_hm3.snplist"

# Output directory
OUT_DIR="ldsc_results"

mkdir -p ${OUT_DIR}

# -----------------------------------------------------------------------------
# GWAS Summary Statistics Files - UPDATE THESE
# -----------------------------------------------------------------------------

# Cleft GWAS
CLEFT_SUMSTATS="cleft_gwas_sumstats.txt"

# Neurodevelopmental/psychiatric GWAS files
ADHD_SUMSTATS="ADHD_gwas.txt"
AUTISM_SUMSTATS="autism_gwas.txt"
SCZ_SUMSTATS="schizophrenia_gwas.txt"
BIPOLAR_SUMSTATS="bipolar_gwas.txt"
DEPRESSION_SUMSTATS="depression_gwas.txt"
ANXIETY_SUMSTATS="anxiety_gwas.txt"
EA_SUMSTATS="educational_attainment_gwas.txt"
INT_SUMSTATS="intelligence_gwas.txt"

# ==============================================================================
# STEP 1: MUNGE SUMMARY STATISTICS
# ==============================================================================

echo "=== Step 1: Munging Summary Statistics ==="

# Munge cleft GWAS
python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${CLEFT_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_cleft \
    --merge-alleles ${HM3_SNPS}

# Munge neurodevelopmental/psychiatric GWAS
python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${ADHD_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_ADHD \
    --merge-alleles ${HM3_SNPS}

python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${AUTISM_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_ASD \
    --merge-alleles ${HM3_SNPS}

python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${SCZ_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_SCZ \
    --merge-alleles ${HM3_SNPS}

python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${BIPOLAR_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_BIP \
    --merge-alleles ${HM3_SNPS}

python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${DEPRESSION_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_DEP \
    --merge-alleles ${HM3_SNPS}

python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${ANXIETY_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_ANX \
    --merge-alleles ${HM3_SNPS}

python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${EA_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_EA \
    --merge-alleles ${HM3_SNPS}

python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats ${INT_SUMSTATS} \
    --out ${OUT_DIR}/MUNGE_INT \
    --merge-alleles ${HM3_SNPS}

# ==============================================================================
# STEP 2: GENETIC CORRELATION ANALYSES
# ==============================================================================

echo "=== Step 2: Running Genetic Correlation Analyses ==="

# Cleft vs ADHD
python ${LDSC_DIR}/ldsc.py \
    --rg ${OUT_DIR}/MUNGE_cleft.sumstats.gz,${OUT_DIR}/MUNGE_ADHD.sumstats.gz \
    --ref-ld-chr ${LD_REF} \
    --w-ld-chr ${LD_REF} \
    --out ${OUT_DIR}/CL_ADHD

# Cleft vs Autism
python ${LDSC_DIR}/ldsc.py \
    --rg ${OUT_DIR}/MUNGE_cleft.sumstats.gz,${OUT_DIR}/MUNGE_ASD.sumstats.gz \
    --ref-ld-chr ${LD_REF} \
    --w-ld-chr ${LD_REF} \
    --out ${OUT_DIR}/CL_ASD

# Cleft vs Schizophrenia
python ${LDSC_DIR}/ldsc.py \
    --rg ${OUT_DIR}/MUNGE_cleft.sumstats.gz,${OUT_DIR}/MUNGE_SCZ.sumstats.gz \
    --ref-ld-chr ${LD_REF} \
    --w-ld-chr ${LD_REF} \
    --out ${OUT_DIR}/CL_SCZ

# Cleft vs Bipolar
python ${LDSC_DIR}/ldsc.py \
    --rg ${OUT_DIR}/MUNGE_cleft.sumstats.gz,${OUT_DIR}/MUNGE_BIP.sumstats.gz \
    --ref-ld-chr ${LD_REF} \
    --w-ld-chr ${LD_REF} \
    --out ${OUT_DIR}/CL_BIP

# Cleft vs Depression
python ${LDSC_DIR}/ldsc.py \
    --rg ${OUT_DIR}/MUNGE_cleft.sumstats.gz,${OUT_DIR}/MUNGE_DEP.sumstats.gz \
    --ref-ld-chr ${LD_REF} \
    --w-ld-chr ${LD_REF} \
    --out ${OUT_DIR}/CL_DEP

# Cleft vs Anxiety
python ${LDSC_DIR}/ldsc.py \
    --rg ${OUT_DIR}/MUNGE_cleft.sumstats.gz,${OUT_DIR}/MUNGE_ANX.sumstats.gz \
    --ref-ld-chr ${LD_REF} \
    --w-ld-chr ${LD_REF} \
    --out ${OUT_DIR}/CL_ANX

# Cleft vs Educational Attainment
python ${LDSC_DIR}/ldsc.py \
    --rg ${OUT_DIR}/MUNGE_cleft.sumstats.gz,${OUT_DIR}/MUNGE_EA.sumstats.gz \
    --ref-ld-chr ${LD_REF} \
    --w-ld-chr ${LD_REF} \
    --out ${OUT_DIR}/CL_EA

# Cleft vs Intelligence
python ${LDSC_DIR}/ldsc.py \
    --rg ${OUT_DIR}/MUNGE_cleft.sumstats.gz,${OUT_DIR}/MUNGE_INT.sumstats.gz \
    --ref-ld-chr ${LD_REF} \
    --w-ld-chr ${LD_REF} \
    --out ${OUT_DIR}/CL_INT

echo "=== LDSC Analysis Complete ==="
echo "Results saved to: ${OUT_DIR}"
echo "Process .log files with 01_ldsc_analysis.R"
