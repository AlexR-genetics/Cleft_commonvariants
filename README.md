# Common Variant Contributions to Neurodevelopmental Risk in Orofacial Clefts


Analysis code accompanying the manuscript:

> **Common Variant Contributions to Neurodevelopmental Risk in Orofacial Clefts**  
> [Author list]  
> [Journal, 2026]

## Overview

This repository contains the analysis scripts used to investigate whether common genetic variants contribute to neurodevelopmental difficulties in children with cleft lip and/or palate (CL/P). The study examines:

1. **Genetic correlations** between CL/P and neurodevelopmental traits using LD Score Regression (LDSC)
2. **Polygenic risk score associations** comparing children with CL/P to population controls
3. **Within-cleft associations** of PRS with behavioural and developmental outcomes
4. **Cleft subtype comparisons** across different cleft phenotypes
5. **ND-CNV carrier analyses** comparing common variant burden between carriers and non-carriers
6. **Mendelian randomization** testing causal effects of genetic liability to CL/P on neurodevelopmental outcomes

## Repository Structure

```
├── scripts/
│   ├── 00a_ldsc_commands.sh        # Server commands for LDSC analysis
│   ├── 00b_prsice_commands.sh      # Server commands for PRS calculation
│   ├── 01_ldsc_analysis.R          # Process LDSC output and visualisation
│   ├── 02_prs_analysis.R           # PRS association analyses
│   └── 03_mendelian_randomization.R # Two-sample MR analysis
├── LICENSE
└── README.md
```

## Analysis Pipeline

### Prerequisites

**Software:**
- [LDSC](https://github.com/bulik/ldsc) (LD Score Regression)
- [PRSice-2](https://www.prsice.info/) (Polygenic Risk Score calculation)
- R (≥4.0) with packages:
  - `tidyverse`, `data.table`, `ggplot2`
  - `TwoSampleMR`, `MRPRESSO` (for MR analysis)
  - `fmsb`, `officer`, `flextable` (for PRS analysis)
  - `pheatmap`, `corrplot` (for LDSC visualisation)

**Data requirements:**
- Quality-controlled genotype data (PLINK format)
- GWAS summary statistics for neurodevelopmental traits
- Phenotype data with cleft classification and developmental measures

### Running the Analyses

1. **LDSC Genetic Correlations**
   ```bash
   # Run on server/HPC
   bash scripts/00a_ldsc_commands.sh
   
   # Process results locally
   Rscript scripts/01_ldsc_analysis.R
   ```

2. **Polygenic Risk Score Calculation**
   ```bash
   # Run on server/HPC
   bash scripts/00b_prsice_commands.sh
   ```

3. **PRS Association Analyses**
   ```r
   # Update file paths in script, then run:
   Rscript scripts/02_prs_analysis.R
   ```

4. **Mendelian Randomization**
   ```r
   # Update file paths in script, then run:
   Rscript scripts/03_mendelian_randomization.R
   ```

## GWAS Summary Statistics

Discovery GWAS summary statistics were obtained from:

| Trait | Reference |
|-------|-----------|
| ADHD | Demontis et al. (2023) |
| Autism spectrum disorder | Grove et al. (2019) |
| Schizophrenia | Trubetskoy et al. (2022) |
| Bipolar disorder | O'Connell et al. (2025) |
| Depression | Als et al. (2023) |
| Anxiety disorders | Purves et al. (2020) |
| Educational attainment | Okbay et al. (2022) |
| Intelligence | Savage et al. (2018) |
| Cleft lip/palate | Dardani et al. (2020) |

## Data Availability

Individual-level data from the Cleft Collective cannot be shared publicly due to ethical restrictions. Access can be requested through the [Cleft Collective](https://www.bristol.ac.uk/cleft-collective/) data access committee.

Control data from the Millennium Cohort Study is available through the [UK Data Service](https://ukdataservice.ac.uk/).

## Citation

If you use this code, please cite:

```
[Citation details to be added upon publication]
```

## Contact

For questions about the analysis code, please open an issue on this repository or contact [vd18986@bristol.ac.uk].

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
