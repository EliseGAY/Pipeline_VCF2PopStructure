# Population Structure Analysis with Whole Genome Sequencing Data

## Overview
These script provides a global analysis of population structure using VCF data. The main analyses performed include:

`Filter_vcf.r : `

- **Depth, genotype and missing data** across all individuals and chromosomes.

`Structure_pop.r : `

- **Principal Component Analysis (PCA)** across all chromosomes
- **F-statistics (FST) computation** for population differentiation
- **sNMF analysis** for ancestry inference
- 
###TODO `Div_indices.r : `
- **Site Frequency Spectrum (SFS) computation**
- **Diversity analysis (DIV)**


## Dependencies
This script requires the following R libraries:

```r
library(LEA)
library(vcfR)
library(hierfstat)
library(pegas)
library(adegenet)
library(SNPRelate)
library(gdsfmt)
library(rlist)
library(ggplot2)
library(withr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(pcadapt)
library(dplyr)
library(stringr)
library(usethis)
library(devtools)
system("git clone https://github.com/EliseGAY/Package_VCF2PopStructure.git")
load_all("../VCF2PopStructure/")```

Ensure these packages are installed before running the script.

---

## Input Data
### 1. Population Metadata
A metadata table containing sample IDs and their associated population labels is required.

**Example format (metadata/Samples_table.txt):**
```plaintext
samples  pop  
sample_1 pop1 
sample_2 pop1 
sample_3 pop1 
sample_4 pop2 
sample_5 pop2 
```

### 2. VCF Files
VCF files containing genotype data should be provided for analysis.

---
