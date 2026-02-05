# Population Structure Analysis with Whole Genome Sequencing Data

## Overview
This script provides a global analysis of population structure using Whole Genome Sequencing (WGS) data. The main analyses performed include:

`Structure_pop.r : `

- **Principal Component Analysis (PCA)** across all chromosomes
- **F-statistics (FST) computation** for population differentiation
- **sNMF analysis** for ancestry inference
  
- **Site Frequency Spectrum (SFS) computation**
- **Diversity analysis (DIV)**

`Div_indices.r : `
# TODO 

## Dependencies
This script requires the following R libraries:

```r
library(LEA)
library(vcfR)
library(rlist)
library(PopGenome)
library(qvalue)
library(pegas)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(withr)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(pcadapt)```

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
Filtered VCF files containing genotype data should be provided for analysis.

---
