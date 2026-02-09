#!/usr/bin/env Rscript

#========================================================================================#
#========================================================================================#
# Pop structure : 
  # - Fstatistics
  # - LD and pruning
  # - PCA
  # - Admixture
 
# Aim : Run general analysis to detect population structure on a gVCF
#========================================================================================#
#========================================================================================#

#============================#
#============================#
# ------ Load libraries ----
#============================#
#============================#
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
load_all("../VCF2PopStructure/")

system("git clone https://github.com/EliseGAY/Package_VCF2PopStructure.git")

#===============================#
#===============================#
# ------ Prepare your data ----
#===============================#
#===============================#

#--------------#
# Metadata pop
#--------------#

# the pop table has to be ordered in the same way as all the VCF header
#'''
#	pop	samples
#	pop1 sample_1
#	pop1 sample_2
#	pop1 sample_3
#	pop2 sample_4
#	pop2 sample_5
#'''

# read metada
metadata=read.table("metadata/metadata.txt", header = TRUE)
metadata=as.data.frame(metadata)
pop=unique(metadata$Population)

#-----------------------------------------------------------#
# Generate Genotype tables needed in different R packages
#-----------------------------------------------------------#

# Read the VCF with vcfR :
VCFR_data=read.vcfR("data/VCF_example.vcf.gz")

# create a pop vector
vec_pop <- metadata$Population[match(colnames(VCFR_data@gt)[-1], metadata$GT_sample)]


# Convert to genind (adegenet and Hierfstat format)
genind_data <- vcfR2genind(VCFR_data) # [110] "1761W1_S352" missing

# Genotype table for hierFstat : 
# Convert genind → hierfstat data.frame. Encode genotypes with 0 (0/0), 1 (0/1), 2 (1/1)
HF_Cata <- genind2hierfstat(genind_data, pop = vec_pop)
HF_Cata_dt = as.data.frame(HF_Cata)
colnames(HF_Cata_dt)[1] = "population"
rownames(HF_Cata_dt) = NULL

# Read the VCF with PEgas (read only 10 000 loci) :
VCFPegas=read.vcf("data/VCF_example.vcf.gz")

# Read the VCF with SNPRealte and convert to GDS
snpgdsVCF2GDS(
  vcf.fn = "data/VCF_example.vcf.gz",
  out.fn = "data/VCF_example.SNP.gds")

genofile_cata <- snpgdsOpen("data/VCF_example.SNP.gds")
read.gdsn(index.gdsn(genofile_cata, "sample.id"))

# get your vector pop in the same order as colnames in vcf

# reorder your metadata the same way : 
metadata <- metadata %>%
  dplyr::slice(match(read.gdsn(index.gdsn(genofile_cata, "sample.id")), metadata$GT_sample))

#===============================================#
#===============================================#
# ---- Global and pariwise FST ----
#===============================================#
#===============================================#

#-----------------------------------------------------------------#
# F-statistics (Weir and Cockerham (1984)) : Pegas, Hierfstat

#   briefly :

#   FST_WC84 = a / (a + b + c)
#
#   a : p variance among populations
#   b : p variance among individuals within populations
#   c : p variance within individuals (heterozygosity)
#-----------------------------------------------------------------#

# ---- With Pegas ---- 
# input : vcf read by pegas and transformed on an object of class "loci"
#         pop : vector of pop name (row)
# Methods : Fst, Fit, Fis computed on each locus (default first 10K loci) 

Fst_pegas = Fst(VCFPegas_data, pop = vec_pop, quiet = FALSE)
summary(Fst_pegas)

# --- Wtih Hierfstat ----
# input : vcf read by vcfR and transformed on an object of class genind
#         pop : vector of pop name (row)

# Compute Fstats
Fstat_HF_Cata = basic.stats(HF_Cata, diploid = 2)

# See what we've got :
Fstat_HF_Cata$overall
#     Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
# 0.0730 0.0749 0.0974 0.0225 0.0999 0.0250 0.2307 0.2501 0.0258 0.0270 

# Get variance component :
vc_cata = varcomp.glob(HF_Cata_dt)
vc_cata$a   # among populations
vc_cata$b   # among individuals within populations
vc_cata$c   # within individuals

# pairwise Fst
pairwisefst = pairwise.WCfst(HF_Cata)
pairwisefst

#         Cesseras  Cucugnan    Termes
# Cesseras        NA 0.3348978 0.3442045
# Cucugnan 0.3348978        NA 0.2151290
# Termes   0.3442045 0.2151290        NA

# --- With SNPRelate----
# input : vcf read by snpgdsVCF2GDS function of SNPRelate and transformed on an object of class GDS
#         pop : vector of pop name (row)

fst_snp <- snpgdsFst(maf = 0.05, 
                     gdsobj = genofile_cata, 
                     autosome.only = F,
                     population = as.factor(vec_pop), 
                     missing.rate = 0.20 ,
                     method = "W&C84")

# select pair of pops manually to compute pairwise Fst
# Create a sub metadata table :
submetadata = metadata[metadata$Population %in% c("Termes", "Cucugnan"),]
submetadata = metadata[metadata$Population %in% c("Termes", "Cesseras"),]
submetadata = metadata[metadata$Population %in% c("Cesseras", "Cucugnan"),]

# get sub gds file
snpgdsCreateGenoSet(src.fn = "data/VCF_example.SNP.gds",
                  dest.fn = "VCF_example.SNP.Cuc-Cess.vcf.gds",
                  sample.id = submetadata$GT_sample,
                  verbose = TRUE)

genofile_TC=snpgdsOpen("VCF_example.SNP.T-C.vcf.gds")
read.gdsn(index.gdsn(genofile_TC, "sample.id"))

genofile_TCess=snpgdsOpen("VCF_example.SNP.T-Cess.vcf.gds")
read.gdsn(index.gdsn(genofile_TCess, "sample.id"))

genofile_Cuc_Cess=snpgdsOpen("VCF_example.SNP.Cuc-Cess.vcf.gds")
read.gdsn(index.gdsn(genofile_Cuc_Cess, "sample.id"))

fst_snp <- snpgdsFst(maf = 0.05, 
                     gdsobj = genofile_Cuc_Cess, 
                     autosome.only = F,
                     population = as.factor(submetadata$Population), 
                     missing.rate = 0.20 ,
                     method = "W&C84")

snpgdsClose(genofile_TC)

# Cucugnan (12), Termes (9) :
# SNPs: 60,571
# $Fst
# [1] 0.1457925

# Cesseras (15), Termes (9)
# SNPs: 40,806
# > fst_snp
# $Fst
# [1] 0.3661034

# Cesseras (15), Cucugnan (12)
# SNPs: 57,223
# > fst_snp
# $Fst
# [1] 0.2628528

#-----------------------------------------------------------------#
# adegenet F - statistics (Reynolds)
#   
#   Reynolds et al. (1983) — drift-based genetic distance
#
#   FST_Reynolds = Σ (p_i − p_j)^2 / Σ [ p̄ (1 − p̄) ]
#
#   p_i, p_j : allele frequencies in populations i and j
#   p̄       : mean allele frequency
#-----------------------------------------------------------------#

# TODO

#-----------------------------------------------------------------#
# Hudson et al. (1992) — diversity-based (π)
#
#   FST_Hudson = (π_between − π_within) / π_between
#
#   π_within  : mean nucleotide diversity within populations
#   π_between : mean nucleotide divergence between populations
#-----------------------------------------------------------------#

# Get genotype table
loci_table = extract.gt(VCFR_data, element = "GT")
loci_table = as.data.frame(loci_table)
colnames(loci_table)

# tranform
loci_table_T_CV = Convert_GT(loci_table)

# If you want allel freq
getAlleleFreqByPop(loci_table = loci_table_T_CV, pop_table = metadata)

# test Fst computation on SNPs subset 
rowrandom = order(round(runif(n = 10000, min = 1, max = 201634), 0))
# Fst by SNPs
Fst_BySnps = getFstBySNP(loci_table_T = as.data.frame(loci_table_T_CV[rowrandom,, drop = F]), 
                       pop_table = metadata)
# Fst Whole SNPs
Fst_Global = getGlobalFst(loci_table_T = loci_table_T_CV[rowrandom, , drop = FALSE], 
                          pop_table = metadata)

# Do it on whole SNPs set
Fst_Global = getGlobalFst(loci_table_T = loci_table_T_CV, 
                          pop_table = metadata)
# results :
# $Cucugnan_Termes
# [1] 0.2813149
# 
# $Cucugnan_Cesseras
# [1] 0.3746937
# 
# $Termes_Cesseras
# [1] 0.3772484

# Permute the Fst to get significance :
# TO DO

#====================================# 
#====================================# 
# ------ LD and Pruning ----
# briefly  : 
# D=P(AB)−p(A)p(B)
# if D = 0 --> NO LD 
# if D != 0 --> LD
# Structure is taking into account in LD computation
# Remove SNP if there are correlated at a rate above "ld.threshold"
#====================================# 
#====================================#

# ---- with SNPRelate ---- 

# Input : a genofile 

# methods : 
#   LD.mat : you may want to compute a LD matrix at some point. Do it on a subset of SNPs
#   pruning : pruning function need a ld.threshold , the more your loose the threshold the more correlated SNPs are kept 

# output : set of 'unlinked' SNPs that you want to keep for some analyses

# prepare subset of SNPs
SNP_vec = read.gdsn(index.gdsn(genofile_cata, "snp.id"))
sub_SNP_vec = sample(SNP_vec, 10000)
rowrandom = order(round(runif(n = 500, min = 1, max = 10000), 0))

# test on whole samples
LD_Mat_Comp = snpgdsLDMat(gdsobj = genofile_cata, slide = -1, method = "composite", snp.id = sub_SNP_vec)
LD_Mat_Corr = snpgdsLDMat(gdsobj = genofile_cata, slide = -1, method = "corr", snp.id = sub_SNP_vec)

# subset to plot :
LD_Mat_Comp_sub = LD_Mat_Comp$LD[rowrandom,rowrandom]
LD_Mat_Corr_sub = LD_Mat_Corr$LD[rowrandom,rowrandom]

# vizualize :
image(t(LD_Mat_Comp_sub^2))
image(t(LD_Mat_Corr_sub^2))

# do pruning on whole data. 
LDcomposite_Pruning <- snpgdsLDpruning(gdsobj = genofile_cata, 
                          method = "composite", 
                          num.thread = 4, autosome.only  = F,
                          missing.rate = 0.2, ld.threshold = 0.8, maf = 0.05)

LDcorr_Pruning <- snpgdsLDpruning(gdsobj = genofile_cata, 
                          method = "corr", 
                          num.thread = 4, 
                          autosome.only  = F,
                          missing.rate = 0.2, ld.threshold = 0.8, maf = 0.05)

# get the SNP subset 
snpgdsCreateGenoSet(src.fn = "data/VCF_example.SNP.gds",
                    dest.fn = "data/VCF_example.SNP.PrunnedCorrMAF05LD01.gds",
                    snp.id = LDcorr_Pruning$chrptg000007l,
                    verbose = TRUE)

snpgdsCreateGenoSet(src.fn = "data/VCF_example.SNP.gds",
                    dest.fn = "data/VCF_example.SNP.PrunnedCompMAF05LD01.gds",
                    snp.id = LDcomposite_Pruning$chrptg000007l,
                    verbose = TRUE)

geno_comp_pruned = openfn.gds("data/VCF_example.SNP.PrunnedCompMAF05.gds")
geno_corr_pruned = openfn.gds("data/VCF_example.SNP.PrunnedCorrMAF05.gds")

#======================# 
#======================# 
# ------ PCA ----
#======================#
#======================#

# INPUT 
#-------#
# VCF with chosen filters : 20% Na and MAF filtered (MAF can also be put directly in the snpgdsPCA function)
# vec_pop = list of populations containing ID in each, created in metadata part

# OUTPUT
#---------#
# PCA plots according to  and population assignation, writen in your current workdir

genofile = genofile_cata # Or your subset samples or prunned SNP

# do the PCA 
pca_all= snpgdsPCA(genofile, num.thread = 4, autosome.only = FALSE )

# SET your PCA variable:
pca = pca_all
samples_list = read.gdsn(index.gdsn(genofile, "sample.id"))
population_list=vec_pop

# getinfo
pca$sample.id
pca$eigenval
pc.percent <- pca$varprop * 100   # percentage of variance explained

# make dataframe to plot
pca_df <- data.frame(PC1= pca$eigenvect[,1],
                     PC2 = pca$eigenvect[,2],
                     PC3 = pca$eigenvect[,3],
                     PC4 = pca$eigenvect[,4],
                     PC5 = pca$eigenvect[,5],
                     PC6 = pca$eigenvect[,6],
                     PC7 = pca$eigenvect[,7],
                     PC8 = pca$eigenvect[,8],
                     PC9 = pca$eigenvect[,9],
                     PC10 = pca$eigenvect[,10],
                     sample.id = pca$sample.id)

# Order PCA to match pop name
pca_df_sorted <- pca_df %>%
  dplyr::slice(match(samples_list, pca_df$sample.id))

pca_df_sorted_meta=cbind(pca_df_sorted, population_list)

# Direct plot and create png
getwd()
nb_axes=4

for (i in 1:nb_axes) {
  
  pc_x <- paste0("PC", i)
  pc_y <- paste0("PC", i + 1)
  
  # File name
  outfile <- paste0("PCA_", i, "_", i + 1, ".png")
  
  # Open PNG device (use cairo if needed)
  png(outfile, width = 1600, height = 1200, res = 150)
  
  print(
    ggplot(pca_df_sorted_meta, 
           aes_string(x = pc_x, y = pc_y, color = "vec_pop")) +
      geom_point(size = 3, alpha = 0.9) +
      geom_text(
         aes(label = samples_list),
         size = 2,
         vjust = -0.5,
         check_overlap = FALSE
       ) +
      labs(
        x = paste0(pc_x, " (", round(pc.percent[i], 2), "%)"),
        y = paste0(pc_y, " (", round(pc.percent[i + 1], 2), "%)"),
        color = "pop"
      ) +
      theme_minimal(base_size = 14)
  )
  
  dev.off()
}

#=====================#
#=====================#
# ---- Admixture ----
#=====================#
#=====================#
# briefly : 
# Test to fit how much your genotype structure fits to K ancestries (K proxi range is given by the PCA you already made).
# Cross entropy measures the -log(likelihood) to seeing your data given the Ki ancestry tested.
# The lower the cross-entropy the better the K

#---- with Snmf ----
# 
# INPUT 
#-------#
# VCFs filtered for Na and DP on R (for NA and DP)

  #	==> MAF at 0.005 added with vcftools with the command line 

# `bcftools view  -i 'F_MISSING<=0.2 && MAF>=0.05' example.SNP.vcf.gz --threads 8 -Oz -o example.SNP.miss02.maf05.vcf.gz
# `bcftools index example.SNP.miss02.maf05.vcf.gz`
#	==> Ped is generated with vcftools with the command line 
# `vcftools --gzvcf example.SNP.miss02.maf05.vcf.gz --out example.ptg000007l.SNP.miss02.maf05 --plink`

  #	==> Samples were ordered by population (formated file for snmf anlaysis)
  #	==> snmf is ran on the ".ped" file

# list_pop = list of populations containing ID in each (list of list type) : create in metadata part

# METHODS
#------------------#
# snmf : "estimates admixture coefficients using sparse Non-Negative Matrix Factorization algorithms"

# OUTPUT
#---------#
# table_all : contains all cross entropy for 20 rep in 100 snmf runs in each K.
# A folder is created  containing the snmf project with 
# For each k :
#   - produces a Q matrice = prop of each ind. in ancestry pops
#   - compute cross entropy n times

# Run sNMF
#----------#
# go in the directory. The absolute path has to be short, ortherwise snmf won't work
setwd("data/")
# create project with nb K + nb repetition chosen
project = snmf("example.SNP.miss02.maf05SUB.ped",
               K=1:8,
               entropy=T,
               repetitions = 20,
               project = "new")

# To re-load the project already created, use:
project = load.snmfProject("example.SNP.miss02.maf05SUB.snmfProject")

# project object :
project@runs
# plot cross-entropy criterion of all runs of the project
plot(project, cex = 2, col = "lightblue", 
     pch = 19, 
     xaxp=c(0,35,35))

# look at the cross entropy for the 20 runs with K=i
cross.entropy(project, K = 3)

# look at the Q matrix (samples = line, k pops = column) :
qmatrix = Q(project, K = 3,run = 5)
nrow(qmatrix)

# Run sNMF with intra and extra repetition
#--------------------------------------------#
# One hundred Runs of snmf on K = 8 , with and with 10 repetition

# initiate table for K = 8
table_all=data.frame("K=2"=numeric(), 
                      "K=3"=numeric(), 
                      "K=4"=numeric(),
                     "K=5"=numeric(), 
                     "K=6"=numeric(), 
                     "K=7"=numeric(),
                     "K=8"=numeric())
table_all

# To choose the best K, with repetition :
# Ten Runs of snmf on K = 8 and with 10 repetition each : 100 run of snmf
for (i in seq(1:10)){
  project = snmf("example.SNP.miss02.maf05.ped",
                 K=1:8,
                 entropy=T,
                 repetitions = 10,
                 project = "new")
  table_i=as.data.frame(cbind(cross.entropy(project, K = 2),
                cross.entropy(project, K = 3),
                cross.entropy(project, K = 4),
                cross.entropy(project, K = 5),
                cross.entropy(project, K = 6),
                cross.entropy(project, K = 7),
                cross.entropy(project, K = 8)))
  table_all=rbind(table_all, table_i)
}

# format table to plot
summary(table_all)

table_all_melt=melt(table_all)

# plot boxplot of cross entropy
ggplot()+
  geom_boxplot(aes(x=table_all_melt$variable, 
              y=table_all_melt$value))

# plot admixture
#----------------#
for (i in c(1,2,3,4)) {

  ce = cross.entropy(project, K = i)
  best = which.min(ce)
  qmatrix = Q(project, K = i,run = best)
  rownames(qmatrix)= metadata$GT_sample # carefull of the samples order in list
  meltedqmatrix = melt(qmatrix)
  
  #One color for one K
  my.colors=c("violet","dodgerblue","goldenrod", rainbow(3))

  
  p=plot(ggplot() +
         geom_bar(data=meltedqmatrix,aes(x=Var1,y=value,fill=Var2),stat="identity",
                  show.legend=F) +
         theme(text = element_blank(), 
               #axis.text = element_text(size = 7, angle = 90),
               axis.text.x=element_blank(),
               panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                               colour = "grey"),
               panel.background = element_rect(fill = "white",
                                               colour = "black",
                                               size = 0.5,
                                               linetype = 'solid'),
               axis.ticks = element_blank()) +
         xlab("Individuals") +
         ylab("Ancestry proportion") +
         ggtitle(paste("K = ",i,", Minimal cross-entropy = ",ce,sep="")) +
         scale_fill_manual(values=my.colors))+
    geom_vline(xintercept = c(47, 71, 105), size=1,linetype="dashed")
	
	name=paste(i, "plot.pdf", sep="")
	ggsave(p, filename = name, width = 10, height = 10, device = "pdf")

}

#---- with SNPRelate ----
RV <- snpgdsEIGMIX(genofile_cata, autosome.only = F, maf = 0.5, missing.rate = 0.2)
RV$eigenval

# make a data.frame
tab <- data.frame(sample.id = read.gdsn(index.gdsn(genofile_cata, "sample.id")), pop = factor(metadata$Population),
                  EV1 = RV$eigenvect[,1],    # the first eigenvector
                  EV2 = RV$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

# draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
     xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:4)


# define groups
samp.id = metadata$GT_sample
groups <- list(Termes = samp.id[metadata$Population == "Termes"],
               Cesseras = samp.id[metadata$Population == "Cesseras"],
               Cucugnan = samp.id[metadata$Population == "Cucugnan"])

prop <- snpgdsAdmixProp(RV, groups=groups)

# draw
snpgdsAdmixPlot(prop, group=metadata$Population)

