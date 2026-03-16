#!/usr/bin/env Rscript

#========================================================================================#
#========================================================================================#
# Pop structure : 
# - Fstatistics
# - PCA
# - Admixture
# - SFS
# - Diversity
# - Imbreeding
# Aim : Run general analysis of population structure on a gVCF
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
library(comprehenr)
library(pheatmap)

# system("git clone https://github.com/EliseGAY/Package_VCF2PopStructure.git")
load_all("../Package_VCF2PopStructure/")

#===============================#
#===============================#
# ------ PREPARE YOUR DATA ----
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
metadata=read.table("metadata/Samples_table.txt", header = TRUE)
metadata=as.data.frame(metadata)
pop=unique(metadata$Population)
head(metadata)

# read chr length
table_chr=read.table("metadata/sc_length.txt", header = T)
# read chr length table
chr_len = table_chr[which(table_chr$scaffold == "ptg000007l"),]$length

# set chr name
chr = "CHR1"

# Read the VCF with vcfR :
VCFR_data=read.vcfR("data/VCF_example_subset.vcf.gz")

# create a pop sorted by VCF colnames
metadata_sorted <- metadata[match(colnames(VCFR_data@gt)[-1], metadata$sample),]
pop_list = split(metadata_sorted$sample, metadata_sorted$Population)

# create a pop vector
vec_pop <- metadata$Population[match(colnames(VCFR_data@gt)[-1], metadata$sample)]

# Convert to genind (adegenet and Hierfstat format)
genind_data <- vcfR2genind(VCFR_data)
# Convert genind → hierfstat data.frame. Encode genotypes with 0 (0/0), 1 (0/1), 2 (1/1)
HF_Cata <- genind2hierfstat(genind_data, pop = metadata_sorted$Population)
HF_Cata_dt = as.data.frame(HF_Cata)
HF_Cata_dt$pop <- as.numeric(as.factor(HF_Cata_dt$pop))
colnames(HF_Cata_dt)[1] = "population"
rownames(HF_Cata_dt) <- NULL

# Read the VCF with SNPRelate and convert to GDS
snpgdsVCF2GDS(vcf.fn = "data/VCF_example_subset.vcf.gz",
              out.fn = "data/VCF_example_subset.gds")

genofile_cata <- snpgdsOpen("data/VCF_example_subset.gds")
read.gdsn(index.gdsn(genofile_cata, "sample.id"))
read.gdsn(index.gdsn(genofile_cata, "snp.id"))

#===============#
#===============#
# ---- FST ----
#===============#
#===============#

#-----------------------------------------------------------------#
# F-statistics (Weir and Cockerham (1984)) : Pegas, Hierfstat

#   briefly :

#   FST_WC84 = a / (a + b + c)
#
#   a : p variance among populations
#   b : p variance among individuals within populations
#   c : p variance within individuals (heterozygosity)
#-----------------------------------------------------------------#

# --- Wtih Hierfstat ----

# input : vcf read by vcfR and transformed on an object of class genind
#         pop : vector of pop name (row)

# Compute Fstats
Fstat_HF_Cata = basic.stats(HF_Cata, diploid = 2)

# See what we've got :
head(Fstat_HF_Cata$perloc)

# HierFstat DOCUMENTATION
# n.ind.samp: A table --with np (number of populations) columns and nl (number of loci) rows-- of genotype counts
# pop.freq: A list containing allele frequencies. Each element of the list is one locus. For each locus, Populations are in columns and alleles in rows
# Ho: A table --with np (number of populations) columns and nl (number of loci) rows-- of observed heterozygosities
# Hs: A table --with np (number of populations) columns and nl (number of loci) rows-- of observed gene diversities
# Fis: A table --with np (number of populations) columns and nl (number of loci) rows--of observed Fis
# perloc: A table --with as many rows as loci-- containing basic statistics Ho, Hs, Ht, Dst, Ht', Dst', Fst, Fst' ,Fis, Dest
# overall: Basic statistics averaged over loci
#----------------------------------------------------------#

Fstat_HF_Cata$overall
#     Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
# 0.0730 0.0749 0.0974 0.0225 0.0999 0.0250 0.2307 0.2501 0.0258 0.0270 

# pairwise Fst
pairwisefst = pairwise.WCfst(HF_Cata)
pairwisefst

#             Pop1      Pop2      Pop3
# Pop1        NA 0.2144688 0.3322326
# Pop2 0.2144688        NA 0.3434954
# Pop3 0.3322326 0.3434954        NA

# --- With SNPRelate----

# input : vcf read by snpgdsVCF2GDS function of SNPRelate and transformed on an object of class GDS
#         pop : vector of pop name (row)

# select pair of pops manually to compute pairwise Fst
# Create a sub metadata table :
submetadata = metadata[metadata$Population %in% c("Pop1", "Pop2"),]

# get sub gds file to compute pairwise. Do it for all pairs
snpgdsCreateGenoSet(src.fn = "data/VCF_example_subset.gds",
                  dest.fn = "data/VCF_example_subset.pop1-2.vcf.gds",
                  sample.id = submetadata$sample,
                  verbose = TRUE)

genofile_TC=snpgdsOpen("data/VCF_example_subset.pop1-2.vcf.gds")
read.gdsn(index.gdsn(genofile_TC, "sample.id"))

# run Fst on subseted dataset
fst_snp <- snpgdsFst(maf = 0.05, 
                     gdsobj = genofile_TC, 
                     autosome.only = F,
                     population = as.factor(submetadata$Population), 
                     missing.rate = 0.20 ,
                    method = "W&C84")

snpgdsClose(genofile_TC)

fst_snp$Fst
#o.36

# Pop2 (9), Pop3 (15) :
# of SNPs: 20,138
# $MeanFst
# [1] 0.2823853


#-----------------------------------------------------------------#
# Hudson et al. (1992) — diversity-based (π)
#
#   FST_Hudson = (π_between − π_within) / π_between
#
#   π_within  : mean nucleotide diversity within populations
#   π_between : mean nucleotide divergence between populations
#-----------------------------------------------------------------#

#---- with Package VCR2PopStructure (home made)----

# Get genotype table
loci_table = extract.gt(VCFR_data, element = "GT")
loci_table = as.data.frame(loci_table)
colnames(loci_table)

# tranform in genotype table with Convert_GT function
loci_table_T_CV = Convert_GT(loci_table)

# Fst Whole SNPs
Fst_GlobalCount = getGlobalFst_Count(loci_table_T = loci_table_T_CV, 
                          pop_table = metadata, Na_rate = 0.20, MAF_threshold = 0.05)


Fst_GlobalCount
# $Pop2_Pop3
# [1] 0.3588299
#
# $Pop2_Pop1
# [1] 0.2480494
#
# $Pop3_Pop1
# [1] 0.3567909

# Permute the Fst to get significance :
# TO DO

#====================================# 
#====================================# 
# ------ LD and PRUNING ----
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

# do pruning on whole data, adapt your methods if you want
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

#----------#
# run ACP
#----------#

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
# ---- ADMIXTURE ----
#=====================#
#=====================#
# briefly : 
# Test to fit how much your genotype structure fits to K ancestries
# Cross entropy measures the -log(likelihood) to seeing your data given the Ki ancestry tested 
# The lower the cross-entropy the better the K

#---- with Snmf ----
 
# INPUT 
#-------#
# VCFs filtered for Na and DP on R (for NA and DP)

  #	==> MAF at 0.005 added with vcftools with the command line 
# `bcftools view  -i 'F_MISSING<=0.2 && MAF>=0.05' VCF_116_samples_TAG.Flowqual_Noindels_Norepeat.ptg000007l.SNP.vcf.gz --threads 8 -Oz -o VCF_116_samples_TAG.Flowqual_Noindels_Norepeat.ptg000007l.SNP.miss02.maf05.vcf.gz
# `bcftools index VCF_116_samples_TAG.Flowqual_Noindels_Norepeat.ptg000007l.SNP.miss02.maf05.vcf.gz`
  #	==> Ped is generated with vcftools with the command line 
# `vcftools --gzvcf VCF_116_samples_TAG.Flowqual_Noindels_Norepeat.ptg000007l.SNP.miss02.maf05.vcf.gz --out VCF_116_samples_TAG.Flowqual_Noindels_Norepeat.ptg000007l.SNP.miss02.maf05 --plink`

  #	==> Samples were ordered by population (formated file for snmf anlaysis)
  #	==> plink is ran on the ".ped" file needed to run snmf.

# list_pop = list of populations containing ID in each (list of list type) : create in metadata part

# METHODS
#------------------#
# snmf : "estimates admixture coefficients using sparse Non-Negative Matrix Factorization algorithms"
# Breifly Snmf : Test for K ancestries and cross validate (by minimizing the loss) the genotype prediction on masked (but known) genotype.

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
project = snmf("example.SNP.miss02.maf05.ped",
               K=1:8,
               entropy=T,
               repetitions = 20,
               project = "new")

# To re-load the project already created, use:
project = load.snmfProject("example.SNP.miss02.maf05.snmfProject")

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
  project = snmf("WGS_21_SUPER_1_TAG_Flowqual_Noindels_Norepeat_SNP_DP10_50_10_200_Na20_MAF005_order_Plink.ped",
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
  rownames(qmatrix)= metadata$sample # carefull of the samples order in list
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
# input : genofile
# Methods : 
# Briefly, the model supposed that ancestral state are hidden in the eigen value computed on the overall SNPs data. It computes the proportion of ancestries eigen value in each samples cooridnates.
# It claims to be unbiases by LD and to be faster than model based admixture.
# check here if you want more detail : "Eigenanalysis of SNP data with an identity by descent interpretation"

RV <- snpgdsEIGMIX(genofile_cata, autosome.only = F, maf = 0.5, missing.rate = 0.2)
RV$eigenval

# make a data.frame
tab <- data.frame(sample.id = read.gdsn(index.gdsn(genofile_cata, "sample.id")), pop = factor(metadata_sorted$Population),
                  EV1 = RV$eigenvect[,1],    # the first eigenvector
                  EV2 = RV$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
# draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
     xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:4)


# define groups
samp.id = metadata_sorted$sample
groups <- list(Pop1 = samp.id[metadata_sorted$Population == "Pop1"],
               Pop2 = samp.id[metadata_sorted$Population == "Pop2"],
               Pop3 = samp.id[metadata_sorted$Population == "Pop3"])

prop <- snpgdsAdmixProp(RV, groups=groups)

# draw
snpgdsAdmixPlot(prop, group=metadata_sorted$Population)




