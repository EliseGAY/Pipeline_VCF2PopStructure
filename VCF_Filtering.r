#!/usr/bin/env Rscript

#==========================================================#
# Filters VCF files by Depth , Missing rate or Heterozygous 
#==========================================================#

#===============#
#===============#
# Load Libraries 
#===============#
#===============#
library(ggplot2)
library(reshape2)
library(gridExtra)
library(stringr)
library(vcfR)
library(usethis)
library(devtools)
system("git clone https://github.com/EliseGAY/Package_VCF2PopStructure.git")

# Load functions to apply filters on genotype and depth
#------------------------------------------------------#
# for filters on sequencing depth 
load_all("Package_VCF2PopStructure/")

#===============================#
#===============================#
# Input VCF Files + Description #
#===============================#
#===============================#

#-----------------------------------------------------------------------------------------------#
# VCF already filters for lowqual, indel and repeat. 
# Made from GATK filters on HPC account
#-----------------------------------------------------------------------------------------------#
VCF=read.vcfR("data/VCF_example.vcf.gz")

# Read metadata 
#----------------#
metadata=read.table("metadata/metadata.txt", header=TRUE)

#------------------------------------------------------------#
#------------------------------------------------------------#
# 1) Sequencing depth rate visualization and filter ----
#------------------------------------------------------------#
#------------------------------------------------------------#

# a. Extract DP table for all position and all individuals ----
#--------------------------------------------------------------#

# get the position sequencing depth
DP<- extract.gt(VCF, element='DP', as.numeric = TRUE) 

# get mean per column (individual)
DPmean_bySample = apply(DP,2,mean,na.rm=T)
DPmean_bySample_tb=as.data.frame(DPmean_bySample)

# get mean per line (position)
DPmean_byPos= apply(DP,1,mean,na.rm=T)
DPmean_byPos_tb=as.data.frame(DPmean_byPos)

# Get the position vector in numeric
str_list=unlist(strsplit(rownames(DP)[1], '_'))
chrname=paste0(paste(str_list[-length(str_list)], collapse = "_"), "_")
position_vector=as.numeric(str_remove(rownames(DP), chrname))

# Add vector of position in the sequencing depth table 
DPmean_byPos_tb=cbind(DPmean_byPos_tb, position_vector)

# b. PLOT detailed DP for each ind or position ----
#-----------------------------------------------------#

# plot DP by positions
ggplot()+
  geom_point(data=DPmean_byPos_tb, 
             aes(y=DPmean_byPos_tb$DPmean_byPos,
                 x=DPmean_byPos_tb$position_vector),
             color="black", pch=20, size=1) +
  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y="mean depth", x = "PB")+
  ggtitle("Average depth sequencing of SNPs")+
  scale_x_continuous(breaks = seq(0, 18e+06, by = 1000000)) +
  # vizualise only range of sequencing depth detph
  scale_y_continuous(breaks = seq(0, 200, by = 10), limits = c(0,200))

# boxplot DP by individuals
# format table 
DP_melt=melt(DP)

ggplot()+
  geom_boxplot(data=DP_melt, aes(x=DP_melt$Var2, 
                                 y=DP_melt$value,
                                 colour="purple"),
               color="black") +
  theme(axis.title.x=element_blank())+
  labs(y="SNP sequencing depth", x = "samples")+
  ggtitle("depth") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # vizualise only range of sequencing depth detph
  scale_y_continuous(breaks = seq(0, 260, by = 20), limits = c(0, 260))

# FILTER FOR DP 
help(filter_vcf_by_dp)
VCF_Cata_DPFilter = filter_vcf_by_dp(vcf = VCF, 
                                     min_dp_pos = 10, max_dp_pos = 100, 
                                     min_dp_ind = 5, max_dp_ind = 200)

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#  2) Vizualisation of chosen genotype frequency or NA ----
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

# a. Extract genotype ----
#----------------------#
geno_table<-extract.gt(VCF_Cata_DPFilter,element="GT",mask=F,as.numeric=F,return.alleles = F,
                 convertNA = F,extract = T)

# b. Raw Count of genotype ----
#------------------------------#
# In rowsum and colsum the number of '.' '1' and '0'and compute genotype count 

Ploidy=2  #or 1
GT = "NA" # Het, NA

if(Ploidy == 1 && GT =="NA") { gt="." 
}else if(Ploidy == 1 && GT =="Het") { gt="1" 
}else if(Ploidy == 2 && GT =="NA") {gt="./."
}else if(Ploidy == 2 && GT =="Het") {gt="0/1"}

# By position
sum_pos_geno_table = as.data.frame(rowSums(geno_table==gt))
# By individuals
sum_ind_geno_table = as.data.frame(colSums(geno_table==gt))

# c. Genotype Frequencies ----
#-----------------------------#

# By position
Prop_Pos_geno=((sum_pos_geno_table/dim(geno_table)[2])*100)
Prop_Pos_geno_table=as.data.frame(Prop_Pos_geno)

# Get the position vector in numeric
str_list=unlist(strsplit(rownames(geno_table)[1], '_'))
chrname=paste0(paste(str_list[-length(str_list)], collapse = "_"), "_")
position_vector=as.numeric(str_remove(rownames(geno_table), chrname))

# Add vector of position in the sequencing depth table 
Prop_Pos_geno_table_pos=cbind(Prop_Pos_geno_table, position_vector)

# By individuals
Prop_Ind_geno=((sum_ind_geno_table/dim(geno_table)[1])*100)
Prop_Ind_geno_table=as.data.frame(Prop_Ind_geno)

#----------------------------------#
# d. PLOT genotype frequencies ----
#----------------------------------#

# By position
ggplot()+
  geom_point(data=Prop_Pos_geno_table_pos, 
             aes(y=Prop_Pos_geno_table_pos[,c(1)] , 
                 x=Prop_Pos_geno_table_pos$position_vector), 
             color="purple",
             size=1)+
  labs(x="pos", y = paste("prop", gt)) +
  ggtitle("Prop gt site by position") +
  scale_x_continuous(breaks = seq(0, position_vector[length(position_vector)], 
                                  by = 1000000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# By individuals
ggplot()+
  geom_point(data=Prop_Ind_geno_table, 
             aes(y=Prop_Ind_geno_table[,c(1)] , 
                 x=row.names(Prop_Ind_geno_table)),
             color="purple",
             size=2) +
  labs(x="Ind", y =  paste("prop", gt)) +
  ggtitle("Prop gt site by individuals") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))

# FILTER FOR NA RATE
help(Filters_na_Ind_Pos)
vcf_Na = Filters_na_Ind_Pos(vcf = VCF_Cata_DPFilter, rate_na_max_POS = 20, rate_na_max_Ind = 20, ploidy = 2)

# FILTER for HET RATE : if needed
vcf_P1_het = Filters_Het_Ind_Pos(vcf = VCF_Cata_116, rate_het_max_POS = 80, rate_het_max_Ind = 10, ploidy = 2)
