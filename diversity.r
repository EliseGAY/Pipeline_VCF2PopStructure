#==============#
#==============#
# ---- SFS ----
#==============#
#==============#

# ---- with VCF2PopStructure ---- 

# INPUT 
#-------#
# VCF with required filters (generally : Na = 0%, Heterozygous rate max = 80%)
# Do this with as following if not done already :
#   `VCF_0Na = Filters_na_Ind_Pos(vcf = VCF_DP, rate_na_max_POS = 0.0, rate_na_max_Ind = 0.0, ploidy = 2)`
#     write.vcf(VCF_0Na, "data/VCF_example.NoNa.SNP.vcf.gz")

# Methods : 
#-------#
# Get a genotype table encoded with 0,1,2 genotypes (position in row and samples in column). See : Convert_GT function
# Compute folded SFS with the function Get_FoldedSFS()

# Extract genotype with vcfR
loci_table = extract.gt(VCF_0Na, element = "GT")
loci_table = as.data.frame(loci_table)
colnames(loci_table)

# Convert loci table
loci_table_convert=Convert_GT(loci_table)
dim(loci_table_convert)

# Get SFS folded
SFS_folded=Get_FoldedSFS(loci_table_Test)

# Normalize by nb of position
SFS_folded_norm=unlist(SFS_folded)/sum(unlist(SFS_folded))
# Plot
ggplot() +
  geom_point(aes(y=SFS_folded_norm, x=seq(1,31)), color="blue")

# ---- with PEGAS ---- 

# METHODS
#------------------#
# PEGAS function : site.spectrum. Create SFS on DNAbin sequences

# OUTPUT
#------------------#
# SFS vector

# STEP1 : Create DNA sequences
DNAbin <-vcfR2DNAbin(VCF_0Na,  extract.indels = TRUE,
                     consensus = FALSE,
                     extract.haps = F,
                     unphased_as_NA = F,
                     asterisk_as_del = FALSE,
                     ref.seq = NULL,
                     start.pos = NULL,
                     verbose = TRUE)

# Compute folded site frequency spectrum with the function site.spectrum (library pegas required)
sfs_folded_GWS<-site.spectrum(DNAbin, folded=T)
sfs_folded_GWS
# [1] 25810 12912  8908  6501  5006  3969  3148  2772  2340  2079  1934  1724  1473  1450  1400  1325  1221  1339  1203  1115  1137  1225  1090  1054  1039  1038  1082   971   962   993   443
# attr(,"class")
# [1] "spectrum"
# attr(,"sample.size")
# [1] 62
# attr(,"folded")
#[1] TRUE
# sum(sfs_folded_GWS)

# Normalized by the spectrum by the number of sites 
# ARGUMENT : Folded SFS
tot_pos =  sum(sfs_folded_GWS)
sfs_folded_norm = (sfs_folded_GWS / tot_pos)
sum(sfs_folded_norm)

# Plot 
ggplot() +
  geom_point(aes(y=sfs_folded_norm, x=seq(1,31)), color="blue")


#==============================#
#==============================#
# ---- DIVERSITY ----
#==============================#
#==============================#
# Briefly 
# θ : Mutation parameters 
# θ = 4*Ne*μ 

# Computation of estimators of Theta (expected and observed)
# π  = nucleotide diversity 
# Under neutral and HW pop π(obs) == E(π)

# Expected π : E(π) = θS (Watterson’s θ) 
#   θS = S /  Σ(1/ 2*nbsample) / nb_total_site
#   Where S = Nb tot segregating site
#         nb_total_site = scale to all the site 
#   
# Estimator π = θπ :
#  θπ =  K / nb_total_site 
#  where K = total allele pairwise diff / nb_pairwise_chr
#  When pop is large : π≈2p(1−p)

# Tajima D compare E(π) and θπ

#==============================#
#==============================#

# ---- with PEGAS  ---- 

# INPUT 
#-------#
# VCF with required filters (generally : Na = 0%, Heterozygous rate max = 80%)

# METHODS PEGAS
#------------------#
# Assume a panmictic population. 
# The number of differentiated pairs are higher when the variant is mid-frequent, it contributes higher to ThetaPi.
# Reversely the singletons contribute lower to thetaPi, if there is a lot of singleton thetaPi will stay low.
# ThetaS uses brutally the nb of difference (whatever the allele frequency in the pop) so it is not weigthed by alt freq.
# So ThetaS >> ThetaPi TD negative means that we have an excess of singleton (which lower the ThetaPi)
# If ThetaS << ThetaPi TD positive means that 

# initiate nb of site
Nb_seg_site = length(seg.sites(DNAbin))
Nb_total_site = 5000000

# Format the alt allele frequencies per site
n_ind <- ncol(loci_table_convert)
p_by_variant <- rowSums(loci_table_convert) / n_ind
theta.obs.variant = theta.h(p_by_variant)
# 0.7465834
# scale for all position (produit en croix) :
theta.obs_by_site = theta.obs.variant*(Nb_seg_site/Nb_total_site)
# 0.01473203

# Expected Theta
thetaS= theta.s(x = length(seg.sites(DNAbin)), n = 2*31)
# 21008.83
thetaS_bysite = thetaS / Nb_total_site
# 0.004201766

# Here the DNABin contain only variant site. Use the real 
TD=tajima.test(DNAbin)
# $D
# [1] -0.5470291

# $Pval.normal
# [1] 0.5843587

# $Pval.beta
# [1] 0.6240548

# ---- with HierFSTAT  ---- 
# Methods (not sure about the difference with PEGAS yet) : 
# It approximate the K by using 2pq formula. 
# By doing that the contribution of intermediate allele to Theta increase 

TajimaD.dosage(loci_table_convert)
# 2.382157
theta.Watt.dosage(loci_table_convert,L=Nb_total_site)
# 4.634872e-07
