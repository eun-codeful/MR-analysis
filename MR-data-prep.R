# loading package 
install.packages("data.table")
install.packages("MendelianRandomization")
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR@0.4.26")

library(data.table) 
library(dplyr)
library(readxl)
library(ggplot2)
library(MendelianRandomization)
library(TwoSampleMR) 

# Read in dementia summary statistics
Dem <- fread("C:/Users/ssani/Downloads/35379992-GCST90027158-MONDO_0004975.h.tsv") 
# Read in GDF15 summary statistics
GDF15<- fread("C:/Users/ssani/OneDrive/Documents/MR analysis/data/4374_45_GDF15_MIC_1.txt") 

colnames(Dem)
colnames(GDF15)

# subset GWAS significant threshold (p< 5*10^(-8))
protein_sig <- GDF15[GDF15$Pval < 5*10^(-8), ]

### To get the largest number of IVs, filter the protein_sig in outcome GWAS before clumping

protein_sig2 <- subset(protein_sig, protein_sig$SNP %in% Dem$SNP)


### Clump the Significant SNPs 
exp <- format_data( protein_sig2,
                    type = "exposure", 
                    phenotype_col = "phenotype",
                    header = TRUE, 
                    snp_col = "hm_rsid", 
                    beta_col = "hm_beta", 
                    se_col = "standard_error", 
                    effect_allele_col = "hm_effect_allele", 
                    eaf_col = "hm_effect_allele_frequency",
                    
                    other_allele_col = "hm_other_allele", 
                    pval_col = "p_value", 
                    samplesize_col = "n_con", 
                    chr_col = "chromosome", 
                    pos_col = "hm_pos", 
)


MR_exposure <- clump_data(
  exp,
  clump_kb = 500,
  clump_r2 = 0.01,
  clump_p1 = 5e-08,
  clump_p2 = 1,
  pop = "EUR"
)


write.csv(MR_exposure, "MR_exposure.csv")### Save your exposure and outcome MR files. 