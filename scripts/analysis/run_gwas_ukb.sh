#!/bin/bash

CHR_START=1
CHR_END=22

imputed_plink="imputed_plink/"
nonimputed_plink="genotype/merged/"

# QC on imputed files
for ((CHR = $CHR_START; CHR <= $CHR_END; CHR++))
do
    # Get list only biallelic variants
    python getBiAllelicVariants.py "$imputed_plink"ukb_imp_chr${CHR}_v3.pvar ukb_imp_chr${CHR}_biallelic_snps.set

    plink2 --pfile "$imputed_plink"ukb_imp_chr${CHR}_v3 \
        --keep samples_caucasian_unrelated.txt \
        --extract bed1 ukb_imp_chr${CHR}_biallelic_snps.set \
        --threads 16 \
        --geno 0.1 \
        --maf 0.001 \
        --max-maf 0.99 \
        --hwe 10e-6 \
        --make-just-pvar \
        --out centralized_imp_snp_qc_chr${CHR}
    
    echo "Finished QC on chromosome"${CHR}
done

# Concatenate all of the postQC-SNP files
head -n 1 centralized_imp_snp_qc_chr1.pvar > centralized_imp_snp_qc_all.pvar
awk 'FNR>1' centralized_imp_snp_qc_chr{1..22}.pvar >> centralized_imp_snp_qc_all.pvar
echo "Finished concatenating post-QC SNP files from all chromosomes"

# Combine post-QC SNPs into a .set file 
# (one that we can use with --extract flag in PLINK2 for filtering SNPs by chromosome and position)
python getSnpFilterByPositionFile.py  centralized_imp_snp_qc_all.pvar  centralized_imp_snp_qc_all.set
echo "Finished creating .set file"

# Extract QC-filtered SNPs from non-imputed dataset and R2 (LD) filter
# VIF threshold settings used: window size=100kb, stepsize=1, R2 < 0.1
plink2 --pfile "$nonimputed_plink"ukb_cal_all \
        --keep samples_caucasian_unrelated.txt \
        --extract bed1 centralized_imp_snp_qc_all.set \
        --threads 16 \
        --indep-pairwise 100kb 1 0.1 \
        --out centralized_nonimp_snps_forpca_R2

# Run PCA on non-imputed fileset with R2-filtered SNPs
plink2 --pfile "$nonimputed_plink"ukb_cal_all \
        --keep samples_caucasian_unrelated.txt \
        --extract centralized_nonimp_snps_forpca_R2.prune.in \
        --threads 16 \
        --pca approx 10 \
        --out centralized_pca_R2_forGWAS 

# Merge PCA results with phenotype, covariate data
# samples_caucasian_unrelated_phenoCov.txt includes ONLY individuals without ANY missing covariates/phenotypes 
python mergePhenoCovWithPCsUkb.py centralized_phenoCov.txt centralized_pca_R2_forGWAS.eigenvec centralized_phenoCovPCs.txt

# Run GWAS with covariates on each chromosome with R2-filtered PCA results
# include sex from .psam file (in --linear flag)
for ((CHR = $CHR_START; CHR <= $CHR_END; CHR++))
do
    plink2 --pfile "$imputed_plink"ukb_imp_chr${CHR}_v3 \
            --keep samples_caucasian_unrelated.txt \
            --extract centralized_imp_snp_qc_chr${CHR}.pvar \
            --threads 16 \
            --pheno centralized_phenoCovPCs.txt \
            --pheno-name BMI \
            --covar centralized_phenoCovPCs.txt \
            --covar-name Age Sex Center Age_Squared Age_Sex Age_Squared_Sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
            --variance-standardize Age Age_Squared Age_Sex Age_Squared_Sex BMI PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
            --linear hide-covar \
            --vif 99999999 \
            --ci 0.95 \
            --out centralized_imp_assoc_pc_R2_chr${CHR}
    echo "Finished GWAS on chromosome"${CHR}
done