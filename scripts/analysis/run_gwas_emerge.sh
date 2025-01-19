#!/bin/bash

CHR_START=1
CHR_END=22

mkdir -p qc
mkdir -p pca

for ((CHR = $CHR_START; CHR <= $CHR_END; CHR++))
do
    # QC on SNPs
    plink2 --pfile centralized_chr${CHR} \
        --threads 16 \
        --geno 0.1 \
        --maf 0.05 \
        --max-maf 0.95 \
        --hwe 10e-6 \
        --make-just-pvar \
        --out qc/centralized_chr${CHR}_postQC 

    # Apply R2 (LD) threshold
    plink2 --pfile centralized_chr${CHR} \
        --extract qc/centralized_chr${CHR}_postQC.pvar \
        --indep-pairwise 100 kb 1 0.1 \
        --threads 16 \
        --out pca/centralized_chr${CHR}_snpsForPcaR2

    # Select SNPs for PCA
    plink2 --pfile centralized_chr${CHR} \
        --extract pca/centralized_chr${CHR}_snpsForPcaR2.prune.in \
        --threads 16 \
        --make-pgen \
        --out pca/centralized_chr${CHR}_pcaInput_r2
done

# Merge chromosomes together
rm -f -- plink_merge_filelist.txt
touch plink_merge_filelist.txt
for ((CHR = $CHR_START; CHR <= $CHR_END; CHR++))
do
  echo pca/centralized_chr${CHR}_pcaInput_r2 >> plink_merge_filelist.txt
done
plink2 --pmerge-list plink_merge_filelist.txt --out pca/centralized_all_pcaInput_r2

# Run PCA 
plink2 --pfile pca/centralized_all_pcaInput_r2 \
        --threads 16 \
        --pca approx 10 \
        --out pca/centralized_pca_r2 

# Merge PCA results with phenotype, covariate data
python mergePhenoCovWithPCsEmerge.py centralized_phenoCov.txt pca/centralized_pca_r2.eigenvec centralized_phenoCovPCs.txt

# Run GWAS with covariates on each chromosome
for ((CHR = $CHR_START; CHR <= $CHR_END; CHR++))
do
    plink2 --pfile centralized_chr${CHR} \
            --extract qc/centralized_chr${CHR}_postQC.pvar \
            --threads 16 \
            --pheno centralized_phenoCovPCs.txt \
            --pheno-name BMI \
            --covar centralized_phenoCovPCs.txt \
            --covar-name AGE SEX AGE_SQUARED CENTER PC1 PC2 PC3 PC4 PC5 \
            --variance-standardize AGE AGE_SQUARED PC1 PC2 PC3 PC4 PC5 BMI \
            --linear hide-covar \
            --vif 99999999 \
            --ci 0.95 \
            --out centralized_assoc_r2_chr${CHR}    
    echo "Finished GWAS on chromosome"${CHR}
done