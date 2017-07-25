# bash script to run pca for within-sample stratification
# napls 'caucasian' subjects; N = 609
# run using these subjs which excludes related individuals
# default is 20 PCs

plink --bfile /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/pca.pruned.data/napls2_merge_pca_samples.menv --keep /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/napls_cauc.txt --make-bed --out /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/pca.pruned.data/pca_napls_cauc
plink --bfile /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/pca.pruned.data/pca_napls_cauc --pca header tabs --out /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/pca.pruned.data/napls_cauc_PCs