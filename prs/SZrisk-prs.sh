# bash script for generating SZ risk scores in NAPLS
# use dosage data
# 0-1 probabilties, 2 columns per subject (format=2)
# list of all the files

# generated SZ risk score lists using 
# /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/preprocessing.R 
# 2014 PGC results

# flags: header; sum score (default = average)

threshold=( sigsnps pTall pT5 pT2 pT1 pT05 pT01 pT001 )
for i in "${threshold[@]}"
do
plink --dosage /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/dosagelist.txt list format=2 --fam dos_scz_scz1_mix_la-qc.hg19.ch.fl.chr9_141_144.out.dosage.fam --score /data/swe_gwas/ABZ/PGC-Results/pTlists/$i-prs.txt header sum --out /data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/prs/napls_SZrisk_$i
done