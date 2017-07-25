# PsychChip data from NAPLS - SNPs
# managing PCs and choosing subjects based on ancestry
# generation of SZ risk scores

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/"
setwd(workdir)

libs <- c("dplyr", "ggplot2")
invisible(lapply(libs, require, character.only = TRUE))


######################## PCs 
#load in data
napls   <- read.csv("napls2_pca/napls2_pca_components_remove_related.txt",
                    header=T, sep=" ", stringsAsFactors = F)
genomes <- read.csv("napls2_pca/napls2_1kg_pca_components_remove_related.txt",
                     header=T, sep=" ", stringsAsFactors = F)

# code 1000 genomes individuals by population; N = 1053
# code napls subjects by case/control; N = 798
# (basically convert FIDs to meaningful groups)
genomes$pop         <- substr(genomes$FID,1,12) %>% as.factor
levels(genomes$pop) <- c("cases","controls","african","american",
                         "asian","european")

# plot PCs
ggplot(genomes, aes(x=C2,y=C1)) + 
  geom_point(aes(colour=pop, shape=pop), size=3) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        panel.grid.minor = element_line(color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  ylab("PC1") +
  xlab("PC2")

# as far as I can tell people just exclude individuals
# based on visually mapping PCs so
# I'm gonna make a call on where the samples divide
# PC1 > -.02 = NOT african
# PC2 > -.03 = NOT asian
# leaves american + european 
# adding within-sample PCs anyway so it's fine right? 
all.cauc   <- subset(genomes,C1 > -.02 & C2 > -.03)
ggplot(all.cauc, aes(x=C2,y=C1)) + 
  geom_point(aes(colour=pop, shape=pop), size=3) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        panel.grid.minor = element_line(color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  ylab("PC1") +
  xlab("PC2")

napls.cauc <- subset(all.cauc, pop == "cases" | pop == "controls") # N = 609
ggplot(napls.cauc, aes(x=C2,y=C1)) + 
  geom_point(aes(colour=pop, shape=pop), size=3) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        panel.grid.minor = element_line(color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  ylab("PC1") +
  xlab("PC2")

write.table(napls.cauc[,1:2],"napls_cauc.txt",col.names=T,
            row.names=F,sep="\t",quote=F)

# run PCA in plink on caucasian subjects
# load in PCs

######################## SZ risk scores
# merge genetics files 
# generate list of file prefixes for the merge-list command
# hard calls 
#             !!!!! MOST OF THIS DATA IS MISSING
#                   (and apparently this is normal)

geno_files <- list.files(path="hard_calls", pattern = "*.fam", all.files = T)
geno_files <- lapply(geno_files, function(x) gsub(".fam","",x)) %>% unlist

write.table(geno_files,"geno_filelist.txt",col.names=F,row.names=F,quote=F,sep="\t")

# dosage data
# seems to be 2 columns for each subject, so format = 2
# 0-1 probabilities (not 0-2)
# need a list of files
dos_files <- list.files(path="dosage", pattern = "*.gz", all.files = T)

write.table(dos_files,"dosagelist.txt",col.names=F,row.names=T,quote=F,sep="\t")

# write input files for SZ risk scores
pgcdir    <- "/data/swe_gwas/ABZ/PGC-Results/"

# 128 best snps
# A1 is reference allele for OR (scz2.readme.pdf)
snps      <- read.table(paste0(pgcdir,"scz2.rep.128.txt"),header=T,stringsAsFactors = F)
snps$A1   <- sapply(snps[,8], function(x) unlist(strsplit(x,""))[1]) # A1
snps$lnOR <- log(snps$or) # natual log of OR (log10 is commmon log)

# all snps
# a1 is reference allele for OR (scz2.readme.pdf)
prs      <- read.table(paste0(pgcdir,"scz2.prs.txt"),header=T,stringsAsFactors = F)
prs$lnOR <- log(prs$or) # natual log of OR

# check overlap of SZ risk vars and .bim file of merged dataset
bim        <- read.table("napls_hardcalls_merged.bim",header=F,stringsAsFactors = F)
names(bim) <- c("chr","snp","pos","bp","a1","a2")
table(prs$snpid %in% bim$snp) # T = 15,593; F = 87,043; total prs snps = 102,636

# write input lists
write.table(snps[,c("snpid","A1","lnOR")],
            paste0(pgcdir,"pTlists/sigsnps-prs.txt"),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(prs[,c("snpid","a1","lnOR")],
            paste0(pgcdir,"pTlists/pTall-prs.txt"),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(prs[prs$p<.5,c("snpid","a1","lnOR")],
            paste0(pgcdir,"pTlists/pT5-prs.txt"),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(prs[prs$p<.2,c("snpid","a1","lnOR")],
            paste0(pgcdir,"pTlists/pT2-prs.txt"),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(prs[prs$p<.1,c("snpid","a1","lnOR")],
            paste0(pgcdir,"pTlists/pT1-prs.txt"),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(prs[prs$p<.05,c("snpid","a1","lnOR")],
            paste0(pgcdir,"pTlists/pT05-prs.txt"),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(prs[prs$p<.01,c("snpid","a1","lnOR")],
            paste0(pgcdir,"pTlists/pT01-prs.txt"),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(prs[prs$p<.001,c("snpid","a1","lnOR")],
            paste0(pgcdir,"pTlists/pT001-prs.txt"),
            col.names=T,row.names=F,quote=F,sep="\t")


