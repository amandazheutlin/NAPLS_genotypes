# PsychChip data from NAPLS - SNPs
# correlations with gray matter longitudinal phenotype

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2", "gdata")
invisible(lapply(libs, require, character.only = TRUE))

#### load in data
# genetics PCs from MDS in plink (see plink_MDS_cauc.sh)
PC.df <- read.table("pca.pruned.data/napls_cauc_PCs.eigenvec",header=T,stringsAsFactors = F)

# SZ risk scores
# read in all the score files (8 pT)
# score input lists from preprocessing.R
# scores generated using prs/SZrisk-prs.sh by plink
score_files   <- list.files(path="prs",pattern="napls_SZrisk_.*\\.profile")
score_files   <- lapply(score_files, function(x) gsub("napls","prs/napls",x)) %>% unlist
scores        <- lapply(score_files, function(x) read.table(x, header=T, stringsAsFactors=F))
names(scores) <- c("pT001","pT01","pT05","pT1","pT2","pT5","pTall","sigsnps")

# just select the score columns
prs.df <- lapply(scores, "[[", 4) %>% as.data.frame
prs.df <- cbind(scores$pT001$FID, scores$pT001$IID, scores$pT001$PHENO, prs.df)
names(prs.df)[1:3] <- c("FID","IID","PHENO")

# convert between IDs
# gene IDs --> imaging / demo IDs
# EDITED XLS: for UNC subjects, Collaborator.Sample.ID_2 seems like 
#             it is site and subject number so manually edited the
#             corresponding columns -- NAPLS2.Site.Number + NAPLS2.Site.ID
# PCs are from caucasian only, so subj loss here
id.conv <- read.xls("NAPLS2_IDconversion_edited.xls", stringsAsFactors = F)
prs.PC  <- merge(PC.df,prs.df,by="IID")

prs.PC$siteID     <- id.conv[match(prs.PC$IID,id.conv$Collaborator.Participant.ID),9]
prs.PC$subjID     <- id.conv[match(prs.PC$IID,id.conv$Collaborator.Participant.ID),10]
prs.PC$SiteSubjID <- paste0(prs.PC$siteID,prs.PC$subjID) %>% as.integer()

# demo + imaging info
demo    <- read.table("../napls_pheno-demo.txt",header=T,stringsAsFactors = F)
demo.df <- demo %>% select(SiteNumber,SubjectNumber,SubjectType,demo_sex,demo_age_ym,
                           demo_racial,demo_hispanic_latino,d_rh_label3_thickness)

demo.df$SiteSubjID  <- paste0(demo.df$SiteNumber,demo.df$SubjectNumber) %>% as.integer()
names(demo.df)[1:8] <- c("site","subj","dx","sex","age","race","ethnicity","gm")

# add dx (ugh)
# conv doc is a list of all converters (N = 94)
conv            <- read.csv("../ConversionAndPopsCriteria_13Jan2016.csv",header=T,stringsAsFactors=F)
conv$SiteSubjID <- paste0(conv$SiteNumber,conv$SubjectNumber) %>% as.integer()
demo.df$dx2     <- demo.df$dx
demo.df$dx2     <- ifelse(demo.df$SiteSubjID %in% conv$SiteSubjID,"Converter",demo.df$dx)

# combined df
allcols  <- merge(prs.PC,demo.df,by="SiteSubjID")
write.table(allcols,"napls-SZrisk-PCs-gm-demo.txt",col.names=T,row.names=F,quote=F,sep="\t")

# caucasian subj with genes + imaging
stats.df <- allcols[complete.cases(allcols$gm),]

#### statistics
# linear regressions
# target = rate of cortical thinning in superior frontal ROI
# predictors = SZ risk vars, age, sex, dx, site, 3 PCs (>2%)
prs.pT <- names(stats.df)[26:33]

prs_corrs <- matrix(nrow=8,ncol=4) # N = 273
for (i in 1:length(prs.pT)){
  f  <- formula(paste0("gm~",prs.pT[i],"+ age + sex + dx2 + site + PC1 + PC2 + PC3"))
  m1 <- lm(f, data=stats.df)
  tval <- summary(m1)$coefficients[2,3]
  prs_corrs[i,1] <- prs.pT[i]
  prs_corrs[i,2] <- tval
  prs_corrs[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$df[2]))
  prs_corrs[i,4] <- summary(m1)$coefficients[2,4]
}

prs_corrs        <- data.frame(prs_corrs, stringsAsFactors=F) 
prs_corrs[,2:4]  <- lapply(prs_corrs[,2:4], as.numeric) %>% as.data.frame
names(prs_corrs) <- c("prs.pT","tvalue","R2","pvalue")
prs_corrs$p.fdr  <- p.adjust(prs_corrs$pvalue, method = "BH")

write.table(prs_corrs,"prs-gm-3PCs.txt",col.names=T,row.names=F,quote=F,sep="\t")

# graph SZ risk x gray matter correlations
ggplot(stats.df, aes(x=pT5,y=gm)) + 
  geom_point() +
  geom_smooth() + 
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  ylab("Rate of Reduction of Gray Matter") +
  xlab("SZ Risk Score")


# just high risk 
stats.df.uhr  <- stats.df %>% dplyr::filter(dx=="Prodromal")

prs_corrs_uhr <- matrix(nrow=8,ncol=4) # N = 188
for (i in 1:length(prs.pT)){
  f  <- formula(paste0("gm~",prs.pT[i],"+ age + sex + site + PC1 + PC2 + PC3"))
  m1 <- lm(f, data=stats.df.uhr)
  tval <- summary(m1)$coefficients[2,3]
  prs_corrs_uhr[i,1] <- prs.pT[i]
  prs_corrs_uhr[i,2] <- tval
  prs_corrs_uhr[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$df[2]))
  prs_corrs_uhr[i,4] <- summary(m1)$coefficients[2,4]
}

prs_corrs_uhr        <- data.frame(prs_corrs_uhr, stringsAsFactors=F) 
prs_corrs_uhr[,2:4]  <- lapply(prs_corrs_uhr[,2:4], as.numeric) %>% as.data.frame
names(prs_corrs_uhr) <- c("prs.pT","tvalue","R2","pvalue")

write.table(prs_corrs,"prs-gm-3PCs-UHR.txt",col.names=T,row.names=F,quote=F,sep="\t")

