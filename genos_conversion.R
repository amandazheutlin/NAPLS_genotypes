# PsychChip data from NAPLS - SNPs
# correlations with conversion and risk

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/"
setwd(workdir)

libs <- c("dplyr", "ggplot2", "gdata", "stringr", "psych", "reshape2", "RColorBrewer")
invisible(lapply(libs, require, character.only = TRUE))

######## load in data
# genetics PCs from MDS in plink (see plink_PCA_cauc.sh)
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
id.conv  <- read.xls("NAPLS2_IDconversion_edited.xls", stringsAsFactors = F)
prs.PC   <- merge(PC.df,prs.df,by="IID")

# add demo data
# add age 
id.conv$SiteSubjID  <- paste0(id.conv$NAPLS2.Site.Number, id.conv$NAPLS2.Site.ID) 

demo_all            <- read.xls("ClientInfo.xlsx", stringsAsFactors = F)
demo_all$SiteSubjID <- paste0(demo_all$SiteNumber,demo_all$SubjectNumber) 

id.conv  <- left_join(id.conv,demo_all[,c("SiteSubjID","demo_age_ym")],by="SiteSubjID")

# sex, site, and dx from ID sheet
stats.df      <- merge(prs.PC,id.conv[,c("Collaborator.Participant.ID",
                                         "Gender","Primary.Disease", 
                                         "NAPLS2.Site.Number", "demo_age_ym")],
                       by.x="IID",by.y="Collaborator.Participant.ID")
names(stats.df)[33:36] <- c("sex","dx","site","age")

stats.df$risk         <- stats.df$dx %>% as.factor
levels(stats.df$risk) <- c("Prodromal", "Prodromal", "Control")

# write.table(stats.df,"napls-SZrisk-PCs-conv.txt",col.names=T,row.names=F,quote=F,sep="\t")

# add IDs for later
# stats.df$subj       <- id.conv[match(stats.df$IID, id.conv$Collaborator.Participant.ID), "NAPLS2.Site.ID"]
# stats.df$SubjSiteID <- id.conv[match(stats.df$IID, id.conv$Collaborator.Participant.ID), "SiteSubjID"]
# 
# stats.df$subj.pad   <- str_pad(stats.df$subj, 4, pad = "0")
# stats.df$naplsid    <- paste0(stats.df$site, ".", stats.df$subj.pad)
# write.table(stats.df,"napls-SZrisk-PCs-conv-IDs.txt",col.names=T,row.names=F,quote=F,sep="\t")


######### stats
# linear regressions for risk
# target = risk (control vs. prodromal)
# predictors = SZ risk vars, sex, site, 3 PCs (>2%)
prs.pT <- names(stats.df)[25:32]

prs_dx <- matrix(nrow=8,ncol=4)   # N = 601
for (i in 1:length(prs.pT)){
  f  <- formula(paste0(prs.pT[i],"~ risk + sex + age + site + PC1 + PC2 + PC3"))
  m1 <- lm(f, data=stats.df, na.action=na.omit)
  tval <- summary(m1)$coefficients[2,3]
  prs_dx[i,1] <- prs.pT[i]
  prs_dx[i,2] <- tval
  prs_dx[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$df[2]))
  prs_dx[i,4] <- summary(m1)$coefficients[2,4]
}

prs_dx        <- data.frame(prs_dx, stringsAsFactors=F) 
prs_dx[,2:4]  <- lapply(prs_dx[,2:4], as.numeric) %>% as.data.frame
names(prs_dx) <- c("prs.pT","tvalue","R2","pvalue")
prs_dx$p.fdr  <- p.adjust(prs_dx$pvalue, method = "BH")

write.table(prs_dx,"prs-dx-3PCs.txt",col.names=T,row.names=F,quote=F,sep="\t")

# linear regressions for conversion
# target = dx2 (converter vs. prodromal)
# predictors = SZ risk vars, age, sex, site, 3 PCs (>2%)
stats.df$dx  <- as.factor(stats.df$dx) 
stats.df.uhr <- stats.df %>% dplyr::filter(risk=="Prodromal") # N = 449

prs_conv <- matrix(nrow=8,ncol=4)
for (i in 1:length(prs.pT)){
  f  <- formula(paste0(prs.pT[i],"~ dx + sex + age + site + PC1 + PC2 + PC3"))
  m1 <- lm(f, data=stats.df.uhr, na.action=na.omit)
  tval <- summary(m1)$coefficients[2,3]
  prs_conv[i,1] <- prs.pT[i]
  prs_conv[i,2] <- tval
  prs_conv[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$df[2]))
  prs_conv[i,4] <- summary(m1)$coefficients[2,4]
}

prs_conv        <- data.frame(prs_conv, stringsAsFactors=F) 
prs_conv[,2:4]  <- lapply(prs_conv[,2:4], as.numeric) %>% as.data.frame
names(prs_conv) <- c("prs.pT","tvalue","R2","pvalue")
prs_conv$p.fdr  <- p.adjust(prs_conv$pvalue, method = "BH")

write.table(prs_conv,"prs-conv-3PCs-UHR.txt",col.names=T,row.names=F,quote=F,sep="\t")


# overall dx
prs_3groups <- matrix(nrow=8,ncol=4)   # N = 601
for (i in 1:length(prs.pT)){
  f  <- formula(paste0(prs.pT[i],"~ dx + sex + age + site + PC1 + PC2 + PC3"))
  m1 <- lm(f, data=stats.df, na.action=na.omit)
  tval <- summary(m1)$coefficients[2,3]
  prs_3groups[i,1] <- prs.pT[i]
  prs_3groups[i,2] <- tval
  prs_3groups[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$df[2]))
  prs_3groups[i,4] <- summary(m1)$coefficients[2,4]
}

# test1 <- lm(pTall ~ dx + sex + age + site + PC1 + PC2 + PC3, data = stats.df, na.action = na.omit)
# anova(test1)
# R2 = F / (F + n-p/p-1); p = # of predictors; n = subjects; F = F-value
# p = 8 (dx is 2); n = 600; F = 13.6377
# 13.6377/(13.6377 + (592/7)) = 0.1388639

prs_3groups        <- data.frame(prs_3groups, stringsAsFactors=F) 
prs_3groups[,2:4]  <- lapply(prs_3groups[,2:4], as.numeric) %>% as.data.frame
names(prs_3groups) <- c("prs.pT","tvalue","R2","pvalue")
prs_3groups$p.fdr  <- p.adjust(prs_3groups$pvalue, method = "BH")

write.table(prs_3groups,"prs-dx3groups-3PCs.txt",col.names=T,row.names=F,quote=F,sep="\t")

# demographics table
# demo info is for all participants included (N = 601)
# add extra vars
id.conv.inc <- id.conv[id.conv$Collaborator.Participant.ID %in% stats.df$IID,]
demo.df     <- demo_all[demo_all$SiteSubjID %in% id.conv.inc$SiteSubjID,]
demo.df$dx  <- id.conv.inc[match(demo.df$SiteSubjID, id.conv.inc$SiteSubjID), "Primary.Disease"]

# sex
table(demo.df$demo_sex,demo.df$dx)
chisq.test(demo.df$demo_sex,demo.df$dx)

# age
describeBy(demo.df$demo_age_ym,group=demo.df$dx)

m1   <- lm(demo_age_ym ~ dx, data = demo.df)
summary(m1)

m1.a <- lm(demo_age_ym ~ dx, data = demo.df[demo.df$dx != "Clinical High Risk-Converted",])
summary(m1.a)
m1.b <- lm(demo_age_ym ~ dx, data = demo.df[demo.df$dx != "Control",])
summary(m1.b)

# SES
# treated as continuous (??)
describeBy(demo.df$demo_income,group=demo.df$dx)
m2 <- lm(demo_income ~ dx, data = demo.df)
summary(m2)

# education
describeBy(demo.df$demo_education_years,group=demo.df$dx)
m3 <- lm(demo_education_years ~ dx, data = demo.df)
summary(m3)

m3.a <- lm(demo_education_years ~ dx, data = demo.df[demo.df$dx != "Clinical High Risk-Converted",])
summary(m3.a)
m3.b <- lm(demo_education_years ~ dx, data = demo.df[demo.df$dx != "Control",])
summary(m3.b)

# SOPS
sops            <- read.table("../SOPS.txt",header=T,stringsAsFactors = F)
sops$SiteSubjID <- paste(sops$SiteNumber,sops$SubjectNumber,sep="") %>% as.integer()
sops$dx         <- demo.df[match(sops$SiteSubjID,demo.df$SiteSubjID),"dx"]
sops.use        <- sops[sops$SiteSubjID %in% demo.df$SiteSubjID,]

sops.vars           <- c("cSOPSPositive","cSOPSNegative","cSOPSDisorganization","cSOPSGeneral")
sops.use$cSOPSTotal <- rowSums(sops.use[,sops.vars])
describeBy(sops.use[,c(sops.vars,"cSOPSTotal")],sops.use$dx)

# s1 <- lm(cSOPSPositive ~ group, data = sops.use[sops.use$group!="UC",])
# s2 <- lm(cSOPSNegative ~ group, data = sops.use[sops.use$group!="UC",])
# s3 <- lm(cSOPSDisorganization ~ group, data = sops.use[sops.use$group!="UC",])
# s4 <- lm(cSOPSGeneral ~ group, data = sops.use[sops.use$group!="UC",])

m4 <- lm(cSOPSTotal ~ dx, data = sops.use)
summary(m4)
m4.a <- lm(cSOPSTotal ~ dx, data = sops.use[sops.use$dx != "Control",])
summary(m4.a)
m4.b <- lm(cSOPSTotal ~ dx, data = sops.use[sops.use$dx != "Clinical High Risk-Converted",])
summary(m4.b)

#### symtpoms x PRS
sops.use$IID   <- id.conv[match(sops.use$SiteSubjID, id.conv$SiteSubjID), "Collaborator.Sample.ID"]
sops.use       <- merge(sops.use, stats.df[,c("IID", "PC1", "PC2", "PC3", prs.pT, "age", "sex")], by = "IID")

# in everyone
prs_sx <- matrix(nrow=8,ncol=4)   # N = 584
for (i in 1:length(prs.pT)){
  f  <- formula(paste0("cSOPSTotal ~", prs.pT[i],"+ sex + age + SiteNumber + PC1 + PC2 + PC3"))
  m1 <- lm(f, data=sops.use, na.action=na.omit)
  tval <- summary(m1)$coefficients[2,3]
  prs_sx[i,1] <- prs.pT[i]
  prs_sx[i,2] <- tval
  prs_sx[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$df[2]))
  prs_sx[i,4] <- summary(m1)$coefficients[2,4]
}

prs_sx        <- data.frame(prs_sx, stringsAsFactors=F) 
prs_sx[,2:4]  <- lapply(prs_sx[,2:4], as.numeric) %>% as.data.frame
names(prs_sx) <- c("prs.pT","tvalue","R2","pvalue")
prs_sx$p.fdr  <- p.adjust(prs_sx$pvalue, method = "BH")

write.table(prs_sx,"prs-sx-3PCs.txt",col.names=T,row.names=F,quote=F,sep="\t")

# in UHR
sops.use.uhr <- sops.use[sops.use$dx!="Control",]
prs_sx_uhr <- matrix(nrow=8,ncol=4)   # N = 434
for (i in 1:length(prs.pT)){
  f  <- formula(paste0("cSOPSTotal ~", prs.pT[i],"+ sex + age + SiteNumber + PC1 + PC2 + PC3"))
  m1 <- lm(f, data=sops.use.uhr, na.action=na.omit)
  tval <- summary(m1)$coefficients[2,3]
  prs_sx_uhr[i,1] <- prs.pT[i]
  prs_sx_uhr[i,2] <- tval
  prs_sx_uhr[i,3] <- sqrt(tval^2/(tval^2+summary(m1)$df[2]))
  prs_sx_uhr[i,4] <- summary(m1)$coefficients[2,4]
}

prs_sx_uhr        <- data.frame(prs_sx_uhr, stringsAsFactors=F) 
prs_sx_uhr[,2:4]  <- lapply(prs_sx_uhr[,2:4], as.numeric) %>% as.data.frame
names(prs_sx_uhr) <- c("prs.pT","tvalue","R2","pvalue")
prs_sx_uhr$p.fdr  <- p.adjust(prs_sx_uhr$pvalue, method = "BH")

write.table(prs_sx_uhr,"prs-sx-3PCs-UHR.txt",col.names=T,row.names=F,quote=F,sep="\t")

# categories within total SOPS scores in UHR
sops.pos <- lm(cSOPSPositive ~ pTall + age + sex + SiteNumber + PC1 + PC2 + PC3, data=sops.use.uhr)
summary(sops.pos)

sops.neg <- lm(cSOPSNegative ~ pTall + age + sex + SiteNumber + PC1 + PC2 + PC3, data=sops.use.uhr)
summary(sops.neg)

sops.org <- lm(cSOPSDisorganization ~ pTall + age + sex + SiteNumber + PC1 + PC2 + PC3, data=sops.use.uhr)
summary(sops.org)

sops.gen <- lm(cSOPSGeneral ~ pTall + age + sex + SiteNumber + PC1 + PC2 + PC3, data=sops.use.uhr)
summary(sops.gen)

######### figure
# mean and SE PRS for each clinical group
# save residual PRS values
prs_resid   <- matrix(nrow=nrow(stats.df),ncol=length(prs.pT)) # N = 601
for (i in 1:length(prs.pT)){
  f  <- formula(paste0(prs.pT[i],"~ sex + age + site + PC1 + PC2 + PC3"))
  m1 <- lm(f, data=stats.df, na.action=na.omit)
  prs_resid[,i] <- residuals(m1)
}

prs_resid <- data.frame(prs_resid, stringsAsFactors=F) 

# figure
# using non-residual values bc they look weird
prs_resid        <- cbind(prs_resid, stats.df[,34])
names(prs_resid) <- c(prs.pT,"dx")
fig1             <- melt(prs_resid[,c(9,8,1:7)], id.var="dx")
# fig1             <- melt(stats.df[,c(32,25:31,34)], id.var="dx")
# fig1$dx          <- factor(fig1$dx)
levels(fig1$dx)  <- c("CHR Converted", "CHR Not Converted", "Control")
levels(fig1$variable) <- c("5x10-8", ".001", ".01", ".05", ".1", ".2", ".5", "1")

fig1.stats    <- fig1 %>% group_by(dx,variable) %>% summarise_each(funs(mean,sd,n()))
fig1.stats$se <- fig1.stats$sd/sqrt(fig1.stats$n)

ggplot(fig1.stats,aes(x=variable, y=mean, ymax = mean + se, ymin = mean - se, fill=dx)) + 
  geom_bar(stat="identity", position = position_dodge(), aes(fill=dx)) + 
  geom_errorbar(position=position_dodge(width=0.9), width=.25) +
  geom_hline(yintercept=0) +
  scale_fill_brewer(palette = "Blues") +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13, color="black"),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size=13, color="black"),
        legend.title = element_blank()) +
  ylab("Mean PRS (Residuals)") +
  xlab("\np-Value Threshold (pT)")



