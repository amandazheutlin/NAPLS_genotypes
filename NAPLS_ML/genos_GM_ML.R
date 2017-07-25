# PsychChip data from NAPLS - SNPs
# gray matter longitudinal phenotype
# try random forest with top SNPs (79/108)

#### housekeeping
workdir <- "/data/swe_gwas/ABZ/NAPLS/PSYCH-CHIP/NAPLS_ML/"
setwd(workdir)

libs <- c("dplyr", "psych", "ggplot2", "gdata", "caret", "randomForest", "doMC")
invisible(lapply(libs, require, character.only = TRUE))

#### extract hard call genotypes + proxies
# 40 / 108 SNPs present in hard call data
gwas  <- read.table("/data/swe_gwas/ABZ/PGC-Results/108-PGC.txt",header=F) 
gwas  <- gwas$V1 %>% as.character # 108 best snps from gwas loci
calls <- read.table("../108_hardcall_napls.map",header=F, stringsAsFactors = F) # 40/108

missing <- gwas[!(gwas %in% calls$V2)] # 68 missing SNPs
# write.table(missing, "missing_108.txt", col.names=F, row.names=F, quote=F, sep="\t")

# SNiPA, .8 r2 proxy search
# 1000 genomes phase 3: European, GRCh37 reference assembly
proxy <- read.csv("proxySearch.results.csv",header=T, sep="\t", stringsAsFactors = F) # using SNiPA, .8 r2
alias <- proxy$RSALIAS %>% as.character # some proxy SNPs have alias rs IDs
alias <- lapply(alias, function(x) strsplit(x, ",")) %>% unlist

proxy.alias <- c(proxy$RSID, alias) 
proxy.alias <- proxy.alias[!is.na(proxy.alias)] # 3528 rs IDs for 59/68 index SNPs
# write.table(proxy.alias, "proxy_NAPLS.txt", col.names=F, row.names=F, quote=F, sep="\t")

p.calls       <- read.table("108_proxy.map",header=F,stringsAsFactors = F)
proxy.present <- proxy[proxy$RSID %in% p.calls$V2,] # proxy SNPs (829/830)
proxy.present <- rbind(proxy.present, proxy[grepl("rs116745399",proxy$RSALIAS)==T,]) # alias rs ID proxy SNP 
index         <- unique(proxy.present$QRSID) # 39 additional SNP proxies (now 79/108)

# select the best proxy SNPs for index SNPs with multiple proxies
# 39 index SNPs; 829 proxy SNPs
# best snps range in r2 from .81 to 1, mean = .95
proxy.man <- proxy.present %>% select(QRSID, RSID, R2) 

best.snps <- data.frame()
for (i in 1:length(index)){
  tmp <- proxy.man[proxy.man$QRSID == index[i], ]
  tmp2 <- tmp[tmp$R2==max(tmp$R2), ]
  best.snps <- rbind(best.snps,tmp2[1,])
}


# missing.still <- missing[!(missing %in% proxy.present$QRSID)] # 29 missing still
# write.table(missing.still, "missing-still_108.txt", col.names=F, row.names=F, quote=F, sep="\t")

# dbSNP-Q, .8 r2 proxy search
# dbSNP IDs are given, so I'm adding 'rs' to them all
# proxy2 <- read.table("query_results.txt", header=T, sep= "\t")
# 
# proxy2$snp_id_rs   <- lapply(proxy2$snp_id, function(x) paste0("rs",x)) %>% unlist
# proxy2$proxy_id_rs <- lapply(proxy2$proxy_id, function(x) paste0("rs",x)) %>% unlist
# write.table(proxy2[,"proxy_id_rs"], "proxy2_NAPLS.txt", col.names=F, row.names=F, quote=F, sep="\t")
# 
# p.calls2 <- read.table("108_proxy2.map",header=F,stringsAsFactors = F)
# proxy.present2 <- proxy2[proxy2$proxy_id_rs %in% p.calls2$V2,] # NO NEW SNPS


#### load in data
demo    <- read.table("../napls-SZrisk-PCs-gm-demo.txt", header=T, stringsAsFactors = F)
genos   <- read.table("../108_hardcall_napls.ped", sep="\t", header=F, stringsAsFactors = F)
p.genos <- read.table("108_proxy.ped", sep="\t", header=F, stringsAsFactors = F)

names(genos)[2]       <- "IID"
names(genos)[7:46]    <- calls$V2
names(p.genos)[2]     <- "IID"
names(p.genos)[7:836] <- p.calls$V2
vars                  <- c("IID", best.snps$RSID)
p.genos               <- p.genos[,vars] # just best SNPs

# recode genotypes to be risk dosage (# of risk alleles)
# 2 = risk, 1 = het, 0 = no risk
alleles     <- read.table("/data/swe_gwas/ABZ/PGC-Results/scz2.rep.128.txt",header=T, stringsAsFactors = F)
alleles$a1  <- substr(alleles$a1a2,1,1)
alleles$a2  <- substr(alleles$a1a2,2,2)
alleles$maj <- ifelse(alleles$frqa>.5,"MAJOR","MINOR") # for proxy coding


# recode call genotypes
genos.rec  <- genos[,c(2,7:46)]
genos.snps <- names(genos.rec)[2:41]

for (i in 1:length(genos.snps)){
  risk <- paste0(alleles[alleles$snpid==genos.snps[i],"a1"], " ", alleles[alleles$snpid==genos.snps[i],"a1"])
  het  <- paste0(alleles[alleles$snpid==genos.snps[i],"a1"], " ", alleles[alleles$snpid==genos.snps[i],"a2"])
  het2 <- paste0(alleles[alleles$snpid==genos.snps[i],"a2"], " ", alleles[alleles$snpid==genos.snps[i],"a1"])
  genos.rec[,genos.snps[i]] <- ifelse(genos.rec[,genos.snps[i]]==risk,2,
                                      ifelse(genos.rec[,genos.snps[i]]==het,1,
                                             ifelse(genos.rec[,genos.snps[i]]==het2,1,0)))
}


# recode proxy genotypes
p.genos.rec  <- p.genos
p.genos.snps <- names(p.genos.rec)[2:40]
  
for (i in 1:length(p.genos.snps)){
  # find index snp for proxy snp in proxy sheet
  indexsnp  <- proxy[proxy$RSID==p.genos.snps[i],"QRSID"]
  # determine if risk is major or minor allele in alleles
  major     <- alleles[alleles$snpid==indexsnp,"maj"]
  # define 'risk' as either major or minor using proxy sheet
  major.not <- ifelse(major=="MAJOR","MINOR","MAJOR")
  risk <- paste0(proxy[proxy$RSID==p.genos.snps[i],major], " ", proxy[proxy$RSID==p.genos.snps[i],major])
  het  <- paste0(proxy[proxy$RSID==p.genos.snps[i],major], " ", proxy[proxy$RSID==p.genos.snps[i],major.not])
  het2 <- paste0(proxy[proxy$RSID==p.genos.snps[i],major.not], " ", proxy[proxy$RSID==p.genos.snps[i],major])
  p.genos.rec[,p.genos.snps[i]] <- ifelse(p.genos.rec[,p.genos.snps[i]]==risk,2,
                                      ifelse(p.genos.rec[,p.genos.snps[i]]==het,1,
                                             ifelse(p.genos.rec[,p.genos.snps[i]]==het2,1,0)))
}


# concatenate SNPs + proxy SNPs
allgenos <- merge(genos.rec, p.genos.rec, by="IID")

# see if PCs correlate with the gm phenotype;
# if yes, then regress out of phenotype (below)
# m1 <- lm(gm ~ PC1, data = df, na.action = na.exclude) # p = .257
# m2 <- lm(gm ~ PC2, data = df, na.action = na.exclude) # p = .707
# m3 <- lm(gm ~ PC3, data = df, na.action = na.exclude) # p = .201

# convert between IDs
# gene IDs --> imaging / demo IDs
# PCs are from caucasian only, so subj loss here
id.conv <- read.xls("../NAPLS2_IDconversion_edited.xls", stringsAsFactors = F)

allgenos$siteID     <- id.conv[match(allgenos$IID,id.conv$Collaborator.Participant.ID),9]
allgenos$subjID     <- id.conv[match(allgenos$IID,id.conv$Collaborator.Participant.ID),10]
allgenos$SiteSubjID <- paste0(allgenos$siteID,allgenos$subjID) %>% as.integer()

# make composite df to verify subjects are matched up across vars
demo.use <- demo %>% select(SiteSubjID, age, sex, site, dx, dx2, gm)
df       <- merge(demo.use, allgenos, by="SiteSubjID")

SNPs <- names(df)[9:87]

# regress out everything to build ML 
# age, sex, site
# PCs don't correlate with phenotype so do not regress out 
res_mod   <- lm(gm ~ age + sex + site, data = df, na.action = na.omit)
df$gm_res <- residuals(res_mod)

#### make matrices 
# target = rate of cortical thinning in superior frontal ROI (residuals)
# predictors = 79 genotypes
# N = 267

target     <- df$gm_res %>% as.data.frame
predictors <- df[,names(df) %in% SNPs] %>% as.matrix

# split sample 
df1 <- df %>% group_by(dx2) %>% sample_frac(.5)
df2 <- df[!(df$SiteSubjID %in% df1$SiteSubjID), ]

target.split     <- df1$gm_res %>% as.data.frame
predictors.split <- df1[,names(df1) %in% SNPs] %>% as.matrix

target.split2     <- df2$gm_res %>% as.data.frame
predictors.split2 <- df2[,names(df2) %in% SNPs] %>% as.matrix

##### ML preparation
registerDoMC(detectCores()-1)
set.seed(1)

## set up cross validation
fitControl <- trainControl(method = "repeatedcv",
                           repeats = 3,
                           number = 10,
                           savePredictions = TRUE,
                           selectionFunction = "oneSE")

## hyperparameter grid search
rfGrid     <- expand.grid(.mtry = c(3,4,10)) 

##### model training

rfPipeline <- function(response, Xmat = predictors.split, grid = rfGrid, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "rf",
                 metric = "Rsquared",
                 tuneGrid = grid,
                 trControl = cvpar,
                 importance = TRUE)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}


rf.loop     <- lapply(target.split, function(i) rfPipeline(i))
rf.models   <- lapply(rf.loop, function(i) i[[1]])
rf.perfs    <- lapply(rf.loop, function(i) i[[2]])

save.image(file = "napls_gm_rf.RData")

##### graph importance of features
# snp-gene table
# 'proxy' column is set of SNPs used in analysis
snp2gene       <- read.csv("/data/ML_genetics/Schizophrenia-Zheutlin/scz2.anneal.108.genes.csv", stringsAsFactors = F)
snp2gene$proxy <- best.snps[match(snp2gene$bestsnp,best.snps$QRSID),"RSID"] # set proxy SNPs
snp2gene$proxy <- ifelse(is.na(snp2gene$proxy), snp2gene$bestsnp, snp2gene$proxy) # fill in index SNPs
table(snp2gene$bestsnp==snp2gene$proxy) ## check that it worked (39 non-matching bc proxy)  

# importance variables for GM
fig_gm          <- caret::varImp(rf.models[[1]], scale=TRUE)[[1]] %>% as.data.frame
fig_gm          <- cbind(rownames(fig_gm), fig_gm) %>% arrange(-Overall)
names(fig_gm)   <- c("SNP", "Importance")
fig_gm$Gene     <- snp2gene[ match(fig_gm$SNP, snp2gene$proxy), "genes"]
fig_gm$SNP      <- factor(fig_gm$SNP, levels = fig_gm$SNP[order(fig_gm$Importance)])
fig_gm$GenePlot <- fig_gm$Gene

top20_genes <- fig_gm[1:20,"Gene"]
write.table(top20_genes, "top20_genes_napls.txt", col.names=F, row.names=F, quote=F, sep="\t")

# top20 are mostly single genes (12/20, 60%)
# how many out of 108 are single genes?
snps.used <- snp2gene[snp2gene$proxy %in% SNPs,]
table(grepl(",",snp2gene$genes)) # 45/108, 42%
table(grepl(",",snps.used$genes)) # 31/79, 39%
chi  <- matrix(c(12, 45, 8, 63), nrow=2)
chi2 <- matrix(c(12, 31, 8, 48), nrow=2)
prop.test(chi, correct = F) # not different, p = .130
prop.test(chi2, correct = F) # not different, p = .094

# manually shorten gene names when relevant
fig_gm[1,"GenePlot"]  <- "ALDOA, ASPHD1, C16orf92*"
fig_gm[7,"GenePlot"] <- "C2orf82, EFHD1, GIGYF2*"
fig_gm[14,"GenePlot"] <- "BTBD18, C11orf31, CLP1*"
fig_gm[18,"GenePlot"] <- "NOSIP, PRR12, PRRG2*"

# figure
ggplot(fig_gm[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(stat="identity",width=.85, fill="#0000FF") +
  geom_text(stat="identity",y=61, hjust=0, aes(label=GenePlot), colour="white", size=4) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust=-.06, size=18, face="bold"),
        plot.margin = unit(c(.75, .75, .75, .75),"cm")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(60,100))


#### test model on baseline imaging data (ind sample)
load("napls_gm_rf.RData")

# baseline data (N = 379 subjects)
# one subject (site = 7, subjID = 12) in longitudinal sample but NOT baseline
bl <- read.table("/data/swe_gwas/ABZ/NAPLS/exome_gwas/NAPLS2_fs_bl_all_new_4_10_2016.csv",
                 header=T,stringsAsFactors = F)

bl$siteID     <- substr(bl$SiteNumber,5,5) %>% as.integer()
bl$subjID     <- substr(bl$imagingID,4,7) %>% as.integer()
bl$SiteSubjID <- paste0(bl$siteID,bl$subjID) %>% as.integer()

bl.only      <- bl[!(bl$SiteSubjID %in% df$SiteSubjID),] # 478 subjects

# exclude caucasian (N = 211 subjects)
cauc            <- read.table("../pca.pruned.data/napls_cauc_PCs.eigenvec",header=T,stringsAsFactors = F)
cauc$siteID     <- id.conv[match(cauc$IID,id.conv$Collaborator.Participant.ID),"NAPLS2.Site.Number"]
cauc$subjID     <- id.conv[match(cauc$IID,id.conv$Collaborator.Participant.ID),"NAPLS2.Site.ID"]
cauc$SiteSubjID <- paste0(cauc$siteID,cauc$subjID) %>% as.integer()

bl.only.cauc    <- bl.only[bl.only$SiteSubjID %in% cauc$SiteSubjID,] # 211 subjects

bl.use  <- bl.only.cauc %>% select(SiteSubjID, age_scan, demo_sex, 
                                   siteID, SubjectType, conversion, 
                                   rh_superiorfrontal_thickness,
                                   lh_superiorfrontal_thickness) 
rep.df  <- merge(bl.use, allgenos, by="SiteSubjID")  # 211 subjects

# regress out age, sex, site
# PCs don't correlate with phenotype so do not regress out 
res_mod_rh.rep      <- lm(rh_superiorfrontal_thickness ~ age_scan + demo_sex + siteID.x, 
                       data = rep.df, na.action = na.omit)
res_mod_lh.rep      <- lm(lh_superiorfrontal_thickness ~ age_scan + demo_sex + siteID.x, 
                       data = rep.df, na.action = na.omit)
rep.df$gm_res_rh    <- residuals(res_mod_rh.rep)
rep.df$gm_res_lh    <- residuals(res_mod_lh.rep)

#### make matrices 
# target = rate of cortical thinning in superior frontal ROI (residuals)
# predictors = 79 genotypes
# N = 211

target.rep.rh         <- rep.df$gm_res_rh %>% as.data.frame
names(target.rep.rh)  <- "gm.rh.res"

# target.rep.rh2        <- rep.df$rh_superiorfrontal_thickness %>% as.data.frame
# names(target.rep.rh2) <- "gm.rh"
# 
# target.rep.lh         <- rep.df$gm_res_lh %>% as.data.frame
# names(target.rep.lh)  <- "gm.lh.res"
# 
# target.rep.lh2        <- rep.df$lh_superiorfrontal_thickness %>% as.data.frame
# names(target.rep.lh2) <- "gm.lh"

predictors.rep        <- rep.df[,names(rep.df) %in% SNPs] %>% as.matrix

## test replication
rf.pred  <- predict(rf.models[[1]], newdata=predictors.rep)
rep.mod  <- cor.test(rf.pred, target.rep.rh$gm.rh.res) # r = -.035, p = .617
# rep.mod2 <- cor.test(rf.pred, target.rep.rh2$gm.rh) # r = -.066, p = .344
# rep.mod3 <- cor.test(rf.pred, target.rep.lh$gm.lh.res) # r = -.044, p = .530
# rep.mod4 <- cor.test(rf.pred, target.rep.lh2$gm.lh) # r = -.076, p = .276

rf.pred.split  <- predict(rf.models[[1]], newdata=predictors.split2)
rep.mod.split  <- cor.test(rf.pred.split, target.split2$.) # r = .059, p = .493

## do the predicted values do anything useful? (no)
rep.df$pred <- rf.pred

# predict risk (166 prodromes; 42 controls)
pred.risk   <- lm(pred ~ SubjectType, data = rep.df, na.action = na.omit)
summary(pred.risk) # p = .121

# predict conversion (24 converters; 142 non-converters)
rep.df.uhr            <- rep.df %>% dplyr::filter(SubjectType=="Prodromal") # N = 166
rep.df.uhr$conversion <- ifelse(is.na(rep.df.uhr$conversion),0,1) 

pred.dx     <- lm(pred ~ conversion, data = rep.df.uhr, na.action = na.omit)
summary(pred.dx) # p = .285


#######################################
#### build model on raw data

# new target; same predictors 
target.raw    <- df$gm %>% as.data.frame

rf.loop.raw   <- lapply(target.raw, function(i) rfPipeline(i))
rf.models.raw <- lapply(rf.loop.raw, function(i) i[[1]])
rf.perfs.raw  <- lapply(rf.loop.raw, function(i) i[[2]])

## test replication
rf.pred.raw  <- predict(rf.models.raw[[1]], newdata=predictors.rep)
rep.mod.raw  <- cor.test(rf.pred.raw, target.rep.rh2$target_gm) # r = -.080, p = .248
rep.mod.raw2 <- cor.test(rf.pred.raw, target.rep.lh2$target_gm) # r = -.040, p = .571


#######################################
#### build model on baseline right superior frontal 

# build model on residual values from baseline right sup front
rfPipeline.bl <- function(response, Xmat = predictors.rep, grid = rfGrid, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "rf",
                 metric = "Rsquared",
                 tuneGrid = grid,
                 trControl = cvpar,
                 importance = TRUE)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}

rf.loop.bl   <- lapply(target.rep.rh, function(i) rfPipeline.bl(i))
rf.models.bl <- lapply(rf.loop.bl, function(i) i[[1]])
rf.perfs.bl  <- lapply(rf.loop.bl, function(i) i[[2]])

## test replication in subjects w two imaging points
bl.trainset <- bl[bl$SiteSubjID %in% df$SiteSubjID,] # 266 subjects
bl.df       <- merge(bl.trainset, allgenos, by="SiteSubjID")  # 266 subjects

# predictors for 266 subjs
predictors.blrep <- bl.df[,names(bl.df) %in% SNPs] %>% as.matrix
rf.pred.bl       <- predict(rf.models.bl[[1]], newdata=predictors.blrep)

# targets
bl.df$gm_long  <- df[match(bl.df$SiteSubjID,df$SiteSubjID),"gm_res"]
target.gm_long <- bl.df$gm_long %>% as.data.frame

res_mod_rh      <- lm(rh_superiorfrontal_thickness ~ age_scan + demo_sex + siteID.x, 
                      data = bl.df, na.action = na.omit)
bl.df$gm_res_rh <- residuals(res_mod_rh)
target.bl_train <- bl.df$gm_res_rh %>% as.data.frame

rep.mod.bl  <- cor.test(rf.pred.bl, target.gm_long[,1]) # r = .004, p = .944 (cort thinning residuals)
rep.mod.bl2 <- cor.test(rf.pred.bl, target.bl_train[,1]) # r = -.027, p = .662 (baseline)




