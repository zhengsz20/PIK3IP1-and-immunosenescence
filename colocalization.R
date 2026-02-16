library(coloc)
library(tidyverse)
library(GWAS.utils)
library(reshape2)
library(susieR)
library(ieugwasr)

out1 <- fread("F:/finngen_R12/finngen_R12_E4_HYTHY_AI_STRICT.gz")
eaf_ref <- fread("E:/eQTLGEN/cis_eqtl/MAF_pos_added.txt.gz")

eaf_ref <- eaf_ref[,c(1,5,9)]
exp_dat <- fread(~/full_cis_eqtl_ENSG00000100100.txt")
exp_dat <- merge(exp_dat,eaf_ref,by="SNP")
exp_dat1 <- subset(exp_dat,exp_dat$AssessedAllele==exp_dat$AlleleB)
exp_dat1$EAF <- exp_dat1$AlleleB_all
exp_dat1 <- exp_dat1[,c(1,18)]
exp_dat2 <- subset(exp_dat,exp_dat$AssessedAllele!=exp_dat$AlleleB)
exp_dat2$EAF <- 1-exp_dat2$AlleleB_all
exp_dat2 <- exp_dat2[,c(1,18)]
exp_dat3 <- rbind(exp_dat1,exp_dat2)
exp_dat <- merge(exp_dat,exp_dat3,by="SNP")
exp_dat$Beta <- exp_dat$Zscore/sqrt(2*exp_dat$EAF*(1-exp_dat$EAF)*(exp_dat$NrSamples+exp_dat$Zscore^2))
exp_dat$SE <- 1/sqrt(2*exp_dat$EAF*(1-exp_dat$EAF)*(exp_dat$NrSamples+exp_dat$Zscore^2))


exp_dat <- exp_dat[!(which(duplicated(exp_dat$SNP))),]
exp_dat <- na.omit(exp_dat)
exp_dat <- subset(exp_dat,exp_dat$SNPPos>=31613555-50000&exp_dat$SNPPos<=31613555+50000)
exp_dat <- exp_dat[which(exp_dat$Beta!=0),]
exp_dat$MAF <- ifelse(exp_dat$EAF>=0.5,1-exp_dat$EAF,exp_dat$EAF)
exp_dat <- exp_dat[,c(1,5,6,19,2,20,13,21)]

merge_snp <- merge(exp_dat,out1,by.x="SNP",by.y ="rsids")
merge_snp <- na.omit(merge_snp)
exp_dat <- merge_snp[,c(1,2,3,4,5,6,7,8)]
colnames(exp_dat) <- c("rsid","expo_a1","expo_a2","expo_eff","expo_pv","expo_se","N","MAF")

out_dat <- merge_snp[,c(1,12,11,16,17,14)]
colnames(out_dat) <- c("rsid","outc_a1","outc_a2","outc_eff", "outc_se", "outc_pv")

tmp <- merge(exp_dat,out_dat,by="rsid")

tmp$expo_varbeta <- (tmp$expo_se)^2
tmp$outc_varbeta <- (tmp$outc_se)^2


bbb <- length(as.vector(which(duplicated(tmp$rsid))))>0
if (bbb>0){
  tmp <- tmp[-which(duplicated(tmp$rsid)),]
}else{
  tmp <- tmp 
}

D1_expo <- list() 
D2_outc <- list()
D1_expo$type <- "quant"
D1_expo$snp <- tmp$rsid
D1_expo$beta <- tmp$expo_eff
D1_expo$varbeta <- tmp$expo_varbeta
D1_expo$N <- tmp$N
D1_expo$MAF <- tmp$MAF

D2_outc$type <- "cc"
D2_outc$snp <- tmp$rsid
D2_outc$beta <- tmp$outc_eff
D2_outc$varbeta <- tmp$outc_varbeta
res <- coloc.abf(dataset1 = D1_expo, dataset2 = D2_outc)
res2 <- t(as.data.frame(res$summary))

###############################################################
library(locuscomparer)
tmp2 <- tmp[,c(1,5,13)]
gwas_fn = tmp2[,c(1,3)]
colnames(gwas_fn) <- c("rsid","pval")
qtl_fn <- tmp2[,c(1,2)]
colnames(qtl_fn) <- c("rsid","pval")

locuscompare(in_fn1 = gwas_fn, in_fn2 = qtl_fn, title1 = 'Hypothyroidism', title2 = 'PIK3IP1',snp = "rs2413035",population = "EUR",combine = F,genome = "hg19")


tmp = tmp %>% filter((expo_a1==outc_a1&expo_a2==outc_a2)|(expo_a1==outc_a2&expo_a2==outc_a1)) 
tmp = tmp %>% mutate(outc_eff = ifelse(expo_a1==outc_a1,outc_eff,-outc_eff))

ld <- readRDS("F:/UKB_LD_MATRIX/22_b37/ukb_b37_0.1_chr22.R_snp.31439918_32664986.RDS")
ld <- na.omit(ld)

ld_variant <- fread("F:/UKB_LD_MATRIX/22_b37/ukb_b37_0.1_chr22.R_snp.31439918_32664986.Rvar")
colnames(ld) <- ld_variant$id
rownames(ld) <- ld_variant$id


snp_id <- as.data.frame(rownames(ld))
colnames(snp_id) <- "rsid"
tmp <- merge(tmp,snp_id,by="rsid")
tmp <- tmp[match(snp_id$rsid,tmp$rsid,),]
tmp <- na.omit(tmp)

snp_id <- as.data.frame(rownames(ld))
colnames(snp_id) <- "rsid"
tmp <- merge(snp_id,tmp,by="rsid")
tmp <- tmp[match(snp_id$rsid,tmp$rsid,),]
tmp <- na.omit(tmp)

z_expo <- tmp$expo_eff/tmp$expo_se
z_outc <- tmp$outc_eff/tmp$outc_se

N_expo <- 31684
N_outc <- 417391

ld <- ld[tmp$rsid, tmp$rsid, drop = FALSE]
#lambda_expo <- estimate_s_rss(z_expo, ld, n = N_expo)
#lambda_outc <- estimate_s_rss(z_outc, ld, n = N_outc)

fit_expo <- susie_rss(z_expo,ld, N_expo)
fit_outc <- susie_rss(z_outc,ld, N_outc)

# coloc susie
res <- coloc.susie(fit_expo, fit_outc)
res
