#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
wd<-"/mnt/sda/gagelo01/Projects/2024/LPL_pathway/"
setwd(wd)
source("Analysis/tosourcedrug.R")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

#####
res_multicis <- fread("Data/Modified/res_multicis_independent.txt")
ctp <- fread("/mnt/sda/gagelo01/Projects/2024/UKB_phewas/Data/Modified/clinical_toxicity_panel.txt")
dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")
dt<-data.table(id=res_multicis[!grepl("^eqtl-|^prot-", id.exposure), unique(id.exposure)])
dt[, c("hgnc", "id_small"):=tstrsplit(id, "_", fixed = TRUE)]
dt_gene_region <- merge(dt_gene_region, dt[,.(id, id_small)], by = "id", all.x = TRUE)
dt_gene_region[is.na(id_small), id_small:=id]
dt_gene_region <- merge(dt_gene_region, df_index[, .(id,sample_size)], by.x = "id_small", by.y = "id", all.x = TRUE)
res_multicis <- merge(res_multicis, dt_gene_region[,.(id, hgnc, sample_size)], by.x = "id.exposure", by.y = "id", all.x = TRUE)
res_multicis <- merge(res_multicis, df_index[, .(id,clean_variable_name)], by.x = "id.outcome", by.y = "id", all.x = TRUE)
res_multicis <- merge(res_multicis, ctp[,.(id, category, subcategory, is_bad_when)], by.x = "id.outcome", by.y = "id", all.x = TRUE)
res_multicis[,c("hgnc", "idscale") := tstrsplit(id.exposure, "_")]
data_tsmr <- res_multicis[method == "Inverse variance weighted"]
data_tsmr[,b:=b*-1]
data_tsmr[,lci:=b-1.96*se]
data_tsmr[,uci:=b+1.96*se]

data_cox <- fread( "Data/Modified/data_cox.txt")
res_cox <- fread("Data/Modified/LPL/res_cox.txt")
res_cox_small<-data.table(res_cox)
res_cox_small[, c("ph_chisq", "ph_df", "ph_p"):=NULL]
res_cox_small<- res_cox_small[cov_inc== "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10", ]
res_cox_small <-res_cox_small[str_count(string = IV, pattern = "\\*")<2,]
colnom <- c("HR", "lci", "uci")
res_cox_small[, (colnom):=lapply(.SD, function(x) 1/x), .SDcols = colnom]

res_continuous <- fread("Data/Modified/LPL/res_continuous.txt")
res_continuous_small <- res_continuous[grepl("PRS", exposure),]
res_continuous_small <- res_continuous_small[cov_inc== "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10", ]
res_continuous_small <- res_continuous_small[str_count(string = IV, pattern = "\\*")<2,]
res_continuous_small[,b:=b*-1]
res_continuous_small[,lci:=b-1.96*se]
res_continuous_small[,uci:=b+1.96*se]

k1 <- res_multicis[grepl("trait-33-", id.outcome)&method%in%c("Inverse variance weighted"), ]
scatter <- dcast(k1, id.outcome ~ hgnc, value.var = c("b", "se"))
colnom<-colnames(scatter)[grepl("^b_", colnames(scatter))]
scatter[,(colnom):=lapply(.SD, function(x) x*-1) , .SDcols = colnom]
######
return_format_data<-function(data) {
  k <- data[, paste0("OR = ", format(round(exp(b), digits = 2), nsmall = 2), ", 95% CI=", format(round(exp(lci), digits = 2), nsmall = 2), " to ",  format(round(exp(uci), digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(id.exposure, "on", outcome)]
  return(k)
}
return_format_data_noexp <-function(data, digits = 2) {
  k<-data[, paste0(format(round(b, digits = digits), nsmall = digits) , " 95% CI=", format(round(lci, digits = digits), nsmall = digits), " to ",  format(round(uci, digits = digits), nsmall = digits), ", p=",pval%>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(id.exposure, "on", outcome)]
  return(k)
}

return_format_data_noexp_HR <-function(data) {
  k<-data[, paste0("HR = ", format(round(HR, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), " to ",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(outcome)]
  return(k)
}

return_format_fstat <-function(data) {
  k <- data[, paste0(N, " SNPs (r2 = ", round(rsq*100, digits =2), "%; F-statistics = ",  round(fstat, digits = 0), ")")]
  names(k) <- data$hgnc
  return(k)
}

summary_logical <- function(x) {
  k1 <- sum(x)
  k2<- length(x)
  paste0(k1,"/", k2, " (", round(k1/k2, digits = 2)*100, "%)")
}

# z-test for comparing two independent regression coefficients.
z_test <- function(b_x, se_x, b_y, se_y) {
  se_diff <- sqrt(se_x^2 + se_y^2)
  Z <- (b_x - b_y) / se_diff
  p_value <- 2 * (1 - pnorm(abs(Z)))
  return(data.table(diff = b_x - b_y, se = se_diff, p_value = p_value))
}

ztest_ondata <-function(data) { #id.exposure, id.outcome, b, se
  k<-data[, unique(id.exposure)] 
  arg <- combn(k, 2,  simplify = TRUE) %>% t(.) %>%
    as.data.frame(.) %>% as.data.table(.)
  
  arg <- expand_grid(arg, id_outcome = data[!grepl("dis-26-|trait-16-", id.outcome), unique(id.outcome)])
  setDT(arg)
  
  b_diff <- map(split(arg, 1:arg[,.N]), function(x) {
    b_x <- data[id.exposure == x$V1 & id.outcome == x$id_outcome, b]
    se_x <- data[id.exposure == x$V1 & id.outcome == x$id_outcome, se]
    b_y <- data[id.exposure == x$V2 & id.outcome == x$id_outcome, b]
    se_y <- data[id.exposure == x$V2 & id.outcome == x$id_outcome, se]
    dt <- z_test(b_x,se_x,b_y,se_y)
    dt[, V1 := x$V1]
    dt[,V2 := x$V2]
    dt[,id.outcome := x$id_outcome]
    return(dt)}) %>% rbindlist(.)
  
  return(b_diff[,.(V1, V2, id.outcome, diff, se, p_value)])
}

#######Abstract ########
data_tsmr[id.outcome %in% ctp$id, length(unique(id.outcome))]
data_tsmr[grepl("dis-26-", id.outcome), length(unique(id.outcome))]
k <- data_tsmr[id.exposure %in% c("ANGPTL4_trait-16-39", "LPL_trait-16-39") & id.outcome == "dis-13-1", ] %>%
ztest_ondata(.) 
k$p_value %>% formatC(., format = "e", digits = 1)
#########Results##########
#para 1 
data_tsmr[grepl("trait-16-39",id.exposure)& grepl("dis-", id.outcome),][!(grepl("dis-26-", id.outcome)),]
data_tsmr[hgnc == "HMGCR" & grepl("dis-", id.outcome),][!(grepl("dis-26-", id.outcome)),]
data_tsmr[hgnc == "PCSK9" & grepl("dis-", id.outcome),][!(grepl("dis-26-", id.outcome)),]
#para 2
data <- data_tsmr[grepl("trait-16-39",id.exposure)& id.outcome %in% c("dis-24-4", "trait-19-5", "dis-13-1"),]
b_diff <- ztest_ondata(data)
b_diff[p_value<0.05,]
data_tsmr[id.exposure == "ANGPTL4_trait-16-39" & id.outcome %in% c("dis-24-4"), ] %>% return_format_data
data_tsmr[id.exposure == "ANGPTL4_trait-16-39" & id.outcome %in% c("trait-19-5"), ] %>% return_format_data_noexp
data_tsmr[id.exposure == "ANGPTL4_trait-16-39" & id.outcome %in% c("dis-13-1"), ] %>% return_format_data
data_tsmr[id.outcome %in% c("trait-10-3", "trait-1-1") & pval < 0.05,]
#para3 
cor(scatter[,.SD,.SDcols = colnom], method = c("pearson"))

###para4
data_ctp <- data_tsmr[id.outcome %in% ctp$id, ]
data_ctp[,pvalfdr:=p.adjust(pval, method = "fdr"), by = "id.exposure"]
data_ctp[is_bad_when == "pos", effect_interpretation := fifelse(b>0, "bad", "good")]
data_ctp[is_bad_when == "neg", effect_interpretation := fifelse(b<0, "bad", "good")]
data_ctp[is_bad_when == "sig", effect_interpretation := fifelse(pvalfdr<0.05, "bad", "neutral")]
data_ctp[, effect_interpretation := fifelse(pvalfdr>0.05, "neutral", effect_interpretation)]

data_ctp[effect_interpretation=="bad" & idscale == "trait-16-39", .(hgnc, clean_variable_name, id.outcome, b, se, pval, pvalfdr)][order(pval)]
data_ctp[effect_interpretation=="good" & idscale == "trait-16-39", .(hgnc,category, clean_variable_name, id.outcome, b, se, pval, pvalfdr)][order(pval)]
data_ctp[effect_interpretation=="good" & idscale == "trait-16-39" & category == "Liver", .(hgnc,category, clean_variable_name, id.outcome, b, se, pval, pvalfdr)][order(pval)]

###para5
ntest<- data_tsmr[grepl("dis-26-", id.outcome), unique(clean_variable_name) %>% length]
phewas <- data_tsmr[grepl("dis-26-", id.outcome) & method == "Inverse variance weighted", ]
phewas[,unique(id.outcome) %>% length]
idtosel <- phewas[idscale=="trait-16-39" & b<0 & pval < 0.05/ntest, .N, by = "id.outcome"][N==3, id.outcome] #super concordant ID
phewas[idscale=="trait-16-39" & id.outcome %in% idtosel, ]
phewas[idscale=="trait-16-39" & id.outcome == "dis-26-5", ] %>% return_format_data

phewas[idscale=="trait-16-39" & b>0 & pval < 0.05/ntest, ]
phewas[hgnc == "LPL" & id.outcome == "dis-26-174", ] %>% return_format_data
phewas[idscale=="trait-16-39" & b>0 & pval < 0.05/ntest, ][, .SD[which.min(pval)], by = "id.exposure"] %>% return_format_data
phewas[idscale=="trait-16-39" & b>0 & pval < 0.05/ntest, ][id.outcome =="dis-26-1077", ] %>% return_format_data

phewas[hgnc=="ANGPTL4" & b>0 & pval < 0.05/10, ][which.min(pval), ] %>% return_format_data
phewas[hgnc=="APOC3" & b>0 & pval < 0.05/10, ][order(pval)][clean_variable_name == "Idiopathic pulmonary fibrosis", ] %>% return_format_data
phewas[hgnc=="PCSK9" & b>0 & pval < 0.05/10, ][which.min(pval), ] %>% return_format_data

inst <-fread("Data/Modified/inst.txt")
k<- inst[id.exposure == "LPL_trait-16-39", ]
gwasvcf::set_bcftools()
gwasvcf::set_plink()
outr<-GagnonMR::extract_outcome_variant(snps = k$SNP, outcomes = "/mnt/sda/gagelo01/Vcffile/Server_vcf/dis-4-1/dis-4-1.vcf.gz", rsq = parameters$proxy_rsq, parameters = parameters)
harm <- TwoSampleMR::harmonise_data(k, outr, 1)
res <- GagnonMR::all_mr_methods(harm)
res[method == "Inverse variance weighted"] %>% return_format_data()
##para6
k <- res_cox_small[!grepl("\\*", IV) & grepl("^PRS", exposure) & pval < 0.05,]
k[outcome == "CAD"]
k[outcome == "T2D"]
k <- res_cox_small[grepl(":", exposure) & pval < 0.05, ]
res_cox_small[IV=="PRS_HMGCR_trait_16_37*PRS_PCSK9_trait_16_37" & outcome == "T2D" & grepl("^PRS", exposure),]
k[exposure=="PRS_HMGCR_trait_16_37:PRS_PCSK9_trait_16_37",] %>% return_format_data_noexp_HR()
res_cox_small[exposure == "PRS_lplpathway" & IV == "PRS_lplpathway", ] %>% return_format_data_noexp_HR
res_cox_small[exposure == "PRS_ldlrpathway" & IV == "PRS_ldlrpathway", ] %>% return_format_data_noexp_HR

k <- res_continuous_small[!grepl("\\*", IV) & grepl("^PRS", exposure) & pval < 0.05,]
k <- res_continuous_small[grepl(":", exposure) & pval < 0.05, ]
k[pval < 0.05/10,.(exposure, outcome, b, se, pval)]

data<- res_cox_small[!grepl("\\*", IV) & grepl("^PRS", exposure),.(exposure, outcome, HR, pval)]
data[,b:=log(HR)]
data[,z := qnorm(pval/2)*sign(b)]
data[,pvaltest := 2*pnorm(q=z, lower.tail=FALSE)]
data[,se:=b/z]
setnames(data, c("exposure", "outcome"), c("id.exposure", "id.outcome"))
b_diff <- ztest_ondata(data)
b_diff[grepl("trait_16_39|lpl", V1) & grepl("trait_16_39|lpl", V2) & p_value<0.05/10,]
b_diff[grepl("trait_16_37|ldlr", V1) & grepl("trait_16_37|ldlr", V2) & p_value<0.05/10,]

res_continuous_small[outcome %in% c("Body_mass_index_BMI", "WHR") & !grepl("\\*", IV) & grepl("^PRS", exposure) & pval < 0.05/5, ] %>% return_format_data_noexp()
res_continuous_small[,id.exposure := exposure]
res_continuous_small[exposure == "PRS_HMGCR_trait_16_2" & IV == "PRS_HMGCR_trait_16_2" & outcome == "Body_mass_index_BMI", ] %>% return_format_data_noexp()
k <- res_continuous_small[outcome%in%c("Body_mass_index_BMI", "WHR") & grepl(":", exposure) & pval < 0.05,]
k <- k[,.(IV, outcome)] %>% unique
sign_interaction<- map(split(k, 1:k[,.N]), function(x) { 
  res_continuous_small[outcome== x$outcome & IV == x$IV,]
}) %>% rbindlist(., fill = TRUE)
sign_interaction

res_continuous_small[grepl("HMGCR", exposure) & grepl("lplpathway", exposure) & outcome == "Body_mass_index_BMI", ] %>% return_format_data_noexp(., digits = 3 )

data<- res_continuous_small[!grepl("\\*", IV) & grepl("^PRS", exposure) & outcome %in% c("WHR", "Body_mass_index_BMI"),.(exposure, outcome, b, se)]
setnames(data, c("exposure", "outcome"), c("id.exposure", "id.outcome"))
b_diff <- ztest_ondata(data)
b_diff[grepl("trait_16_4|lpl", V1) & grepl("trait_16_4|lpl", V2) & p_value<0.05,]
b_diff[grepl("trait_16_2|ldlr", V1) & grepl("trait_16_2|ldlr", V2) & p_value<0.05,]

##########Methods#############
k<-df_index[id %in% res_multicis[!grepl("dis-26-|trait-33-", id.outcome), id.outcome], ]
k[,.(id, clean_variable_name, ncase, ncontrol, sample_size)]
df_index[id%in%"dis-24-4",]
res_multicis[grepl("dis-26-", id.outcome), unique(id.outcome)]
df_index[grepl("dis-26-", id), max(as.numeric(sample_size))]
df_index[id%in%"trait-19-5",]

#######Discussion##########
# prar 1
data_tsmr[id.outcome %in% ctp$id, length(unique(id.outcome))]
data_tsmr[grepl("dis-26-", id.outcome), length(unique(id.outcome))]
data_tsmr[grepl("trait-33-", id.outcome), length(unique(id.outcome))]


#para3
cor(scatter[,.SD,.SDcols = colnom], method = c("pearson"))^2

data_tsmr[hgnc == "ANGPTL4" & id.outcome == "dis-13-1", ] %>% return_format_data()
data_tsmr[hgnc == "LPL" & id.outcome == "dis-13-1", ] %>% return_format_data()
ztest_ondata(data_tsmr[hgnc %in% c("ANGPTL4", "LPL") & id.outcome == "dis-13-1", ])[,p_value] %>% formatC(., format = "e", digits = 1)

