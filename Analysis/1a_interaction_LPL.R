#!/usr/bin/env Rscript
library("bigsnpr")
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)
library(survival)
library(lubridate)
library(survminer)

wd<-"/mnt/sda/gagelo01/Projects/2024/LPL_pathway/"
setwd(wd)

dictionnary <- fread("Data/Modified/dictionnary.txt")
dat <- fread("Data/Modified/data_cox.txt")

########
k <- colnames(dat)
k <- k[grepl("^PRS", k)] %>% unique
k<-factor(k, levels = k[order(k)])
k<-k[order(k)]

k<-expand_grid(prs1= k, prs2= k)
setDT(k)
k <- k[!apply(k, 1, is.unsorted), ]
k[,hgnc1:=str_match(prs1, "^PRS_\\s*(.*?)\\s*_")[,2]]
k[,hgnc2:=str_match(prs2, "^PRS_\\s*(.*?)\\s*_")[,2]]
k <- k[hgnc1!=hgnc2,]
k[,IV_vec := paste(as.character(prs1), as.character(prs2), sep = "*")]
##########
prs_list<-list(PRS_lplpathway = c("PRS_ANGPTL4_trait_16_39", "PRS_APOC3_trait_16_39", "PRS_LPL_trait_16_39"),
PRS_ldlrpathway = c("PRS_HMGCR_trait_16_37", "PRS_PCSK9_trait_16_37"))
map(1:length(prs_list), function(i) { 
dat[, (names(prs_list)[i]) := apply(.SD, 1, function(x) sum(x, na.rm = TRUE)),.SDcols = prs_list[[i]]]
dat[,(names(prs_list)[i]):=lapply(.SD, function(x) (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)) ,.SDcols = names(prs_list)[i]]
})
dat[,WHR:=Waist_circumference/Hip_circumference]
##########
IV_sum <- c(paste(names(prs_list), collapse = "*"), names(prs_list), "PRS_lplpathway*PRS_HMGCR_trait_16_37", "PRS_lplpathway*PRS_PCSK9_trait_16_37",
            paste(prs_list$PRS_lplpathway, "PRS_ldlrpathway", sep =  "*"))
IV_vec <- c(unique(as.character(c(k$prs1,k$prs2))), k$IV_vec,IV_sum)
DV_vec <- c("CAD", "T2D")
cov_inc_vec <-c("+ age_enrollment + sex", " ", paste(c("+ age_enrollment", "sex", paste0("PCA", 1:10)), collapse =  " + "))

arguments <- tidyr::expand_grid(IV =IV_vec, DV=DV_vec, cov_inc = cov_inc_vec) %>% unique(.)

options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)


run_coxph_wrapper_safely<-safely(GagnonMR::run_coxph_wrapper)

res_cox <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  run_coxph_wrapper_safely(dat = dat[get(paste0(x$DV, "_toinclude"))==1,], IV = x$IV, DV = x$DV, cov_inc = x$cov_inc)
}, .options = furrr_options(seed = TRUE))

res_cox <- map(res_cox, function(x) x$result)%>%rbindlist(.,fill = TRUE)
fwrite(res_cox, "Data/Modified/LPL/res_cox.txt")

setDT(arguments)
res_logit <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  k<-glm(paste0(x$DV, "_censored ~", x$IV, x$cov_inc),data =  dat[get(paste0(x$DV, "_toinclude"))==1,], family = "binomial")
  k <- broom::tidy(k) %>% as.data.table(.)
  colnames(k)<-c("exposure", "b", "se", "todump", "pval")
  k[,todump:=NULL]
  k$IV <- x$IV
  k$outcome <- x$DV
  k$cov_inc <- x$cov_inc
  return(k)}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.)


####continous outcomes#####
DV_vec <- c(dictionnary[var_name%in%colnames(dat) & ValueType == "Continuous",]$var_name, "WHR")
arguments <- tidyr::expand_grid(IV =IV_vec, DV=DV_vec, cov_inc = cov_inc_vec) 
res_continuous <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  k<-glm(paste0(x$DV, "~", x$IV, x$cov_inc),data =  dat[toinclude==1,], family = "gaussian")
  k <- broom::tidy(k) %>% as.data.table(.)
  colnames(k)<-c("exposure", "b", "se", "todump", "pval")
  k[,todump:=NULL]
  k$IV <- x$IV
  k$outcome <- x$DV
  k$cov_inc <- x$cov_inc
  return(k)}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.)


fwrite(res_continuous, "Data/Modified/LPL/res_continuous.txt")
fwrite(res_logit, "Data/Modified/LPL/res_logit.txt")
fwrite(dat, "Data/Modified/LPL/data_cox.txt")
message("This script finished without errors")

