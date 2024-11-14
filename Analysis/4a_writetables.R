#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(openxlsx)

wd<-"/mnt/sda/gagelo01/Projects/2024/LPL_pathway/"
setwd(wd)
res_data <- fread("Data/Modified/res_multicis_independent.txt")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
df_index[id=="trait-10-3", clean_variable_name := "Waist-to-hip ratio"]
dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")
ctp <- fread("/mnt/sda/gagelo01/Projects/2024/UKB_phewas/Data/Modified/clinical_toxicity_panel.txt")
source("/mnt/sda/gagelo01/Projects/Pipelines/plot_and_tables/Analysis/plot_figures_function.R")
source("/mnt/sda/gagelo01/Projects/Pipelines/plot_and_tables/Analysis/create_tables_function.R")
inst <- fread("Data/Modified/inst.txt" )

res_data <- prepare_res_data(res_data = res_data, should_inverse_association = TRUE)
res_data_list <- res_data_tolist(res_data, ctp = ctp)
inst <- format_inst(inst)
dataset <- format_df_index(df_index, res_data)

list_supdat <- c(list(dataset = dataset, Instrument = inst), res_data_list)
names(list_supdat)<-paste0("ST", 1:length(list_supdat))
k <- map(list_supdat[3:length(list_supdat)], function(x) x[1,.(type_outcome, type_method)]) %>% rbindlist()
dt_title <- data.table(title = paste0("ST", 1:length(list_supdat)),
                       caption = c( "GWAS data sets used","Instruments", paste0("Mendelian randomization estimates using ", k$type_method, " for ", k$type_outcome)))


col_description<- vector(mode = "list", length = length(list_supdat))
col_description[[1]] <- return_col_description(description = "dataset")
col_description[[2]] <- return_col_description(description = "inst")
generic<-return_col_description(description = "generic")

col_description[c(3,5,7,9)] <- map(1:4, function(x) rbind(generic, data.table(x="pvalfdr", y="The false discovery rate corrected p-value for the number of tests")))
col_description[c(4,6,8,10)] <- map(1:4, function(x) generic)

goodorder<-c(1,2,5,6,7,8,9,10)
dt_title <- dt_title[goodorder,]
col_description <- col_description[goodorder]
list_supdat <- list_supdat[goodorder]
##############
res_multicis <- fread("Data/Modified/res_multicis_independent.txt")
res_multicis[,c("hgnc", "idscale") := tstrsplit(id.exposure, "_")]

df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")
dt<-data.table(id=res_multicis[!grepl("^eqtl-|^prot-", id.exposure), unique(id.exposure)])
dt[, c("hgnc", "id_small"):=tstrsplit(id, "_", fixed = TRUE)]
dt_gene_region <- merge(dt_gene_region, dt[,.(id, id_small)], by = "id", all.x = TRUE)
dt_gene_region[is.na(id_small), id_small:=id]
dt_gene_region <- merge(dt_gene_region, df_index[, .(id,sample_size)], by.x = "id_small", by.y = "id", all.x = TRUE)
#cor
k1 <- res_multicis[grepl("trait-33-", id.outcome)&method%in%c("Inverse variance weighted"), ]
scatter <- dcast(k1, id.outcome ~ hgnc, value.var = c("b", "se"))
colnom<-colnames(scatter)[grepl("^b_", colnames(scatter))]
scatter[,(colnom):=lapply(.SD, function(x) x*-1) , .SDcols = colnom]
pearson_cor <- cor(scatter[,.SD,.SDcols = colnom], method = c("pearson"))
dt_pearson_cor <- pearson_cor %>% as.data.table(., row.names = TRUE)
dt_pearson_cor[,rownames:=rownames(pearson_cor)]
dt_pearson_cor <- dt_pearson_cor[, .SD,.SDcols = c("rownames",rownames(pearson_cor))]
#cox
res_logit<-fread("Data/Modified/LPL/res_logit.txt")
res_continuous<-fread("Data/Modified/LPL/res_continuous.txt")

res_logit <- res_logit[grepl("PRS", exposure),]
res_logit[, c("ph_chisq", "ph_df", "ph_p"):=NULL]
res_logit <- res_logit[cov_inc== "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10", ]
res_logit <- res_logit[str_count(string = IV, pattern = "\\*")<2,]
res_logit <- res_logit[!grepl("PRS_ldlrpathway|PRS_lplpathway", IV),]
res_logit <- res_logit[,c("exposure1", "exposure2") := tstrsplit(IV, "\\*")]
res_logit <- rbind(res_logit[is.na(exposure2)],res_logit[grepl("trait_16_39", exposure1) & grepl("trait_16_37", exposure2) | 
            grepl("trait_16_37", exposure1) & grepl("trait_16_39", exposure2), ])
res_logit[,c("exposure1", "exposure2") := NULL]
res_logit <-  res_logit[order(outcome, exposure)]


res_continuous_small <- res_continuous[grepl("PRS", exposure) & outcome %in% c("LDL_direct", "Triglycerides"),]
res_continuous_small <- res_continuous_small[cov_inc== "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10", ]
res_continuous_small <- res_continuous_small[str_count(string = IV, pattern = "\\*")==0,]
res_continuous_small <- res_continuous_small[!(exposure %in% c("PRS_ldlrpathway","PRS_lplpathway")),]
res_continuous_small <- res_continuous_small[order(outcome, exposure)]


res_coloc <- fread("Data/Modified/dt_coloc.txt")
res_coloc<- merge(res_coloc, df_index[, .(id,clean_variable_name)], by.x = "id.outcome", by.y = "id", all.x = TRUE)
res_coloc <- merge(res_coloc, dt_gene_region[,.(id, hgnc, sample_size)], by.x = "id.exposure", by.y = "id", all.x = TRUE)
res_coloc[,exposure := hgnc]
res_coloc[,outcome := clean_variable_name]
res_coloc[,c("hgnc", "clean_variable_name", "sample_size"):=NULL]
#########
list_supdat <- c(list_supdat[1:4], list(res_coloc), list_supdat[5:6], list(pearson_cor), list(res_continuous_small), list(res_logit), list_supdat[7:length(list_supdat)])
names(list_supdat)<-paste0("ST",1:length(list_supdat))
#########
dt_title2 <- data.table(title = paste0("ST", c(5,8,9,10)),
                        caption = c( "Colocalisation results",
                                     "Correlation matrix of the effect of different lipid targets on metabolites",                                   
                                     "Results of the linear regression analysis testing the association between GRS of lipid lowering targer and their interaction on LDL-c and TG.",
                                     "Results of the logistic regression analysis testing the association between GRS of lipid lowering target and their interaction on CAD and T2D."
                                     ))

dt_title <- rbindlist(list(dt_title[1:4,], dt_title2[1,], dt_title[5:6,], dt_title2[2,], dt_title2[3:.N,], dt_title[7:.N]))
dt_title[,title := paste0("ST", c(1:.N))]

col_description2<- vector(mode = "list", length = dt_title2[,.N])
col_description2[[1]] <- tribble(
  ~x, ~y,
  "id.exposure", "id of the exposure",
  "id.outcome", "id of the outcome",
  "exposure", "the name of the exposure",
  "outcome", "the name of the outcome",
  "posprob_coloc_PPH0", "posterior probability of no causal variants in tthe two traits",
  "posprob_coloc_PPH1", "posterior probability of only trait one with causal variant",
  "posprob_coloc_PPH2", "posterior probability of only trait two with causal variant",
  "posprob_coloc_PPH3", "posterior probability of both traits with distinct causal variant",
  "posprob_coloc_PPH4", "posterior probability of both traits with same causal variant",
  "posprob_colocH4.SNP", "the SNP prioritised as shared causal variant",
  "posprob_colocH4.SNPexplained_var", "the variance explained of the priotirtised shared causal variant",
) %>% as.data.table(.)


col_description2[[2]] <- tribble(
  ~x, ~y,
  "row", "lipid lowering targets",
  "column", "lipid lowering targets",
) %>% as.data.table(.)

col_description2[[3]] <- tribble(
  ~x, ~y,
  "exposure", "the exposure from which the statistics are reported.",
  "b", "the effect size",
  "pval", "p-value",
  "lci", "lower 95% confidence interval",
  "uci", "upper 95% confidence interval",
  "IV", "the independent variable (including the interaction term)",
  "outcome", "the outcome",
  "cov_inc", "the covariate that were included in the model",
) %>% as.data.table(.)

col_description2[[4]] <- tribble(
  ~x, ~y,
  "exposure", "the exposure from which the statistics are reported",
  "b", "log(OR)",
  "pval", "p-value",
  "lci", "lower 95% confidence interval",
  "uci", "upper 95% confidence interval",
  "IV", "the independent variable (including the interaction term)",
  "outcome", "the outcome",
  "cov_inc", "the covariate that were included in the model",
) %>% as.data.table(.)




########
col_description <- c(col_description[1:4], col_description2[1], col_description[5:6], col_description2[2], col_description2[3:length(col_description2)], col_description[7:length(col_description)])

save_openxlsx(list_supdat=list_supdat, col_description=col_description, dt_title, file_name = "Results/supplementary_tables.xlsx")

message("this script finished without errors")
