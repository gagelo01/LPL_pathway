#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(openxlsx)
wd<-"/mnt/sda/gagelo01/Projects/2024/LPL_pathway/"
setwd(wd)
system(paste0("/mnt/sda/gagelo01/Projects/Pipelines/plot_and_tables/Analysis/2a_create_table.R", " ", wd, " ", "res_multicis_independent"))
list_supdat <- readRDS("Data/Modified/suptable/list_supdat.Rdata")
col_description <- readRDS("Data/Modified/suptable/col_description.Rdata")
dt_title <- fread("Data/Modified//suptable/dt_title.txt")
goodorder<-c(1,2,5,6,7,8,3,4,9,10)
dt_title <- dt_title[goodorder,]
col_description <- col_description[goodorder]
list_supdat <- list_supdat[goodorder]

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

#cox
res_cox<-fread("Data/Modified/LPL/res_cox.txt")
res_continuous<-fread("Data/Modified/LPL/res_continuous.txt")

res_cox_small <- res_cox[grepl("PRS", exposure),]
res_cox_small[, c("ph_chisq", "ph_df", "ph_p"):=NULL]
res_cox_small <- res_cox_small[cov_inc== "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10", ]
res_cox_small <- res_cox_small[str_count(string = IV, pattern = "\\*")<2,]
res_cox_small <-  res_cox_small[order(outcome, exposure)]
colnom <- c("HR", "lci", "uci")
res_cox_small[, (colnom):=lapply(.SD, function(x) 1/x), .SDcols = colnom]
setnames(res_cox_small, c("lci", "uci"), c("uci", "lci"))

res_continuous_small <- res_continuous[grepl("PRS", exposure),]
res_continuous_small <- res_continuous_small[cov_inc== "+ age_enrollment + sex + PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + PCA9 + PCA10", ]
res_continuous_small <- res_continuous_small[str_count(string = IV, pattern = "\\*")<2,]
res_continuous_small <-  res_continuous_small[order(outcome, exposure)]
res_continuous_small[,b:=b*-1]
res_continuous_small[,lci:=b-1.96*se]
res_continuous_small[,uci:=b+1.96*se]

res_coloc <- fread("Data/Modified/dt_coloc.txt")
res_coloc<- merge(res_coloc, df_index[, .(id,clean_variable_name)], by.x = "id.outcome", by.y = "id", all.x = TRUE)
res_coloc <- merge(res_coloc, dt_gene_region[,.(id, hgnc, sample_size)], by.x = "id.exposure", by.y = "id", all.x = TRUE)
res_coloc[,exposure := hgnc]
res_coloc[,outcome := clean_variable_name]
res_coloc[,c("hgnc", "clean_variable_name", "sample_size"):=NULL]
#########
list_supdat <- c(list_supdat[1:4], list(res_coloc), list_supdat[5:6], list(pearson_cor), list_supdat[7:length(list_supdat)], list(res_continuous_small, res_cox_small))

#########
dt_title2 <- data.table(title = paste0("ST", c(5,8,13,14)),
                        caption = c( "Colocalisation results",
                                     "Correlation matrix of the effect of different lipid targets on metabolites",
                                     "Results of cox regression analysis testing the association between GRS of lipid lowering targer and their interaction on CAD and T2D.",
                                     "Results of the linear regression analysis testing the association between GRS of lipid lowering targer and their interaction on BMI and lipids."))

dt_title <- rbindlist(list(dt_title[1:4,], dt_title2[1,], dt_title[5:6,], dt_title2[2,], dt_title[7:.N], dt_title2[3:.N,]))
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
  "exposure", "the exposure from which the statistics are reported, the 'sum_PRS' sum the APOC3, ANGPTL4, and LPL PRS",
  "HR", "Hazard ratio",
  "pval", "p-value",
  "lci", "lower 95% confidence interval",
  "uci", "upper 95% confidence interval",
  "IV", "the independent variable (including the interaction term)",
  "outcome", "the outcome",
  "cov_inc", "the covariate that were included in the model",
) %>% as.data.table(.)

col_description2[[4]] <- tribble(
  ~x, ~y,
  "exposure", "the exposure from which the statistics are reported. the 'sum_PRS' sum the APOC3, ANGPTL4, and LPL PRS.",
  "b", "the effect size",
  "pval", "p-value",
  "lci", "lower 95% confidence interval",
  "uci", "upper 95% confidence interval",
  "IV", "the independent variable (including the interaction term)",
  "outcome", "the outcome",
  "cov_inc", "the covariate that were included in the model",
) %>% as.data.table(.)


col_description <- c(col_description[1:4], col_description2[1], col_description[5:6], col_description2[2], col_description[7:length(col_description)], col_description2[3:length(col_description2)])


bold_st <- createStyle(textDecoration = "Bold")
wb <- createWorkbook()
for(i in 1:length(list_supdat)) {
  print(i)
  addWorksheet(wb, sheetName =  dt_title[i, title])
  
  title <- dt_title[i,paste0(title, " : ", caption)]
  writeData(wb, sheet = i, x = title, startCol = 1, startRow = 1)
  addStyle(wb, sheet = i, style = bold_st, rows = 1, cols = 1:2)
  writeData(wb, sheet = i, x = col_description[[i]], startCol = 1, startRow = 2)
  addStyle(wb, sheet = i, style = bold_st, rows = 2:col_description[[i]][,.N+2], cols = 1)
  deleteData(wb, sheet = i, rows = 2, cols = 1:2, gridExpand = TRUE)
  writeData(wb, sheet = i, x = list_supdat[[i]], startCol = 1, startRow = col_description[[i]][,.N+4])
  addStyle(wb, sheet = i, style = bold_st, rows = col_description[[i]][,.N+4], cols = 1:length(colnames(list_supdat[[i]])), gridExpand = TRUE, stack = TRUE)
}
saveWorkbook(wb, file = "Results/supplementary_tables.xlsx", overwrite = TRUE)


########Table 1######
##########Table 1###########
table1 <- res_cox_small[grepl("PRS_ldlrpathway", IV) & grepl("PRS", exposure) & exposure != "PRS_ldlrpathway", ]
table1[,exposure := gsub("_trait_16_4", "", exposure)] 
format_data_noexp_HR <-function(data, x= "HR", digits = 2) {
  k<-data[, paste0(format(round(get(x), digits = digits), nsmall = digits), " (", format(round(lci, digits = digits), nsmall = digits), " to ",  format(round(uci, digits = digits), nsmall = digits), "), p=",pval %>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(outcome)]
  return(k)
}
table1$result <- format_data_noexp_HR(table1)

####second try#####
k<-data.table(table1)
k[,id := gsub(":PRS_ldlrpathway", "", exposure)]
k1 <- k[grepl(":", exposure),.(id, outcome, pval)]
k1[,pval:=pval %>% formatC(., format = "e", digits = 1)]
setnames(k1, "pval", "pval_interaction")
k <- merge(k[!grepl(":", exposure)], k1, by = c("id", "outcome"))
k<- dcast(k, id ~ outcome, value.var = c("result", "pval_interaction"))

top <- data.table(id  = c("", "Exposure"),   result_CAD  = c("CAD", "HR (95%CI)"), pval_interaction_CAD = c("CAD","Interaction p-value"),  result_T2D = c("T2D", "HR (95%CI)"), pval_interaction_T2D = c("T2D", "Interaction p-value"))
k <- rbind(top, k)


bold_st <- createStyle(textDecoration = "Bold", halign = "center", valign = "center")
row2_st <- createStyle(textDecoration = "Bold", border = "bottom",  halign = "center", valign = "center")

wb <- createWorkbook()
addWorksheet(wb, sheetName =  "Table1")
writeData(wb, sheet = 1, x = k, startCol = 1, startRow = 1, colNames = FALSE, borders = openxlsx_getOp("borders"))
mergeCells(wb, sheet = 1, cols = c(2,3), rows = 1)
mergeCells(wb, sheet = 1, cols = 4:5, rows = 1)
addStyle(wb, sheet = 1, style = bold_st, rows = 1, cols = 1:ncol(k))
addStyle(wb, sheet = 1, style = row2_st, rows = 2, cols = 1:ncol(k))
setColWidths(wb, sheet = 1, cols = 1:ncol(k), widths = "auto")
saveWorkbook(wb, "Results/Table1.xlsx", overwrite = TRUE)

res_coloc <- fread("Data/Modified/dt_coloc.txt")
write.xlsx(res_coloc, "Data/Modified/res_coloc.xlsx")

message("This script finished without errors")
