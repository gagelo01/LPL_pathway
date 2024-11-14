#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

wd<-"/mnt/sda/gagelo01/Projects/2024/LPL_pathway/"
setwd(wd)
gwasvcf::set_bcftools()
###chose instrument
gwasvcf::set_plink()
ldref<-"/home/couchr02/Mendel_Commun/Christian/LDlocal/big_EUR_rs"
x <- paste0("awk -F ' ' '{print $2}' ","/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", ".bim")
dt_ref_snp <- fread(cmd = x, header = FALSE) #Those are SNP in 1000G with MAF>0.01
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
if(file.exists("Data/Modified/dt_gene_region.txt")){dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")}
ncores<-40
options(future.globals.maxSize= 5e9, future.globals.onReference = "warning")
plan(multicore, workers = ncores, gc = TRUE)

#####select inst#####
should_select_pan<-FALSE
######change parameters#########
######
parameters <- GagnonMR::default_param()
parameters$path <- c(GagnonMR::default_param()$path, paste0(wd, "/Data/Modified/Gtex_vcf/"))
parameters$uni_cis_minpval <- 1
parameters$ldref <- ldref
parameters$snp_bim <- dt_ref_snp$V1 #I only include SNP in 1000G with maf > 0.01
##########Choose the gene you wish to include##########
gene_toinclude <- c("LPL", "ANGPTL4", "APOC3", "ANGPTL3", "PCSK9", "HMGCR")
should_QTL_instrument_selection<-FALSE
thewindow<-1e5
typeof_sel<- "multicis_independent_clumping"
############Choose the outcome you wish to include  ###########
ctp <- fread("/mnt/sda/gagelo01/Projects/2024/UKB_phewas/Data/Modified/clinical_toxicity_panel.txt")

ID_server_out <- c( "trait-1-1", "dis-13-1", "dis-24-4", "trait-10-3", "trait-19-5",
                    ctp$id,
                    df_index[grepl("trait-33-", id), id],
                    df_index[grepl("trait-16", id) & grepl("^HDL|^LDL|^logTG", trait) & population %in% c("Mixed") & sex == "Males and Females"& consortium == "GLGC", id],
                    df_index[grepl("dis-26-", id) & ncase>1000,id]) %>% unique(.)

# df_index[grepl("trait-33-", id), id], #I wnat to have the impact on metabolites to do a genetic mimicry analysis comparision of LPL with ANGPTL4 and APOC3.
# # df_index[grepl("dis-21-", id) & ncase>1000,id],
# "trait-19-2", "trait-19-5",
# "trait-7-2", "trait-6-1","dis-14-7", "dis-18-1", "dis-19-1",
# "dis-2-2", "trait-14-8", "trait-14-10", paste0("trait-27-", 1:2), "trait-2-2", "trait-2-4", "dis-23-2",
# "trait-12-2", "dis-7-1", "trait-25-16", "trait-25-1", "trait-29-20", "dis-12-1",
# "trait-13-1", "trait-13-2", "dis-20-1", "dis-1-1", "dis-16-1","dis-4-1", "trait-7-1",
# df_index[grepl("dis-14-", id) & population == "European",id]

out_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_out, "/", ID_server_out, ".vcf.gz")
ID_mrbase_out <- NULL
out_mrbase<-NULL

###### #####
k <- fread("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/Tribal/Data/Modified/dt_gene_region.txt")
k <- k[hgnc%in%gene_toinclude, ]
k[, gene_region := paste0(chr, ":", (start-thewindow)%>%ifelse(.<1, 1, .), "-", (end+thewindow))]
k[hgnc %in% c("LPL", "ANGPTL4", "APOC3", "ANGPTL3"), ID:= df_index[grepl("glgc", tolower(consortium)) & clean_variable_name == "Triglyceride" & sex == "Males and Females" & population == "Mixed"& consortium == "GLGC", ]$id]
k[hgnc %in% c("PCSK9", "HMGCR"), ID:= df_index[grepl("glgc", tolower(consortium)) & clean_variable_name == "LDL cholesterol" & sex == "Males and Females" & population == "Mixed"& consortium == "GLGC", ]$id]
arguments_exposure_proxies <- unique(k[,.(hgnc, ID, gene_region)])
arguments_exposure_proxies[, vcffile_inst := paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID, "/", ID, ".vcf.gz")]

######split_large_outcome#####
split_outcome<-TRUE
all_mr_methods_short = FALSE
all_mr_methods_skip_presso = FALSE
tsmr_robust<-c("mr_weighted_median",  "mr_egger_regression")

####Dummy######
study_to_selectinst<-list()

make.dir <- function(fp) {
  if(!file.exists(fp)) {  # If the folder does not exist, create a new one
    dir.create(fp)
  } else {   # If it existed, delete and replace with a new one  
    unlink(fp, recursive = TRUE)
    dir.create(fp)
  }
}
