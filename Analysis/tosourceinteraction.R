#!/usr/bin/env Rscript
library("bigsnpr")
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)
library(lubridate)

wd<-"/mnt/sda/gagelo01/Projects/2024/LPL_pathway/"
setwd(wd)

#########The instruments#########
instrument_path<-"Data/Modified/inst.txt"

#########The outcomes ####
ukbpath<-"/home/couchr02/Mendel_UKB/Source/Phenotype/February_2023_Update/ukb671338.tab"
pca_path<-"/home/couchr02/Mendel_UKB/Source/Phenotype/October_30_2019_Refresh/ukb38266.tab"
fieldid<- c("21001", #BMI
            "30760", #HDL-C
            "30870", #TG
            "30780", "48", "30640", "49")
######The diseases definition#####
arg_cad <- data.table(code_inclusion_icd10 = list(c("I21", "I22", "I23", "I241", "I252")),
                      code_inclusion_opcs = list(c(paste0("K40", 1:4), paste0("K41", 1:4), paste0("K45", 1:5), paste0("K49", c(1:2,8:9), "K502", paste0("K75", c(1:4, 8, 9))))),
                      code_exclusion_icd10 = list(NULL),
                      code_exclusion_opcs = list(NULL),
                      outcome_name = "CAD")

arg_t2d <- data.table(code_inclusion_icd10 = list(paste0("E1",1:4)),
                      code_inclusion_opcs = list(NULL),
                      code_exclusion_icd10 = list(NULL),
                      code_exclusion_opcs = list(NULL),
                      outcome_name = "T2D")


disease_definition <- rbindlist(list(arg_cad,arg_t2d))
disease_definition[,nrows:=Inf]


########date_end_followup######
date_end_followup  <- "2019-12-31"
date_end_followup <- ymd(date_end_followup) %>% as.IDate(.)
####see derivation in ("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/Pancreatitis_UKB_Baass/Analysis/2d_grs_pheno.R")