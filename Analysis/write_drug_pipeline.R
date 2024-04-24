#!/usr/bin/env Rscript
library(GagnonMR)

setwd("/mnt/sda/gagelo01/Projects/2024/LPL_pathway/")
myArray <- GagnonMR::list_masterscript("/mnt/sda/gagelo01/Projects/Pipelines/Robust_efficient_drug_pipeline/Analysis/")
myArray <- myArray[c(2:5,7, 8)] #c(2:5,7, 8)
GagnonMR::write_masterscript(myArray = myArray, mysource = "tosourcedrug.R", output_name = "Analysis/drug_pipeline.sh")
message("This script finished without errors")