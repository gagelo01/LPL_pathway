#!/usr/bin/env Rscript
library(GagnonMR)

setwd("/mnt/sda/gagelo01/Projects/2024/LPL_pathway/")
myArray <- GagnonMR::list_masterscript("/mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/")
GagnonMR::write_masterscript(myArray = myArray[1:5], mysource = "tosourceinteraction.R", output_name = "Analysis/interaction_phewas.sh")
message("This script finished without errors")