#!/bin/bash
start_time=$(date +%s)
myArray=(/mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2a_extract_snp_individual_info.R /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2b_computeGRS.R /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2c_ukb_pheno.R /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2d_UKBdisease.R /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2e_getcoxdat.R)
for str in ${myArray[@]}; do
chmod u+x $str
done
echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2a_extract_snp_individual_info.R' && /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2a_extract_snp_individual_info.R tosourceinteraction.R && echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2b_computeGRS.R' && /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2b_computeGRS.R tosourceinteraction.R && echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2c_ukb_pheno.R' && /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2c_ukb_pheno.R tosourceinteraction.R && echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2d_UKBdisease.R' && /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2d_UKBdisease.R tosourceinteraction.R && echo 'Initializing /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2e_getcoxdat.R' && /mnt/sda/gagelo01/Projects/Pipelines/Interaction_phewas/Analysis/2e_getcoxdat.R tosourceinteraction.R && echo 'The master script finished without errors'
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "The script took $elapsed_time seconds to run."
