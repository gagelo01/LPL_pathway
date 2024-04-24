#!/bin/bash
start_time=$(date +%s)
myArray=(/mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/drug_pipeline.sh /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/interaction_phewas.sh /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/1a_interaction_LPL.R /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/3c_plot_all_figures.R /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/4a_writetables.R /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/5a_write_text.R)
for str in ${myArray[@]}; do
chmod u+x $str
done
echo 'Initializing /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/drug_pipeline.sh' && /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/drug_pipeline.sh  && echo 'Initializing /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/interaction_phewas.sh' && /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/interaction_phewas.sh  && echo 'Initializing /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/1a_interaction_LPL.R' && /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/1a_interaction_LPL.R  && echo 'Initializing /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/3c_plot_all_figures.R' && /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/3c_plot_all_figures.R  && echo 'Initializing /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/4a_writetables.R' && /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/4a_writetables.R  && echo 'Initializing /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/5a_write_text.R' && /mnt/sda/gagelo01/Projects/2024/LPL_pathway/Analysis/5a_write_text.R  && echo 'The master script finished without errors'
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "The script took $elapsed_time seconds to run."
