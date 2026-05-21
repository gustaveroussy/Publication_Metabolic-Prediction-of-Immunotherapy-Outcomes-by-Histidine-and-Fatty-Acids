#!/bin/bash

########################################################################
## using: sbatch /mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/script/01_launcher.sh
########################################################################

#SBATCH --job-name=B26026_LAZI_01
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --partition=longq

module load singularity/3.10.5

singularity exec --no-home -B /mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/ \
/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/script/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/R_SuissaD_2026.simg \
Rscript /mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/script/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/01_main_workflow_reviewed.R

echo "Finish!"
