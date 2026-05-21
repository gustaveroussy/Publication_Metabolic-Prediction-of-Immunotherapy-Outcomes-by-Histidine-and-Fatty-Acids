#!/bin/bash

########################################################################
## using: sbatch /mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/script/02_launcher.sh
########################################################################

#SBATCH --job-name=B26026_LAZI_01
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=shortq

module load singularity/3.10.5

mkdir -p /mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/data_output/report/

singularity exec --contain -B /mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/ \
/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/script/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/Rapport_Rmd.simg \
Rscript -e 'rmarkdown::render("/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/script/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/02_report.Rmd", output_file = "/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/data_output/report/B26026_LAZI_01_report.html")'
