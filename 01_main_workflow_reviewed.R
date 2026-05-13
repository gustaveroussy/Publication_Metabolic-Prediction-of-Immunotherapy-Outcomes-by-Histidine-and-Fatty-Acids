library(readxl)
library(DynForest)
library(timeROC)
library(nftbart)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(dplyr)
library(tidyverse)
library(pec)
library(prodlim)
library(foreach)
library(doParallel)
library(stringr)

project_path <- "/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/"

# Load main functions
source(paste0(project_path,"script/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/main_functions.r"))
source(paste0(project_path,"script/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/evaluation_performances_functions.r"))


# Load data (replace comma by point and convert into numeric automatically)
df <- read_excel(paste0(project_path,"com/input/Table_Repo_Hist_Group_cohorts_VF2.xlsx"), sheet = 2, skip = 1)

# Convert progression into character, ID into numeric and rename columns
df$id_numeric <- as.numeric(gsub("Pt", "", df$ID))
df$distDebImmuno_months <- df$`distDebImmuno (Months)`
colnames(df) <- gsub("\\(", "", colnames(df))
colnames(df) <- gsub("\\)", "", colnames(df))
df <- df %>%
  rename_with(
    ~ .x %>%
      str_remove(" AreaQCCorrLog2") %>%  # enlever le suffixe
      str_replace_all(" ", "_") %>%            # espaces -> _
      paste0("Metabolite_", .),                # ajouter préfixe
    .cols = contains("AreaQCCorrLog2")
  )
colnames(df) <- gsub("[^a-zA-Z0-9_]", "_", colnames(df))
df$PFS <- df$PFS_Months
df$PD <- df$Progression

#convert character into factor
#char_cols <- sapply(df, is.character)
#df[char_cols] <- lapply(df[char_cols], as.factor)

# les IDs doivent être des entiers consécutifs à partir de 1 :
#id_mapping <- data.frame(
#  id_original = sort(unique(df$id_numeric)),
#  id_new      = seq_along(unique(df$id_numeric))
#)
#df$id_numeric <- id_mapping$id_new[match(df$id_numeric, id_mapping$id_original)]

# Check data formating
head(df)
str(df)
summary(df)

# Run
metabolites_names <- grep("Metabolite", colnames(df), value = TRUE)
results_ablation_iteration_SABR_IML1 <- run_dynforest_ablation_iteration(data = df, 
                                                                         time_chr = "PFS_Months", event_chr = "Progression",
                                                                         columns_metabolites_init = metabolites_names,
                                                                         fixed_vars_init = c("Treatment_type","Gender", "Age", "BMI", "Tumor", "Stage"), 
                                                                         output_folder = "/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/data_output/result_training/")

save(results_ablation_iteration_SABR_IML1, file = paste0(project_path,"data_output/result_training/results_ablation_iteration_SABR_IML1.RData"))

#results_ablation_iteration_SABR_IML1 <- readRDS("/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/data_output/result_training/ablation_partial_results.rds")







