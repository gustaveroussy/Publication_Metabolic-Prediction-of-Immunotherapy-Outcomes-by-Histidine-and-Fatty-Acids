# This branch
This branch was created as a dedicated environment for code review and reproducibility assessment of the original analysis pipeline.  
The goal was to independently rerun the analyses, validate the reported results, and evaluate the robustness of the identified predictors.  
This separation from the main development repository ensures traceability of the review process and allows all modifications, tests, and validation steps performed during the review to be documented independently from the original workflow.  
The original code is under [the main branch](https://github.com/gustaveroussy/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/tree/main) of this repository.  

This analysis was performed by M. Aglave from BIGR platform (bigr@gustaveroussy.fr) in May 2026 based on the code and initial data, provided by L. Zitvogel's team.  

The review process was carried out in two stages:  
- regeneration of the results using the script 01_main_workflow_reviewed.R.  
- comparison of the results obtained with those presented in the publication in the notebook using the 02_report.Rmd script.  

To work in a controlled environment, the calculations were performed using singularity images on the Gustave Roussy HPC. They are available on Zenodo [here](https://zenodo.org/records/20324581).  
The 0X_launcher.sh script allows to run computation througth singularity images on the Gustave Roussy HPC.  

The report of the reviewing is [here](https://github.com/gustaveroussy/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/blob/reviewed/B26026_LAZI_01_report.html).  

---

# Overview

This repository provides a complete R implementation of an iterative ablation workflow for DynForest survival models applied to longitudinal metabolomics and clinical data.

The workflow was developed starting from the original DynForest codebase (Devaux et al., 2023) and extended to enable stepwise variable ablation: at each iteration, the predictor (metabolite or fixed covariate) with the lowest variable importance (VIMP) is removed, and the model is retrained.

The repository contains three main R scripts:

## File 1 – Main workflow (run_DynForest_ablation_workflow.R)
Orchestrates the iterative ablation pipeline, including:

- Running ablation loops across iterations.

- Computing internal performances (Out-of-Bag Integrated Brier Score, AUC).

- Performing external validations on independent cohorts.

- Generating visualizations (IBS plots, AUC curves, VIMP grids, alluvial plots).

## File 2 – Core functions (main_functions.R)
Defines the training and ablation functions:

- train_DynForest_model() → fits a single DynForest model with longitudinal and fixed predictors.

- run_dynforest_ablation_iteration() → main iterative ablation loop that removes the worst predictor at each step and stores results.

- extract_predictors_over_iterations() → tracks predictor sets across iterations.

## File 3 – Evaluation functions (evaluation_performances_functions.R)
Provides performance assessment utilities:

- get_AUC_model_based_on_pred() → internal AUC with confidence intervals from timeROC.

- predict_DynForest_model() → wrapper for generating patient-level predictions.

- extract_auc_df() → aggregates AUCs across iterations.

- eval_dynforest_on_external_cohort() and eval_dynforest_cindex_on_external_cohort() → external validation using AUC and C-index.

- eval_dynforest_on_internal_cohort_with_bootstrap() → internal validation via bootstrap resampling.

- plot_auc_evolution() → visualization of AUC trends across iterations.

# Workflow Summary

1. Training phase – Build DynForest survival models using longitudinal metabolite trajectories and fixed clinical covariates.

2. Ablation iterations – Iteratively drop the least informative predictor (lowest VIMP) until no variables remain or performance deteriorates.

3. Performance evaluation –

- Internal: OOB IBS, AUC at clinical landmarks (e.g., 6–12 months).
- External: Validation on independent cohorts
- Bootstrap: Robust estimates of internal predictive performance.

4. Visualization – Generate IBS curves, ROC/AUC plots, VIMP grids, and alluvial diagrams showing predictor elimination paths.

# Requirements

- R (≥ 4.2)
- Main dependencies:  
DynForest  
timeROC  
nftbart  
ggplot2, patchwork, ggalluvial  
dplyr, tidyverse  
pec, prodlim  
foreach, doParallel  

# Usage

1. Prepare your dataset:

- Longitudinal variable: distDebImmuno_months (time from ICB start to sampling).

- Patient ID: id_numeric (numeric).

- Outcome: survival object (time_chr = PFS/OS, event_chr = PD/Death).

- Predictors: metabolites + fixed covariates.

2. Run the main workflow (File 1) with your dataset:

results <- run_dynforest_ablation_iteration(
  data = my_data,
  time_chr = "PFS",
  event_chr = "PD",
  columns_metabolites_init = my_metabolites,
  fixed_vars_init = c("Age", "Gender", "BMI", "Tumor", "Stage"),
  output_folder = "path/to/results"
)


3. Evaluate and visualize performances using functions from File 3.

# Citation

If you use this workflow, please cite:

Devaux, A., Helmer, C., Genuer, R., & Proust-Lima, C. (2023). Random survival forests with multivariate longitudinal endogenous covariates. Statistical methods in medical research, 32(12), 2331–2346. https://doi.org/10.1177/09622802231206477

Suissa D, et al. (in preparation). Metabolic Prediction of Immunotherapy Outcomes by Histidine and Fatty Acids.
