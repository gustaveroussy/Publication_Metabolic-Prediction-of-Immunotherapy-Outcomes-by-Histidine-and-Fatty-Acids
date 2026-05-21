# =====================================================
# --- SET ENVIRONMENT ---
# =====================================================
# Load libraries
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

# Set output path
project_path <- "/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/"
dir.create(paste0(project_path,"data_output/"), showWarnings = FALSE, recursive = TRUE)

# Load main functions
source(paste0(project_path,"script/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/main_functions.r"))
source(paste0(project_path,"script/Publication_Metabolic-Prediction-of-Immunotherapy-Outcomes-by-Histidine-and-Fatty-Acids/evaluation_performances_functions.r"))

# =====================================================
# --- LOAD INITAL DATA ---
# =====================================================

# Load data (replace comma by point and convert into numeric automatically)
data_for_ablation_iteration <- read_excel(paste0(project_path,"com/input/Table_Repo_Hist_Group_cohorts_VF2.xlsx"), sheet = 2, skip = 1)

# Convert ID into numeric and rename columns
data_for_ablation_iteration$id_numeric <- as.numeric(gsub("Pt", "", data_for_ablation_iteration$ID))
data_for_ablation_iteration$distDebImmuno_months <- data_for_ablation_iteration$`distDebImmuno (Months)`
colnames(data_for_ablation_iteration) <- gsub("\\(", "", colnames(data_for_ablation_iteration))
colnames(data_for_ablation_iteration) <- gsub("\\)", "", colnames(data_for_ablation_iteration))
data_for_ablation_iteration <- data_for_ablation_iteration %>%
  rename_with(
    ~ .x %>%
      str_remove(" AreaQCCorrLog2") %>%  # remove suffix
      str_replace_all(" ", "_") %>%      # spaces -> _
      paste0("Metabolite_", .),          # add préfix
    .cols = contains("AreaQCCorrLog2")
  )
colnames(data_for_ablation_iteration) <- gsub("[^a-zA-Z0-9_]", "_", colnames(data_for_ablation_iteration))
data_for_ablation_iteration$PFS <- data_for_ablation_iteration$PFS_Months
data_for_ablation_iteration$PD <- data_for_ablation_iteration$Progression

# Check data formatting
head(data_for_ablation_iteration)
str(data_for_ablation_iteration)
summary(data_for_ablation_iteration)

# =====================================================
# --- RUN MAIN FUNCTION: MODEL BUILDING ---
# =====================================================
metabolites_names <- grep("Metabolite", colnames(data_for_ablation_iteration), value = TRUE)
results_ablation_iteration_SABR_IML1 <- run_dynforest_ablation_iteration(data = data_for_ablation_iteration, 
                                                                         time_chr = "PFS_Months", event_chr = "Progression",
                                                                         columns_metabolites_init = metabolites_names,
                                                                         fixed_vars_init = c("Treatment_type","Gender", "Age", "BMI", "Tumor", "Stage"), 
                                                                         output_folder = paste0(project_path,"data_output/result_training/"))

#save(results_ablation_iteration_SABR_IML1, file = paste0(project_path,"data_output/result_training/results_ablation_iteration_SABR_IML1.RData")) #should be the same as "/mnt/beegfs01/scratch/bioinfo_core/B26026_LAZI_01/data_output/result_training/ablation_partial_results.rds"
results_ablation_iteration_SABR_IML1 <- readRDS(paste0(project_path,"data_output/result_training/ablation_partial_results.rds"))

# =================================================================
# --- CALCULATE PERFORMANCES EVALUATION based on Out-Of-Bag IBS ---
# =================================================================

res_dyn_OOB_all <- list()
for (iteration_num in 145:158) {
  cat(paste0(">>> Traitement de iteration_", iteration_num, "...\n"))
  dynforest_obj = results_ablation_iteration_SABR_IML1[[paste0("iteration_", iteration_num)]]$model
  IBS.min <-0
  IBS.max <- NULL
  cause <- 1
  ncores <- 2
  if (!methods::is(dynforest_obj,"dynforest")){
    cli_abort(c(
      "{.var dynforest_obj} must be a dynforest object",
      "x" = "You've supplied a {.cls {class(dynforest_obj)}} object"
    ))
  }
  if (dynforest_obj$type=="surv"){
    if (is.null(IBS.max)){
      IBS.max <- max(dynforest_obj$data$Y$Y[,1])
    }
  }
  rf <- dynforest_obj
  Longitudinal <- rf$data$Longitudinal
  Numeric <- rf$data$Numeric
  Factor <- rf$data$Factor
  timeVar <- rf$timeVar
  Y <- rf$data$Y
  ntree <- ncol(rf$rf)
  # ncores
  if (is.null(ncores)==TRUE){
    ncores <- parallel::detectCores()-1
  }
  Inputs <- c(Longitudinal$type, Numeric$type, Factor$type)
  err <- rep(NA,length(unique(Y$id)))
  i <- NULL
  Longitudinal_courant <- NULL
  Numeric_courant <- NULL
  Factor_courant <- NULL
  allTimes <- sort(unique(c(0,Y$Y[,1])))
  if (is.null(IBS.max)){IBS.max <- max(allTimes)}
  Y.surv <- data.frame(id = Y$id,time.event = Y$Y[,1],event = ifelse(Y$Y[,2]==cause,1,0))
  Y.surv <- Y.surv[order(Y.surv$time.event, -Y.surv$event),]
  allTimes_IBS <- allTimes[which(allTimes>=IBS.min&allTimes<=IBS.max)]
  Y.surv <- Y.surv[which(Y.surv$time.event>=IBS.min),]
  # KM estimate of the survival function for the censoring
  G <- pec::ipcw(formula = Surv(time.event, event) ~ 1,data = Y.surv,method = "marginal", times=allTimes_IBS, subjectTimes=allTimes_IBS)
  OOB_IBS <- sort(Y$id[which(Y$Y[, 1] >= IBS.min)])
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  library(foreach)
  library(doParallel)
  pck <- .packages()
  dir0 <- find.package(pck)
  dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
  parallel::clusterExport(cl,list("pck","dir"),envir=environment())
  parallel::clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))
  res.oob <- foreach(i=1:length(OOB_IBS), .combine='c', .multicombine = TRUE
                     #,.packages = c("pec", "prodlim")
  ) %dopar%
    {
      # for (i in 1:length(OOB_IBS)){
      indiv <- OOB_IBS[i]
      w_y <- which(Y$id==indiv)
      if (is.element("Longitudinal",Inputs)==TRUE){
        if (IBS.min==0){
          w_XLongitudinal <- which(Longitudinal$id==indiv) # all measurements until IBS.min
        }else{
          w_XLongitudinal <- which(Longitudinal$id==indiv&Longitudinal$time<=IBS.min) # only measurements until IBS.min
        }
        Longitudinal_courant <- list(type="Longitudinal", X=Longitudinal$X[w_XLongitudinal,, drop=FALSE], id=Longitudinal$id[w_XLongitudinal], time=Longitudinal$time[w_XLongitudinal],
                                     model=Longitudinal$model)
      }
      if (is.element("Numeric",Inputs)==TRUE){
        w_XNumeric <- which(Numeric$id==indiv)
        Numeric_courant <- list(type="Numeric", X=Numeric$X[w_XNumeric,, drop=FALSE], id=Numeric$id[w_XNumeric])
      }
      if (is.element("Factor",Inputs)==TRUE){
        w_XFactor <- which(Factor$id==indiv)
        Factor_courant <- list(type="Factor", X=Factor$X[w_XFactor,, drop=FALSE], id=Factor$id[w_XFactor])
      }
      pred.mat <- matrix(NA, nrow = ncol(rf$rf), ncol = length(allTimes_IBS))
      for (t in 1:ncol(rf$rf)){
        BOOT <- rf$rf[,t]$boot
        oob <- setdiff(unique(Y$id), BOOT)
        if (is.element(indiv, oob)== TRUE){
          pred_node <- tryCatch(DynForest:::pred.MMT(rf$rf[,t], Longitudinal = Longitudinal_courant,Numeric = Numeric_courant, Factor = Factor_courant,timeVar = timeVar),
                                error = function(e) return(NA))
          pred_node_chr <- as.character(pred_node)
          if (is.na(pred_node_chr)){
            pred.mat[t,] <- NA
          }else{
            if (IBS.min == 0){
              pi_t <- rf$rf[,t]$Y_pred[[pred_node_chr]][[as.character(cause)]]$traj[which(allTimes%in%allTimes_IBS)]
              pred.mat[t,] <- pi_t
            }else{
              pi_t <- rf$rf[,t]$Y_pred[[pred_node_chr]][[as.character(cause)]]$traj[which(allTimes%in%allTimes_IBS)] # pi(t)
              pi_s <- rf$rf[,t]$Y_pred[[pred_node_chr]][[as.character(cause)]]$traj[sum(allTimes<IBS.min)] # pi(s)
              s_s <- 1 - sum(unlist(lapply(rf$rf[,t]$Y_pred[[pred_node_chr]], FUN = function(x){
                return(x$traj[sum(allTimes<IBS.min)])
              }))) # s(s)
              pred.mat[t,] <- (pi_t - pi_s)/s_s # P(S<T<S+t|T>S)
            } } } }
      oob.pred <- apply(pred.mat, 2, mean, na.rm = TRUE)
      if (all(is.nan(oob.pred))) {
        return(NA)
      }
      # IPCW
      Wi_event <- (ifelse(Y$Y[w_y,1] <= allTimes_IBS, 1, 0)*ifelse(Y$Y[w_y,2]!=0,1,0))/(G$IPCW.subjectTimes[which(Y.surv$id==indiv)])
      Wi_censored <- ifelse(Y$Y[w_y,1] > allTimes_IBS, 1, 0)/(G$IPCW.times)
      Wi <- Wi_event + Wi_censored
      # Individual Brier Score
      Di <- ifelse(Y$Y[w_y,1] <= allTimes_IBS, 1, 0)*ifelse(Y$Y[w_y,2]==cause,1,0) # D(t) = 1(s<Ti<s+t, event = cause)
      pec.res <- list()
      pec.res$AppErr$matrix <- Wi*(Di-oob.pred)^2 # BS(t)
      pec.res$models <- "matrix"
      pec.res$time <- allTimes_IBS
      class(pec.res) <- "pec"
      err <- pec::ibs(pec.res, start = IBS.min, times = max(allTimes_IBS))[1] # IBS
      print(err)
    }
  parallel::stopCluster(cl)
  
  cat(paste0("IBS médian iteration_", iteration_num, ": ", round(median(res.oob, na.rm=TRUE), 4), "\n"))
  res_dyn_OOB_all[[paste0("iteration_", iteration_num)]] <- res.oob
  cat(paste0(">>> Résultats stockés pour iteration_", iteration_num, "\n"))
}

# Print or save results
print(res_dyn_OOB_all)
# Create data.frame based on results
oob_df <- data.frame(
  iteration = as.numeric(gsub("iteration_", "", names(res_dyn_OOB_all))),
  IBS = as.numeric(unlist(res_dyn_OOB_all)))
# Increasing order of iteration
oob_df <- oob_df %>% arrange(iteration)
# Filter
oob_df <- oob_df %>%filter(iteration >=145)
# Save IBS summary table
oob_summary <- data.frame(
  iteration = names(res_dyn_OOB_all),
  IBS_median = sapply(res_dyn_OOB_all, function(x) median(x, na.rm = TRUE)),
  IBS_mean   = sapply(res_dyn_OOB_all, function(x) mean(x, na.rm = TRUE)),
  IBS_sd     = sapply(res_dyn_OOB_all, function(x) sd(x, na.rm = TRUE))
) %>% arrange(iteration_num)
print(oob_summary)
dir.create(paste0(project_path,"data_output/performances_evaluation/"), showWarnings = FALSE, recursive = TRUE)
write.csv(oob_summary, file = paste0(project_path, "data_output/performances_evaluation/IBS_summary_per_iteration.csv"), row.names = FALSE)
# Plot
p <- ggplot(oob_df, aes(x = iteration, y = IBS)) +
  geom_point(size = 4, color = "steelblue") +
  geom_line(size = 1, color = "steelblue") +
  geom_text(aes(label = round(IBS, 3)), vjust = -0.8, size = 3.5, show.legend = FALSE) +
  labs( title = "Out-Of-Bag Integrated Brier Score per iteration",x = "iteration Number",y = "Integrated Brier Score (IBS)" ) +
  theme_minimal(base_size = 18) +
  theme(  plot.title = element_text(hjust = 0.5, face = "bold", size = 18),axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"),panel.grid.major = element_line(color = "gray80"), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = oob_df$iteration) +
  ylim(0, max(oob_df$IBS) * 1.1)
ggsave(paste0(project_path,"data_output/performances_evaluation/IBS_ablation_model_SABR_IML1.pdf"), p, width = 8, height = 5)



# --- PERFORMANCES EVALUATION: Calculation of internal AUC --- # 
results_external_IML1_SABR_ablation_iteration_SABR_IML1 <- list()
last_valid_iteration <- NA
for (iteration_num in c(145:158)) {
  cat(paste0("\n===== Running external validation on iteration ", iteration_num, " =====\n"))
  result <- NULL
  tryCatch({
    result <- eval_dynforest_on_external_cohort(
      results_list = results_ablation_iteration_SABR_IML1,
      iteration_number = iteration_num,
      external_data = data_for_ablation_iteration,
      study_name = c("IML1_SABR"),
      time_chr = "PFS",
      event_chr = "PD",
      clinical_times = c(3, 6, 12, 18, 24),
      plot_path = paste0(project_path,"data_output/performances_evaluation/AUC_external_IML1_model_IML1_SABR", iteration_num, ".pdf"))
  }, error = function(e) { 
      cat("Error iteration", iteration_num, ":", e$message, "\n")
      return(NULL)  
  })
  if (is.null(result)) {
    cat(paste0("🛑 Early stopping: iteration ", iteration_num, " fails, go to the next.\n"))
    next
  }
  print(result$AUC_table)
  auc_tbl <- result$AUC_table
  auc_tbl_clean <- auc_tbl %>% filter(!is.na(AUC) & !is.na(CI_lower))
  valid <- any(auc_tbl_clean$AUC > 70 & auc_tbl_clean$CI_lower > 50)
  if (!valid) {
    cat(paste0("🛑 Early stopping: iteration ", iteration_num, " fails AUC/CI thresholds.\n"))
  }
  # If valid, download
  results_external_IML1_SABR_ablation_iteration_SABR_IML1[[paste0("iteration_", iteration_num)]] <- result
  last_valid_iteration <- iteration_num
  print(result$AUC_table)
}
if (!is.na(last_valid_iteration)) {
  cat(paste0("✅ Last valid iteration: ", last_valid_iteration, "\n"))
} else {
  cat("❌ No valid iteration found.\n")
}

# Transform into a dataframe
df_auc <- imap_dfr(results_external_IML1_SABR_ablation_iteration_SABR_IML1, ~{
  .x$AUC_table %>%
    mutate(iteration = .y)
})

# Ensure right types
df_auc <- df_auc %>%
  mutate(
    Time_months = as.numeric(Time_months),
    AUC = as.numeric(AUC),
    CI_lower = as.numeric(CI_lower),
    CI_upper = as.numeric(CI_upper))
# Save table
write.csv(df_auc[df_auc$Time_months %in% c(6,12),], file = paste0(project_path, "data_output/performances_evaluation/validation_on_IML1_SABR_model_IML1_SABR_AUC_t6_t12.csv"), row.names = FALSE)
# Plot
p <- ggplot(df_auc, aes(x = Time_months, y = AUC, color = iteration)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.6, size = 0.7) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = "External Validation of DynForest Ablation iterations on PANDORE",
    x = "Time (months)",
    y = "AUC (%)",
    color = "iteration"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
ggsave(paste0(project_path,"data_output/performances_evaluation/validation_on_IML1_SABR_model_IML1_SABR.pdf"), p, width = 10, height = 6)



# --- PERFORMANCES EVALUATION: Calculation of external AUC --- # 
# Load external data
data_for_ablation_iteration_PANDORE <- read_excel(paste0(project_path,"com/input/Table_Repo_Hist_Group_cohorts_VF2.xlsx"), sheet = 3, skip = 1)

# Convert ID into numeric and rename columns
data_for_ablation_iteration_PANDORE$id_numeric <- as.numeric(gsub("Pt", "", data_for_ablation_iteration_PANDORE$ID))
data_for_ablation_iteration_PANDORE$distDebImmuno_months <- data_for_ablation_iteration_PANDORE$`distDebImmuno (Months)`
colnames(data_for_ablation_iteration_PANDORE) <- gsub("\\(", "", colnames(data_for_ablation_iteration_PANDORE))
colnames(data_for_ablation_iteration_PANDORE) <- gsub("\\)", "", colnames(data_for_ablation_iteration_PANDORE))
data_for_ablation_iteration_PANDORE <- data_for_ablation_iteration_PANDORE %>%
  rename_with(
    ~ .x %>%
      str_remove(" AreaQCCorrLog2") %>%  # remove suffix
      str_replace_all(" ", "_") %>%      # spaces -> _
      paste0("Metabolite_", .),          # add préfix
    .cols = contains("AreaQCCorrLog2")
  )
colnames(data_for_ablation_iteration_PANDORE) <- gsub("[^a-zA-Z0-9_]", "_", colnames(data_for_ablation_iteration_PANDORE))
data_for_ablation_iteration_PANDORE$PFS <- data_for_ablation_iteration_PANDORE$PFS_Months
data_for_ablation_iteration_PANDORE$PD <- data_for_ablation_iteration_PANDORE$Progression

# Compute performance
results_external_PANDORE_ablation_iteration_SABR_IML1 <- list()
last_valid_iteration <- NA
for (iteration_num in c(145:158)) {
  cat(paste0("\n===== Running external validation on iteration ", iteration_num, " =====\n"))
  result <- NULL
  tryCatch({
    result <- eval_dynforest_on_external_cohort(
      results_list = results_ablation_iteration_SABR_IML1,
      iteration_number = iteration_num,
      external_data = data_for_ablation_iteration_PANDORE,
      study_name = c("PANDORE"),
      time_chr = "PFS",
      event_chr = "PD",
      clinical_times = c(3, 6, 12, 18, 24),
      plot_path = paste0(project_path,"data_output/performances_evaluation/AUC_external_PANDORE_model_IML1_SABR", iteration_num, ".pdf"))
  }, error = function(e) { 
    cat("Error iteration", iteration_num, ":", e$message, "\n")
    return(NULL)  
  })
  if (is.null(result)) {
    cat(paste0("🛑 Early stopping: iteration ", iteration_num, " fails, go to the next.\n"))
    next
  }
  print(result$AUC_table)
  auc_tbl <- result$AUC_table
  auc_tbl_clean <- auc_tbl %>% filter(!is.na(AUC) & !is.na(CI_lower))
  valid <- any(auc_tbl_clean$AUC > 70 & auc_tbl_clean$CI_lower > 50)
  if (!valid) {
    cat(paste0("🛑 Early stopping: iteration ", iteration_num, " fails AUC/CI thresholds.\n"))
  }
  # If valide, download
  results_external_PANDORE_ablation_iteration_SABR_IML1[[paste0("iteration_", iteration_num)]] <- result
  last_valid_iteration <- iteration_num
  print(result$AUC_table)
}
if (!is.na(last_valid_iteration)) {
  cat(paste0("✅ Last valid iteration: ", last_valid_iteration, "\n"))
} else {
  cat("❌ No valid iteration found.\n")
}
# Transform into a dataframe
df_auc <- imap_dfr(results_external_PANDORE_ablation_iteration_SABR_IML1, ~{
  .x$AUC_table %>%
    mutate(iteration = .y)
})

# Ensure right types
df_auc <- df_auc %>%
  mutate(
    Time_months = as.numeric(Time_months),
    AUC = as.numeric(AUC),
    CI_lower = as.numeric(CI_lower),
    CI_upper = as.numeric(CI_upper))
# Save table
write.csv(df_auc[df_auc$Time_months %in% c(6,12),], file = paste0(project_path, "data_output/performances_evaluation/validation_on_PANDORE_model_IML1_SABR_AUC_t6_t12.csv"), row.names = FALSE)
# Plot
p <- ggplot(df_auc, aes(x = Time_months, y = AUC, color = iteration)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.6, size = 0.7) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = "External Validation of DynForest Ablation iterations on PANDORE",
    x = "Time (months)",
    y = "AUC (%)",
    color = "iteration"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
ggsave(paste0(project_path,"data_output/performances_evaluation/validation_on_PANDORE_model_IML1_SABR.pdf"), p, width = 10, height = 6)

# ===========================
# --- Produce VIMP Plots ---
# ===========================

library(patchwork)
vimp_plots <- list()
# Loop on iterations 145–158 ## depending on your number of predictors
for (iteration_num in 145:158) {
  cat(paste0(">>> Generate VIMP plot for iteration_", iteration_num, "...\n"))
  model <- results_ablation_iteration_SABR_IML1[[paste0("iteration_", iteration_num)]]$model
  vimp_result <- DynForest::compute_vimp(model)
  p <- plot(x = vimp_result, PCT = TRUE) +
    ggtitle(paste0("iteration ", iteration_num)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) 
  vimp_plots[[paste0("iteration_", iteration_num)]] <- p
}
# Combine plots into a grid
grid_plot <- wrap_plots(vimp_plots, ncol = 4)
# Sauvegarder sur Desktop
dir.create(paste0(project_path,"data_output/VIMP/"), showWarnings = FALSE, recursive = TRUE)
ggsave(paste0(project_path,"data_output/VIMP/vimp_grid_145_158.pdf"), grid_plot, width = 20, height = 18)






# ===========================
# --- Produce Alluvial Plots ---
# ===========================
library(tidyverse)
library(ggalluvial)
df_mets_results_ablation_iteration_SABR_IML1 <- extract_predictors_over_iterations(results_ablation_iteration_SABR_IML1)
df_alluvial_results_ablation_iteration_SABR_IML1 <- df_mets_results_ablation_iteration_SABR_IML1 %>%
  mutate(Predictor = factor(Predictor),
         iteration_num = parse_number(iteration)) # Extract iteration number
# Identify last 6 dernières iterations
#last_iterations_num <- sort(unique(df_alluvial_results_ablation_iteration_SABR_IML1$iteration_num), decreasing = TRUE)[1:6] %>% sort()
last_iterations_num <- sort(unique(df_alluvial_results_ablation_iteration_SABR_IML1$iteration_num), decreasing = TRUE)[1:31] %>% sort()
last_iterations <- paste0("iteration_", last_iterations_num)
df_alluvial_results_ablation_iteration_SABR_IML1 <- df_alluvial_results_ablation_iteration_SABR_IML1 %>%
  mutate(iteration = factor(iteration, levels = last_iterations),
         alluvium = Predictor) %>% filter(iteration %in% last_iterations)

df_alluvial_results_ablation_iteration_SABR_IML1$Predictor <- gsub("Metabolite_", "", df_alluvial_results_ablation_iteration_SABR_IML1$Predictor)
df_alluvial_results_ablation_iteration_SABR_IML1$alluvium <- gsub("Metabolite_", "", df_alluvial_results_ablation_iteration_SABR_IML1$alluvium)


p <- ggplot(df_alluvial_results_ablation_iteration_SABR_IML1,
            aes(x = iteration, stratum = Predictor, alluvium = alluvium,
                fill = Predictor, label = Predictor)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.8) +
  geom_stratum(width = 0.6, height = 1.3) +
  geom_text(stat = "stratum", size = 2.3, color = "#443232") +
  theme_minimal() +
  guides(fill = FALSE) +
  labs(title = "Predictor flow across last 10 ablation iterations")
dir.create(paste0(project_path,"data_output/Alluvial/"), showWarnings = FALSE, recursive = TRUE)
ggsave(paste0(project_path,"data_output/Alluvial/df_alluvial_results_ablation_iteration_SABR_IML1.pdf"), p, width = 30, height = 6)


