# ------------------------------
# --- CORE FUNCTIONS ---
# ------------------------------

# Training function for a single iteration of the ablation moodel
train_DynForest_model <- function(data, time_chr, event_chr, columns_metabolites, fixed_vars) {
  timeData_train <- data[, c("id_numeric", "distDebImmuno_months", columns_metabolites)]
  timeVarModel_train <- lapply(columns_metabolites, function(var) {
    list(fixed = as.formula(paste0("`", var, "` ~ distDebImmuno_months")), random = ~ distDebImmuno_months)
  })
  names(timeVarModel_train) <- columns_metabolites
  timeData_train <- as.data.frame(timeData_train)
  
  fixedData_train <- unique(data[, c("id_numeric", fixed_vars),drop=F])
  fixedData_train <- if (!is.null(fixed_vars) && length(fixed_vars) > 0) {
                        as.data.frame(unique(data[, c("id_numeric", fixed_vars), drop = FALSE]))
                     } else NULL
  df <- as.data.frame(unique(data[, c("id_numeric", time_chr, event_chr)]))
  Y <- list(type = "surv", Y = df)
  res_dyn_train <- dynforest(
    timeData = timeData_train,
    fixedData = fixedData_train,
    timeVar = "distDebImmuno_months",
    idVar = "id_numeric",
    timeVarModel = timeVarModel_train,
    Y = Y,
    ntree = 50, nodesize = 5, minsplit = 5,
    cause = 1, ncores = 10, seed = 1234
  )
  return(res_dyn_train)
}

# Main function with removal of the worst predictor (metabolite OR fixed)
run_dynforest_ablation_iteration <- function(data, time_chr, event_chr, columns_metabolites_init, fixed_vars_init, output_folder, t_landmarks = c(6, 12)) {
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  all_results <- list()
  current_metabolites <- columns_metabolites_init
  current_fixed_vars <- fixed_vars_init
  stop_condition <- FALSE
  iteration_index <- 1

  while (!stop_condition) {
    saveRDS(all_results, file = file.path(output_folder, "ablation_partial_results.rds"))

    cat(paste0("\n===== Starting iteration ", iteration_index, " with ", length(current_metabolites), " metabolites and ", length(current_fixed_vars), " fixed vars =====\n"))

        # 🛑 Check stop condition immediately
      if (length(current_metabolites) == 0) {
        cat("[!] No metabolites remaining, stopping.\n")
        break
      }
    model <- train_DynForest_model(data, time_chr, event_chr, current_metabolites, current_fixed_vars)
    vimp <- DynForest::compute_vimp(model)
    vimp_data <- plot(vimp, PCT = TRUE)$data
        # Predictions + AUC
    predictions <- predict_DynForest_model(model, time_chr, event_chr, data, landmark = 0, columns_metabolites_here = current_metabolites, fixed_vars = current_fixed_vars)
    aucs <- list()
    for (t_lm in t_landmarks) {
      plot_path <- file.path(output_folder, paste0("ROC_train_iteration", iteration_index, "_t", t_lm, "months.pdf"))
      auc_info <- get_AUC_model_based_on_pred(predictions, data, time_chr, event_chr, landmark = t_lm, plot_path)
      aucs[[paste0("t", t_lm)]] <- auc_info
    }
    auc_vals <- sapply(aucs, function(x) as.numeric(x$AUC))
    lower_ci <- sapply(aucs, function(x) as.numeric(x$lower_CI))
    all_results[[paste0("iteration_", iteration_index)]] <- list(
      model = model,
      metabolites = current_metabolites,
      fixed_vars = current_fixed_vars,
      aucs = aucs
    )
    if (length(current_metabolites)==0) {
      cat("[!] Too few predictors remaining, stopping.\n")
      stop_condition <- TRUE
    } else {
      iteration_index <- iteration_index + 1
    }
    # Identify the worst global predictor (metabolite or fixed var)
    worst_pred <- vimp_data$var[which.min(vimp_data$vimp)]
    cat(paste0("Removed worst predictor: ", worst_pred, "\n"))
    if (worst_pred %in% current_metabolites) {
      current_metabolites <- setdiff(current_metabolites, worst_pred)
    } else if (worst_pred %in% current_fixed_vars) {
      if (length(current_fixed_vars)==1){
        current_fixed_vars <- NULL
      }else{
        current_fixed_vars <- setdiff(current_fixed_vars, worst_pred)
      }
    } else {
      cat("[!] Predictor not in known lists, skipping update\n")
    }
  }
  return(all_results)
}

extract_predictors_over_iterations <- function(results_ablation) {
  all_rows <- list()
  for (iteration_name in names(results_ablation)) {
    # Retrieve current predictors
    mets <- results_ablation[[iteration_name]]$metabolites
    fixed <- results_ablation[[iteration_name]]$fixed_vars
    df_mets <- data.frame(
      iteration = iteration_name,
      Predictor = mets,
      Class = "Metabolite",
      stringsAsFactors = FALSE)
    if (!is.null(fixed) && length(fixed) > 0) {
      df_fixed <- data.frame(
        iteration = iteration_name,
        Predictor = fixed,
        Class     = "Fixed",
        stringsAsFactors = FALSE)
      all_rows[[iteration_name]] <- rbind(df_mets, df_fixed)
    } else {
      all_rows[[iteration_name]] <- df_mets
    }
  }
  final_df <- do.call(rbind, all_rows)
  return(final_df)
}
