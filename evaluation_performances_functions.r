# =======================================================================
# --- PERFORMANCE EVALUATION: Main AUC based on internal predictions ---
# =======================================================================

# time_chr can be either "OS" or "PFS", "PFS" was used for the paper, hence "PD" for event_chr
get_AUC_model_based_on_pred <- function(pred, data, time_chr, event_chr, landmark, output_pdf_path){
    # Transform predictions to scale 0-100
    pred_indiv_100 <- pred$pred_indiv * 100
    # Clean data based on time variable
    data <- data %>% filter(!is.na(.data[[time_chr]]), !is.na(.data[[event_chr]]))
    # Keep only individuals present in pred_indiv_100
    data <- subset(data, as.character(id_numeric) %in% rownames(pred_indiv_100))
    # Select the first measure per individual
    data <- data %>% dplyr::group_by(id_numeric) %>%dplyr::slice(which.min(distDebImmuno_months)) %>%dplyr::ungroup()
    # Prepare structures to store results
    AUC_s <- c()
    lower <- c()
    upper <- c()
    list_indexes <- numeric()
    # Identify indexes of corresponding times
    for (k in seq(landmark, 60, 1)) { # change landmark if needed, in months here
        list_indexes[k] <- which.min(abs(pred$times - k))
    }
    # Open PDF file to save plots
    pdf(output_pdf_path, width = 8, height = 6)
    a <- 1
    for (j in seq_along(na.omit(list_indexes))) {
        tryCatch({
            # Compute AUC with timeROC
            ROC <- timeROC(
                T = data[[time_chr]],
                delta = data[[event_chr]],
                marker = pred_indiv_100[, as.numeric(list_indexes)[j]] / 100,
                cause = 1,  
                weighting = "marginal", 
                times = a,  
                ROC = TRUE,
                iid = TRUE)
            # Extract AUC and CI values
            auc_value <- ROC["AUC"]$AUC[2][[1]]
            lower_ci <- confint(ROC)$CI_AUC[[1]]
            upper_ci <- confint(ROC)$CI_AUC[[2]]
            AUC_s <- c(AUC_s, auc_value)
            lower <- c(lower, lower_ci)
            upper <- c(upper, upper_ci)
            plot(ROC, time = a, main = paste0("ROC Curve at Time", a," months"))
            legend("bottomright",  legend = paste0("AUC: ", round(auc_value, 3), "\n95% CI: [", round(lower_ci, 3), "-", round(upper_ci, 3), "]"), bty = "n")
        }, error = function(e) {
            message("Error in the ROC calculation at t=", a, ": ", e$message)
        })
        a <- a + 1
    }
    dev.off()
    return(list(AUC = AUC_s, lower_CI = lower, upper_CI = upper))
}

predict_DynForest_model <- function(DynModel, time_chr, event_chr, data, landmark, columns_metabolites_here, fixed_vars) {
  # Select patients to predict
  data <- as.data.frame(data)
  id_pred <- as.numeric(unique(data$id_numeric))
  data_pred_IML <- data[data$id_numeric %in% id_pred, ]
  # Longitudinal data
  timeData_pred_IML <- data_pred_IML[, c("id_numeric", "distDebImmuno_months", columns_metabolites_here)]
  # Fixed data
  fixedData_pred_IML <- if (!is.null(fixed_vars) && length(fixed_vars) > 0) {
                            unique(data_pred_IML[, c("id_numeric", fixed_vars)])
                        } else NULL
  # Prediction
  pred_dyn <- predict(
    object = DynModel,
    timeData = timeData_pred_IML,
    fixedData = fixedData_pred_IML,
    idVar = "id_numeric",
    timeVar = "distDebImmuno_months",
    t0 = landmark)
  return(pred_dyn)
}

extract_auc_df <- function(results_list) {
  df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
    iteration_name <- names(results_list)[i]
    aucs <- results_list[[i]]$aucs
    metabolites <- paste(results_list[[i]]$metabolites, collapse = ";")
    fixed_vars <- paste(results_list[[i]]$fixed_vars, collapse = ";")
    do.call(rbind, lapply(names(aucs), function(tname) {
      auc_val <- aucs[[tname]]$AUC
      lower_ci <- aucs[[tname]]$lower_CI
      upper_ci <- aucs[[tname]]$upper_CI
      # Ensure lengths match (otherwise error)
      len <- length(auc_val)
      data.frame(
        iteration = rep(iteration_name, len),
        Time = rep(tname, len),
        AUC = auc_val,
        CI_lower = lower_ci,
        CI_upper = upper_ci,
        Metabolites = rep(metabolites, len),
        Fixed_Vars = rep(fixed_vars, len),
        stringsAsFactors = FALSE)}))}))
  return(df)
}

eval_dynforest_on_external_cohort <- function(
    results_list,
    iteration_number,
    external_data,
    study_name,
    time_chr,
    event_chr,
    clinical_times = c(3, 6, 12, 18, 24),
    plot_path = NULL
) {

    external_data <- as.data.frame(external_data)
  # --- Extract model, metabolites and fixed_vars from this iteration ---
  model <- results_list[[paste0("iteration_", iteration_number)]]$model
  metabolites <- results_list[[paste0("iteration_", iteration_number)]]$metabolites
  fixed_vars <- results_list[[paste0("iteration_", iteration_number)]]$fixed_vars
  print(fixed_vars)
  # --- Subset external cohort ---
  # --- DynForest predictions ---
  predictions <- predict_DynForest_model(
    DynModel = model,
    time_chr = time_chr,
    event_chr = event_chr,
    data = external_data,
    landmark = 6,
    columns_metabolites_here = metabolites,
    fixed_vars = fixed_vars)

  # --- Initialize AUC vectors ---
  AUC_vec <- numeric(length(clinical_times))
  CI_lower_vec <- numeric(length(clinical_times))
  CI_upper_vec <- numeric(length(clinical_times))

  for (i in seq_along(clinical_times)) {
    t_interest <- clinical_times[i]
    idx_time <- which.min(abs(predictions$times - t_interest))
    cat(paste0("  → Time: ", t_interest, " months | Pred idx: ", idx_time, "\n"))

    pred_ids <- as.numeric(rownames(predictions$pred_indiv[, idx_time, drop = FALSE]))
    marker_vals <- predictions$pred_indiv[, idx_time]

    data_ext_distinct <- distinct(external_data, id_numeric, .keep_all = TRUE)
    data_ext_distinct <- as.data.frame(data_ext_distinct)

    true_times <- data_ext_distinct[which(data_ext_distinct$id_numeric %in% pred_ids), time_chr]
    true_events <- data_ext_distinct[which(data_ext_distinct$id_numeric %in% pred_ids), event_chr]

    roc_obj <- timeROC(
      T = true_times,
      delta = true_events,
      marker = marker_vals,
      cause = 1,
      times = t_interest,
      weighting = "marginal",
      ROC = TRUE,
      iid = TRUE
    )
    AUC_vec[i] <- roc_obj$AUC[2] * 100
    CI <- confint(roc_obj)$CI_AUC
    CI_lower_vec[i] <- CI[1]
    CI_upper_vec[i] <- CI[2]
  }
  auc_df <- data.frame(
    Time_months = clinical_times,
    AUC = AUC_vec,
    CI_lower = CI_lower_vec,
    CI_upper = CI_upper_vec)
  # --- Plot AUC curve ---
  if (!is.null(plot_path)) {
    p <- ggplot(auc_df, aes(x = Time_months, y = AUC)) +
      geom_line(color = "steelblue", size = 1.2) +
      geom_point(size = 3) +
      geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "steelblue", alpha = 0.2) +
      scale_y_continuous(limits = c(-30, 120), breaks = seq(-30, 120, by = 10)) +
      labs(
        title = paste0("External Validation - AUC over Time (", paste(study_name, collapse = "_"), ", iteration ", iteration_number, ")"),
        x = "Months",
        y = "AUC (95% CI)"
      ) +
      theme_minimal(base_size = 14)
        ggsave(plot_path, p, width = 12, height = 6)}
  return(list(predictions = predictions, AUC_table = auc_df))
}

eval_dynforest_cindex_on_external_cohort <- function(
    results_list,
    iteration_number,
    external_data,
    study_name,
    time_chr,
    event_chr,
    plot_path = NULL
) {
  library(nftbart)
  external_data <- as.data.frame(external_data)

  # --- Extract model, metabolites and fixed_vars from this iteration ---
  model <- results_list[[paste0("iteration_", iteration_number)]]$model
  metabolites <- results_list[[paste0("iteration_", iteration_number)]]$metabolites
  fixed_vars <- results_list[[paste0("iteration_", iteration_number)]]$fixed_vars
  print(fixed_vars)

  # --- DynForest predictions ---
  predictions <- predict_DynForest_model(
    DynModel = model,
    time_chr = time_chr,
    event_chr = event_chr,
    data = external_data,
    landmark = 6,
    columns_metabolites_here = metabolites,
    fixed_vars = fixed_vars
  )

  # Use predictions at landmark 6 (or adjust if needed)
  idx_time <- which.min(abs(predictions$times - 6))
  pred_ids <- as.numeric(rownames(predictions$pred_indiv[, idx_time, drop = FALSE]))
  risk_scores <- predictions$pred_indiv[, idx_time]

  data_ext_distinct <- distinct(external_data, id_numeric, .keep_all = TRUE) %>% as.data.frame()
  matched_idx <- which(data_ext_distinct$id_numeric %in% pred_ids)

  true_times <- data_ext_distinct[matched_idx, time_chr]
  true_events <- data_ext_distinct[matched_idx, event_chr]

  # Compute C-index
  c_index <- Cindex(risk_scores, true_times, true_events)

  cat(paste0("Computed C-index: ", round(c_index, 4), "\n"))
  return(list(
    c_index = c_index,
    predictions = predictions
  ))
}

eval_dynforest_on_internal_cohort_with_bootstrap <- function(
  results_list,
  iteration_number,
  external_data,
  study_name,
  time_chr,
  event_chr,
  clinical_times = c(6, 12), # We focus on a 6 and 12 months predictions
  n_bootstrap = 50, # 50 bootstrap resamplings
  seed = 123,
  plot_path = NULL
) {
  set.seed(seed + iteration_number)
  
  metabolites <- results_list[[paste0("iteration_", iteration_number)]]$metabolites
  fixed_vars <- results_list[[paste0("iteration_", iteration_number)]]$fixed_vars
  
  results_by_boot <- list()
  b <- 1  # compteur externe pour suivre combien de bootstraps valides on a
  
  while (b <= n_bootstrap) {
    cat(paste0("Bootstrap iteration ", b, "/", n_bootstrap, "...\n"))
    
    ids <- unique(external_data$id_numeric)
    sampled_ids <- sample(ids, replace = TRUE)
    oob_ids <- setdiff(ids, sampled_ids)
    if (length(oob_ids) == 0) {
      warning(paste0("Bootstrap ", b, " → no OOB samples, skipping."))
      next
    }
    
    data_boot <- external_data %>% filter(id_numeric %in% sampled_ids)
    data_oob <- external_data %>% filter(id_numeric %in% oob_ids)
    
    model <- train_DynForest_model(data_boot, time_chr, event_chr, metabolites, fixed_vars)
    predictions <- predict_DynForest_model(
      DynModel = model,
      time_chr = time_chr,
      event_chr = event_chr,
      data = data_oob,
      landmark = 6, #Landmark at 6 months, informations up until this date is used for predictions
      columns_metabolites_here = metabolites,
      fixed_vars = fixed_vars
    )
    
    # Verify if landmark is too far
    check_times <- sapply(clinical_times, function(t_interest) {
      idx_time <- which.min(abs(predictions$times - t_interest))
      deviation <- abs(predictions$times - t_interest)[idx_time]
      deviation >= 0.5
    })
    
    if (any(check_times)) {
      cat("⚠ Deviation ≥ 0.5 detected → skipping this bootstrap and retrying...\n") # Skip if landmark is to far from prediction - 1/2 month max is allowed
      next  # recommence la boucle sans incrémenter b
    }
    
    AUC_vec <- numeric(length(clinical_times))
    CI_lower_vec <- numeric(length(clinical_times))
    CI_upper_vec <- numeric(length(clinical_times))
    for (i in seq_along(clinical_times)) {
      t_interest <- clinical_times[i]
      idx_time <- which.min(abs(predictions$times - t_interest))
      pred_ids <- as.numeric(rownames(predictions$pred_indiv[, idx_time, drop = FALSE]))
      marker_vals <- predictions$pred_indiv[, idx_time]
      data_oob_distinct <- distinct(data_oob, id_numeric, .keep_all = TRUE) %>% as.data.frame()
      true_times <- data_oob_distinct[which(data_oob_distinct$id_numeric %in% pred_ids), time_chr]
      true_events <- data_oob_distinct[which(data_oob_distinct$id_numeric %in% pred_ids), event_chr]
      tryCatch({
        roc_obj <- timeROC( T = true_times, delta = true_events,  marker = marker_vals, cause = 1, times = t_interest, weighting = "marginal",ROC = TRUE,  iid = TRUE )
        auc_val <- roc_obj$AUC[2]
        ci_vals <- confint(roc_obj)$CI_AUC
        if (!is.na(auc_val) && !any(is.na(ci_vals))) {
          AUC_vec[i] <- auc_val * 100
          CI_lower_vec[i] <- ci_vals[1]
          CI_upper_vec[i] <- ci_vals[2]
        } else {
          warning(paste0("NA AUC or CI at t=", t_interest, " (bootstrap ", b, ")"))
          AUC_vec[i] <- NA
          CI_lower_vec[i] <- NA
          CI_upper_vec[i] <- NA
        }
      }, error = function(e) {
        warning(paste0("Error in timeROC at t=", t_interest, " (bootstrap ", b, "): ", e$message))
        AUC_vec[i] <- NA
        CI_lower_vec[i] <- NA
        CI_upper_vec[i] <- NA
      }) }
    auc_df <- data.frame(  Bootstrap = b,  Time_months = clinical_times,AUC = AUC_vec,CI_lower = CI_lower_vec,CI_upper = CI_upper_vec )
    results_by_boot[[paste0("boot", b)]] <- auc_df
    b <- b + 1  # increment only after a valid boostrap
  }
  return(results_by_boot)
}

# Plot AUC over time for a given iteration of the model 
plot_auc_evolution <- function(auc_df, save_path = NULL) {
  p <- ggplot(auc_df, aes(x = iteration, y = AUC, group = Time, color = Time)) +
    geom_line() +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    labs(title = "Evolution of AUC over ablation iterations",
         y = "AUC (95% CI)", x = "Ablation iteration") +
    theme_minimal(base_size = 14)
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 6)
  }
  return(p)
}