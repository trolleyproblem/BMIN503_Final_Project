# Functions to take in regression_data tibble, perform
# 10-fold CV, create ROC plot, and report difference in Brier scores

compare_SPINC_performance <- function(regression_data,title) {
  
  # Assume that regression_data contains columns
  # with PFT measurements, binary outcome called "outcome"
  # and an arbitrary number of mu variables
  
  # Ensure that outcome is one-hot encoded
  regression_data <- 
    regression_data %>%
    mutate(outcome = if_else(outcome == first(outcome), 1, 0)) %>%
    mutate(outcome = factor(outcome))
  
  # Define the columns you want to include if present
  potential_predictors <- c("fev1", "fvc", "f_ratio", "pef", "fef257")
  
  # Filter out only columns that are actually present in your tibble
  available_predictors <- intersect(names(regression_data), potential_predictors)
  
  # Create the formula by combining "outcome" with the available predictors
  formula_base <- as.formula(paste("outcome ~", paste(available_predictors, collapse = " + ")))
  
  #Perform 10-fold cross validation on the training data
  set.seed(1234)
  df_folds <- vfold_cv(regression_data, v = 10)
  
  lr_spec <- 
    logistic_reg() |> 
    set_engine("glm")
  
  #Create a workflow() for fitting the glm
  glm_wf <- workflow() |>
    add_model(lr_spec) |>
    add_formula(formula_base)
  
  glm_fit_cv <- 
    glm_wf |>
    fit_resamples(df_folds, control = control_resamples(save_pred = TRUE))
  
  #Collect predictions out of folds into one tibble
  glm_cv_preds <- collect_predictions(glm_fit_cv)
  
  print(paste0("Brier base model: ",
               brier_class(glm_cv_preds, outcome, .pred_0)))
  
  glm_base_auc <-
    scales::percent(roc_auc(glm_cv_preds, outcome, .pred_0)$.estimate,
                    accuracy = 0.001)
  
  r0 <- 
    roc_curve(glm_cv_preds, outcome, .pred_0) %>%
    mutate(group = "base")
  
  # Select columns starting with "mu_" and store their names
  mu_predictors <- names(regression_data) %>% 
    .[grepl("^mu_", .)]
  
  # Construct the formula using the available "mu_" predictors
  formula_mu <- as.formula(paste("outcome ~", 
                                 paste(c(available_predictors, mu_predictors), collapse = " + ")))
  
  #Create a workflow() for fitting the glm
  glm_wf_mu <- workflow() |>
    add_model(lr_spec) |>
    add_formula(formula_mu)
  
  glm_fit_cv_mu <- 
    glm_wf_mu |>
    fit_resamples(df_folds, control = control_resamples(save_pred = TRUE))
  
  #Collect predictions out of folds into one tibble
  glm_cv_preds_mu <- collect_predictions(glm_fit_cv_mu)
  
  print(paste0("Brier augmented model: ",
               brier_class(glm_cv_preds_mu, outcome, .pred_0)))
  
  glm_mu_auc <-
    scales::percent(roc_auc(glm_cv_preds_mu, outcome, .pred_0)$.estimate,
                    accuracy = 0.001)
  
  r1 <- 
    roc_curve(glm_cv_preds_mu, outcome, .pred_0) %>%
    mutate(group = "mu")
  
  
  # Generate combined plot
  
  auc_labels <- tibble(
    group = c("mu","base"),
    auc = c(glm_mu_auc, glm_base_auc),
    x = c(0.2, 0.2),
    y = c(0.9, 0.8)
  )
  
  auc_plot <- bind_rows(r0, r1) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(aes(color=group)) +
    geom_abline(lty = 3) +
    geom_text(data = auc_labels, aes(label = auc, x = x, y = y, color=group)) +
    coord_equal() +
    theme_bw() +
    ggtitle(title)
  
  # Bootstrap for differences in Brier
  
  # Prepare data by binding predictions based on row alignment
  boot_data <- 
    tibble(
      outcome = as.numeric(glm_cv_preds$outcome)-1,
      .pred_0_base = glm_cv_preds$.pred_0,
      .pred_0_mu = glm_cv_preds_mu$.pred_0
    ) %>%
    drop_na()
  
  # Function to calculate the Brier Score
  calculate_brier_score <- function(predictions, truth) {
    mean((predictions - truth)^2)
  }
  
  # Function to calculate the difference in Brier Scores
  brier_diff <- function(data, indices) {
    sampled_data <- data[indices, ]
    
    brier_base <- calculate_brier_score(sampled_data$.pred_0_base, as.numeric(sampled_data$outcome))
    brier_mu <- calculate_brier_score(sampled_data$.pred_0_mu, as.numeric(sampled_data$outcome))
    
    diff <- brier_base - brier_mu
    if (is.na(diff)) diff <- 0  # Handle NA gracefully
    
    return(diff)
    
  }
  
  # Bootstrap the Brier Score difference with 1000 replications
  set.seed(42)
  boot_result <- boot(boot_data, statistic = brier_diff, R = 1000)
  
  # Calculate confidence interval for the difference
  boot_ci <- boot.ci(boot_result, type = "basic")
  
  # Calculate p-value
  observed_diff <- brier_diff(boot_data, 1:nrow(boot_data))
  if (observed_diff > 0) {
    # For a positive observed difference, calculate proportion of bootstraps with diff >= observed
    p_value <- mean(boot_result$t >= observed_diff) * 2
  } else {
    # For a negative observed difference, calculate proportion of bootstraps with diff <= observed
    p_value <- mean(boot_result$t <= observed_diff) * 2
  }
  
  # Print results
  print(boot_result)
  print(boot_ci)
  cat("p-value:", p_value, "\n")
  
  return(auc_plot)
  
}
