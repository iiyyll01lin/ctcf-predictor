# Bayesian Information Criterion for PWM model comparison
# Author: PWM Improvement Team

library(Biostrings)

# Calculate BIC for PWM model
calculate_bic <- function(pwm, sequences) {
  # Calculate log-likelihood of sequences under PWM model
  log_likelihood <- calculate_log_likelihood(pwm, sequences)
  
  # Number of parameters (3 free parameters per position due to sum constraint)
  n_params <- 3 * ncol(pwm)
  
  # Number of observations
  n_obs <- length(sequences)
  
  # BIC = -2 × log(L) + k × log(n)
  bic <- -2 * log_likelihood + n_params * log(n_obs)
  
  return(list(
    log_likelihood = log_likelihood,
    n_params = n_params,
    n_obs = n_obs,
    bic = bic
  ))
}

# Calculate log-likelihood of sequences under PWM
calculate_log_likelihood <- function(pwm, sequences) {
  total_ll <- 0
  
  for (seq in sequences) {
    bases <- strsplit(seq, "")[[1]]
    seq_ll <- 0
    
    for (pos in seq_along(bases)) {
      if (pos <= ncol(pwm)) {
        base <- bases[pos]
        if (base %in% rownames(pwm)) {
          prob <- pwm[base, pos]
          if (prob > 0) {
            seq_ll <- seq_ll + log(prob)
          } else {
            seq_ll <- seq_ll + log(1e-10)  # Small pseudocount
          }
        } else {
          # Handle unknown bases (like N)
          seq_ll <- seq_ll + log(0.25)  # Uniform probability
        }
      }
    }
    
    total_ll <- total_ll + seq_ll
  }
  
  return(total_ll)
}

# Compare multiple PWM models using BIC
compare_pwm_models_bic <- function(pwm_files, sequences) {
  cat("Comparing PWM models using Bayesian Information Criterion...\n")
  
  bic_results <- list()
  
  for (file in pwm_files) {
    pwm_name <- basename(file)
    
    # Load PWM
    pwm_data <- readRDS(file)
    if (is.list(pwm_data) && "pwm" %in% names(pwm_data)) {
      pwm <- pwm_data$pwm
    } else if (is.matrix(pwm_data)) {
      pwm <- pwm_data
    } else {
      cat("Warning: Cannot extract PWM from", file, "\n")
      next
    }
    
    # Calculate BIC
    bic_result <- calculate_bic(pwm, sequences)
    bic_results[[pwm_name]] <- bic_result
    
    cat("Model:", pwm_name, "BIC:", round(bic_result$bic, 2), 
        "Log-likelihood:", round(bic_result$log_likelihood, 2), "\n")
  }
  
  # Rank models by BIC (lower is better)
  bic_values <- sapply(bic_results, function(x) x$bic)
  ranked_models <- names(sort(bic_values))
  
  cat("\nModel ranking (best to worst):\n")
  for (i in seq_along(ranked_models)) {
    model <- ranked_models[i]
    cat(i, ".", model, "- BIC:", round(bic_values[model], 2), "\n")
  }
  
  return(list(
    bic_results = bic_results,
    ranked_models = ranked_models,
    bic_values = bic_values
  ))
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  results_dir <- if (length(args) >= 1) args[1] else "results"
  sequences_file <- if (length(args) >= 2) args[2] else "data/training_sequences.fasta"
  output_file <- if (length(args) >= 3) args[3] else "results/bic_comparison.rds"
  
  cat("BIC Model Comparison for PWMs\n")
  cat("=============================\n")
  cat("Results directory:", results_dir, "\n")
  cat("Sequences file:", sequences_file, "\n")
  cat("Output file:", output_file, "\n\n")
  
  # Validate inputs
  if (!dir.exists(results_dir)) {
    stop("Results directory not found: ", results_dir)
  }
  if (!file.exists(sequences_file)) {
    stop("Sequences file not found: ", sequences_file)
  }
  
  # Find PWM files
  pwm_files <- list.files(results_dir, pattern = ".*pwm.*\\.rds$", full.names = TRUE)
  pwm_files <- pwm_files[!grepl("null", pwm_files)]
  
  if (length(pwm_files) == 0) {
    stop("No PWM files found in ", results_dir)
  }
  
  cat("Found", length(pwm_files), "PWM files:\n")
  for (file in pwm_files) {
    cat(" -", basename(file), "\n")
  }
  
  # Load sequences
  sequences <- readDNAStringSet(sequences_file)
  sequences <- as.character(sequences)
  
  cat("Loaded", length(sequences), "sequences\n\n")
  
  # Compare models
  comparison_results <- compare_pwm_models_bic(pwm_files, sequences)
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save results
  saveRDS(comparison_results, output_file)
  
  cat("\nBIC comparison results saved to:", output_file, "\n")
}
