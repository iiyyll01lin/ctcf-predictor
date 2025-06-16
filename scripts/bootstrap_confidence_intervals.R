# Bootstrap confidence intervals for PWM validation
# Author: PWM Improvement Team

library(Biostrings)

# Calculate bootstrap confidence intervals for PWM information content
bootstrap_confidence_interval <- function(sequences, n_bootstrap = 1000, alpha = 0.05) {
  cat("Calculating bootstrap confidence intervals...\n")
  
  # Load PWM building function
  source("scripts/build_pwm.R")
  
  bootstrap_ics <- replicate(n_bootstrap, {
    # Sample with replacement
    boot_indices <- sample(length(sequences), replace = TRUE)
    boot_sequences <- sequences[boot_indices]
    
    # Build PWM and calculate IC
    boot_pwm <- build_simple_pwm(boot_sequences)
    
    # Calculate information content
    info_content <- apply(boot_pwm, 2, function(x) {
      x[x == 0] <- 1e-10
      2 + sum(x * log2(x))
    })
    
    return(sum(info_content))
  })
  
  # Calculate confidence interval
  lower <- quantile(bootstrap_ics, alpha/2)
  upper <- quantile(bootstrap_ics, 1 - alpha/2)
  
  return(list(
    mean = mean(bootstrap_ics),
    lower = lower,
    upper = upper,
    bootstrap_distribution = bootstrap_ics
  ))
}

# Build simple PWM function (if not available from build_pwm.R)
build_simple_pwm <- function(sequences, pseudocount = 0.01) {
  # Convert to matrix
  seq_matrix <- do.call(rbind, strsplit(sequences, ""))
  
  # Count bases at each position
  count_matrix <- matrix(0, nrow = 4, ncol = ncol(seq_matrix))
  rownames(count_matrix) <- c("A", "C", "G", "T")
  
  for (i in 1:ncol(seq_matrix)) {
    counts <- table(factor(seq_matrix[, i], levels = c("A", "C", "G", "T")))
    count_matrix[, i] <- as.numeric(counts)
  }
  
  # Add pseudocounts and normalize
  pwm <- (count_matrix + pseudocount) / (nrow(seq_matrix) + 4 * pseudocount)
  
  return(pwm)
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  sequences_file <- if (length(args) >= 1) args[1] else "data/training_sequences.fasta"
  output_file <- if (length(args) >= 2) args[2] else "results/bootstrap_ci.rds"
  n_bootstrap <- if (length(args) >= 3) as.numeric(args[3]) else 1000
  
  cat("Bootstrap Confidence Intervals for PWM\n")
  cat("=====================================\n")
  cat("Input file:", sequences_file, "\n")
  cat("Output file:", output_file, "\n")
  cat("Bootstrap replicates:", n_bootstrap, "\n\n")
  
  # Validate inputs
  if (!file.exists(sequences_file)) {
    stop("Input file not found: ", sequences_file)
  }
  
  # Load sequences
  sequences <- readDNAStringSet(sequences_file)
  sequences <- as.character(sequences)
  
  cat("Loaded", length(sequences), "sequences\n")
  
  # Calculate bootstrap CI
  ci_results <- bootstrap_confidence_interval(sequences, n_bootstrap = n_bootstrap)
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save results
  saveRDS(ci_results, output_file)
  
  cat("\nBootstrap Results:\n")
  cat("Mean IC:", round(ci_results$mean, 3), "bits\n")
  cat("95% CI:", round(ci_results$lower, 3), "-", round(ci_results$upper, 3), "bits\n")
  cat("Results saved to:", output_file, "\n")
}
