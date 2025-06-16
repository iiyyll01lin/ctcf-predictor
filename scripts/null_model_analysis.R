#!/usr/bin/env Rscript

# CTCF PWM Testing Pipeline - Null Model Analysis
# Background comparison and null model testing for PWM validation
# Author: CTCF Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(Biostrings)
  library(parallel)
})

# Logging function
log_message <- function(level, message, script = "null_model_analysis.R") {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(sprintf("%s %s [%s] %s\n", timestamp, level, script, message))
}

# Main null model analysis function
analyze_null_models <- function(pwm_file, sequences_file, 
                               output_dir = "results/null_analysis",
                               n_null_models = 100, n_cores = 4) {
  
  log_message("INFO", "Starting null model analysis")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    log_message("INFO", paste("Created output directory:", output_dir))
  }
  
  # Load PWM and sequences
  log_message("INFO", paste("Loading PWM from:", pwm_file))
  pwm_data <- load_pwm_data(pwm_file)
  
  log_message("INFO", paste("Loading sequences from:", sequences_file))
  sequences <- readDNAStringSet(sequences_file)
  log_message("INFO", paste("Loaded", length(sequences), "sequences"))
  
  # Perform null model analyses
  null_results <- list()
  
  # 1. Random Background Model
  log_message("INFO", "Analyzing random background model")
  null_results$random_background <- analyze_random_background_model(pwm_data, sequences, n_null_models)
  
  # 2. Shuffled Sequence Model
  log_message("INFO", "Analyzing shuffled sequence model")
  null_results$shuffled_sequences <- analyze_shuffled_sequence_model(pwm_data, sequences, n_null_models)
  
  # 3. Dinucleotide Shuffled Model
  log_message("INFO", "Analyzing dinucleotide shuffled model")
  null_results$dinucleotide_shuffled <- analyze_dinucleotide_model(pwm_data, sequences, n_null_models)
  
  # 4. Uniform PWM Model
  log_message("INFO", "Analyzing uniform PWM model")
  null_results$uniform_pwm <- analyze_uniform_pwm_model(pwm_data, sequences, n_null_models)
  
  # 5. Position-Shuffled PWM Model
  log_message("INFO", "Analyzing position-shuffled PWM model")
  null_results$position_shuffled <- analyze_position_shuffled_model(pwm_data, sequences, n_null_models)
  
  # Compare with observed PWM
  log_message("INFO", "Comparing observed PWM with null models")
  comparison_results <- compare_with_null_models(pwm_data, sequences, null_results)
  
  # Generate comprehensive report
  generate_null_analysis_report(pwm_data, null_results, comparison_results, output_dir)
  
  # Save results
  results_file <- file.path(output_dir, "null_model_results.rds")
  saveRDS(list(
    pwm_data = pwm_data,
    null_results = null_results,
    comparison_results = comparison_results
  ), results_file)
  
  log_message("INFO", paste("Null model analysis completed. Results saved to:", results_file))
  return(comparison_results)
}

# Load PWM data from file
load_pwm_data <- function(pwm_file) {
  if (grepl("\\.rds$", pwm_file, ignore.case = TRUE)) {
    pwm_data <- readRDS(pwm_file)
    if (is.list(pwm_data)) {
      return(pwm_data)
    } else if (is.matrix(pwm_data)) {
      return(list(pwm = pwm_data))
    }
  } else {
    # Load as text matrix
    pwm_matrix <- as.matrix(read.table(pwm_file, header = TRUE, row.names = 1))
    return(list(pwm = pwm_matrix))
  }
}

# Analyze random background model
analyze_random_background_model <- function(pwm_data, sequences, n_models = 100) {
  log_message("DEBUG", "Generating random background models")
  
  pwm <- pwm_data$pwm
  motif_length <- ncol(pwm)
  
  # Calculate background nucleotide frequencies
  all_bases <- unlist(strsplit(as.character(sequences), ""))
  background_freq <- table(factor(all_bases, levels = c("A", "C", "G", "T")))
  background_freq <- background_freq / sum(background_freq)
  
  # Generate random PWMs based on background frequencies
  null_metrics <- list(
    information_content = numeric(n_models),
    max_position_ic = numeric(n_models),
    conservation_score = numeric(n_models)
  )
  
  for (i in 1:n_models) {
    # Generate random PWM with background frequencies
    random_pwm <- matrix(0, nrow = 4, ncol = motif_length)
    rownames(random_pwm) <- c("A", "C", "G", "T")
    
    for (pos in 1:motif_length) {
      # Add small random variation to background frequencies
      random_freqs <- background_freq + runif(4, -0.05, 0.05)
      random_freqs <- pmax(random_freqs, 0.01)  # Ensure positive
      random_freqs <- random_freqs / sum(random_freqs)  # Normalize
      random_pwm[, pos] <- random_freqs
    }
    
    # Calculate metrics
    null_metrics$information_content[i] <- calculate_total_ic(random_pwm)
    null_metrics$max_position_ic[i] <- max(calculate_position_ic(random_pwm))
    null_metrics$conservation_score[i] <- calculate_conservation_score(random_pwm)
  }
  
  return(null_metrics)
}

# Analyze shuffled sequence model
analyze_shuffled_sequence_model <- function(pwm_data, sequences, n_models = 100) {
  log_message("DEBUG", "Generating shuffled sequence models")
  
  pwm <- pwm_data$pwm
  motif_length <- ncol(pwm)
  
  null_metrics <- list(
    information_content = numeric(n_models),
    max_position_ic = numeric(n_models),
    conservation_score = numeric(n_models)
  )
  
  for (i in 1:n_models) {
    # Shuffle sequences
    shuffled_sequences <- shuffle_sequences_composition(sequences)
    
    # Build PWM from shuffled sequences
    shuffled_pwm <- build_pwm_from_sequences(shuffled_sequences, motif_length)
    
    if (!is.null(shuffled_pwm)) {
      null_metrics$information_content[i] <- calculate_total_ic(shuffled_pwm)
      null_metrics$max_position_ic[i] <- max(calculate_position_ic(shuffled_pwm))
      null_metrics$conservation_score[i] <- calculate_conservation_score(shuffled_pwm)
    }
  }
  
  return(null_metrics)
}

# Analyze dinucleotide shuffled model
analyze_dinucleotide_model <- function(pwm_data, sequences, n_models = 100) {
  log_message("DEBUG", "Generating dinucleotide shuffled models")
  
  pwm <- pwm_data$pwm
  motif_length <- ncol(pwm)
  
  null_metrics <- list(
    information_content = numeric(n_models),
    max_position_ic = numeric(n_models),
    conservation_score = numeric(n_models)
  )
  
  for (i in 1:n_models) {
    # Shuffle sequences preserving dinucleotide composition
    dinuc_shuffled <- shuffle_sequences_dinucleotide(sequences)
    
    # Build PWM from shuffled sequences
    shuffled_pwm <- build_pwm_from_sequences(dinuc_shuffled, motif_length)
    
    if (!is.null(shuffled_pwm)) {
      null_metrics$information_content[i] <- calculate_total_ic(shuffled_pwm)
      null_metrics$max_position_ic[i] <- max(calculate_position_ic(shuffled_pwm))
      null_metrics$conservation_score[i] <- calculate_conservation_score(shuffled_pwm)
    }
  }
  
  return(null_metrics)
}

# Analyze uniform PWM model
analyze_uniform_pwm_model <- function(pwm_data, sequences, n_models = 100) {
  log_message("DEBUG", "Generating uniform PWM models")
  
  pwm <- pwm_data$pwm
  motif_length <- ncol(pwm)
  
  null_metrics <- list(
    information_content = numeric(n_models),
    max_position_ic = numeric(n_models),
    conservation_score = numeric(n_models)
  )
  
  for (i in 1:n_models) {
    # Create uniform PWM with small random variations
    uniform_pwm <- matrix(0.25, nrow = 4, ncol = motif_length)
    rownames(uniform_pwm) <- c("A", "C", "G", "T")
    
    # Add small random noise
    for (pos in 1:motif_length) {
      noise <- runif(4, -0.02, 0.02)
      uniform_pwm[, pos] <- uniform_pwm[, pos] + noise
      uniform_pwm[, pos] <- pmax(uniform_pwm[, pos], 0.01)
      uniform_pwm[, pos] <- uniform_pwm[, pos] / sum(uniform_pwm[, pos])
    }
    
    null_metrics$information_content[i] <- calculate_total_ic(uniform_pwm)
    null_metrics$max_position_ic[i] <- max(calculate_position_ic(uniform_pwm))
    null_metrics$conservation_score[i] <- calculate_conservation_score(uniform_pwm)
  }
  
  return(null_metrics)
}

# Analyze position-shuffled PWM model
analyze_position_shuffled_model <- function(pwm_data, sequences, n_models = 100) {
  log_message("DEBUG", "Generating position-shuffled PWM models")
  
  pwm <- pwm_data$pwm
  motif_length <- ncol(pwm)
  
  null_metrics <- list(
    information_content = numeric(n_models),
    max_position_ic = numeric(n_models),
    conservation_score = numeric(n_models)
  )
  
  for (i in 1:n_models) {
    # Shuffle PWM positions
    shuffled_pwm <- pwm[, sample(ncol(pwm))]
    
    null_metrics$information_content[i] <- calculate_total_ic(shuffled_pwm)
    null_metrics$max_position_ic[i] <- max(calculate_position_ic(shuffled_pwm))
    null_metrics$conservation_score[i] <- calculate_conservation_score(shuffled_pwm)
  }
  
  return(null_metrics)
}

# Helper functions
calculate_total_ic <- function(pwm) {
  sum(apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col))
  }))
}

calculate_position_ic <- function(pwm) {
  apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col))
  })
}

calculate_conservation_score <- function(pwm) {
  mean(apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col))
  }))
}

shuffle_sequences_composition <- function(sequences) {
  DNAStringSet(lapply(sequences, function(seq) {
    seq_str <- as.character(seq)
    bases <- strsplit(seq_str, "")[[1]]
    shuffled_bases <- sample(bases)
    DNAString(paste(shuffled_bases, collapse = ""))
  }))
}

shuffle_sequences_dinucleotide <- function(sequences) {
  # Simplified dinucleotide shuffling
  DNAStringSet(lapply(sequences, function(seq) {
    seq_str <- as.character(seq)
    
    # Extract dinucleotides
    dinucs <- character(nchar(seq_str) - 1)
    for (i in 1:(nchar(seq_str) - 1)) {
      dinucs[i] <- substr(seq_str, i, i + 1)
    }
    
    # Shuffle dinucleotides and reconstruct
    shuffled_dinucs <- sample(dinucs)
    
    # Reconstruct sequence (simplified)
    if (length(shuffled_dinucs) > 0) {
      result <- substr(shuffled_dinucs[1], 1, 1)
      for (i in 1:length(shuffled_dinucs)) {
        result <- paste0(result, substr(shuffled_dinucs[i], 2, 2))
      }
      return(DNAString(result))
    } else {
      return(seq)
    }
  }))
}

build_pwm_from_sequences <- function(sequences, motif_length) {
  if (length(sequences) == 0) return(NULL)
  
  # Sample sequences and extract motifs
  sample_size <- min(1000, length(sequences))
  sample_seqs <- sample(sequences, sample_size)
  
  # Extract center regions
  motifs <- character(length(sample_seqs))
  for (i in 1:length(sample_seqs)) {
    seq_str <- as.character(sample_seqs[i])
    seq_len <- nchar(seq_str)
    
    if (seq_len >= motif_length) {
      center_start <- max(1, floor((seq_len - motif_length) / 2) + 1)
      center_end <- center_start + motif_length - 1
      motifs[i] <- substr(seq_str, center_start, center_end)
    } else {
      motifs[i] <- seq_str
    }
  }
  
  # Build PWM
  valid_motifs <- motifs[nchar(motifs) == motif_length]
  if (length(valid_motifs) == 0) return(NULL)
  
  pwm <- matrix(0, nrow = 4, ncol = motif_length)
  rownames(pwm) <- c("A", "C", "G", "T")
  
  for (pos in 1:motif_length) {
    base_counts <- table(factor(substr(valid_motifs, pos, pos), levels = c("A", "C", "G", "T")))
    base_counts <- base_counts + 0.1  # Pseudocount
    pwm[, pos] <- base_counts / sum(base_counts)
  }
  
  return(pwm)
}

# Compare observed PWM with null models
compare_with_null_models <- function(pwm_data, sequences, null_results) {
  log_message("DEBUG", "Comparing observed PWM with null models")
  
  pwm <- pwm_data$pwm
  
  # Calculate observed metrics
  observed_ic <- calculate_total_ic(pwm)
  observed_max_ic <- max(calculate_position_ic(pwm))
  observed_conservation <- calculate_conservation_score(pwm)
  
  # Compare with each null model
  comparisons <- list()
  
  for (null_model in names(null_results)) {
    null_data <- null_results[[null_model]]
    
    # Calculate p-values (empirical)
    ic_p_value <- sum(null_data$information_content >= observed_ic) / length(null_data$information_content)
    max_ic_p_value <- sum(null_data$max_position_ic >= observed_max_ic) / length(null_data$max_position_ic)
    conservation_p_value <- sum(null_data$conservation_score >= observed_conservation) / length(null_data$conservation_score)
    
    # Effect sizes
    ic_effect_size <- (observed_ic - mean(null_data$information_content)) / sd(null_data$information_content)
    max_ic_effect_size <- (observed_max_ic - mean(null_data$max_position_ic)) / sd(null_data$max_position_ic)
    conservation_effect_size <- (observed_conservation - mean(null_data$conservation_score)) / sd(null_data$conservation_score)
    
    comparisons[[null_model]] <- list(
      null_model = null_model,
      ic_p_value = ic_p_value,
      max_ic_p_value = max_ic_p_value,
      conservation_p_value = conservation_p_value,
      ic_effect_size = ic_effect_size,
      max_ic_effect_size = max_ic_effect_size,
      conservation_effect_size = conservation_effect_size,
      observed_ic = observed_ic,
      null_ic_mean = mean(null_data$information_content),
      null_ic_sd = sd(null_data$information_content)
    )
  }
  
  return(list(
    observed = list(
      information_content = observed_ic,
      max_position_ic = observed_max_ic,
      conservation_score = observed_conservation
    ),
    comparisons = comparisons
  ))
}

# Generate null analysis report
generate_null_analysis_report <- function(pwm_data, null_results, comparison_results, output_dir) {
  report_file <- file.path(output_dir, "null_model_analysis_report.txt")
  
  sink(report_file)
  cat("CTCF PWM Null Model Analysis Report\n")
  cat("==================================\n\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # PWM Information
  cat("PWM Information:\n")
  cat("  Dimensions:", nrow(pwm_data$pwm), "x", ncol(pwm_data$pwm), "\n")
  cat("  Observed Information Content:", round(comparison_results$observed$information_content, 3), "bits\n")
  cat("  Observed Max Position IC:", round(comparison_results$observed$max_position_ic, 3), "bits\n")
  cat("  Observed Conservation Score:", round(comparison_results$observed$conservation_score, 3), "\n\n")
  
  # Null Model Comparisons
  cat("Null Model Comparisons:\n")
  cat("======================\n\n")
  
  for (null_model in names(comparison_results$comparisons)) {
    comp <- comparison_results$comparisons[[null_model]]
    
    cat(toupper(gsub("_", " ", null_model)), "\n")
    cat(paste(rep("-", nchar(null_model) + 5), collapse = ""), "\n")
    
    cat("  Information Content:\n")
    cat("    Observed:", round(comp$observed_ic, 3), "bits\n")
    cat("    Null mean:", round(comp$null_ic_mean, 3), "± ", round(comp$null_ic_sd, 3), "bits\n")
    cat("    P-value:", formatC(comp$ic_p_value, format = "e", digits = 3), "\n")
    cat("    Effect size:", round(comp$ic_effect_size, 3), "\n")
    
    cat("  Max Position IC:\n")
    cat("    P-value:", formatC(comp$max_ic_p_value, format = "e", digits = 3), "\n")
    cat("    Effect size:", round(comp$max_ic_effect_size, 3), "\n")
    
    cat("  Conservation Score:\n")
    cat("    P-value:", formatC(comp$conservation_p_value, format = "e", digits = 3), "\n")
    cat("    Effect size:", round(comp$conservation_effect_size, 3), "\n")
    
    # Overall significance
    min_p_value <- min(comp$ic_p_value, comp$max_ic_p_value, comp$conservation_p_value)
    if (min_p_value < 0.001) {
      significance <- "HIGHLY SIGNIFICANT"
    } else if (min_p_value < 0.01) {
      significance <- "SIGNIFICANT"
    } else if (min_p_value < 0.05) {
      significance <- "MARGINALLY SIGNIFICANT"
    } else {
      significance <- "NOT SIGNIFICANT"
    }
    
    cat("  Overall Assessment:", significance, "\n\n")
  }
  
  # Summary
  cat("Summary:\n")
  cat("========\n")
  all_p_values <- unlist(lapply(comparison_results$comparisons, function(x) {
    c(x$ic_p_value, x$max_ic_p_value, x$conservation_p_value)
  }))
  
  highly_sig <- sum(all_p_values < 0.001)
  sig <- sum(all_p_values < 0.01)
  margin_sig <- sum(all_p_values < 0.05)
  
  cat("  Highly significant tests (p < 0.001):", highly_sig, "\n")
  cat("  Significant tests (p < 0.01):", sig, "\n")
  cat("  Marginally significant tests (p < 0.05):", margin_sig, "\n")
  cat("  Total tests:", length(all_p_values), "\n")
  
  if (highly_sig >= 6) {
    cat("  PWM Assessment: EXCELLENT - Highly significant across all null models\n")
  } else if (sig >= 9) {
    cat("  PWM Assessment: GOOD - Significant improvement over null models\n")
  } else if (margin_sig >= 6) {
    cat("  PWM Assessment: FAIR - Some improvement over null models\n")
  } else {
    cat("  PWM Assessment: POOR - Little improvement over null models\n")
  }
  
  sink()
  
  log_message("INFO", paste("Null model analysis report written to:", report_file))
}

# Command line interface
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript null_model_analysis.R <pwm_file> <sequences_file> [output_dir] [n_null_models]\n")
    cat("Example: Rscript null_model_analysis.R results/pwm_robust.rds data/filtered_sequences.fasta results/null_analysis 100\n")
    quit(status = 1)
  }
  
  pwm_file <- args[1]
  sequences_file <- args[2]
  output_dir <- ifelse(length(args) >= 3, args[3], "results/null_analysis")
  n_null_models <- ifelse(length(args) >= 4, as.numeric(args[4]), 100)
  
  tryCatch({
    result <- analyze_null_models(pwm_file, sequences_file, output_dir, n_null_models)
    cat("\nNull model analysis completed successfully!\n")
  }, error = function(e) {
    log_message("ERROR", paste("Null model analysis failed:", e$message))
    quit(status = 1)
  })
}

# Run if called directly
if (!interactive()) {
  main()
}
