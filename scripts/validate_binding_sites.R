#!/usr/bin/env Rscript

# Validate Binding Sites
# Validates PWM predictions against experimental binding data
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

#' Main function for binding site validation
#' @param pwm_file Path to PWM results file
#' @param test_sequences_file Path to test sequences file
#' @param true_labels_file Path to true binding labels (optional)
#' @param output_file Output file for validation results
#' @param score_thresholds Score thresholds to test
#' @param validation_mode Validation mode (roc, precision_recall, confusion)
#' @param verbose Enable verbose output
validate_binding_sites <- function(pwm_file, test_sequences_file,
                                  true_labels_file = NULL,
                                  output_file = NULL,
                                  score_thresholds = seq(0.1, 0.9, 0.1),
                                  validation_mode = "roc",
                                  verbose = FALSE) {
  
  if (verbose) cat("Validating binding site predictions...\n")
  
  # Load PWM model
  if (verbose) cat("Loading PWM model...\n")
  pwm_result <- load_pwm_model(pwm_file, verbose)
  
  # Load test sequences
  if (verbose) cat("Loading test sequences...\n")
  test_data <- load_test_data(test_sequences_file, true_labels_file, verbose)
  
  # Score sequences with PWM
  if (verbose) cat("Scoring sequences with PWM...\n")
  sequence_scores <- score_test_sequences(test_data$sequences, pwm_result, verbose)
  
  # Perform validation analysis
  if (verbose) cat("Performing validation analysis...\n")
  validation_results <- perform_validation_analysis(sequence_scores, test_data$labels,
                                                   score_thresholds, validation_mode, verbose)
  
  # Calculate performance metrics
  if (verbose) cat("Calculating performance metrics...\n")
  performance_metrics <- calculate_performance_metrics(validation_results, verbose)
  
  # Create comprehensive results
  final_results <- list(
    pwm_file = pwm_file,
    test_sequences_file = test_sequences_file,
    true_labels_file = true_labels_file,
    n_test_sequences = length(test_data$sequences),
    n_positive_labels = sum(test_data$labels, na.rm = TRUE),
    validation_mode = validation_mode,
    score_thresholds = score_thresholds,
    sequence_scores = sequence_scores,
    validation_results = validation_results,
    performance_metrics = performance_metrics,
    timestamp = Sys.time()
  )
  
  # Save results
  if (!is.null(output_file)) {
    save_validation_results(final_results, output_file, verbose)
  }
  
  # Display summary
  display_validation_summary(performance_metrics, verbose)
  
  return(final_results)
}

#' Load PWM model from file
load_pwm_model <- function(pwm_file, verbose) {
  
  if (!file.exists(pwm_file)) {
    stop("PWM file not found: ", pwm_file)
  }
  
  pwm_result <- readRDS(pwm_file)
  
  # Extract PWM matrix
  if ("pwm" %in% names(pwm_result)) {
    pwm_matrix <- pwm_result$pwm
  } else if ("prob_matrix" %in% names(pwm_result)) {
    pwm_matrix <- pwm_result$prob_matrix
  } else {
    stop("No PWM matrix found in results file")
  }
  
  # Validate matrix format
  if (nrow(pwm_matrix) != 4) {
    if (ncol(pwm_matrix) == 4) {
      pwm_matrix <- t(pwm_matrix)
    } else {
      stop("PWM matrix must have 4 rows (bases)")
    }
  }
  
  if (is.null(rownames(pwm_matrix))) {
    rownames(pwm_matrix) <- c("A", "C", "G", "T")
  }
  
  pwm_result$pwm <- pwm_matrix
  
  if (verbose) cat("  PWM loaded:", nrow(pwm_matrix), "bases x", ncol(pwm_matrix), "positions\n")
  
  return(pwm_result)
}

#' Load test data (sequences and labels)
load_test_data <- function(test_sequences_file, true_labels_file, verbose) {
  
  # Load test sequences
  if (!file.exists(test_sequences_file)) {
    stop("Test sequences file not found: ", test_sequences_file)
  }
  
  sequences <- readDNAStringSet(test_sequences_file)
  
  if (length(sequences) == 0) {
    stop("No sequences found in test file")
  }
  
  # Load true labels if provided
  labels <- NULL
  if (!is.null(true_labels_file)) {
    if (file.exists(true_labels_file)) {
      
      # Try different formats
      if (grepl("\\.txt$|\\.csv$", true_labels_file)) {
        # Text/CSV format
        label_data <- read.table(true_labels_file, header = FALSE, stringsAsFactors = FALSE)
        labels <- as.logical(label_data[, 1])
      } else if (grepl("\\.rds$", true_labels_file)) {
        # RDS format
        labels <- readRDS(true_labels_file)
      } else {
        warning("Unknown label file format. Assuming text format.")
        label_data <- read.table(true_labels_file, header = FALSE, stringsAsFactors = FALSE)
        labels <- as.logical(label_data[, 1])
      }
      
      # Validate labels
      if (length(labels) != length(sequences)) {
        warning("Number of labels (", length(labels), ") does not match number of sequences (", 
                length(sequences), "). Truncating to minimum.")
        min_length <- min(length(labels), length(sequences))
        labels <- labels[1:min_length]
        sequences <- sequences[1:min_length]
      }
      
    } else {
      warning("True labels file not found: ", true_labels_file)
    }
  }
  
  # If no labels provided, generate synthetic labels based on sequence names
  if (is.null(labels)) {
    if (verbose) cat("  No true labels provided. Generating synthetic labels based on sequence names.\n")
    
    seq_names <- names(sequences)
    if (!is.null(seq_names)) {
      # Look for positive indicators in sequence names
      positive_patterns <- c("positive", "pos", "true", "binding", "bound", "ctcf")
      negative_patterns <- c("negative", "neg", "false", "nonbinding", "unbound", "control")
      
      labels <- rep(NA, length(sequences))
      
      for (i in seq_along(seq_names)) {
        name_lower <- tolower(seq_names[i])
        if (any(sapply(positive_patterns, function(p) grepl(p, name_lower)))) {
          labels[i] <- TRUE
        } else if (any(sapply(negative_patterns, function(p) grepl(p, name_lower)))) {
          labels[i] <- FALSE
        }
      }
      
      # Remove sequences without clear labels
      valid_labels <- !is.na(labels)
      if (sum(valid_labels) > 0) {
        sequences <- sequences[valid_labels]
        labels <- labels[valid_labels]
      } else {
        # Default: assume first half positive, second half negative
        n_seq <- length(sequences)
        labels <- c(rep(TRUE, n_seq %/% 2), rep(FALSE, n_seq - n_seq %/% 2))
        warning("No clear labels found. Using default split: first half positive, second half negative.")
      }
    } else {
      # No sequence names - default split
      n_seq <- length(sequences)
      labels <- c(rep(TRUE, n_seq %/% 2), rep(FALSE, n_seq - n_seq %/% 2))
      warning("No sequence names or labels. Using default split: first half positive, second half negative.")
    }
  }
  
  if (verbose) {
    cat("  Test sequences loaded:", length(sequences), "\n")
    if (!is.null(labels)) {
      cat("  Positive labels:", sum(labels, na.rm = TRUE), "\n")
      cat("  Negative labels:", sum(!labels, na.rm = TRUE), "\n")
    }
  }
  
  return(list(sequences = sequences, labels = labels))
}

#' Score test sequences with PWM
score_test_sequences <- function(sequences, pwm_result, verbose) {
  
  pwm_matrix <- pwm_result$pwm
  pwm_length <- ncol(pwm_matrix)
  
  sequence_scores <- data.frame(
    sequence_id = 1:length(sequences),
    sequence_name = names(sequences) %||% paste0("seq_", 1:length(sequences)),
    sequence_length = width(sequences),
    best_score = numeric(length(sequences)),
    best_position = integer(length(sequences)),
    best_strand = character(length(sequences)),
    stringsAsFactors = FALSE
  )
  
  # Score each sequence
  for (i in seq_along(sequences)) {
    
    sequence <- sequences[[i]]
    seq_length <- length(sequence)
    
    if (seq_length < pwm_length) {
      # Sequence too short
      sequence_scores$best_score[i] <- -Inf
      sequence_scores$best_position[i] <- NA
      sequence_scores$best_strand[i] <- "N/A"
      next
    }
    
    # Score forward strand
    forward_scores <- score_sequence_with_pwm(sequence, pwm_matrix)
    
    # Score reverse strand
    reverse_scores <- score_sequence_with_pwm(reverseComplement(sequence), pwm_matrix)
    
    # Find best score
    best_forward <- if (length(forward_scores) > 0) max(forward_scores, na.rm = TRUE) else -Inf
    best_reverse <- if (length(reverse_scores) > 0) max(reverse_scores, na.rm = TRUE) else -Inf
    
    if (best_forward >= best_reverse) {
      sequence_scores$best_score[i] <- best_forward
      sequence_scores$best_position[i] <- which.max(forward_scores)
      sequence_scores$best_strand[i] <- "+"
    } else {
      sequence_scores$best_score[i] <- best_reverse
      sequence_scores$best_position[i] <- which.max(reverse_scores)
      sequence_scores$best_strand[i] <- "-"
    }
  }
  
  if (verbose) cat("  Scored", nrow(sequence_scores), "sequences\n")
  
  return(sequence_scores)
}

#' Score single sequence with PWM
score_sequence_with_pwm <- function(sequence, pwm_matrix) {
  
  seq_length <- length(sequence)
  pwm_length <- ncol(pwm_matrix)
  
  if (seq_length < pwm_length) {
    return(numeric(0))
  }
  
  # Convert to character
  seq_chars <- as.character(sequence)
  
  scores <- numeric(seq_length - pwm_length + 1)
  
  for (pos in 1:(seq_length - pwm_length + 1)) {
    
    # Extract subsequence
    subseq <- substr(seq_chars, pos, pos + pwm_length - 1)
    subseq_chars <- strsplit(subseq, "")[[1]]
    
    # Calculate log-odds score
    score <- 0
    valid_position <- TRUE
    
    for (i in 1:pwm_length) {
      base <- subseq_chars[i]
      if (base %in% c("A", "C", "G", "T")) {
        # Log-odds score: log(PWM / background)
        pwm_prob <- pwm_matrix[base, i]
        background_prob <- 0.25
        score <- score + log2((pwm_prob + 1e-10) / background_prob)
      } else {
        # Handle ambiguous bases (N, etc.)
        valid_position <- FALSE
        break
      }
    }
    
    scores[pos] <- if (valid_position) score else -Inf
  }
  
  return(scores)
}

#' Perform validation analysis
perform_validation_analysis <- function(sequence_scores, true_labels, score_thresholds, 
                                       validation_mode, verbose) {
  
  if (is.null(true_labels)) {
    stop("True labels are required for validation analysis")
  }
  
  # Remove sequences with invalid scores
  valid_scores <- is.finite(sequence_scores$best_score)
  sequence_scores <- sequence_scores[valid_scores, ]
  true_labels <- true_labels[valid_scores]
  
  if (nrow(sequence_scores) == 0) {
    stop("No valid sequence scores for validation")
  }
  
  validation_results <- list()
  
  if (validation_mode %in% c("roc", "all")) {
    # ROC analysis
    validation_results$roc <- calculate_roc_curve(sequence_scores$best_score, true_labels, verbose)
  }
  
  if (validation_mode %in% c("precision_recall", "all")) {
    # Precision-Recall analysis
    validation_results$precision_recall <- calculate_precision_recall_curve(
      sequence_scores$best_score, true_labels, verbose)
  }
  
  if (validation_mode %in% c("confusion", "all")) {
    # Confusion matrix for each threshold
    validation_results$confusion_matrices <- calculate_confusion_matrices(
      sequence_scores$best_score, true_labels, score_thresholds, verbose)
  }
  
  # Threshold optimization
  validation_results$optimal_thresholds <- find_optimal_thresholds(
    sequence_scores$best_score, true_labels, verbose)
  
  return(validation_results)
}

#' Calculate ROC curve
calculate_roc_curve <- function(scores, labels, verbose) {
  
  # Sort by score (descending)
  order_idx <- order(scores, decreasing = TRUE)
  sorted_scores <- scores[order_idx]
  sorted_labels <- labels[order_idx]
  
  # Calculate cumulative true positives and false positives
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  
  tpr <- numeric(length(scores) + 1)  # True Positive Rate
  fpr <- numeric(length(scores) + 1)  # False Positive Rate
  
  tpr[1] <- 0  # Start at (0,0)
  fpr[1] <- 0
  
  tp <- 0
  fp <- 0
  
  for (i in seq_along(sorted_labels)) {
    if (sorted_labels[i]) {
      tp <- tp + 1
    } else {
      fp <- fp + 1
    }
    
    tpr[i + 1] <- tp / n_pos
    fpr[i + 1] <- fp / n_neg
  }
  
  # Calculate AUC using trapezoidal rule
  auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
  
  roc_result <- list(
    fpr = fpr,
    tpr = tpr,
    auc = auc,
    thresholds = c(Inf, sorted_scores)
  )
  
  if (verbose) cat("  ROC AUC:", round(auc, 4), "\n")
  
  return(roc_result)
}

#' Calculate Precision-Recall curve
calculate_precision_recall_curve <- function(scores, labels, verbose) {
  
  # Sort by score (descending)
  order_idx <- order(scores, decreasing = TRUE)
  sorted_scores <- scores[order_idx]
  sorted_labels <- labels[order_idx]
  
  # Calculate cumulative precision and recall
  n_pos <- sum(labels)
  
  precision <- numeric(length(scores) + 1)
  recall <- numeric(length(scores) + 1)
  
  precision[1] <- 1  # Start at (0,1) - perfect precision at recall=0
  recall[1] <- 0
  
  tp <- 0
  total_predictions <- 0
  
  for (i in seq_along(sorted_labels)) {
    total_predictions <- i
    
    if (sorted_labels[i]) {
      tp <- tp + 1
    }
    
    precision[i + 1] <- tp / total_predictions
    recall[i + 1] <- tp / n_pos
  }
  
  # Calculate AUC-PR using trapezoidal rule
  auc_pr <- sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)
  
  pr_result <- list(
    recall = recall,
    precision = precision,
    auc_pr = auc_pr,
    thresholds = c(Inf, sorted_scores)
  )
  
  if (verbose) cat("  Precision-Recall AUC:", round(auc_pr, 4), "\n")
  
  return(pr_result)
}

#' Calculate confusion matrices for different thresholds
calculate_confusion_matrices <- function(scores, labels, thresholds, verbose) {
  
  confusion_results <- list()
  
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    
    # Make predictions
    predictions <- scores >= threshold
    
    # Calculate confusion matrix
    tp <- sum(predictions & labels)
    fp <- sum(predictions & !labels)
    tn <- sum(!predictions & !labels)
    fn <- sum(!predictions & labels)
    
    # Calculate metrics
    sensitivity <- tp / (tp + fn)  # True Positive Rate
    specificity <- tn / (tn + fp)  # True Negative Rate
    precision <- if (tp + fp > 0) tp / (tp + fp) else 0
    f1_score <- if (precision + sensitivity > 0) 2 * (precision * sensitivity) / (precision + sensitivity) else 0
    accuracy <- (tp + tn) / (tp + fp + tn + fn)
    
    confusion_results[[i]] <- list(
      threshold = threshold,
      tp = tp,
      fp = fp,
      tn = tn,
      fn = fn,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      f1_score = f1_score,
      accuracy = accuracy
    )
  }
  
  return(confusion_results)
}

#' Find optimal thresholds using different criteria
find_optimal_thresholds <- function(scores, labels, verbose) {
  
  # Test a range of thresholds
  test_thresholds <- seq(min(scores, na.rm = TRUE), max(scores, na.rm = TRUE), length.out = 100)
  
  best_f1 <- 0
  best_f1_threshold <- NA
  best_youden <- 0
  best_youden_threshold <- NA
  
  for (threshold in test_thresholds) {
    predictions <- scores >= threshold
    
    tp <- sum(predictions & labels)
    fp <- sum(predictions & !labels)
    tn <- sum(!predictions & !labels)
    fn <- sum(!predictions & labels)
    
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    precision <- if (tp + fp > 0) tp / (tp + fp) else 0
    
    # F1 score optimization
    f1_score <- if (precision + sensitivity > 0) 2 * (precision * sensitivity) / (precision + sensitivity) else 0
    
    if (f1_score > best_f1) {
      best_f1 <- f1_score
      best_f1_threshold <- threshold
    }
    
    # Youden's J statistic optimization (sensitivity + specificity - 1)
    youden_j <- sensitivity + specificity - 1
    
    if (youden_j > best_youden) {
      best_youden <- youden_j
      best_youden_threshold <- threshold
    }
  }
  
  optimal_thresholds <- list(
    f1_optimal = list(threshold = best_f1_threshold, f1_score = best_f1),
    youden_optimal = list(threshold = best_youden_threshold, youden_j = best_youden)
  )
  
  if (verbose) {
    cat("  Optimal F1 threshold:", round(best_f1_threshold, 4), "(F1 =", round(best_f1, 4), ")\n")
    cat("  Optimal Youden threshold:", round(best_youden_threshold, 4), "(J =", round(best_youden, 4), ")\n")
  }
  
  return(optimal_thresholds)
}

#' Calculate comprehensive performance metrics
calculate_performance_metrics <- function(validation_results, verbose) {
  
  metrics <- list()
  
  # ROC metrics
  if ("roc" %in% names(validation_results)) {
    metrics$roc_auc <- validation_results$roc$auc
  }
  
  # Precision-Recall metrics
  if ("precision_recall" %in% names(validation_results)) {
    metrics$pr_auc <- validation_results$precision_recall$auc_pr
  }
  
  # Optimal threshold metrics
  if ("optimal_thresholds" %in% names(validation_results)) {
    opt_thresh <- validation_results$optimal_thresholds
    metrics$optimal_f1_threshold <- opt_thresh$f1_optimal$threshold
    metrics$optimal_f1_score <- opt_thresh$f1_optimal$f1_score
    metrics$optimal_youden_threshold <- opt_thresh$youden_optimal$threshold
    metrics$optimal_youden_j <- opt_thresh$youden_optimal$youden_j
  }
  
  # Performance grade
  if (!is.null(metrics$roc_auc)) {
    if (metrics$roc_auc >= 0.9) {
      metrics$performance_grade <- "Excellent"
    } else if (metrics$roc_auc >= 0.8) {
      metrics$performance_grade <- "Very Good"
    } else if (metrics$roc_auc >= 0.7) {
      metrics$performance_grade <- "Good"
    } else if (metrics$roc_auc >= 0.6) {
      metrics$performance_grade <- "Fair"
    } else {
      metrics$performance_grade <- "Poor"
    }
  } else {
    metrics$performance_grade <- "Unknown"
  }
  
  return(metrics)
}

#' Save validation results
save_validation_results <- function(results, output_file, verbose) {
  
  if (grepl("\\.rds$", output_file)) {
    saveRDS(results, output_file)
  } else if (grepl("\\.json$", output_file)) {
    writeLines(toJSON(results, pretty = TRUE), output_file)
  } else {
    # Default to RDS
    saveRDS(results, paste0(output_file, ".rds"))
  }
  
  if (verbose) cat("Validation results saved to:", output_file, "\n")
}

#' Display validation summary
display_validation_summary <- function(performance_metrics, verbose) {
  
  cat("\n=== Binding Site Validation Summary ===\n")
  
  if (!is.null(performance_metrics$roc_auc)) {
    cat("ROC AUC:", round(performance_metrics$roc_auc, 4), "\n")
  }
  
  if (!is.null(performance_metrics$pr_auc)) {
    cat("Precision-Recall AUC:", round(performance_metrics$pr_auc, 4), "\n")
  }
  
  if (!is.null(performance_metrics$optimal_f1_score)) {
    cat("Optimal F1 Score:", round(performance_metrics$optimal_f1_score, 4), 
        "at threshold", round(performance_metrics$optimal_f1_threshold, 4), "\n")
  }
  
  if (!is.null(performance_metrics$optimal_youden_j)) {
    cat("Optimal Youden's J:", round(performance_metrics$optimal_youden_j, 4), 
        "at threshold", round(performance_metrics$optimal_youden_threshold, 4), "\n")
  }
  
  cat("Performance Grade:", performance_metrics$performance_grade %||% "Unknown", "\n")
  
  # Performance recommendation
  if (!is.null(performance_metrics$roc_auc)) {
    if (performance_metrics$roc_auc >= 0.8) {
      cat("\n✅ RECOMMENDED: High-performance PWM suitable for binding site prediction\n")
    } else if (performance_metrics$roc_auc >= 0.6) {
      cat("\n⚠️  ACCEPTABLE: Moderate performance, may need improvement\n")
    } else {
      cat("\n❌ NOT RECOMMENDED: Poor performance, significant improvement needed\n")
    }
  }
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-p", "--pwm"), type = "character", default = NULL,
                help = "PWM results file", metavar = "character"),
    make_option(c("-t", "--test"), type = "character", default = NULL,
                help = "Test sequences file", metavar = "character"),
    make_option(c("-l", "--labels"), type = "character", default = NULL,
                help = "True binding labels file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file for validation results", metavar = "character"),
    make_option(c("-m", "--mode"), type = "character", default = "roc",
                help = "Validation mode (roc, precision_recall, confusion, all) [default: %default]", 
                metavar = "character"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Validate PWM Binding Site Predictions")
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$pwm) || is.null(opt$test)) {
    print_help(opt_parser)
    stop("PWM file and test sequences file are required.", call. = FALSE)
  }
  
  if (!file.exists(opt$pwm)) {
    stop("PWM file does not exist: ", opt$pwm, call. = FALSE)
  }
  
  if (!file.exists(opt$test)) {
    stop("Test sequences file does not exist: ", opt$test, call. = FALSE)
  }
  
  # Generate default output filename if not provided
  if (is.null(opt$output)) {
    base_name <- gsub("\\.[^.]*$", "", basename(opt$pwm))
    opt$output <- file.path("results", paste0(base_name, "_binding_validation.rds"))
  }
  
  # Validate binding sites
  tryCatch({
    results <- validate_binding_sites(
      pwm_file = opt$pwm,
      test_sequences_file = opt$test,
      true_labels_file = opt$labels,
      output_file = opt$output,
      validation_mode = opt$mode,
      verbose = opt$verbose
    )
    
    cat("Binding site validation completed successfully.\n")
    
    # Return appropriate exit code based on performance
    if (!is.null(results$performance_metrics$roc_auc)) {
      if (results$performance_metrics$roc_auc >= 0.7) {
        quit(status = 0)
      } else {
        quit(status = 1)
      }
    } else {
      quit(status = 0)
    }
    
  }, error = function(e) {
    cat("Error in binding site validation:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
