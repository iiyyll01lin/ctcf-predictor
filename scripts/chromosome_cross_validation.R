# Chromosome-based cross-validation for PWM validation
# Author: PWM Improvement Team

library(Biostrings)

# Perform chromosome-based k-fold cross-validation
chromosome_cross_validation <- function(sequences, chromosomes, k = 5) {
  cat("Performing chromosome-based", k, "-fold cross-validation...\n")
  
  # Load PWM building functions
  if (file.exists("scripts/build_pwm.R")) {
    source("scripts/build_pwm.R")
  }
  
  unique_chrs <- unique(chromosomes)
  chr_folds <- split(unique_chrs, cut(seq_along(unique_chrs), k, labels = FALSE))
  
  cv_results <- lapply(seq_along(chr_folds), function(fold_idx) {
    test_chrs <- chr_folds[[fold_idx]]
    
    # Training set: sequences NOT on test chromosomes
    train_idx <- !chromosomes %in% test_chrs
    train_sequences <- sequences[train_idx]
    
    # Test set: sequences on test chromosomes  
    test_idx <- chromosomes %in% test_chrs
    test_sequences <- sequences[test_idx]
    
    cat("Fold", fold_idx, "- Training:", length(train_sequences), 
        "Testing:", length(test_sequences), "\n")
    
    # Build PWM on training data
    train_pwm <- build_simple_pwm(train_sequences)
    
    # Calculate training IC
    train_info <- apply(train_pwm, 2, function(x) {
      x[x == 0] <- 1e-10
      2 + sum(x * log2(x))
    })
    train_ic <- sum(train_info)
    
    # Build PWM on test data for comparison
    test_pwm <- build_simple_pwm(test_sequences)
    test_info <- apply(test_pwm, 2, function(x) {
      x[x == 0] <- 1e-10
      2 + sum(x * log2(x))
    })
    test_ic <- sum(test_info)
    
    return(list(
      fold = fold_idx,
      train_ic = train_ic,
      test_ic = test_ic,
      train_size = length(train_sequences),
      test_size = length(test_sequences),
      test_chromosomes = test_chrs
    ))
  })
  
  # Summarize results
  train_ics <- sapply(cv_results, function(x) x$train_ic)
  test_ics <- sapply(cv_results, function(x) x$test_ic)
  
  summary <- list(
    cv_results = cv_results,
    mean_train_ic = mean(train_ics),
    sd_train_ic = sd(train_ics),
    mean_test_ic = mean(test_ics),
    sd_test_ic = sd(test_ics),
    correlation = cor(train_ics, test_ics)
  )
  
  return(summary)
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

# Extract chromosome information from sequence names
extract_chromosome <- function(seq_name) {
  # Extract chromosome from sequence name (format: chr1:123-456)
  if (grepl("chr[0-9XY]+", seq_name)) {
    chr_match <- regexpr("chr[0-9XY]+", seq_name)
    return(substr(seq_name, chr_match, chr_match + attr(chr_match, "match.length") - 1))
  } else {
    return("unknown")
  }
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  sequences_file <- if (length(args) >= 1) args[1] else "data/training_sequences.fasta"
  output_file <- if (length(args) >= 2) args[2] else "results/cv_results.rds"
  k_folds <- if (length(args) >= 3) as.numeric(args[3]) else 5
  
  cat("Chromosome-based Cross-validation for PWM\n")
  cat("=========================================\n")
  cat("Input file:", sequences_file, "\n")
  cat("Output file:", output_file, "\n")
  cat("K-folds:", k_folds, "\n\n")
  
  # Validate inputs
  if (!file.exists(sequences_file)) {
    stop("Input file not found: ", sequences_file)
  }
  
  # Load data
  sequences <- readDNAStringSet(sequences_file)
  sequences <- as.character(sequences)
  
  # Extract chromosomes from sequence names
  chromosomes <- sapply(names(sequences), extract_chromosome)
  
  cat("Loaded", length(sequences), "sequences\n")
  cat("Unique chromosomes:", length(unique(chromosomes)), "\n")
  cat("Chromosome distribution:\n")
  print(table(chromosomes))
  
  # Perform cross-validation
  cv_results <- chromosome_cross_validation(sequences, chromosomes, k = k_folds)
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save results
  saveRDS(cv_results, output_file)
  
  cat("\nCross-validation Results:\n")
  cat("Mean training IC:", round(cv_results$mean_train_ic, 3), "bits\n")
  cat("Mean test IC:", round(cv_results$mean_test_ic, 3), "bits\n")
  cat("Correlation:", round(cv_results$correlation, 3), "\n")
  cat("Results saved to:", output_file, "\n")
}
