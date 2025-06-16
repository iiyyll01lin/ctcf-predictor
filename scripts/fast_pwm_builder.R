#!/usr/bin/env Rscript

# Fast PWM Builder
# Optimized PWM building for speed with large datasets
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
})

#' Fast PWM Builder - optimized for speed
#' @param sequences_file Path to input sequences file
#' @param output_file Output file for PWM results
#' @param max_sequences Maximum sequences to use (for speed)
#' @param pseudocount Pseudocount for matrix smoothing
#' @param subsample_method Method for subsampling (random, quality, length)
#' @param quality_threshold Quality threshold for filtering
#' @param verbose Enable verbose output
fast_pwm_builder <- function(sequences_file, output_file = NULL,
                            max_sequences = 10000, pseudocount = 0.1,
                            subsample_method = "random", quality_threshold = 0.8,
                            verbose = FALSE) {
  
  if (verbose) cat("Fast PWM Builder starting...\n")
  
  start_time <- Sys.time()
  
  # Load sequences efficiently
  if (verbose) cat("Loading sequences...\n")
  sequences <- load_sequences_fast(sequences_file, max_sequences, subsample_method, 
                                  quality_threshold, verbose)
  
  # Quick quality check
  if (length(sequences) < 100) {
    warning("Very few sequences (", length(sequences), "). Results may be unreliable.")
  }
  
  # Fast PWM building
  if (verbose) cat("Building PWM matrix...\n")
  pwm_result <- build_pwm_fast(sequences, pseudocount, verbose)
  
  # Quick quality assessment
  if (verbose) cat("Assessing PWM quality...\n")
  quality_assessment <- assess_pwm_quality_fast(pwm_result, verbose)
  
  # Combine results
  final_result <- list(
    pwm = pwm_result$pwm,
    freq_matrix = pwm_result$freq_matrix,
    info_content = pwm_result$info_content,
    total_info = pwm_result$total_info,
    consensus = pwm_result$consensus,
    n_sequences = length(sequences),
    max_sequences = max_sequences,
    pseudocount = pseudocount,
    subsample_method = subsample_method,
    quality_assessment = quality_assessment,
    method = "fast_pwm_builder",
    timestamp = Sys.time()
  )
  
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  final_result$processing_time <- processing_time
  
  # Save results
  if (!is.null(output_file)) {
    save_pwm_results_fast(final_result, output_file, verbose)
  }
  
  # Display summary
  display_fast_pwm_summary(final_result, verbose)
  
  if (verbose) cat("Fast PWM building completed in", round(processing_time, 2), "seconds\n")
  
  return(final_result)
}

#' Load sequences with fast methods
load_sequences_fast <- function(sequences_file, max_sequences, subsample_method, 
                               quality_threshold, verbose) {
  
  if (!file.exists(sequences_file)) {
    stop("Sequences file not found: ", sequences_file)
  }
  
  # Load sequences
  sequences <- readDNAStringSet(sequences_file)
  
  if (verbose) cat("  Loaded", length(sequences), "sequences\n")
  
  # Quick subsample if needed
  if (length(sequences) > max_sequences) {
    if (verbose) cat("  Subsampling to", max_sequences, "sequences using", subsample_method, "method\n")
    
    if (subsample_method == "random") {
      # Random sampling
      sample_indices <- sample(length(sequences), max_sequences)
      sequences <- sequences[sample_indices]
      
    } else if (subsample_method == "quality") {
      # Quality-based sampling
      sequences <- subsample_by_quality_fast(sequences, max_sequences, quality_threshold, verbose)
      
    } else if (subsample_method == "length") {
      # Length-based sampling (prefer median length)
      sequences <- subsample_by_length_fast(sequences, max_sequences, verbose)
      
    } else {
      # Default to random
      sample_indices <- sample(length(sequences), max_sequences)
      sequences <- sequences[sample_indices]
    }
  }
  
  # Quick length filtering - keep most common length
  seq_lengths <- width(sequences)
  common_length <- as.numeric(names(sort(table(seq_lengths), decreasing = TRUE))[1])
  sequences <- sequences[seq_lengths == common_length]
  
  if (verbose) cat("  Using", length(sequences), "sequences of length", common_length, "\n")
  
  return(sequences)
}

#' Fast quality-based subsampling
subsample_by_quality_fast <- function(sequences, max_sequences, quality_threshold, verbose) {
  
  # Calculate simple quality metrics
  seq_lengths <- width(sequences)
  n_content <- letterFrequency(sequences, "N") / seq_lengths
  
  # Quality score: low N content
  quality_scores <- 1 - n_content
  
  # Select high-quality sequences first
  high_quality <- which(quality_scores >= quality_threshold)
  
  if (length(high_quality) >= max_sequences) {
    # Enough high-quality sequences
    selected <- sample(high_quality, max_sequences)
  } else {
    # Include all high-quality, then fill with best remaining
    remaining <- setdiff(1:length(sequences), high_quality)
    remaining_scores <- quality_scores[remaining]
    remaining_sorted <- remaining[order(remaining_scores, decreasing = TRUE)]
    
    n_additional <- max_sequences - length(high_quality)
    additional <- head(remaining_sorted, n_additional)
    
    selected <- c(high_quality, additional)
  }
  
  return(sequences[selected])
}

#' Fast length-based subsampling
subsample_by_length_fast <- function(sequences, max_sequences, verbose) {
  
  seq_lengths <- width(sequences)
  median_length <- median(seq_lengths)
  
  # Calculate distance from median length
  length_distances <- abs(seq_lengths - median_length)
  
  # Select sequences closest to median length
  selected_indices <- order(length_distances)[1:min(max_sequences, length(sequences))]
  
  return(sequences[selected_indices])
}

#' Fast PWM building
build_pwm_fast <- function(sequences, pseudocount, verbose) {
  
  if (length(sequences) == 0) {
    stop("No sequences provided for PWM building")
  }
  
  # Ensure all sequences same length
  seq_lengths <- width(sequences)
  if (length(unique(seq_lengths)) > 1) {
    stop("All sequences must have the same length for fast PWM building")
  }
  
  seq_length <- seq_lengths[1]
  n_sequences <- length(sequences)
  
  # Initialize frequency matrix
  bases <- c("A", "C", "G", "T")
  freq_matrix <- matrix(0, nrow = 4, ncol = seq_length, 
                       dimnames = list(bases, paste0("pos", 1:seq_length)))
  
  # Fast counting using vectorized operations
  if (verbose) cat("  Counting nucleotides at each position...\n")
  
  # Convert sequences to matrix for fast processing
  seq_matrix <- as.matrix(sequences)
  
  # Count frequencies for each position
  for (pos in 1:seq_length) {
    pos_bases <- seq_matrix[, pos]
    
    # Count each base
    freq_matrix["A", pos] <- sum(pos_bases == "A")
    freq_matrix["C", pos] <- sum(pos_bases == "C")
    freq_matrix["G", pos] <- sum(pos_bases == "G")
    freq_matrix["T", pos] <- sum(pos_bases == "T")
  }
  
  # Add pseudocounts
  freq_matrix <- freq_matrix + pseudocount
  
  # Convert to probabilities
  prob_matrix <- sweep(freq_matrix, 2, colSums(freq_matrix), FUN = "/")
  
  # Calculate information content quickly
  if (verbose) cat("  Calculating information content...\n")
  
  background_prob <- 0.25
  ic_per_pos <- apply(prob_matrix, 2, function(col) {
    # Avoid log(0) issues
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col / background_prob))
  })
  
  total_ic <- sum(ic_per_pos)
  
  # Generate consensus sequence
  consensus <- paste(bases[apply(prob_matrix, 2, which.max)], collapse = "")
  
  # Return PWM result
  pwm_result <- list(
    pwm = prob_matrix,
    freq_matrix = freq_matrix,
    info_content = ic_per_pos,
    total_info = total_ic,
    consensus = consensus
  )
  
  return(pwm_result)
}

#' Fast PWM quality assessment
assess_pwm_quality_fast <- function(pwm_result, verbose) {
  
  total_ic <- pwm_result$total_info
  ic_per_pos <- pwm_result$info_content
  
  # Basic quality metrics
  mean_ic <- mean(ic_per_pos)
  max_ic <- max(ic_per_pos)
  n_high_ic <- sum(ic_per_pos > 1.5)
  n_conserved <- sum(ic_per_pos > 1.0)
  
  # Quality grade
  if (total_ic >= 20) {
    quality_grade <- "Excellent"
  } else if (total_ic >= 15) {
    quality_grade <- "Very Good"
  } else if (total_ic >= 10) {
    quality_grade <- "Good"
  } else if (total_ic >= 5) {
    quality_grade <- "Fair"
  } else {
    quality_grade <- "Poor"
  }
  
  # Pattern strength
  pattern_strength <- n_high_ic / length(ic_per_pos)
  
  # Conservation ratio
  conservation_ratio <- n_conserved / length(ic_per_pos)
  
  quality_assessment <- list(
    total_information_content = total_ic,
    mean_ic = mean_ic,
    max_ic = max_ic,
    n_high_ic_positions = n_high_ic,
    n_conserved_positions = n_conserved,
    quality_grade = quality_grade,
    pattern_strength = pattern_strength,
    conservation_ratio = conservation_ratio,
    pwm_length = length(ic_per_pos)
  )
  
  return(quality_assessment)
}

#' Save PWM results efficiently
save_pwm_results_fast <- function(result, output_file, verbose) {
  
  # Default to RDS for speed
  if (!grepl("\\.[^.]*$", output_file)) {
    output_file <- paste0(output_file, ".rds")
  }
  
  if (grepl("\\.rds$", output_file)) {
    saveRDS(result, output_file)
  } else if (grepl("\\.json$", output_file)) {
    # JSON is slower but sometimes needed
    writeLines(jsonlite::toJSON(result, pretty = TRUE), output_file)
  } else {
    # Default to RDS
    saveRDS(result, paste0(output_file, ".rds"))
  }
  
  if (verbose) cat("  Results saved to:", output_file, "\n")
}

#' Display fast PWM summary
display_fast_pwm_summary <- function(result, verbose) {
  
  cat("\n=== Fast PWM Builder Summary ===\n")
  
  cat("Sequences used:", result$n_sequences, "\n")
  cat("Processing time:", round(result$processing_time, 2), "seconds\n")
  cat("PWM length:", length(result$info_content), "positions\n")
  
  cat("\nInformation Content:\n")
  cat("  Total:", round(result$total_info, 3), "bits\n")
  cat("  Mean per position:", round(mean(result$info_content), 3), "bits\n")
  cat("  Max per position:", round(max(result$info_content), 3), "bits\n")
  
  cat("\nQuality Assessment:\n")
  qa <- result$quality_assessment
  cat("  Grade:", qa$quality_grade, "\n")
  cat("  High IC positions:", qa$n_high_ic_positions, "\n")
  cat("  Conserved positions:", qa$n_conserved_positions, "\n")
  cat("  Pattern strength:", round(qa$pattern_strength, 3), "\n")
  cat("  Conservation ratio:", round(qa$conservation_ratio, 3), "\n")
  
  cat("\nConsensus sequence:", result$consensus, "\n")
  
  # Performance metrics
  sequences_per_second <- result$n_sequences / result$processing_time
  cat("Processing rate:", round(sequences_per_second, 1), "sequences/second\n")
  
  # Quality recommendation
  if (qa$quality_grade %in% c("Excellent", "Very Good")) {
    cat("\n✅ RECOMMENDED: High-quality PWM suitable for applications\n")
  } else if (qa$quality_grade == "Good") {
    cat("\n⚠️  ACCEPTABLE: Good quality PWM, suitable for most applications\n")
  } else {
    cat("\n❌ NOT RECOMMENDED: Poor quality PWM, consider parameter adjustment\n")
  }
}

# Performance optimization settings
options(stringsAsFactors = FALSE)

# Command-line interface
if (!interactive()) {
  # Simple argument parsing for speed
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript fast_pwm_builder.R <sequences_file> [output_file] [max_sequences] [pseudocount]\n")
    cat("Example: Rscript fast_pwm_builder.R data/sequences.fasta results/fast_pwm.rds 5000 0.1\n")
    quit(status = 1)
  }
  
  sequences_file <- args[1]
  output_file <- if (length(args) >= 2) args[2] else NULL
  max_sequences <- if (length(args) >= 3) as.numeric(args[3]) else 10000
  pseudocount <- if (length(args) >= 4) as.numeric(args[4]) else 0.1
  
  # Basic validation
  if (!file.exists(sequences_file)) {
    cat("Error: Sequences file not found:", sequences_file, "\n")
    quit(status = 1)
  }
  
  if (is.null(output_file)) {
    base_name <- gsub("\\.[^.]*$", "", basename(sequences_file))
    output_file <- paste0("results/", base_name, "_fast_pwm.rds")
  }
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Run fast PWM builder
  tryCatch({
    cat("Starting fast PWM building...\n")
    
    result <- fast_pwm_builder(
      sequences_file = sequences_file,
      output_file = output_file,
      max_sequences = max_sequences,
      pseudocount = pseudocount,
      subsample_method = "random",
      quality_threshold = 0.8,
      verbose = TRUE
    )
    
    cat("Fast PWM building completed successfully.\n")
    
    # Exit with appropriate status
    if (result$quality_assessment$quality_grade %in% c("Excellent", "Very Good", "Good")) {
      quit(status = 0)
    } else {
      quit(status = 1)
    }
    
  }, error = function(e) {
    cat("Error in fast PWM building:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
