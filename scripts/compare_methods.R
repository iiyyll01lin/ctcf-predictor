#!/usr/bin/env Rscript

# Compare Alignment Methods and Select Best
# Analyzes results from parallel alignment runs and selects optimal method
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

#' Main function for comparing alignment methods
#' @param results_dir Directory containing alignment results
#' @param output_file Output file for comparison report
#' @param selection_criteria Criteria for method selection
#' @param config_file Configuration file
#' @param verbose Enable verbose output
compare_alignment_methods <- function(results_dir, output_file = NULL,
                                     selection_criteria = "quality",
                                     config_file = NULL, verbose = FALSE) {
  
  if (verbose) cat("Comparing alignment methods...\n")
  
  # Load configuration
  config <- load_comparison_config(config_file)
  
  # Find alignment method directories
  alignments_dir <- file.path(results_dir, "alignments")
  
  if (!dir.exists(alignments_dir)) {
    stop("Alignments directory not found: ", alignments_dir)
  }
  
  method_dirs <- list.dirs(alignments_dir, full.names = FALSE, recursive = FALSE)
  
  if (length(method_dirs) == 0) {
    stop("No alignment method directories found")
  }
  
  if (verbose) cat("Found methods:", paste(method_dirs, collapse = ", "), "\n")
  
  # Analyze each method
  method_results <- list()
  
  for (method in method_dirs) {
    if (verbose) cat("Analyzing method:", method, "\n")
    
    method_dir <- file.path(alignments_dir, method)
    analysis <- analyze_alignment_method(method, method_dir, config, verbose)
    
    if (!is.null(analysis)) {
      method_results[[method]] <- analysis
    }
  }
  
  if (length(method_results) == 0) {
    stop("No valid alignment results found")
  }
  
  # Compare methods and select best
  comparison <- compare_methods(method_results, selection_criteria, config, verbose)
  
  # Generate comprehensive report
  report <- generate_comparison_report(method_results, comparison, config)
  
  # Save results
  if (!is.null(output_file)) {
    save_comparison_results(report, output_file, verbose)
  }
  
  # Display summary
  display_comparison_summary(comparison, verbose)
  
  return(list(
    comparison = comparison,
    method_results = method_results,
    report = report
  ))
}

#' Load comparison configuration
load_comparison_config <- function(config_file) {
  default_config <- list(
    weights = list(
      information_content = 0.4,
      conservation = 0.3,
      alignment_quality = 0.2,
      processing_time = 0.1
    ),
    thresholds = list(
      min_information_content = 8.0,
      min_conservation_ratio = 0.6,
      min_sequences = 100,
      max_processing_time = 3600  # 1 hour
    ),
    selection_criteria = "weighted_score",
    quality_metrics = c("information_content", "conservation", "alignment_stability")
  )
  
  if (!is.null(config_file) && file.exists(config_file)) {
    if (grepl("\\.json$", config_file)) {
      user_config <- fromJSON(config_file)
    } else if (grepl("\\.yml$|\\.yaml$", config_file)) {
      user_config <- yaml::read_yaml(config_file)
    } else {
      warning("Unsupported config format. Using defaults.")
      return(default_config)
    }
    
    # Deep merge configs
    default_config <- merge_configs(default_config, user_config)
  }
  
  return(default_config)
}

#' Merge configuration objects recursively
merge_configs <- function(default, user) {
  for (name in names(user)) {
    if (is.list(user[[name]]) && is.list(default[[name]])) {
      default[[name]] <- merge_configs(default[[name]], user[[name]])
    } else {
      default[[name]] <- user[[name]]
    }
  }
  return(default)
}

#' Analyze a single alignment method
analyze_alignment_method <- function(method, method_dir, config, verbose) {
  
  # Check for required files
  aligned_file <- file.path(method_dir, "aligned_sequences.fasta")
  metadata_file <- file.path(method_dir, "metadata.json")
  
  if (!file.exists(aligned_file)) {
    if (verbose) cat("  No aligned sequences found for", method, "\n")
    return(NULL)
  }
  
  # Load metadata if available
  metadata <- list()
  if (file.exists(metadata_file)) {
    metadata <- fromJSON(metadata_file)
  }
  
  # Check method status
  if (!is.null(metadata$status) && metadata$status != "success") {
    if (verbose) cat("  Method", method, "failed:", metadata$status, "\n")
    return(NULL)
  }
  
  # Read aligned sequences
  tryCatch({
    sequences <- readDNAStringSet(aligned_file)
    
    if (length(sequences) < config$thresholds$min_sequences) {
      if (verbose) cat("  Too few sequences for", method, ":", length(sequences), "\n")
      return(NULL)
    }
    
    # Calculate quality metrics
    quality_metrics <- calculate_method_quality(sequences, config)
    
    # Processing time from metadata
    processing_time <- metadata$duration_seconds %||% NA
    
    # Combine all metrics
    analysis <- list(
      method = method,
      n_sequences = length(sequences),
      sequence_length = unique(width(sequences)),
      processing_time = processing_time,
      metadata = metadata,
      quality = quality_metrics
    )
    
    if (verbose) {
      cat("  Method", method, "analysis:\n")
      cat("    Sequences:", analysis$n_sequences, "\n")
      cat("    Information Content:", round(quality_metrics$information_content, 3), "\n")
      cat("    Conservation Ratio:", round(quality_metrics$conservation_ratio, 3), "\n")
      if (!is.na(processing_time)) {
        cat("    Processing Time:", processing_time, "seconds\n")
      }
    }
    
    return(analysis)
    
  }, error = function(e) {
    if (verbose) cat("  Error analyzing", method, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}

#' Calculate quality metrics for alignment method
calculate_method_quality <- function(sequences, config) {
  
  # Ensure sequences are same length
  seq_lengths <- width(sequences)
  if (length(unique(seq_lengths)) > 1) {
    return(list(
      information_content = 0,
      conservation_ratio = 0,
      alignment_quality = 0,
      error = "Variable sequence lengths"
    ))
  }
  
  seq_length <- seq_lengths[1]
  
  # Build frequency matrix
  bases <- c("A", "C", "G", "T")
  freq_matrix <- matrix(0, nrow = 4, ncol = seq_length, 
                       dimnames = list(bases, paste0("pos", 1:seq_length)))
  
  # Count nucleotides at each position
  for (pos in 1:seq_length) {
    pos_chars <- as.character(subseq(sequences, pos, pos))
    for (base in bases) {
      freq_matrix[base, pos] <- sum(pos_chars == base)
    }
  }
  
  # Convert to probabilities
  prob_matrix <- sweep(freq_matrix, 2, colSums(freq_matrix), FUN = "/")
  
  # Calculate information content
  background_prob <- 0.25
  ic_per_pos <- apply(prob_matrix, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col / background_prob))
  })
  
  total_ic <- sum(ic_per_pos)
  
  # Conservation metrics
  high_ic_positions <- sum(ic_per_pos > 1.5)
  conservation_ratio <- high_ic_positions / seq_length
  
  # Alignment stability (consistency across positions)
  max_prob_per_pos <- apply(prob_matrix, 2, max)
  alignment_stability <- mean(max_prob_per_pos)
  
  # Pattern recognition quality
  pattern_score <- calculate_pattern_score(prob_matrix, ic_per_pos)
  
  # Overall alignment quality
  alignment_quality <- (alignment_stability + pattern_score) / 2
  
  return(list(
    information_content = total_ic,
    conservation_ratio = conservation_ratio,
    alignment_stability = alignment_stability,
    alignment_quality = alignment_quality,
    pattern_score = pattern_score,
    position_ic = ic_per_pos,
    max_prob_per_pos = max_prob_per_pos,
    prob_matrix = prob_matrix
  ))
}

#' Calculate pattern recognition score
calculate_pattern_score <- function(prob_matrix, ic_per_pos) {
  
  # Look for CTCF-like patterns
  ctcf_patterns <- list(
    c("C", "C"),      # CC dinucleotide
    c("G", "C"),      # GC dinucleotide  
    c("G", "G"),      # GG dinucleotide
    c("C", "A", "G"), # CAG trinucleotide
    c("G", "G", "C")  # GGC trinucleotide
  )
  
  pattern_scores <- numeric()
  
  # Check each pattern
  for (pattern in ctcf_patterns) {
    pattern_len <- length(pattern)
    best_score <- 0
    
    # Slide pattern across positions
    for (start_pos in 1:(ncol(prob_matrix) - pattern_len + 1)) {
      positions <- start_pos:(start_pos + pattern_len - 1)
      
      # Calculate pattern match score
      match_score <- 1
      for (i in seq_along(pattern)) {
        pos <- positions[i]
        base <- pattern[i]
        prob <- prob_matrix[base, pos]
        match_score <- match_score * prob
      }
      
      # Weight by information content
      ic_weight <- mean(ic_per_pos[positions])
      weighted_score <- match_score * ic_weight
      
      best_score <- max(best_score, weighted_score)
    }
    
    pattern_scores <- c(pattern_scores, best_score)
  }
  
  # Return average pattern score
  return(mean(pattern_scores))
}

#' Compare methods and select best
compare_methods <- function(method_results, selection_criteria, config, verbose) {
  
  if (verbose) cat("Comparing methods using criteria:", selection_criteria, "\n")
  
  # Calculate scores for each method
  method_scores <- list()
  
  for (method_name in names(method_results)) {
    result <- method_results[[method_name]]
    
    # Calculate weighted score
    score <- calculate_weighted_score(result, config)
    
    method_scores[[method_name]] <- list(
      method = method_name,
      weighted_score = score$weighted_score,
      component_scores = score$components,
      meets_thresholds = check_thresholds(result, config),
      result = result
    )
  }
  
  # Filter methods that meet thresholds
  valid_methods <- method_scores[sapply(method_scores, function(x) x$meets_thresholds)]
  
  if (length(valid_methods) == 0) {
    warning("No methods meet quality thresholds. Selecting best available.")
    valid_methods <- method_scores
  }
  
  # Select best method
  best_method <- select_best_method(valid_methods, selection_criteria, verbose)
  
  # Rank all methods
  all_scores <- sapply(method_scores, function(x) x$weighted_score)
  ranking <- names(sort(all_scores, decreasing = TRUE))
  
  return(list(
    best_method = best_method,
    method_scores = method_scores,
    ranking = ranking,
    selection_criteria = selection_criteria,
    valid_methods = names(valid_methods)
  ))
}

#' Calculate weighted score for a method
calculate_weighted_score <- function(result, config) {
  
  weights <- config$weights
  quality <- result$quality
  
  # Normalize scores to 0-1 range
  ic_score <- min(1.0, quality$information_content / 20.0)  # Assume max IC of 20
  conservation_score <- quality$conservation_ratio
  alignment_score <- quality$alignment_quality
  
  # Time score (faster is better, normalized to 0-1)
  time_score <- if (!is.na(result$processing_time)) {
    max(0, 1 - (result$processing_time / config$thresholds$max_processing_time))
  } else {
    0.5  # Default if time unknown
  }
  
  # Calculate weighted score
  weighted_score <- (
    weights$information_content * ic_score +
    weights$conservation * conservation_score +
    weights$alignment_quality * alignment_score +
    weights$processing_time * time_score
  )
  
  return(list(
    weighted_score = weighted_score,
    components = list(
      information_content = ic_score,
      conservation = conservation_score,
      alignment_quality = alignment_score,
      processing_time = time_score
    )
  ))
}

#' Check if method meets quality thresholds
check_thresholds <- function(result, config) {
  
  thresholds <- config$thresholds
  quality <- result$quality
  
  checks <- list(
    min_information_content = quality$information_content >= thresholds$min_information_content,
    min_conservation_ratio = quality$conservation_ratio >= thresholds$min_conservation_ratio,
    min_sequences = result$n_sequences >= thresholds$min_sequences,
    max_processing_time = is.na(result$processing_time) || 
                         result$processing_time <= thresholds$max_processing_time
  )
  
  return(all(unlist(checks)))
}

#' Select best method based on criteria
select_best_method <- function(valid_methods, criteria, verbose) {
  
  if (length(valid_methods) == 0) {
    return(NULL)
  }
  
  if (length(valid_methods) == 1) {
    return(valid_methods[[1]])
  }
  
  # Select based on criteria
  if (criteria == "quality" || criteria == "weighted_score") {
    scores <- sapply(valid_methods, function(x) x$weighted_score)
    best_idx <- which.max(scores)
  } else if (criteria == "information_content") {
    scores <- sapply(valid_methods, function(x) x$result$quality$information_content)
    best_idx <- which.max(scores)
  } else if (criteria == "speed") {
    times <- sapply(valid_methods, function(x) x$result$processing_time)
    best_idx <- which.min(times[!is.na(times)])
  } else {
    # Default to weighted score
    scores <- sapply(valid_methods, function(x) x$weighted_score)
    best_idx <- which.max(scores)
  }
  
  best_method <- valid_methods[[best_idx]]
  
  if (verbose) {
    cat("Selected best method:", best_method$method, "\n")
    cat("  Weighted score:", round(best_method$weighted_score, 3), "\n")
  }
  
  return(best_method)
}

#' Generate comprehensive comparison report
generate_comparison_report <- function(method_results, comparison, config) {
  
  report <- list(
    summary = list(
      timestamp = Sys.time(),
      n_methods_tested = length(method_results),
      n_methods_valid = length(comparison$valid_methods),
      best_method = comparison$best_method$method,
      selection_criteria = comparison$selection_criteria
    ),
    method_comparison = list(),
    quality_analysis = list(),
    recommendations = list()
  )
  
  # Method comparison details
  for (method_name in names(comparison$method_scores)) {
    method_score <- comparison$method_scores[[method_name]]
    method_result <- method_score$result
    
    report$method_comparison[[method_name]] <- list(
      method = method_name,
      weighted_score = round(method_score$weighted_score, 4),
      meets_thresholds = method_score$meets_thresholds,
      n_sequences = method_result$n_sequences,
      information_content = round(method_result$quality$information_content, 3),
      conservation_ratio = round(method_result$quality$conservation_ratio, 3),
      alignment_quality = round(method_result$quality$alignment_quality, 3),
      processing_time = method_result$processing_time,
      rank = which(comparison$ranking == method_name)
    )
  }
  
  # Quality analysis
  ic_values <- sapply(method_results, function(x) x$quality$information_content)
  conservation_values <- sapply(method_results, function(x) x$quality$conservation_ratio)
  
  report$quality_analysis <- list(
    information_content = list(
      min = min(ic_values),
      max = max(ic_values),
      mean = mean(ic_values),
      sd = sd(ic_values)
    ),
    conservation_ratio = list(
      min = min(conservation_values),
      max = max(conservation_values),
      mean = mean(conservation_values),
      sd = sd(conservation_values)
    )
  )
  
  # Generate recommendations
  report$recommendations <- generate_recommendations(comparison, config)
  
  return(report)
}

#' Generate recommendations based on comparison results
generate_recommendations <- function(comparison, config) {
  
  recommendations <- list()
  
  # Best method recommendation
  if (!is.null(comparison$best_method)) {
    best_score <- comparison$best_method$weighted_score
    
    if (best_score > 0.8) {
      recommendations$best_method <- "Excellent: Use recommended method for production"
    } else if (best_score > 0.6) {
      recommendations$best_method <- "Good: Method acceptable but consider parameter tuning"
    } else {
      recommendations$best_method <- "Poor: Consider alternative approaches or data quality review"
    }
  }
  
  # Method-specific recommendations
  method_scores <- comparison$method_scores
  
  # Check for processing time issues
  slow_methods <- names(method_scores)[sapply(method_scores, function(x) {
    !is.na(x$result$processing_time) && x$result$processing_time > 1800  # 30 minutes
  })]
  
  if (length(slow_methods) > 0) {
    recommendations$performance <- paste(
      "Consider optimizing slow methods:", paste(slow_methods, collapse = ", ")
    )
  }
  
  # Check for quality consistency
  scores <- sapply(method_scores, function(x) x$weighted_score)
  if (sd(scores) > 0.2) {
    recommendations$consistency <- "High variability in method performance. Review input data quality."
  }
  
  # Threshold warnings
  valid_methods <- comparison$valid_methods
  if (length(valid_methods) < length(method_scores)) {
    failed_methods <- setdiff(names(method_scores), valid_methods)
    recommendations$thresholds <- paste(
      "Methods failed quality thresholds:", paste(failed_methods, collapse = ", ")
    )
  }
  
  return(recommendations)
}

#' Save comparison results
save_comparison_results <- function(report, output_file, verbose) {
  
  if (grepl("\\.json$", output_file)) {
    writeLines(toJSON(report, pretty = TRUE), output_file)
  } else if (grepl("\\.rds$", output_file)) {
    saveRDS(report, output_file)
  } else {
    # Default to JSON
    writeLines(toJSON(report, pretty = TRUE), paste0(output_file, ".json"))
  }
  
  if (verbose) cat("Comparison report saved to:", output_file, "\n")
}

#' Display comparison summary
display_comparison_summary <- function(comparison, verbose) {
  
  cat("\n=== Alignment Method Comparison Results ===\n")
  
  if (!is.null(comparison$best_method)) {
    cat("Best Method:", comparison$best_method$method, "\n")
    cat("Weighted Score:", round(comparison$best_method$weighted_score, 3), "\n")
    
    result <- comparison$best_method$result
    cat("Information Content:", round(result$quality$information_content, 3), "\n")
    cat("Conservation Ratio:", round(result$quality$conservation_ratio, 3), "\n")
    cat("Sequences:", result$n_sequences, "\n")
    
    if (!is.na(result$processing_time)) {
      cat("Processing Time:", result$processing_time, "seconds\n")
    }
  }
  
  cat("\nMethod Ranking:\n")
  for (i in seq_along(comparison$ranking)) {
    method <- comparison$ranking[i]
    score <- comparison$method_scores[[method]]$weighted_score
    valid <- method %in% comparison$valid_methods
    status <- if (valid) "✓" else "✗"
    cat(sprintf("%d. %s %s (%.3f)\n", i, status, method, score))
  }
  
  cat("\nValid Methods:", length(comparison$valid_methods), "of", length(comparison$method_scores), "\n")
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-d", "--dir"), type = "character", default = NULL,
                help = "Results directory from parallel alignment", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file for comparison report", metavar = "character"),
    make_option(c("-s", "--selection"), type = "character", default = "quality",
                help = "Selection criteria [default: %default]", metavar = "character"),
    make_option(c("-c", "--config"), type = "character", default = NULL,
                help = "Configuration file", metavar = "character"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Compare Alignment Methods and Select Best")
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$dir)) {
    print_help(opt_parser)
    stop("Results directory is required.", call. = FALSE)
  }
  
  if (!dir.exists(opt$dir)) {
    stop("Results directory does not exist: ", opt$dir, call. = FALSE)
  }
  
  # Run method comparison
  tryCatch({
    result <- compare_alignment_methods(
      results_dir = opt$dir,
      output_file = opt$output,
      selection_criteria = opt$selection,
      config_file = opt$config,
      verbose = opt$verbose
    )
    
    if (!is.null(result$comparison$best_method)) {
      cat("\n✓ Method comparison completed successfully\n")
      cat("Recommended method:", result$comparison$best_method$method, "\n")
      quit(status = 0)
    } else {
      cat("\n✗ No suitable method found\n")
      quit(status = 1)
    }
    
  }, error = function(e) {
    cat("Error in method comparison:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
