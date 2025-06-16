#!/usr/bin/env Rscript

# Test Alignment Quality Script
# Tests and compares different alignment methods for PWM construction
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
  library(parallel)
})

# Command line options
option_list <- list(
  make_option(c("--sequences", "-s"), type="character", default=NULL,
              help="Input sequences file (FASTA format) [required]", metavar="character"),
  make_option(c("--methods", "-m"), type="character", default="center,consensus,integrated",
              help="Alignment methods to test (comma-separated) [default: %default]", metavar="character"),
  make_option(c("--output", "-o"), type="character", default="results/alignment_validation.html",
              help="Output validation report file [default: %default]", metavar="character"),
  make_option(c("--n-cores"), type="integer", default=1,
              help="Number of cores for parallel processing [default: %default]", metavar="integer"),
  make_option(c("--max-sequences"), type="integer", default=5000,
              help="Maximum sequences to use for testing [default: %default]", metavar="integer"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Test and compare alignment methods for PWM construction")
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$sequences)) {
  cat("Error: Input sequences file is required\n")
  print_help(opt_parser)
  quit(status=1)
}

if (!file.exists(opt$sequences)) {
  cat("Error: Input sequences file does not exist:", opt$sequences, "\n")
  quit(status=1)
}

# Create output directory
output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Parse methods
methods <- trimws(strsplit(opt$methods, ",")[[1]])
valid_methods <- c("center", "consensus", "integrated", "muscle", "clustal")
invalid_methods <- setdiff(methods, valid_methods)
if (length(invalid_methods) > 0) {
  cat("Error: Invalid alignment methods:", paste(invalid_methods, collapse=", "), "\n")
  cat("Valid methods:", paste(valid_methods, collapse=", "), "\n")
  quit(status=1)
}

if (opt$verbose) {
  cat("Starting alignment quality testing...\n")
  cat("Input sequences:", opt$sequences, "\n")
  cat("Methods to test:", paste(methods, collapse=", "), "\n")
  cat("Output file:", opt$output, "\n")
  cat("Max sequences:", opt$`max-sequences`, "\n")
}

#' Test alignment quality for different methods
#' @param sequences_file Path to input sequences
#' @param methods Vector of alignment methods to test
#' @param max_sequences Maximum number of sequences to use
#' @param n_cores Number of cores for parallel processing
#' @param verbose Enable verbose output
test_alignment_quality <- function(sequences_file, methods = c("center", "consensus", "integrated"),
                                 max_sequences = 5000, n_cores = 1, verbose = FALSE) {
  
  # Load sequences
  if (verbose) cat("Loading sequences from:", sequences_file, "\n")
  sequences <- readDNAStringSet(sequences_file)
  
  # Subsample if needed
  if (length(sequences) > max_sequences) {
    if (verbose) cat("Subsampling to", max_sequences, "sequences\n")
    sequences <- sequences[sample(length(sequences), max_sequences)]
  }
  
  if (verbose) cat("Using", length(sequences), "sequences for testing\n")
  
  # Test each method
  results <- list()
  
  for (method in methods) {
    if (verbose) cat("Testing alignment method:", method, "\n")
    
    method_result <- tryCatch({
      test_single_alignment_method(sequences, method, verbose)
    }, error = function(e) {
      list(
        method = method,
        success = FALSE,
        error = e$message,
        metrics = NULL
      )
    })
    
    results[[method]] <- method_result
  }
  
  # Compare methods
  comparison <- compare_alignment_methods(results, verbose)
  
  return(list(
    individual_results = results,
    comparison = comparison,
    best_method = comparison$best_method,
    summary = generate_alignment_summary(results, comparison)
  ))
}

#' Test a single alignment method
test_single_alignment_method <- function(sequences, method, verbose = FALSE) {
  start_time <- Sys.time()
  
  # Calculate pre-alignment information content
  pre_alignment_ic <- calculate_unaligned_ic(sequences)
  
  # Perform alignment
  aligned_sequences <- perform_alignment(sequences, method, verbose)
  
  if (is.null(aligned_sequences)) {
    return(list(
      method = method,
      success = FALSE,
      error = "Alignment failed",
      metrics = NULL
    ))
  }
  
  # Calculate post-alignment metrics
  post_alignment_ic <- calculate_aligned_ic(aligned_sequences)
  
  # Calculate alignment quality metrics
  metrics <- calculate_alignment_metrics(sequences, aligned_sequences, method)
  metrics$pre_alignment_ic <- pre_alignment_ic
  metrics$post_alignment_ic <- post_alignment_ic
  metrics$ic_improvement <- post_alignment_ic - pre_alignment_ic
  metrics$ic_fold_change <- post_alignment_ic / pre_alignment_ic
  metrics$execution_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  return(list(
    method = method,
    success = TRUE,
    aligned_sequences = aligned_sequences,
    metrics = metrics
  ))
}

#' Perform sequence alignment using specified method
perform_alignment <- function(sequences, method, verbose = FALSE) {
  switch(method,
    "center" = align_by_center(sequences, verbose),
    "consensus" = align_by_consensus(sequences, verbose),
    "integrated" = align_integrated(sequences, verbose),
    "muscle" = align_with_muscle(sequences, verbose),
    "clustal" = align_with_clustal(sequences, verbose),
    {
      if (verbose) cat("Unknown method:", method, "\n")
      NULL
    }
  )
}

#' Align sequences by centering
align_by_center <- function(sequences, verbose = FALSE) {
  if (verbose) cat("  Aligning by center method\n")
  
  # Find the median length
  lengths <- width(sequences)
  target_length <- median(lengths)
  
  # Center each sequence
  aligned <- DNAStringSet(sapply(sequences, function(seq) {
    seq_len <- width(seq)
    if (seq_len >= target_length) {
      # Trim from both ends
      start_pos <- floor((seq_len - target_length) / 2) + 1
      end_pos <- start_pos + target_length - 1
      subseq(seq, start_pos, end_pos)
    } else {
      # Pad with Ns
      pad_total <- target_length - seq_len
      pad_left <- floor(pad_total / 2)
      pad_right <- pad_total - pad_left
      paste0(paste(rep("N", pad_left), collapse=""), as.character(seq), paste(rep("N", pad_right), collapse=""))
    }
  }))
  
  return(aligned)
}

#' Align sequences by consensus
align_by_consensus <- function(sequences, verbose = FALSE) {
  if (verbose) cat("  Aligning by consensus method\n")
  
  # Simple consensus-based alignment
  # This is a placeholder - in practice would use more sophisticated methods
  
  # Find the most common length
  lengths <- width(sequences)
  common_length <- as.numeric(names(sort(table(lengths), decreasing = TRUE))[1])
  
  # Filter to sequences of common length
  common_seqs <- sequences[width(sequences) == common_length]
  
  if (length(common_seqs) < length(sequences) * 0.1) {
    # If too few sequences of common length, use median length approach
    return(align_by_center(sequences, verbose))
  }
  
  # Build consensus from common length sequences
  consensus <- build_simple_consensus(common_seqs)
  
  # Align all sequences to consensus (simplified)
  aligned <- align_to_consensus(sequences, consensus, verbose)
  
  return(aligned)
}

#' Integrated alignment method
align_integrated <- function(sequences, verbose = FALSE) {
  if (verbose) cat("  Aligning using integrated method\n")
  
  # Combine center and consensus approaches
  # First pass: center alignment
  center_aligned <- align_by_center(sequences, verbose)
  
  # Second pass: consensus refinement
  consensus <- build_simple_consensus(center_aligned)
  final_aligned <- align_to_consensus(center_aligned, consensus, verbose)
  
  return(final_aligned)
}

#' Align with MUSCLE (external tool)
align_with_muscle <- function(sequences, verbose = FALSE) {
  if (verbose) cat("  Aligning with MUSCLE (placeholder)\n")
  
  # This would call external MUSCLE alignment tool
  # For now, return center alignment as fallback
  return(align_by_center(sequences, verbose))
}

#' Align with ClustalW (external tool)
align_with_clustal <- function(sequences, verbose = FALSE) {
  if (verbose) cat("  Aligning with ClustalW (placeholder)\n")
  
  # This would call external ClustalW alignment tool
  # For now, return center alignment as fallback
  return(align_by_center(sequences, verbose))
}

#' Build simple consensus sequence
build_simple_consensus <- function(sequences) {
  if (length(sequences) == 0) return(NULL)
  
  # Convert to matrix
  seq_matrix <- do.call(rbind, strsplit(as.character(sequences), ""))
  
  # Find most common base at each position
  consensus_chars <- apply(seq_matrix, 2, function(col) {
    tab <- table(col)
    names(tab)[which.max(tab)]
  })
  
  return(paste(consensus_chars, collapse = ""))
}

#' Align sequences to consensus
align_to_consensus <- function(sequences, consensus, verbose = FALSE) {
  # Simplified alignment to consensus
  # In practice, would use dynamic programming alignment
  
  target_length <- nchar(consensus)
  
  aligned <- DNAStringSet(sapply(sequences, function(seq) {
    seq_str <- as.character(seq)
    seq_len <- nchar(seq_str)
    
    if (seq_len >= target_length) {
      # Simple approach: take middle portion
      start_pos <- floor((seq_len - target_length) / 2) + 1
      substr(seq_str, start_pos, start_pos + target_length - 1)
    } else {
      # Pad sequence
      pad_total <- target_length - seq_len
      pad_left <- floor(pad_total / 2)
      pad_right <- pad_total - pad_left
      paste0(paste(rep("N", pad_left), collapse=""), seq_str, paste(rep("N", pad_right), collapse=""))
    }
  }))
  
  return(aligned)
}

#' Calculate information content for unaligned sequences
calculate_unaligned_ic <- function(sequences) {
  # For unaligned sequences, calculate average per-position IC
  if (length(sequences) == 0) return(0)
  
  # Use the most common length
  lengths <- width(sequences)
  common_length <- as.numeric(names(sort(table(lengths), decreasing = TRUE))[1])
  
  # Filter to common length and calculate IC
  common_seqs <- sequences[width(sequences) == common_length]
  if (length(common_seqs) < 10) return(0)
  
  return(calculate_aligned_ic(common_seqs))
}

#' Calculate information content for aligned sequences
calculate_aligned_ic <- function(sequences) {
  if (length(sequences) == 0) return(0)
  
  # Convert to matrix
  seq_matrix <- do.call(rbind, strsplit(as.character(sequences), ""))
  
  # Calculate information content for each position
  ic_values <- apply(seq_matrix, 2, function(col) {
    # Remove N's for IC calculation
    valid_bases <- col[col %in% c("A", "T", "G", "C")]
    if (length(valid_bases) < 5) return(0)
    
    # Calculate frequencies
    freqs <- table(valid_bases)
    probs <- freqs / sum(freqs)
    
    # Calculate information content
    # IC = 2 + sum(p * log2(p)) where p are probabilities
    ic <- 2 + sum(probs * log2(probs + 1e-10))  # Add small value to avoid log(0)
    return(max(0, ic))  # IC should be non-negative
  })
  
  return(sum(ic_values))
}

#' Calculate alignment quality metrics
calculate_alignment_metrics <- function(original_sequences, aligned_sequences, method) {
  metrics <- list(method = method)
  
  # Basic metrics
  metrics$n_sequences <- length(aligned_sequences)
  metrics$alignment_length <- unique(width(aligned_sequences))[1]
  metrics$original_length_range <- range(width(original_sequences))
  
  # Conservation metrics
  seq_matrix <- do.call(rbind, strsplit(as.character(aligned_sequences), ""))
  
  # Calculate conservation score for each position
  conservation_scores <- apply(seq_matrix, 2, function(col) {
    valid_bases <- col[col %in% c("A", "T", "G", "C")]
    if (length(valid_bases) < 5) return(0)
    
    freqs <- table(valid_bases)
    max_freq <- max(freqs)
    return(max_freq / length(valid_bases))
  })
  
  metrics$mean_conservation <- mean(conservation_scores)
  metrics$max_conservation <- max(conservation_scores)
  metrics$conserved_positions <- sum(conservation_scores > 0.7)
  
  # Gap analysis
  n_count <- sum(seq_matrix == "N")
  total_bases <- length(seq_matrix)
  metrics$gap_percentage <- (n_count / total_bases) * 100
  
  return(metrics)
}

#' Compare alignment methods
compare_alignment_methods <- function(results, verbose = FALSE) {
  if (verbose) cat("Comparing alignment methods...\n")
  
  # Extract successful results
  successful_results <- results[sapply(results, function(x) x$success)]
  
  if (length(successful_results) == 0) {
    return(list(
      best_method = NULL,
      ranking = NULL,
      comparison_table = NULL
    ))
  }
  
  # Create comparison table
  comparison_data <- do.call(rbind, lapply(successful_results, function(result) {
    metrics <- result$metrics
    data.frame(
      method = metrics$method,
      ic_improvement = metrics$ic_improvement,
      ic_fold_change = metrics$ic_fold_change,
      mean_conservation = metrics$mean_conservation,
      conserved_positions = metrics$conserved_positions,
      gap_percentage = metrics$gap_percentage,
      execution_time = metrics$execution_time,
      stringsAsFactors = FALSE
    )
  }))
  
  # Calculate overall score (weighted combination of metrics)
  comparison_data$overall_score <- (
    0.4 * scale_metric(comparison_data$ic_improvement) +
    0.3 * scale_metric(comparison_data$mean_conservation) +
    0.2 * scale_metric(comparison_data$conserved_positions) +
    0.1 * scale_metric(-comparison_data$gap_percentage)  # Negative because lower is better
  )
  
  # Rank methods
  comparison_data <- comparison_data[order(comparison_data$overall_score, decreasing = TRUE), ]
  comparison_data$rank <- 1:nrow(comparison_data)
  
  best_method <- comparison_data$method[1]
  
  return(list(
    best_method = best_method,
    ranking = comparison_data$method,
    comparison_table = comparison_data
  ))
}

#' Scale metric to 0-1 range
scale_metric <- function(x) {
  if (length(x) <= 1 || var(x) == 0) return(rep(0.5, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

#' Generate alignment summary
generate_alignment_summary <- function(results, comparison) {
  summary <- list(
    n_methods_tested = length(results),
    n_successful = sum(sapply(results, function(x) x$success)),
    best_method = comparison$best_method
  )
  
  if (!is.null(comparison$comparison_table)) {
    best_result <- comparison$comparison_table[1, ]
    summary$best_ic_improvement <- best_result$ic_improvement
    summary$best_conservation <- best_result$mean_conservation
    summary$best_execution_time <- best_result$execution_time
  }
  
  return(summary)
}

#' Generate alignment validation report
generate_alignment_report <- function(test_results, output_file) {
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>CTCF Alignment Quality Validation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #f0f8ff; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }
        .best { background-color: #e8f5e8; }
        .poor { background-color: #ffeaea; }
        table { border-collapse: collapse; width: 100%; margin: 10px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .metric { margin: 10px 0; }
        .success { color: green; font-weight: bold; }
        .failure { color: red; font-weight: bold; }
    </style>
</head>
<body>
    <div class='header'>
        <h1>CTCF Alignment Quality Validation Report</h1>
        <p>Generated: ", Sys.time(), "</p>
        <p>Methods Tested: ", test_results$summary$n_methods_tested, "</p>
        <p>Successful Methods: ", test_results$summary$n_successful, "</p>
        <p>Best Method: <strong>", ifelse(is.null(test_results$best_method), "None", test_results$best_method), "</strong></p>
    </div>")
  
  # Individual method results
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Individual Method Results</h2>")
  
  for (method_name in names(test_results$individual_results)) {
    result <- test_results$individual_results[[method_name]]
    
    section_class <- if (result$success) {
      if (method_name == test_results$best_method) "best" else ""
    } else {
      "poor"
    }
    
    html_content <- paste0(html_content, "
        <div class='section ", section_class, "'>
            <h3>", method_name, " ", if (result$success) "<span class='success'>✓</span>" else "<span class='failure'>✗</span>", "</h3>")
    
    if (result$success) {
      metrics <- result$metrics
      html_content <- paste0(html_content, "
            <div class='metric'>IC Improvement: ", round(metrics$ic_improvement, 2), " bits</div>
            <div class='metric'>IC Fold Change: ", round(metrics$ic_fold_change, 2), "x</div>
            <div class='metric'>Mean Conservation: ", round(metrics$mean_conservation, 3), "</div>
            <div class='metric'>Conserved Positions: ", metrics$conserved_positions, "</div>
            <div class='metric'>Gap Percentage: ", round(metrics$gap_percentage, 1), "%</div>
            <div class='metric'>Execution Time: ", round(metrics$execution_time, 1), " seconds</div>")
    } else {
      html_content <- paste0(html_content, "
            <div class='metric'><span class='failure'>Error: ", result$error, "</span></div>")
    }
    
    html_content <- paste0(html_content, "
        </div>")
  }
  
  html_content <- paste0(html_content, "
    </div>")
  
  # Comparison table
  if (!is.null(test_results$comparison$comparison_table)) {
    html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Method Comparison</h2>
        <table>
            <tr>
                <th>Rank</th><th>Method</th><th>IC Improvement</th><th>IC Fold Change</th>
                <th>Mean Conservation</th><th>Conserved Positions</th><th>Gap %</th>
                <th>Time (s)</th><th>Overall Score</th>
            </tr>")
    
    for (i in 1:nrow(test_results$comparison$comparison_table)) {
      row <- test_results$comparison$comparison_table[i, ]
      row_class <- if (i == 1) "best" else ""
      
      html_content <- paste0(html_content, "
            <tr class='", row_class, "'>
                <td>", row$rank, "</td>
                <td>", row$method, "</td>
                <td>", round(row$ic_improvement, 2), "</td>
                <td>", round(row$ic_fold_change, 2), "</td>
                <td>", round(row$mean_conservation, 3), "</td>
                <td>", row$conserved_positions, "</td>
                <td>", round(row$gap_percentage, 1), "</td>
                <td>", round(row$execution_time, 1), "</td>
                <td>", round(row$overall_score, 3), "</td>
            </tr>")
    }
    
    html_content <- paste0(html_content, "
        </table>
    </div>")
  }
  
  html_content <- paste0(html_content, "
</body>
</html>")
  
  writeLines(html_content, output_file)
}

# Main execution
test_results <- test_alignment_quality(
  sequences_file = opt$sequences,
  methods = methods,
  max_sequences = opt$`max-sequences`,
  n_cores = opt$`n-cores`,
  verbose = opt$verbose
)

if (opt$verbose) cat("Generating alignment validation report...\n")

generate_alignment_report(test_results, opt$output)

if (opt$verbose) {
  cat("Alignment quality testing complete!\n")
  if (!is.null(test_results$best_method)) {
    cat("Best method:", test_results$best_method, "\n")
  }
  cat("Report saved to:", opt$output, "\n")
}

# Exit with appropriate status
if (test_results$summary$n_successful > 0) {
  quit(status = 0)
} else {
  quit(status = 1)
}
