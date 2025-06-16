#!/usr/bin/env Rscript

# Test PWM Construction Script
# Tests robustness and quality of PWM construction methods
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
              help="Input aligned sequences file (FASTA format) [required]", metavar="character"),
  make_option(c("--bootstrap-runs"), type="integer", default=100,
              help="Number of bootstrap runs [default: %default]", metavar="integer"),
  make_option(c("--output", "-o"), type="character", default="results/pwm_validation.html",
              help="Output validation report file [default: %default]", metavar="character"),
  make_option(c("--pseudocounts"), type="character", default="0.01,0.1,0.5",
              help="Pseudocount values to test (comma-separated) [default: %default]", metavar="character"),
  make_option(c("--n-cores"), type="integer", default=1,
              help="Number of cores for parallel processing [default: %default]", metavar="integer"),
  make_option(c("--confidence-level"), type="double", default=0.95,
              help="Confidence level for intervals [default: %default]", metavar="double"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Test PWM construction robustness and quality")
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

# Parse pseudocounts
pseudocounts <- as.numeric(trimws(strsplit(opt$pseudocounts, ",")[[1]]))
if (any(is.na(pseudocounts)) || any(pseudocounts <= 0)) {
  cat("Error: Invalid pseudocount values. Must be positive numbers.\n")
  quit(status=1)
}

if (opt$verbose) {
  cat("Starting PWM construction validation...\n")
  cat("Input sequences:", opt$sequences, "\n")
  cat("Bootstrap runs:", opt$`bootstrap-runs`, "\n")
  cat("Pseudocounts:", paste(pseudocounts, collapse=", "), "\n")
  cat("Output file:", opt$output, "\n")
  cat("Confidence level:", opt$`confidence-level`, "\n")
}

#' Test PWM construction robustness
#' @param sequences_file Path to aligned sequences file
#' @param bootstrap_runs Number of bootstrap replicates
#' @param pseudocounts Vector of pseudocount values to test
#' @param confidence_level Confidence level for intervals
#' @param n_cores Number of cores for parallel processing
#' @param verbose Enable verbose output
test_pwm_construction <- function(sequences_file, bootstrap_runs = 100, 
                                pseudocounts = c(0.01, 0.1, 0.5),
                                confidence_level = 0.95, n_cores = 1, verbose = FALSE) {
  
  # Load sequences
  if (verbose) cat("Loading aligned sequences from:", sequences_file, "\n")
  sequences <- readDNAStringSet(sequences_file)
  
  if (verbose) cat("Loaded", length(sequences), "aligned sequences\n")
  
  # Test different pseudocount values
  pseudocount_results <- list()
  
  for (pc in pseudocounts) {
    if (verbose) cat("Testing pseudocount:", pc, "\n")
    
    pc_result <- test_pseudocount_robustness(
      sequences, pc, bootstrap_runs, confidence_level, n_cores, verbose
    )
    
    pseudocount_results[[as.character(pc)]] <- pc_result
  }
  
  # Test construction stability
  if (verbose) cat("Testing construction stability...\n")
  stability_result <- test_construction_stability(sequences, bootstrap_runs, verbose)
  
  # Test parameter sensitivity
  if (verbose) cat("Testing parameter sensitivity...\n")
  sensitivity_result <- test_parameter_sensitivity(sequences, verbose)
  
  # Compare methods and select best
  comparison <- compare_construction_methods(pseudocount_results, verbose)
  
  return(list(
    pseudocount_results = pseudocount_results,
    stability_result = stability_result,
    sensitivity_result = sensitivity_result,
    comparison = comparison,
    best_pseudocount = comparison$best_pseudocount,
    summary = generate_construction_summary(pseudocount_results, stability_result, sensitivity_result)
  ))
}

#' Test PWM construction with different pseudocounts
test_pseudocount_robustness <- function(sequences, pseudocount, bootstrap_runs, 
                                      confidence_level, n_cores, verbose = FALSE) {
  
  if (verbose) cat("  Running", bootstrap_runs, "bootstrap replicates\n")
  
  # Bootstrap replicates
  if (n_cores > 1) {
    bootstrap_results <- mclapply(1:bootstrap_runs, function(i) {
      boot_sequences <- sequences[sample(length(sequences), replace = TRUE)]
      build_pwm_with_metrics(boot_sequences, pseudocount)
    }, mc.cores = n_cores)
  } else {
    bootstrap_results <- lapply(1:bootstrap_runs, function(i) {
      boot_sequences <- sequences[sample(length(sequences), replace = TRUE)]
      build_pwm_with_metrics(boot_sequences, pseudocount)
    })
  }
  
  # Filter successful results
  successful_results <- bootstrap_results[sapply(bootstrap_results, function(x) !is.null(x))]
  
  if (length(successful_results) == 0) {
    return(list(
      pseudocount = pseudocount,
      success = FALSE,
      error = "All bootstrap runs failed"
    ))
  }
  
  # Calculate original PWM
  original_pwm <- build_pwm_with_metrics(sequences, pseudocount)
  
  # Extract metrics from bootstrap results
  ic_values <- sapply(successful_results, function(x) x$total_ic)
  conservation_scores <- sapply(successful_results, function(x) x$mean_conservation)
  entropy_values <- sapply(successful_results, function(x) x$mean_entropy)
  
  # Calculate confidence intervals
  alpha <- 1 - confidence_level
  ic_ci <- quantile(ic_values, c(alpha/2, 1-alpha/2))
  conservation_ci <- quantile(conservation_scores, c(alpha/2, 1-alpha/2))
  entropy_ci <- quantile(entropy_values, c(alpha/2, 1-alpha/2))
  
  # Calculate stability metrics
  ic_cv <- sd(ic_values) / mean(ic_values)  # Coefficient of variation
  conservation_cv <- sd(conservation_scores) / mean(conservation_scores)
  
  return(list(
    pseudocount = pseudocount,
    success = TRUE,
    n_successful = length(successful_results),
    n_total = bootstrap_runs,
    original_pwm = original_pwm,
    bootstrap_stats = list(
      ic = list(
        mean = mean(ic_values),
        sd = sd(ic_values),
        cv = ic_cv,
        ci = ic_ci,
        original = original_pwm$total_ic
      ),
      conservation = list(
        mean = mean(conservation_scores),
        sd = sd(conservation_scores),
        cv = conservation_cv,
        ci = conservation_ci,
        original = original_pwm$mean_conservation
      ),
      entropy = list(
        mean = mean(entropy_values),
        sd = sd(entropy_values),
        ci = entropy_ci
      )
    ),
    stability_score = calculate_stability_score(ic_cv, conservation_cv)
  ))
}

#' Build PWM and calculate metrics
build_pwm_with_metrics <- function(sequences, pseudocount = 0.1) {
  tryCatch({
    # Convert sequences to matrix
    seq_matrix <- do.call(rbind, strsplit(as.character(sequences), ""))
    
    # Remove positions with too many Ns
    valid_positions <- apply(seq_matrix, 2, function(col) {
      sum(col %in% c("A", "T", "G", "C")) >= length(col) * 0.5
    })
    
    if (sum(valid_positions) < 5) {
      return(NULL)  # Not enough valid positions
    }
    
    seq_matrix <- seq_matrix[, valid_positions, drop = FALSE]
    
    # Build PWM
    pwm <- build_pwm_matrix(seq_matrix, pseudocount)
    
    # Calculate metrics
    ic_values <- calculate_position_ic(pwm)
    total_ic <- sum(ic_values)
    
    conservation_scores <- apply(seq_matrix, 2, function(col) {
      valid_bases <- col[col %in% c("A", "T", "G", "C")]
      if (length(valid_bases) < 5) return(0)
      freqs <- table(valid_bases)
      max(freqs) / length(valid_bases)
    })
    
    entropy_values <- apply(pwm, 2, function(col) {
      -sum(col * log2(col + 1e-10))
    })
    
    return(list(
      pwm = pwm,
      ic_values = ic_values,
      total_ic = total_ic,
      conservation_scores = conservation_scores,
      mean_conservation = mean(conservation_scores),
      entropy_values = entropy_values,
      mean_entropy = mean(entropy_values),
      n_positions = ncol(pwm),
      pseudocount = pseudocount
    ))
    
  }, error = function(e) {
    return(NULL)
  })
}

#' Build PWM matrix from sequence matrix
build_pwm_matrix <- function(seq_matrix, pseudocount = 0.1) {
  bases <- c("A", "C", "G", "T")
  n_positions <- ncol(seq_matrix)
  pwm <- matrix(0, nrow = 4, ncol = n_positions)
  rownames(pwm) <- bases
  
  for (pos in 1:n_positions) {
    col <- seq_matrix[, pos]
    valid_bases <- col[col %in% bases]
    
    if (length(valid_bases) > 0) {
      counts <- table(factor(valid_bases, levels = bases))
      freqs <- (counts + pseudocount) / (sum(counts) + 4 * pseudocount)
      pwm[, pos] <- freqs
    } else {
      # Equal frequencies if no valid bases
      pwm[, pos] <- rep(0.25, 4)
    }
  }
  
  return(pwm)
}

#' Calculate information content for each position
calculate_position_ic <- function(pwm) {
  apply(pwm, 2, function(col) {
    # IC = 2 + sum(p * log2(p))
    2 + sum(col * log2(col + 1e-10))
  })
}

#' Test construction stability across different subsets
test_construction_stability <- function(sequences, n_tests = 50, verbose = FALSE) {
  
  subset_sizes <- c(0.5, 0.7, 0.8, 0.9)  # Proportion of sequences to use
  stability_results <- list()
  
  for (subset_size in subset_sizes) {
    if (verbose) cat("  Testing subset size:", subset_size, "\n")
    
    n_seqs <- round(length(sequences) * subset_size)
    subset_ics <- replicate(n_tests, {
      subset_seqs <- sequences[sample(length(sequences), n_seqs)]
      pwm_result <- build_pwm_with_metrics(subset_seqs)
      if (is.null(pwm_result)) NA else pwm_result$total_ic
    })
    
    subset_ics <- subset_ics[!is.na(subset_ics)]
    
    if (length(subset_ics) > 5) {
      stability_results[[as.character(subset_size)]] <- list(
        subset_size = subset_size,
        mean_ic = mean(subset_ics),
        sd_ic = sd(subset_ics),
        cv_ic = sd(subset_ics) / mean(subset_ics),
        n_successful = length(subset_ics)
      )
    }
  }
  
  return(stability_results)
}

#' Test parameter sensitivity
test_parameter_sensitivity <- function(sequences, verbose = FALSE) {
  
  # Test different minimum sequence thresholds
  min_seqs_thresholds <- c(10, 50, 100, 200)
  min_seq_results <- list()
  
  for (min_seqs in min_seqs_thresholds) {
    if (length(sequences) >= min_seqs) {
      subset_seqs <- sequences[1:min_seqs]
      pwm_result <- build_pwm_with_metrics(subset_seqs)
      
      if (!is.null(pwm_result)) {
        min_seq_results[[as.character(min_seqs)]] <- list(
          min_sequences = min_seqs,
          total_ic = pwm_result$total_ic,
          mean_conservation = pwm_result$mean_conservation
        )
      }
    }
  }
  
  return(list(
    min_sequences_test = min_seq_results
  ))
}

#' Calculate stability score
calculate_stability_score <- function(ic_cv, conservation_cv) {
  # Lower CV indicates higher stability
  # Convert to 0-1 scale where 1 is most stable
  ic_score <- exp(-ic_cv * 10)  # Exponential decay
  conservation_score <- exp(-conservation_cv * 10)
  
  # Combined score
  (ic_score + conservation_score) / 2
}

#' Compare construction methods
compare_construction_methods <- function(pseudocount_results, verbose = FALSE) {
  
  successful_results <- pseudocount_results[sapply(pseudocount_results, function(x) x$success)]
  
  if (length(successful_results) == 0) {
    return(list(
      best_pseudocount = NULL,
      comparison_table = NULL,
      ranking = NULL
    ))
  }
  
  # Create comparison table
  comparison_data <- do.call(rbind, lapply(names(successful_results), function(pc) {
    result <- successful_results[[pc]]
    data.frame(
      pseudocount = as.numeric(pc),
      mean_ic = result$bootstrap_stats$ic$mean,
      ic_cv = result$bootstrap_stats$ic$cv,
      mean_conservation = result$bootstrap_stats$conservation$mean,
      conservation_cv = result$bootstrap_stats$conservation$cv,
      stability_score = result$stability_score,
      success_rate = result$n_successful / result$n_total,
      stringsAsFactors = FALSE
    )
  }))
  
  # Calculate overall score
  comparison_data$overall_score <- (
    0.3 * scale_metric(comparison_data$mean_ic) +
    0.2 * scale_metric(-comparison_data$ic_cv) +  # Negative because lower CV is better
    0.2 * scale_metric(comparison_data$mean_conservation) +
    0.2 * scale_metric(comparison_data$stability_score) +
    0.1 * scale_metric(comparison_data$success_rate)
  )
  
  # Rank methods
  comparison_data <- comparison_data[order(comparison_data$overall_score, decreasing = TRUE), ]
  comparison_data$rank <- 1:nrow(comparison_data)
  
  best_pseudocount <- comparison_data$pseudocount[1]
  
  return(list(
    best_pseudocount = best_pseudocount,
    comparison_table = comparison_data,
    ranking = comparison_data$pseudocount
  ))
}

#' Scale metric to 0-1 range
scale_metric <- function(x) {
  if (length(x) <= 1 || var(x) == 0) return(rep(0.5, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

#' Generate construction summary
generate_construction_summary <- function(pseudocount_results, stability_result, sensitivity_result) {
  successful_pseudocounts <- sum(sapply(pseudocount_results, function(x) x$success))
  
  summary <- list(
    n_pseudocounts_tested = length(pseudocount_results),
    n_successful = successful_pseudocounts,
    stability_tests_completed = length(stability_result) > 0,
    sensitivity_tests_completed = length(sensitivity_result$min_sequences_test) > 0
  )
  
  # Add best metrics if available
  successful_results <- pseudocount_results[sapply(pseudocount_results, function(x) x$success)]
  if (length(successful_results) > 0) {
    best_ic <- max(sapply(successful_results, function(x) x$bootstrap_stats$ic$mean))
    best_stability <- max(sapply(successful_results, function(x) x$stability_score))
    
    summary$best_ic <- best_ic
    summary$best_stability <- best_stability
  }
  
  return(summary)
}

#' Generate PWM construction validation report
generate_construction_report <- function(test_results, output_file) {
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>CTCF PWM Construction Validation Report</title>
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
        .ci { font-style: italic; color: #666; }
    </style>
</head>
<body>
    <div class='header'>
        <h1>CTCF PWM Construction Validation Report</h1>
        <p>Generated: ", Sys.time(), "</p>
        <p>Pseudocounts Tested: ", test_results$summary$n_pseudocounts_tested, "</p>
        <p>Successful Tests: ", test_results$summary$n_successful, "</p>
        <p>Best Pseudocount: <strong>", ifelse(is.null(test_results$best_pseudocount), "None", test_results$best_pseudocount), "</strong></p>
    </div>")
  
  # Pseudocount comparison
  if (!is.null(test_results$comparison$comparison_table)) {
    html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Pseudocount Comparison</h2>
        <table>
            <tr>
                <th>Rank</th><th>Pseudocount</th><th>Mean IC</th><th>IC CV</th>
                <th>Mean Conservation</th><th>Stability Score</th><th>Success Rate</th>
                <th>Overall Score</th>
            </tr>")
    
    for (i in 1:nrow(test_results$comparison$comparison_table)) {
      row <- test_results$comparison$comparison_table[i, ]
      row_class <- if (i == 1) "best" else ""
      
      html_content <- paste0(html_content, "
            <tr class='", row_class, "'>
                <td>", row$rank, "</td>
                <td>", row$pseudocount, "</td>
                <td>", round(row$mean_ic, 2), "</td>
                <td>", round(row$ic_cv, 3), "</td>
                <td>", round(row$mean_conservation, 3), "</td>
                <td>", round(row$stability_score, 3), "</td>
                <td>", round(row$success_rate, 3), "</td>
                <td>", round(row$overall_score, 3), "</td>
            </tr>")
    }
    
    html_content <- paste0(html_content, "
        </table>
    </div>")
  }
  
  # Individual pseudocount results
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Detailed Bootstrap Results</h2>")
  
  for (pc_name in names(test_results$pseudocount_results)) {
    result <- test_results$pseudocount_results[[pc_name]]
    
    section_class <- if (result$success) {
      if (as.numeric(pc_name) == test_results$best_pseudocount) "best" else ""
    } else {
      "poor"
    }
    
    html_content <- paste0(html_content, "
        <div class='section ", section_class, "'>
            <h3>Pseudocount: ", pc_name, " ", if (result$success) "<span class='success'>✓</span>" else "<span class='failure'>✗</span>", "</h3>")
    
    if (result$success) {
      ic_stats <- result$bootstrap_stats$ic
      cons_stats <- result$bootstrap_stats$conservation
      
      html_content <- paste0(html_content, "
            <div class='metric'>Bootstrap Success Rate: ", round(result$n_successful / result$n_total * 100, 1), "%</div>
            <div class='metric'>Information Content: ", round(ic_stats$mean, 2), " ± ", round(ic_stats$sd, 2), " bits</div>
            <div class='metric'>IC Confidence Interval: <span class='ci'>[", round(ic_stats$ci[1], 2), ", ", round(ic_stats$ci[2], 2), "]</span></div>
            <div class='metric'>IC Coefficient of Variation: ", round(ic_stats$cv, 3), "</div>
            <div class='metric'>Mean Conservation: ", round(cons_stats$mean, 3), " ± ", round(cons_stats$sd, 3), "</div>
            <div class='metric'>Conservation CI: <span class='ci'>[", round(cons_stats$ci[1], 3), ", ", round(cons_stats$ci[2], 3), "]</span></div>
            <div class='metric'>Stability Score: ", round(result$stability_score, 3), "</div>")
    } else {
      html_content <- paste0(html_content, "
            <div class='metric'><span class='failure'>Error: ", result$error, "</span></div>")
    }
    
    html_content <- paste0(html_content, "
        </div>")
  }
  
  html_content <- paste0(html_content, "
    </div>")
  
  # Stability results
  if (length(test_results$stability_result) > 0) {
    html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Construction Stability Analysis</h2>
        <table>
            <tr><th>Subset Size</th><th>Mean IC</th><th>IC Standard Deviation</th><th>IC CV</th><th>Successful Runs</th></tr>")
    
    for (subset_name in names(test_results$stability_result)) {
      stability <- test_results$stability_result[[subset_name]]
      html_content <- paste0(html_content, "
            <tr>
                <td>", round(stability$subset_size * 100, 0), "%</td>
                <td>", round(stability$mean_ic, 2), "</td>
                <td>", round(stability$sd_ic, 2), "</td>
                <td>", round(stability$cv_ic, 3), "</td>
                <td>", stability$n_successful, "</td>
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
test_results <- test_pwm_construction(
  sequences_file = opt$sequences,
  bootstrap_runs = opt$`bootstrap-runs`,
  pseudocounts = pseudocounts,
  confidence_level = opt$`confidence-level`,
  n_cores = opt$`n-cores`,
  verbose = opt$verbose
)

if (opt$verbose) cat("Generating PWM construction validation report...\n")

generate_construction_report(test_results, opt$output)

if (opt$verbose) {
  cat("PWM construction validation complete!\n")
  if (!is.null(test_results$best_pseudocount)) {
    cat("Best pseudocount:", test_results$best_pseudocount, "\n")
  }
  cat("Report saved to:", opt$output, "\n")
}

# Exit with appropriate status
if (test_results$summary$n_successful > 0) {
  quit(status = 0)
} else {
  quit(status = 1)
}
