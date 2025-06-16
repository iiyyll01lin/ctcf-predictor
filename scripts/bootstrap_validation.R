#!/usr/bin/env Rscript

# Bootstrap Validation Script
# Provides bootstrap confidence intervals for PWM validation
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
  make_option(c("--bootstrap-samples"), type="integer", default=1000,
              help="Number of bootstrap samples [default: %default]", metavar="integer"),
  make_option(c("--confidence-level"), type="double", default=0.95,
              help="Confidence level for intervals [default: %default]", metavar="double"),
  make_option(c("--output", "-o"), type="character", default="results/bootstrap_validation.html",
              help="Output validation report file [default: %default]", metavar="character"),
  make_option(c("--pseudocount"), type="double", default=0.1,
              help="Pseudocount for PWM construction [default: %default]", metavar="double"),
  make_option(c("--n-cores"), type="integer", default=1,
              help="Number of cores for parallel processing [default: %default]", metavar="integer"),
  make_option(c("--metrics"), type="character", default="ic,conservation,entropy",
              help="Metrics to bootstrap (comma-separated) [default: %default]", metavar="character"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Generate bootstrap confidence intervals for PWM validation")
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

if (opt$`confidence-level` <= 0 || opt$`confidence-level` >= 1) {
  cat("Error: Confidence level must be between 0 and 1\n")
  quit(status=1)
}

# Create output directory
output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Parse metrics
metrics <- trimws(strsplit(opt$metrics, ",")[[1]])
valid_metrics <- c("ic", "conservation", "entropy", "gc_content", "length")
invalid_metrics <- setdiff(metrics, valid_metrics)
if (length(invalid_metrics) > 0) {
  cat("Error: Invalid metrics:", paste(invalid_metrics, collapse=", "), "\n")
  cat("Valid metrics:", paste(valid_metrics, collapse=", "), "\n")
  quit(status=1)
}

if (opt$verbose) {
  cat("Starting bootstrap validation...\n")
  cat("Input sequences:", opt$sequences, "\n")
  cat("Bootstrap samples:", opt$`bootstrap-samples`, "\n")
  cat("Confidence level:", opt$`confidence-level`, "\n")
  cat("Metrics:", paste(metrics, collapse=", "), "\n")
  cat("Output file:", opt$output, "\n")
}

#' Bootstrap confidence intervals for PWM validation
#' @param sequences_file Path to aligned sequences file
#' @param n_bootstrap Number of bootstrap samples
#' @param confidence_level Confidence level for intervals
#' @param pseudocount Pseudocount for PWM construction
#' @param metrics Vector of metrics to calculate
#' @param n_cores Number of cores for parallel processing
#' @param verbose Enable verbose output
bootstrap_validation <- function(sequences_file, n_bootstrap = 1000, confidence_level = 0.95,
                                pseudocount = 0.1, metrics = c("ic", "conservation", "entropy"),
                                n_cores = 1, verbose = FALSE) {
  
  # Load sequences
  if (verbose) cat("Loading aligned sequences from:", sequences_file, "\n")
  sequences <- readDNAStringSet(sequences_file)
  
  if (verbose) cat("Loaded", length(sequences), "aligned sequences\n")
  
  # Build original PWM
  if (verbose) cat("Building original PWM...\n")
  original_pwm <- build_pwm_with_bootstrap_metrics(sequences, pseudocount, metrics, verbose)
  
  if (is.null(original_pwm)) {
    stop("Failed to build original PWM")
  }
  
  # Perform bootstrap sampling
  if (verbose) cat("Performing", n_bootstrap, "bootstrap samples...\n")
  
  if (n_cores > 1) {
    bootstrap_results <- mclapply(1:n_bootstrap, function(i) {
      if (verbose && i %% 100 == 0) cat("  Bootstrap sample", i, "\n")
      
      # Sample with replacement
      boot_indices <- sample(length(sequences), replace = TRUE)
      boot_sequences <- sequences[boot_indices]
      
      # Build PWM and calculate metrics
      build_pwm_with_bootstrap_metrics(boot_sequences, pseudocount, metrics, verbose = FALSE)
    }, mc.cores = n_cores)
  } else {
    bootstrap_results <- lapply(1:n_bootstrap, function(i) {
      if (verbose && i %% 100 == 0) cat("  Bootstrap sample", i, "\n")
      
      # Sample with replacement
      boot_indices <- sample(length(sequences), replace = TRUE)
      boot_sequences <- sequences[boot_indices]
      
      # Build PWM and calculate metrics
      build_pwm_with_bootstrap_metrics(boot_sequences, pseudocount, metrics, verbose = FALSE)
    })
  }
  
  # Filter successful results
  successful_results <- bootstrap_results[!sapply(bootstrap_results, is.null)]
  
  if (length(successful_results) == 0) {
    stop("All bootstrap samples failed")
  }
  
  if (verbose) cat("Successfully completed", length(successful_results), "out of", n_bootstrap, "bootstrap samples\n")
  
  # Calculate confidence intervals
  confidence_intervals <- calculate_bootstrap_ci(successful_results, original_pwm, confidence_level, metrics, verbose)
  
  # Calculate bootstrap statistics
  bootstrap_stats <- calculate_bootstrap_statistics(successful_results, original_pwm, metrics, verbose)
  
  # Assess bootstrap quality
  bootstrap_quality <- assess_bootstrap_quality(successful_results, original_pwm, confidence_intervals, verbose)
  
  return(list(
    original_pwm = original_pwm,
    bootstrap_results = successful_results,
    confidence_intervals = confidence_intervals,
    bootstrap_statistics = bootstrap_stats,
    bootstrap_quality = bootstrap_quality,
    n_successful = length(successful_results),
    n_total = n_bootstrap,
    success_rate = length(successful_results) / n_bootstrap,
    confidence_level = confidence_level,
    metrics = metrics
  ))
}

#' Build PWM with bootstrap metrics
build_pwm_with_bootstrap_metrics <- function(sequences, pseudocount = 0.1, 
                                            metrics = c("ic", "conservation", "entropy"), 
                                            verbose = FALSE) {
  
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
    
    # Calculate requested metrics
    result <- list(
      pwm = pwm,
      n_sequences = length(sequences),
      n_positions = ncol(pwm),
      pseudocount = pseudocount
    )
    
    # Information content
    if ("ic" %in% metrics) {
      ic_values <- calculate_position_ic(pwm)
      result$ic_values <- ic_values
      result$total_ic <- sum(ic_values)
      result$mean_ic <- mean(ic_values)
      result$max_ic <- max(ic_values)
    }
    
    # Conservation scores
    if ("conservation" %in% metrics) {
      conservation_scores <- apply(seq_matrix, 2, function(col) {
        valid_bases <- col[col %in% c("A", "T", "G", "C")]
        if (length(valid_bases) < 5) return(0)
        freqs <- table(valid_bases)
        max(freqs) / length(valid_bases)
      })
      result$conservation_scores <- conservation_scores
      result$mean_conservation <- mean(conservation_scores)
      result$max_conservation <- max(conservation_scores)
      result$conserved_positions <- sum(conservation_scores > 0.7)
    }
    
    # Entropy
    if ("entropy" %in% metrics) {
      entropy_values <- apply(pwm, 2, function(col) {
        -sum(col * log2(col + 1e-10))
      })
      result$entropy_values <- entropy_values
      result$mean_entropy <- mean(entropy_values)
      result$min_entropy <- min(entropy_values)
    }
    
    # GC content
    if ("gc_content" %in% metrics) {
      gc_content <- calculate_sequence_gc_content(sequences)
      result$mean_gc_content <- mean(gc_content)
      result$sd_gc_content <- sd(gc_content)
    }
    
    # Sequence length
    if ("length" %in% metrics) {
      lengths <- width(sequences)
      result$mean_length <- mean(lengths)
      result$sd_length <- sd(lengths)
    }
    
    return(result)
    
  }, error = function(e) {
    if (verbose) cat("Error in PWM construction:", e$message, "\n")
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
    max(0, 2 + sum(col * log2(col + 1e-10)))  # Ensure non-negative
  })
}

#' Calculate GC content for sequences
calculate_sequence_gc_content <- function(sequences) {
  letterFrequency(sequences, letters = "GC", as.prob = TRUE)[,1]
}

#' Calculate bootstrap confidence intervals
calculate_bootstrap_ci <- function(bootstrap_results, original_pwm, confidence_level, metrics, verbose = FALSE) {
  
  if (verbose) cat("Calculating confidence intervals...\n")
  
  alpha <- 1 - confidence_level
  ci_list <- list()
  
  # Information content CIs
  if ("ic" %in% metrics) {
    total_ics <- sapply(bootstrap_results, function(x) x$total_ic)
    mean_ics <- sapply(bootstrap_results, function(x) x$mean_ic)
    max_ics <- sapply(bootstrap_results, function(x) x$max_ic)
    
    ci_list$ic <- list(
      total_ic = quantile(total_ics, c(alpha/2, 1-alpha/2), na.rm = TRUE),
      mean_ic = quantile(mean_ics, c(alpha/2, 1-alpha/2), na.rm = TRUE),
      max_ic = quantile(max_ics, c(alpha/2, 1-alpha/2), na.rm = TRUE)
    )
  }
  
  # Conservation CIs
  if ("conservation" %in% metrics) {
    mean_conservations <- sapply(bootstrap_results, function(x) x$mean_conservation)
    max_conservations <- sapply(bootstrap_results, function(x) x$max_conservation)
    conserved_positions <- sapply(bootstrap_results, function(x) x$conserved_positions)
    
    ci_list$conservation <- list(
      mean_conservation = quantile(mean_conservations, c(alpha/2, 1-alpha/2), na.rm = TRUE),
      max_conservation = quantile(max_conservations, c(alpha/2, 1-alpha/2), na.rm = TRUE),
      conserved_positions = quantile(conserved_positions, c(alpha/2, 1-alpha/2), na.rm = TRUE)
    )
  }
  
  # Entropy CIs
  if ("entropy" %in% metrics) {
    mean_entropies <- sapply(bootstrap_results, function(x) x$mean_entropy)
    min_entropies <- sapply(bootstrap_results, function(x) x$min_entropy)
    
    ci_list$entropy <- list(
      mean_entropy = quantile(mean_entropies, c(alpha/2, 1-alpha/2), na.rm = TRUE),
      min_entropy = quantile(min_entropies, c(alpha/2, 1-alpha/2), na.rm = TRUE)
    )
  }
  
  # GC content CIs
  if ("gc_content" %in% metrics) {
    mean_gc_contents <- sapply(bootstrap_results, function(x) x$mean_gc_content)
    
    ci_list$gc_content <- list(
      mean_gc_content = quantile(mean_gc_contents, c(alpha/2, 1-alpha/2), na.rm = TRUE)
    )
  }
  
  # Length CIs
  if ("length" %in% metrics) {
    mean_lengths <- sapply(bootstrap_results, function(x) x$mean_length)
    
    ci_list$length <- list(
      mean_length = quantile(mean_lengths, c(alpha/2, 1-alpha/2), na.rm = TRUE)
    )
  }
  
  return(ci_list)
}

#' Calculate bootstrap statistics
calculate_bootstrap_statistics <- function(bootstrap_results, original_pwm, metrics, verbose = FALSE) {
  
  if (verbose) cat("Calculating bootstrap statistics...\n")
  
  stats_list <- list()
  
  # Information content statistics
  if ("ic" %in% metrics) {
    total_ics <- sapply(bootstrap_results, function(x) x$total_ic)
    
    stats_list$ic <- list(
      original_total_ic = original_pwm$total_ic,
      bootstrap_mean = mean(total_ics, na.rm = TRUE),
      bootstrap_sd = sd(total_ics, na.rm = TRUE),
      bootstrap_cv = sd(total_ics, na.rm = TRUE) / mean(total_ics, na.rm = TRUE),
      bias = mean(total_ics, na.rm = TRUE) - original_pwm$total_ic,
      bias_corrected = 2 * original_pwm$total_ic - mean(total_ics, na.rm = TRUE)
    )
  }
  
  # Conservation statistics
  if ("conservation" %in% metrics) {
    mean_conservations <- sapply(bootstrap_results, function(x) x$mean_conservation)
    
    stats_list$conservation <- list(
      original_mean_conservation = original_pwm$mean_conservation,
      bootstrap_mean = mean(mean_conservations, na.rm = TRUE),
      bootstrap_sd = sd(mean_conservations, na.rm = TRUE),
      bootstrap_cv = sd(mean_conservations, na.rm = TRUE) / mean(mean_conservations, na.rm = TRUE),
      bias = mean(mean_conservations, na.rm = TRUE) - original_pwm$mean_conservation
    )
  }
  
  # Entropy statistics
  if ("entropy" %in% metrics) {
    mean_entropies <- sapply(bootstrap_results, function(x) x$mean_entropy)
    
    stats_list$entropy <- list(
      original_mean_entropy = original_pwm$mean_entropy,
      bootstrap_mean = mean(mean_entropies, na.rm = TRUE),
      bootstrap_sd = sd(mean_entropies, na.rm = TRUE),
      bootstrap_cv = sd(mean_entropies, na.rm = TRUE) / mean(mean_entropies, na.rm = TRUE)
    )
  }
  
  return(stats_list)
}

#' Assess bootstrap quality
assess_bootstrap_quality <- function(bootstrap_results, original_pwm, confidence_intervals, verbose = FALSE) {
  
  if (verbose) cat("Assessing bootstrap quality...\n")
  
  quality_assessment <- list()
  
  # Check if original values fall within confidence intervals
  quality_assessment$original_in_ci <- list()
  
  if ("ic" %in% names(confidence_intervals)) {
    total_ic_ci <- confidence_intervals$ic$total_ic
    quality_assessment$original_in_ci$total_ic <- 
      original_pwm$total_ic >= total_ic_ci[1] && original_pwm$total_ic <= total_ic_ci[2]
  }
  
  if ("conservation" %in% names(confidence_intervals)) {
    conservation_ci <- confidence_intervals$conservation$mean_conservation
    quality_assessment$original_in_ci$mean_conservation <- 
      original_pwm$mean_conservation >= conservation_ci[1] && original_pwm$mean_conservation <= conservation_ci[2]
  }
  
  # Check coefficient of variation (stability)
  quality_assessment$stability <- list()
  
  if ("ic" %in% names(confidence_intervals)) {
    total_ics <- sapply(bootstrap_results, function(x) x$total_ic)
    cv <- sd(total_ics, na.rm = TRUE) / mean(total_ics, na.rm = TRUE)
    quality_assessment$stability$ic_cv <- cv
    quality_assessment$stability$ic_stable <- cv < 0.1  # CV < 10% indicates good stability
  }
  
  # Overall quality grade
  ci_checks <- unlist(quality_assessment$original_in_ci)
  stability_checks <- unlist(quality_assessment$stability[grepl("stable", names(quality_assessment$stability))])
  
  all_checks <- c(ci_checks, stability_checks)
  pass_rate <- sum(all_checks, na.rm = TRUE) / length(all_checks)
  
  if (pass_rate >= 0.8) {
    quality_assessment$overall_grade <- "EXCELLENT"
  } else if (pass_rate >= 0.6) {
    quality_assessment$overall_grade <- "GOOD"
  } else if (pass_rate >= 0.4) {
    quality_assessment$overall_grade <- "ACCEPTABLE"
  } else {
    quality_assessment$overall_grade <- "POOR"
  }
  
  quality_assessment$pass_rate <- pass_rate
  
  return(quality_assessment)
}

#' Generate bootstrap validation report
generate_bootstrap_report <- function(validation_results, output_file) {
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>CTCF Bootstrap Validation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #f0f8ff; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }
        .excellent { background-color: #e8f5e8; }
        .good { background-color: #e8f0ff; }
        .acceptable { background-color: #fff8e8; }
        .poor { background-color: #ffeaea; }
        table { border-collapse: collapse; width: 100%; margin: 10px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .metric { margin: 10px 0; }
        .ci { font-weight: bold; color: #0066cc; }
        .original { font-weight: bold; color: #ff6600; }
    </style>
</head>
<body>
    <div class='header'>
        <h1>CTCF Bootstrap Validation Report</h1>
        <p>Generated: ", Sys.time(), "</p>
        <p>Bootstrap Samples: ", validation_results$n_total, " (", validation_results$n_successful, " successful)</p>
        <p>Success Rate: ", round(validation_results$success_rate * 100, 1), "%</p>
        <p>Confidence Level: ", validation_results$confidence_level * 100, "%</p>
        <p>Overall Quality: <strong>", validation_results$bootstrap_quality$overall_grade, "</strong></p>
    </div>")
  
  # Bootstrap statistics summary
  html_content <- paste0(html_content, "
    <div class='section ", tolower(validation_results$bootstrap_quality$overall_grade), "'>
        <h2>Bootstrap Statistics Summary</h2>")
  
  # Information content results
  if ("ic" %in% validation_results$metrics) {
    ic_stats <- validation_results$bootstrap_statistics$ic
    ic_ci <- validation_results$confidence_intervals$ic$total_ic
    
    html_content <- paste0(html_content, "
        <h3>Information Content</h3>
        <div class='metric'>Original Total IC: <span class='original'>", round(ic_stats$original_total_ic, 2), " bits</span></div>
        <div class='metric'>Bootstrap Mean: ", round(ic_stats$bootstrap_mean, 2), " ± ", round(ic_stats$bootstrap_sd, 2), " bits</div>
        <div class='metric'>Confidence Interval: <span class='ci'>[", round(ic_ci[1], 2), ", ", round(ic_ci[2], 2), "] bits</span></div>
        <div class='metric'>Coefficient of Variation: ", round(ic_stats$bootstrap_cv, 3), "</div>
        <div class='metric'>Bias: ", round(ic_stats$bias, 3), " bits</div>
        <div class='metric'>Bias-Corrected Estimate: ", round(ic_stats$bias_corrected, 2), " bits</div>")
  }
  
  # Conservation results
  if ("conservation" %in% validation_results$metrics) {
    cons_stats <- validation_results$bootstrap_statistics$conservation
    cons_ci <- validation_results$confidence_intervals$conservation$mean_conservation
    
    html_content <- paste0(html_content, "
        <h3>Conservation</h3>
        <div class='metric'>Original Mean Conservation: <span class='original'>", round(cons_stats$original_mean_conservation, 3), "</span></div>
        <div class='metric'>Bootstrap Mean: ", round(cons_stats$bootstrap_mean, 3), " ± ", round(cons_stats$bootstrap_sd, 3), "</div>
        <div class='metric'>Confidence Interval: <span class='ci'>[", round(cons_ci[1], 3), ", ", round(cons_ci[2], 3), "]</span></div>
        <div class='metric'>Coefficient of Variation: ", round(cons_stats$bootstrap_cv, 3), "</div>")
  }
  
  # Entropy results
  if ("entropy" %in% validation_results$metrics) {
    ent_stats <- validation_results$bootstrap_statistics$entropy
    ent_ci <- validation_results$confidence_intervals$entropy$mean_entropy
    
    html_content <- paste0(html_content, "
        <h3>Entropy</h3>
        <div class='metric'>Original Mean Entropy: <span class='original'>", round(ent_stats$original_mean_entropy, 3), "</span></div>
        <div class='metric'>Bootstrap Mean: ", round(ent_stats$bootstrap_mean, 3), " ± ", round(ent_stats$bootstrap_sd, 3), "</div>
        <div class='metric'>Confidence Interval: <span class='ci'>[", round(ent_ci[1], 3), ", ", round(ent_ci[2], 3), "]</span></div>")
  }
  
  html_content <- paste0(html_content, "
    </div>")
  
  # Quality assessment
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Quality Assessment</h2>
        <table>
            <tr><th>Check</th><th>Result</th><th>Details</th></tr>")
  
  # Original in CI checks
  for (check_name in names(validation_results$bootstrap_quality$original_in_ci)) {
    check_result <- validation_results$bootstrap_quality$original_in_ci[[check_name]]
    html_content <- paste0(html_content, "
            <tr>
                <td>", gsub("_", " ", check_name), " in CI</td>
                <td>", if (check_result) "PASS" else "FAIL", "</td>
                <td>Original value ", if (check_result) "within" else "outside", " confidence interval</td>
            </tr>")
  }
  
  # Stability checks
  for (check_name in names(validation_results$bootstrap_quality$stability)) {
    if (grepl("stable", check_name)) {
      check_result <- validation_results$bootstrap_quality$stability[[check_name]]
      html_content <- paste0(html_content, "
            <tr>
                <td>", gsub("_", " ", check_name), "</td>
                <td>", if (check_result) "PASS" else "FAIL", "</td>
                <td>Coefficient of variation ", if (check_result) "< 10%" else ">= 10%", "</td>
            </tr>")
    }
  }
  
  html_content <- paste0(html_content, "
        </table>
        <div class='metric'>Overall Pass Rate: ", round(validation_results$bootstrap_quality$pass_rate * 100, 1), "%</div>
    </div>")
  
  html_content <- paste0(html_content, "
</body>
</html>")
  
  writeLines(html_content, output_file)
}

# Main execution
validation_results <- bootstrap_validation(
  sequences_file = opt$sequences,
  n_bootstrap = opt$`bootstrap-samples`,
  confidence_level = opt$`confidence-level`,
  pseudocount = opt$pseudocount,
  metrics = metrics,
  n_cores = opt$`n-cores`,
  verbose = opt$verbose
)

if (opt$verbose) cat("Generating bootstrap validation report...\n")

generate_bootstrap_report(validation_results, opt$output)

if (opt$verbose) {
  cat("Bootstrap validation complete!\n")
  cat("Overall quality:", validation_results$bootstrap_quality$overall_grade, "\n")
  cat("Success rate:", round(validation_results$success_rate * 100, 1), "%\n")
  cat("Report saved to:", opt$output, "\n")
}

# Exit with appropriate status
if (validation_results$bootstrap_quality$overall_grade %in% c("EXCELLENT", "GOOD")) {
  quit(status = 0)
} else {
  quit(status = 1)
}
