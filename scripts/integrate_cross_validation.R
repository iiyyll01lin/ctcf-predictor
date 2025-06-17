#!/usr/bin/env Rscript

# Cross-Validation Integration Script
# Integrates existing cross-validation scripts into main pipeline
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(optparse)
  library(parallel)
  library(jsonlite)
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input sequences file (FASTA format)", metavar="file"),
  make_option(c("-o", "--output-dir"), type="character", default="results/validation_results",
              help="Output directory for validation results [default: %default]"),
  make_option(c("-m", "--methods"), type="character", default="bootstrap,chromosome,kfold",
              help="Validation methods (comma-separated) [default: %default]"),
  make_option(c("-k", "--k-folds"), type="integer", default=5,
              help="Number of folds for k-fold cross-validation [default: %default]"),
  make_option(c("-b", "--bootstrap-reps"), type="integer", default=1000,
              help="Number of bootstrap replicates [default: %default]"),
  make_option(c("-c", "--confidence-level"), type="numeric", default=0.95,
              help="Confidence level for intervals [default: %default]"),
  make_option(c("--test-chromosomes"), type="character", default="chr21,chr22",
              help="Test chromosomes for chromosome split [default: %default]"),
  make_option(c("-n", "--n-cores"), type="integer", default=1,
              help="Number of cores for parallel processing [default: %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Integrate cross-validation methods for PWM validation")
opt <- parse_args(opt_parser)

# Source required functions
source_validation_scripts <- function(verbose = FALSE) {
  if (verbose) cat("Sourcing validation scripts...\n")
  
  required_scripts <- c(
    "scripts/bootstrap_validation.R",
    "scripts/validate_chromosome_split.R",
    "scripts/chromosome_cross_validation.R"
  )
  
  for (script in required_scripts) {
    if (file.exists(script)) {
      source(script)
      if (verbose) cat("  Sourced:", script, "\n")
    } else {
      warning("Script not found:", script)
    }
  }
}

#' Run bootstrap validation
#' @param sequences_file Path to sequences file
#' @param output_dir Output directory
#' @param n_bootstrap Number of bootstrap replicates
#' @param confidence_level Confidence level
#' @param verbose Enable verbose output
run_bootstrap_validation <- function(sequences_file, output_dir, n_bootstrap = 1000,
                                   confidence_level = 0.95, verbose = FALSE) {
  if (verbose) cat("Running bootstrap validation...\n")
  
  # Create bootstrap output directory
  bootstrap_dir <- file.path(output_dir, "bootstrap")
  if (!dir.exists(bootstrap_dir)) {
    dir.create(bootstrap_dir, recursive = TRUE)
  }
  
  # Run bootstrap validation using existing script
  bootstrap_output <- file.path(bootstrap_dir, "bootstrap_results.json")
  
  # Check if bootstrap_validation.R has the required function
  if (exists("bootstrap_confidence_interval")) {
    # Load sequences
    sequences <- readDNAStringSet(sequences_file)
    
    # Run bootstrap analysis
    alpha <- 1 - confidence_level
    bootstrap_results <- bootstrap_confidence_interval(
      sequences, n_bootstrap = n_bootstrap, alpha = alpha
    )
    
    # Save results
    writeLines(toJSON(bootstrap_results, pretty = TRUE), bootstrap_output)
    
    if (verbose) {
      cat("  Bootstrap validation complete\n")
      cat("  Results saved to:", bootstrap_output, "\n")
    }
    
    return(bootstrap_results)
  } else {
    warning("bootstrap_confidence_interval function not available")
    return(NULL)
  }
}

#' Run chromosome split validation
#' @param sequences_file Path to sequences file
#' @param output_dir Output directory
#' @param test_chromosomes Vector of test chromosomes
#' @param verbose Enable verbose output
run_chromosome_split_validation <- function(sequences_file, output_dir, 
                                          test_chromosomes = c("chr21", "chr22"),
                                          verbose = FALSE) {
  if (verbose) cat("Running chromosome split validation...\n")
  
  # Create chromosome output directory
  chromosome_dir <- file.path(output_dir, "chromosome_split")
  if (!dir.exists(chromosome_dir)) {
    dir.create(chromosome_dir, recursive = TRUE)
  }
  
  # Run chromosome split validation
  chromosome_output <- file.path(chromosome_dir, "chromosome_split_results.json")
  
  # Use existing chromosome split validation function if available
  if (exists("chromosome_split_validation")) {
    chromosome_results <- chromosome_split_validation(
      input_file = sequences_file,
      output_dir = chromosome_dir,
      test_chromosomes = paste(test_chromosomes, collapse = ","),
      verbose = verbose
    )
    
    if (verbose) {
      cat("  Chromosome split validation complete\n")
      cat("  Results saved to:", chromosome_output, "\n")
    }
    
    return(chromosome_results)
  } else {
    warning("chromosome_split_validation function not available")
    return(NULL)
  }
}

#' Run k-fold cross-validation
#' @param sequences_file Path to sequences file
#' @param output_dir Output directory
#' @param k_folds Number of folds
#' @param verbose Enable verbose output
run_kfold_validation <- function(sequences_file, output_dir, k_folds = 5,
                                verbose = FALSE) {
  if (verbose) cat("Running k-fold cross-validation...\n")
  
  # Create k-fold output directory
  kfold_dir <- file.path(output_dir, "kfold")
  if (!dir.exists(kfold_dir)) {
    dir.create(kfold_dir, recursive = TRUE)
  }
  
  # Load sequences
  sequences <- readDNAStringSet(sequences_file)
  
  # Extract chromosome information from sequence names
  seq_names <- names(sequences)
  chromosomes <- sapply(seq_names, function(name) {
    # Try to extract chromosome from sequence name
    chr_match <- regmatches(name, regexpr("chr[0-9XY]+", name))
    if (length(chr_match) > 0) {
      return(chr_match)
    } else {
      return(paste0("seq_", sample(1:k_folds, 1)))  # Random assignment if no chr info
    }
  })
  
  # Run chromosome-based cross-validation if function exists
  if (exists("chromosome_cross_validation")) {
    kfold_results <- chromosome_cross_validation(sequences, chromosomes, k = k_folds)
    
    # Save results
    kfold_output <- file.path(kfold_dir, "kfold_results.json")
    writeLines(toJSON(kfold_results, pretty = TRUE), kfold_output)
    
    if (verbose) {
      cat("  K-fold cross-validation complete\n")
      cat("  Results saved to:", kfold_output, "\n")
    }
    
    return(kfold_results)
  } else {
    warning("chromosome_cross_validation function not available")
    return(NULL)
  }
}

#' Generate comprehensive validation summary
#' @param validation_results List of validation results
#' @param output_dir Output directory
#' @param verbose Enable verbose output
generate_validation_summary <- function(validation_results, output_dir, verbose = FALSE) {
  if (verbose) cat("Generating validation summary...\n")
  
  summary_data <- list(
    timestamp = Sys.time(),
    methods_run = names(validation_results)[!sapply(validation_results, is.null)],
    results = validation_results
  )
  
  # Calculate summary statistics
  summary_stats <- list()
  
  # Bootstrap summary
  if (!is.null(validation_results$bootstrap)) {
    bootstrap_res <- validation_results$bootstrap
    summary_stats$bootstrap <- list(
      mean_ic = bootstrap_res$mean_ic %||% NA,
      confidence_interval = bootstrap_res$confidence_interval %||% c(NA, NA),
      cv_coefficient = bootstrap_res$cv_coefficient %||% NA
    )
  }
  
  # Chromosome split summary
  if (!is.null(validation_results$chromosome_split)) {
    chr_res <- validation_results$chromosome_split
    summary_stats$chromosome_split <- list(
      training_performance = chr_res$training_performance %||% NA,
      testing_performance = chr_res$testing_performance %||% NA,
      generalization_gap = chr_res$generalization_gap %||% NA
    )
  }
  
  # K-fold summary
  if (!is.null(validation_results$kfold)) {
    kfold_res <- validation_results$kfold
    summary_stats$kfold <- list(
      mean_performance = kfold_res$mean_performance %||% NA,
      std_deviation = kfold_res$std_deviation %||% NA,
      consistency_grade = kfold_res$consistency_grade %||% "UNKNOWN"
    )
  }
  
  summary_data$summary_statistics <- summary_stats
  
  # Overall validation grade
  validation_grade <- determine_validation_grade(summary_stats)
  summary_data$overall_grade <- validation_grade
  
  # Save comprehensive summary
  summary_file <- file.path(output_dir, "validation_summary.json")
  writeLines(toJSON(summary_data, pretty = TRUE), summary_file)
  
  # Generate HTML report
  html_report <- generate_validation_html_report(summary_data, output_dir)
  
  if (verbose) {
    cat("  Validation summary generated\n")
    cat("  JSON summary:", summary_file, "\n")
    cat("  HTML report:", html_report, "\n")
    cat("  Overall validation grade:", validation_grade, "\n")
  }
  
  return(summary_data)
}

#' Determine overall validation grade
#' @param summary_stats Summary statistics
#' @return Validation grade
determine_validation_grade <- function(summary_stats) {
  grades <- c()
  
  # Bootstrap grade
  if (!is.null(summary_stats$bootstrap)) {
    cv_coef <- summary_stats$bootstrap$cv_coefficient
    if (!is.null(cv_coef) && !is.na(cv_coef)) {
      if (cv_coef < 0.05) {
        grades <- c(grades, "EXCELLENT")
      } else if (cv_coef < 0.10) {
        grades <- c(grades, "GOOD")
      } else {
        grades <- c(grades, "ACCEPTABLE")
      }
    }
  }
  
  # Chromosome split grade
  if (!is.null(summary_stats$chromosome_split)) {
    gap <- summary_stats$chromosome_split$generalization_gap
    if (!is.null(gap) && !is.na(gap)) {
      if (gap < 1.0) {
        grades <- c(grades, "EXCELLENT")
      } else if (gap < 2.0) {
        grades <- c(grades, "GOOD")
      } else {
        grades <- c(grades, "ACCEPTABLE")
      }
    }
  }
  
  # K-fold grade
  if (!is.null(summary_stats$kfold)) {
    consistency <- summary_stats$kfold$consistency_grade
    if (!is.null(consistency) && consistency != "UNKNOWN") {
      grades <- c(grades, consistency)
    }
  }
  
  # Determine overall grade
  if (length(grades) == 0) {
    return("UNKNOWN")
  }
  
  grade_scores <- sapply(grades, function(g) {
    switch(g,
           "EXCELLENT" = 4,
           "GOOD" = 3,
           "ACCEPTABLE" = 2,
           "POOR" = 1,
           0)
  })
  
  avg_score <- mean(grade_scores)
  
  if (avg_score >= 3.5) {
    return("EXCELLENT")
  } else if (avg_score >= 2.5) {
    return("GOOD")
  } else if (avg_score >= 1.5) {
    return("ACCEPTABLE")
  } else {
    return("POOR")
  }
}

#' Generate HTML validation report
#' @param summary_data Summary data
#' @param output_dir Output directory
#' @return HTML report file path
generate_validation_html_report <- function(summary_data, output_dir) {
  html_file <- file.path(output_dir, "validation_report.html")
  
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>PWM Cross-Validation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        h2 { color: #3498db; }
        .grade { font-size: 1.2em; font-weight: bold; padding: 10px; border-radius: 5px; }
        .excellent { background-color: #d5f4e6; color: #27ae60; }
        .good { background-color: #fef9e7; color: #f39c12; }
        .acceptable { background-color: #ebf3fd; color: #3498db; }
        .poor { background-color: #fadbd8; color: #e74c3c; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h1>PWM Cross-Validation Report</h1>
    <p><strong>Generated:</strong> ", format(summary_data$timestamp, "%Y-%m-%d %H:%M:%S"), "</p>
    
    <h2>Overall Validation Grade</h2>
    <div class='grade ", tolower(summary_data$overall_grade), "'>", 
    summary_data$overall_grade, "</div>
    
    <h2>Methods Run</h2>
    <ul>")
  
  for (method in summary_data$methods_run) {
    html_content <- paste0(html_content, "<li>", method, "</li>")
  }
  
  html_content <- paste0(html_content, "
    </ul>
    
    <h2>Summary Statistics</h2>
    <table>
        <tr><th>Method</th><th>Metric</th><th>Value</th></tr>")
  
  # Add bootstrap statistics
  if (!is.null(summary_data$summary_statistics$bootstrap)) {
    bs <- summary_data$summary_statistics$bootstrap
    html_content <- paste0(html_content, "
        <tr><td rowspan='3'>Bootstrap</td><td>Mean IC</td><td>", 
        round(bs$mean_ic %||% 0, 3), " bits</td></tr>
        <tr><td>95% CI</td><td>[", 
        round(bs$confidence_interval[1] %||% 0, 3), ", ",
        round(bs$confidence_interval[2] %||% 0, 3), "]</td></tr>
        <tr><td>CV Coefficient</td><td>", 
        round(bs$cv_coefficient %||% 0, 3), "</td></tr>")
  }
  
  # Add chromosome split statistics
  if (!is.null(summary_data$summary_statistics$chromosome_split)) {
    cs <- summary_data$summary_statistics$chromosome_split
    html_content <- paste0(html_content, "
        <tr><td rowspan='3'>Chromosome Split</td><td>Training Performance</td><td>", 
        round(cs$training_performance %||% 0, 3), " bits</td></tr>
        <tr><td>Testing Performance</td><td>", 
        round(cs$testing_performance %||% 0, 3), " bits</td></tr>
        <tr><td>Generalization Gap</td><td>", 
        round(cs$generalization_gap %||% 0, 3), " bits</td></tr>")
  }
  
  # Add k-fold statistics
  if (!is.null(summary_data$summary_statistics$kfold)) {
    kf <- summary_data$summary_statistics$kfold
    html_content <- paste0(html_content, "
        <tr><td rowspan='3'>K-Fold CV</td><td>Mean Performance</td><td>", 
        round(kf$mean_performance %||% 0, 3), " bits</td></tr>
        <tr><td>Standard Deviation</td><td>", 
        round(kf$std_deviation %||% 0, 3), "</td></tr>
        <tr><td>Consistency Grade</td><td>", 
        kf$consistency_grade %||% "UNKNOWN", "</td></tr>")
  }
  
  html_content <- paste0(html_content, "
    </table>
</body>
</html>")
  
  writeLines(html_content, html_file)
  return(html_file)
}

# Helper function for null-coalescing operator
`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

# Main execution
if (!interactive()) {
  # Validate inputs
  if (is.null(opt$input)) {
    cat("Error: Input sequences file is required\n")
    quit(status = 1)
  }
  
  if (!file.exists(opt$input)) {
    cat("Error: Input file does not exist:", opt$input, "\n")
    quit(status = 1)
  }
  
  if (opt$verbose) {
    cat("Cross-Validation Integration\n")
    cat("===========================\n")
    cat("Input file:", opt$input, "\n")
    cat("Output directory:", opt$`output-dir`, "\n")
    cat("Methods:", opt$methods, "\n\n")
  }
  
  # Create output directory
  if (!dir.exists(opt$`output-dir`)) {
    dir.create(opt$`output-dir`, recursive = TRUE)
    if (opt$verbose) cat("Created output directory:", opt$`output-dir`, "\n")
  }
  
  # Source validation scripts
  source_validation_scripts(opt$verbose)
  
  # Parse methods
  methods <- trimws(strsplit(opt$methods, ",")[[1]])
  validation_results <- list()
  
  # Run requested validation methods
  for (method in methods) {
    if (method == "bootstrap") {
      validation_results$bootstrap <- run_bootstrap_validation(
        opt$input, opt$`output-dir`, opt$`bootstrap-reps`, 
        opt$`confidence-level`, opt$verbose
      )
    } else if (method == "chromosome") {
      test_chrs <- trimws(strsplit(opt$`test-chromosomes`, ",")[[1]])
      validation_results$chromosome_split <- run_chromosome_split_validation(
        opt$input, opt$`output-dir`, test_chrs, opt$verbose
      )
    } else if (method == "kfold") {
      validation_results$kfold <- run_kfold_validation(
        opt$input, opt$`output-dir`, opt$`k-folds`, opt$verbose
      )
    } else {
      warning("Unknown validation method:", method)
    }
  }
  
  # Generate comprehensive summary
  summary_data <- generate_validation_summary(validation_results, opt$`output-dir`, opt$verbose)
  
  if (opt$verbose) {
    cat("\nCross-validation integration complete!\n")
    cat("Overall validation grade:", summary_data$overall_grade, "\n")
    cat("Methods completed:", length(validation_results), "/", length(methods), "\n")
  }
}
