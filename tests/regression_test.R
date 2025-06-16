#!/usr/bin/env Rscript

#' Regression Testing Script for CTCF PWM Testing Pipeline
#' 
#' This script performs regression testing by comparing current results against
#' known good reference results to detect any regressions in the pipeline.
#' 
#' Usage:
#'   Rscript tests/regression_test.R \
#'     --reference-results tests/reference_outputs/ \
#'     --current-results results/ \
#'     --tolerance 0.01

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

#' Load and parse PWM result files
#' @param file_path Path to PWM result file (RDS, JSON, or CSV)
#' @return Parsed PWM result
load_pwm_result <- function(file_path) {
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  ext <- tools::file_ext(file_path)
  
  tryCatch({
    switch(ext,
      "rds" = readRDS(file_path),
      "json" = fromJSON(file_path),
      "csv" = read.csv(file_path, stringsAsFactors = FALSE),
      {
        warning("Unsupported file format:", ext)
        NULL
      }
    )
  }, error = function(e) {
    warning("Error loading", file_path, ":", e$message)
    NULL
  })
}

#' Compare two PWM results
#' @param reference Reference PWM result
#' @param current Current PWM result
#' @param tolerance Numerical tolerance for comparisons
#' @return List with comparison results
compare_pwm_results <- function(reference, current, tolerance = 0.01) {
  if (is.null(reference) || is.null(current)) {
    return(list(
      status = "ERROR",
      message = "Missing reference or current result",
      passed = FALSE
    ))
  }
  
  results <- list(
    status = "PASS",
    message = "All tests passed",
    passed = TRUE,
    detailed_results = list()
  )
  
  # Test 1: Information Content Comparison
  if ("information_content" %in% names(reference) && 
      "information_content" %in% names(current)) {
    
    ic_ref <- reference$information_content
    ic_cur <- current$information_content
    
    if (length(ic_ref) != length(ic_cur)) {
      results$detailed_results$ic_length <- list(
        test = "Information Content Length",
        passed = FALSE,
        message = sprintf("Length mismatch: reference=%d, current=%d", 
                         length(ic_ref), length(ic_cur))
      )
      results$passed <- FALSE
    } else {
      # Compare IC values
      ic_diff <- abs(ic_ref - ic_cur)
      max_diff <- max(ic_diff, na.rm = TRUE)
      
      if (max_diff > tolerance) {
        results$detailed_results$ic_values <- list(
          test = "Information Content Values",
          passed = FALSE,
          message = sprintf("Max difference %.6f > tolerance %.6f", max_diff, tolerance),
          max_difference = max_diff
        )
        results$passed <- FALSE
      } else {
        results$detailed_results$ic_values <- list(
          test = "Information Content Values",
          passed = TRUE,
          message = sprintf("Max difference %.6f <= tolerance %.6f", max_diff, tolerance),
          max_difference = max_diff
        )
      }
    }
  }
  
  # Test 2: Total IC Comparison
  if ("total_ic" %in% names(reference) && "total_ic" %in% names(current)) {
    ic_diff <- abs(reference$total_ic - current$total_ic)
    
    if (ic_diff > tolerance) {
      results$detailed_results$total_ic <- list(
        test = "Total Information Content",
        passed = FALSE,
        message = sprintf("Difference %.6f > tolerance %.6f", ic_diff, tolerance),
        reference_value = reference$total_ic,
        current_value = current$total_ic,
        difference = ic_diff
      )
      results$passed <- FALSE
    } else {
      results$detailed_results$total_ic <- list(
        test = "Total Information Content",
        passed = TRUE,
        message = sprintf("Difference %.6f <= tolerance %.6f", ic_diff, tolerance),
        reference_value = reference$total_ic,
        current_value = current$total_ic,
        difference = ic_diff
      )
    }
  }
  
  # Test 3: PWM Matrix Comparison
  if ("pwm" %in% names(reference) && "pwm" %in% names(current)) {
    pwm_ref <- reference$pwm
    pwm_cur <- current$pwm
    
    if (is.matrix(pwm_ref) && is.matrix(pwm_cur)) {
      if (!identical(dim(pwm_ref), dim(pwm_cur))) {
        results$detailed_results$pwm_dimensions <- list(
          test = "PWM Matrix Dimensions",
          passed = FALSE,
          message = sprintf("Dimension mismatch: reference=%s, current=%s",
                           paste(dim(pwm_ref), collapse="x"),
                           paste(dim(pwm_cur), collapse="x"))
        )
        results$passed <- FALSE
      } else {
        # Compare matrix values
        pwm_diff <- abs(pwm_ref - pwm_cur)
        max_diff <- max(pwm_diff, na.rm = TRUE)
        
        if (max_diff > tolerance) {
          results$detailed_results$pwm_values <- list(
            test = "PWM Matrix Values",
            passed = FALSE,
            message = sprintf("Max difference %.6f > tolerance %.6f", max_diff, tolerance),
            max_difference = max_diff
          )
          results$passed <- FALSE
        } else {
          results$detailed_results$pwm_values <- list(
            test = "PWM Matrix Values",
            passed = TRUE,
            message = sprintf("Max difference %.6f <= tolerance %.6f", max_diff, tolerance),
            max_difference = max_diff
          )
        }
      }
    }
  }
  
  # Test 4: Quality Grade Comparison
  if ("quality_grade" %in% names(reference) && "quality_grade" %in% names(current)) {
    if (reference$quality_grade != current$quality_grade) {
      results$detailed_results$quality_grade <- list(
        test = "Quality Grade",
        passed = FALSE,
        message = sprintf("Grade mismatch: reference='%s', current='%s'",
                         reference$quality_grade, current$quality_grade),
        reference_value = reference$quality_grade,
        current_value = current$quality_grade
      )
      results$passed <- FALSE
    } else {
      results$detailed_results$quality_grade <- list(
        test = "Quality Grade",
        passed = TRUE,
        message = "Quality grades match",
        value = reference$quality_grade
      )
    }
  }
  
  # Update overall status
  if (!results$passed) {
    results$status <- "FAIL"
    results$message <- "Some tests failed"
  }
  
  return(results)
}

#' Find matching reference file for a current result file
#' @param current_file Current result file path
#' @param reference_dir Reference directory
#' @return Path to matching reference file or NULL
find_reference_file <- function(current_file, reference_dir) {
  basename_file <- basename(current_file)
  
  # Try exact match first
  ref_file <- file.path(reference_dir, basename_file)
  if (file.exists(ref_file)) {
    return(ref_file)
  }
  
  # Try with different extensions
  base_name <- tools::file_path_sans_ext(basename_file)
  extensions <- c("rds", "json", "csv")
  
  for (ext in extensions) {
    ref_file <- file.path(reference_dir, paste0(base_name, ".", ext))
    if (file.exists(ref_file)) {
      return(ref_file)
    }
  }
  
  return(NULL)
}

#' Run regression tests
#' @param reference_dir Directory containing reference results
#' @param current_dir Directory containing current results
#' @param tolerance Numerical tolerance for comparisons
#' @return List with test results
run_regression_tests <- function(reference_dir, current_dir, tolerance = 0.01) {
  
  if (!dir.exists(reference_dir)) {
    stop("Reference directory does not exist: ", reference_dir)
  }
  
  if (!dir.exists(current_dir)) {
    stop("Current results directory does not exist: ", current_dir)
  }
  
  # Find all result files in current directory
  result_files <- list.files(current_dir, 
                            pattern = "\\.(rds|json|csv)$", 
                            full.names = TRUE, 
                            recursive = TRUE)
  
  if (length(result_files) == 0) {
    warning("No result files found in current directory: ", current_dir)
    return(list(
      overall_status = "ERROR",
      message = "No result files found",
      total_tests = 0,
      passed_tests = 0,
      failed_tests = 0,
      test_results = list()
    ))
  }
  
  cat("Found", length(result_files), "result files to test\n")
  
  test_results <- list()
  passed_count <- 0
  failed_count <- 0
  
  for (current_file in result_files) {
    cat("Testing:", basename(current_file), "... ")
    
    # Find matching reference file
    reference_file <- find_reference_file(current_file, reference_dir)
    
    if (is.null(reference_file)) {
      cat("SKIP (no reference)\n")
      test_results[[basename(current_file)]] <- list(
        status = "SKIP",
        message = "No matching reference file found",
        passed = NA
      )
      next
    }
    
    # Load results
    reference_result <- load_pwm_result(reference_file)
    current_result <- load_pwm_result(current_file)
    
    # Compare results
    comparison <- compare_pwm_results(reference_result, current_result, tolerance)
    test_results[[basename(current_file)]] <- comparison
    
    if (comparison$passed) {
      cat("PASS\n")
      passed_count <- passed_count + 1
    } else {
      cat("FAIL\n")
      failed_count <- failed_count + 1
    }
  }
  
  # Overall status
  total_tests <- passed_count + failed_count
  overall_status <- if (failed_count == 0) "PASS" else "FAIL"
  
  return(list(
    overall_status = overall_status,
    message = sprintf("Passed: %d, Failed: %d, Total: %d", 
                     passed_count, failed_count, total_tests),
    total_tests = total_tests,
    passed_tests = passed_count,
    failed_tests = failed_count,
    tolerance = tolerance,
    test_results = test_results
  ))
}

#' Generate regression test report
#' @param test_results Test results from run_regression_tests
#' @param output_file Output report file path
generate_regression_report <- function(test_results, output_file = NULL) {
  
  # Generate text report
  report_lines <- c(
    "CTCF PWM Pipeline Regression Test Report",
    paste("Generated:", Sys.time()),
    "",
    "=== SUMMARY ===",
    paste("Overall Status:", test_results$overall_status),
    paste("Total Tests:", test_results$total_tests),
    paste("Passed:", test_results$passed_tests),
    paste("Failed:", test_results$failed_tests),
    paste("Tolerance:", test_results$tolerance),
    "",
    "=== DETAILED RESULTS ==="
  )
  
  for (test_name in names(test_results$test_results)) {
    result <- test_results$test_results[[test_name]]
    
    report_lines <- c(
      report_lines,
      "",
      paste("Test:", test_name),
      paste("Status:", result$status),
      paste("Message:", result$message)
    )
    
    if ("detailed_results" %in% names(result)) {
      for (detail_name in names(result$detailed_results)) {
        detail <- result$detailed_results[[detail_name]]
        report_lines <- c(
          report_lines,
          paste("  -", detail$test, ":", 
                if(detail$passed) "PASS" else "FAIL"),
          paste("   ", detail$message)
        )
      }
    }
  }
  
  # Print to console
  cat(paste(report_lines, collapse = "\n"), "\n")
  
  # Save to file if specified
  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    writeLines(report_lines, output_file)
    cat("Report saved to:", output_file, "\n")
  }
  
  return(invisible(report_lines))
}

#' Main execution function
main <- function() {
  # Parse command line arguments
  option_list <- list(
    make_option(c("--reference-results"), type = "character", 
                default = "tests/reference_outputs/",
                help = "Directory containing reference results [default: %default]"),
    make_option(c("--current-results"), type = "character", 
                default = "results/",
                help = "Directory containing current results [default: %default]"),
    make_option(c("--tolerance"), type = "double", default = 0.01,
                help = "Numerical tolerance for comparisons [default: %default]"),
    make_option(c("--output-report"), type = "character", default = NULL,
                help = "Output report file path [default: console only]"),
    make_option(c("--help"), action = "store_true", default = FALSE,
                help = "Show this help message and exit")
  )
  
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  opt <- parse_args(parser)
  
  if (opt$help) {
    print_help(parser)
    return()
  }
  
  # Validate parameters
  if (opt$tolerance <= 0) {
    stop("Tolerance must be positive")
  }
  
  # Run regression tests
  cat("Starting regression tests...\n")
  cat("Reference directory:", opt$`reference-results`, "\n")
  cat("Current directory:", opt$`current-results`, "\n")
  cat("Tolerance:", opt$tolerance, "\n\n")
  
  test_results <- run_regression_tests(
    opt$`reference-results`, 
    opt$`current-results`, 
    opt$tolerance
  )
  
  # Generate report
  generate_regression_report(test_results, opt$`output-report`)
  
  # Exit with appropriate code
  if (test_results$overall_status == "PASS") {
    cat("\n✓ All regression tests passed!\n")
    quit(status = 0)
  } else {
    cat("\n✗ Some regression tests failed!\n")
    quit(status = 1)
  }
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
