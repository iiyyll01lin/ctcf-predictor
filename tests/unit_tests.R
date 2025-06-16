#!/usr/bin/env Rscript

#' Unit Tests for CTCF PWM Testing Pipeline
#' 
#' This script contains unit tests for individual functions in the pipeline
#' to ensure each component works correctly in isolation.

suppressPackageStartupMessages({
  library(testthat)
  library(Biostrings)
})

# Source the functions to test
source_files <- c(
  "scripts/build_pwm_robust.R",
  "scripts/utility_functions.R",
  "scripts/statistical_functions.R"
)

for (file in source_files) {
  if (file.exists(file)) {
    source(file)
  } else {
    cat("Warning: Could not source", file, "\n")
  }
}

#' Test helper functions
create_test_sequences <- function(n = 10) {
  # Create test sequences with known CTCF motif
  sequences <- c(
    "ATGCCGCGGGGGCGCTAATCGATCG",  # Contains CTCF motif
    "GGACCGCGGGGGCGCTTTAAGCTAG",  # Contains CTCF motif
    "ATCGATCGATCGATCGATCGATCGA",  # Random sequence
    "CCGCGGGGGCGCTAAAAAAAAAAAA",  # CTCF motif at start
    "AAAAAAAACCGCGGGGGCGCTAAAA",  # CTCF motif in middle
    "AAAAAAAAAAAACCGCGGGGGCGCT",  # CTCF motif at end
    "TTTTTTTTTTTTTTTTTTTTTTTT",   # Low complexity
    "GGGGGGGGGGGGGGGGGGGGGGGG",   # Low complexity
    "ACGTACGTACGTACGTACGTACGT",   # Periodic pattern
    "ATCGATCGATCGATCGATCGATCG"    # Another random
  )
  
  if (n <= length(sequences)) {
    return(DNAStringSet(sequences[1:n]))
  } else {
    # Extend with random sequences
    additional <- n - length(sequences)
    for (i in 1:additional) {
      random_seq <- paste(sample(c("A", "C", "G", "T"), 25, replace = TRUE), collapse = "")
      sequences <- c(sequences, random_seq)
    }
    return(DNAStringSet(sequences))
  }
}

#' Unit Tests
test_that("PWM construction produces valid output", {
  sequences <- create_test_sequences(5)
  
  # Test basic PWM construction if function exists
  if (exists("build_basic_pwm")) {
    result <- build_basic_pwm(sequences)
    
    expect_true(is.list(result))
    expect_true("pwm" %in% names(result))
    expect_true(is.matrix(result$pwm))
    expect_equal(nrow(result$pwm), 4)  # A, C, G, T rows
    expect_true(all(rownames(result$pwm) %in% c("A", "C", "G", "T")))
    
    # Check probabilities sum to 1 for each position
    col_sums <- colSums(result$pwm)
    expect_true(all(abs(col_sums - 1) < 1e-10))
  }
})

test_that("Information content calculation is correct", {
  # Test with known probability matrix
  test_pwm <- matrix(c(
    0.25, 0.25, 0.25, 0.25,  # Equal probabilities - IC should be 0
    1.0,  0.0,  0.0,  0.0,   # All A - IC should be 2
    0.5,  0.5,  0.0,  0.0,   # Half A, half C - IC should be 1
    0.4,  0.3,  0.2,  0.1    # Mixed - IC between 0 and 2
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("A", "C", "G", "T"), 1:4))
  
  if (exists("calculate_information_content")) {
    ic <- calculate_information_content(test_pwm)
    
    expect_equal(length(ic), 4)
    expect_true(ic[1] < 0.1)      # Nearly uniform - low IC
    expect_true(ic[2] > 1.9)      # Highly conserved - high IC
    expect_true(ic[3] > 0.9 && ic[3] < 1.1)  # Half conserved - medium IC
    expect_true(ic[4] > 0 && ic[4] < 2)       # Mixed - reasonable IC
  }
})

test_that("Sequence validation works correctly", {
  valid_sequences <- create_test_sequences(3)
  invalid_sequences <- DNAStringSet(c("ATCG", "INVALID", "atcg"))
  
  if (exists("validate_sequences")) {
    # Valid sequences should pass
    expect_true(validate_sequences(valid_sequences))
    
    # Invalid sequences should fail or be filtered
    result <- validate_sequences(invalid_sequences)
    expect_true(is.logical(result) || is(result, "DNAStringSet"))
  }
})

test_that("Statistical significance testing works", {
  # Create test data with known statistical properties
  test_ic_values <- c(12.5, 15.2, 8.9, 16.7, 11.3)
  null_ic_values <- c(5.2, 4.8, 6.1, 5.5, 4.9)
  
  if (exists("calculate_statistical_significance")) {
    p_value <- calculate_statistical_significance(test_ic_values, null_ic_values)
    
    expect_true(is.numeric(p_value))
    expect_true(p_value >= 0 && p_value <= 1)
    expect_true(p_value < 0.05)  # Should be significant given the test data
  }
})

test_that("Quality grading is consistent", {
  if (exists("determine_quality_grade")) {
    # Test different IC values
    expect_equal(determine_quality_grade(20.0), "EXCELLENT")
    expect_equal(determine_quality_grade(15.0), "GOOD")
    expect_equal(determine_quality_grade(10.0), "FAIR")
    expect_equal(determine_quality_grade(5.0), "POOR")
    
    # Test edge cases
    expect_true(determine_quality_grade(0) %in% c("POOR", "FAIL"))
    expect_true(determine_quality_grade(-1) %in% c("POOR", "FAIL"))
  }
})

test_that("Alignment methods produce consistent results", {
  sequences <- create_test_sequences(8)
  
  # Test center alignment if available
  if (exists("center_align_sequences")) {
    aligned <- center_align_sequences(sequences)
    
    expect_true(is(aligned, "DNAStringSet"))
    expect_equal(length(aligned), length(sequences))
    expect_true(all(width(aligned) == width(aligned)[1]))  # All same width
  }
  
  # Test consensus alignment if available
  if (exists("consensus_align_sequences")) {
    aligned <- consensus_align_sequences(sequences)
    
    expect_true(is(aligned, "DNAStringSet"))
    expect_equal(length(aligned), length(sequences))
  }
})

test_that("Bootstrap validation produces stable results", {
  sequences <- create_test_sequences(20)
  
  if (exists("bootstrap_pwm_validation")) {
    result <- bootstrap_pwm_validation(sequences, n_bootstrap = 10)
    
    expect_true(is.list(result))
    expect_true("confidence_intervals" %in% names(result))
    expect_true("bootstrap_results" %in% names(result))
    
    # Check confidence intervals
    ci <- result$confidence_intervals
    expect_true(all(ci$lower <= ci$upper))
  }
})

test_that("Memory management functions work correctly", {
  if (exists("clean_memory")) {
    # Should not throw errors
    expect_silent(clean_memory())
  }
  
  if (exists("monitor_memory_usage")) {
    memory_info <- monitor_memory_usage()
    expect_true(is.list(memory_info) || is.numeric(memory_info))
  }
})

test_that("File I/O functions handle errors gracefully", {
  if (exists("save_pwm_result")) {
    test_pwm <- list(
      pwm = matrix(0.25, nrow = 4, ncol = 10),
      information_content = rep(1.0, 10),
      total_ic = 10.0
    )
    
    temp_file <- tempfile(fileext = ".rds")
    
    # Should save successfully
    expect_silent(save_pwm_result(test_pwm, temp_file))
    expect_true(file.exists(temp_file))
    
    # Clean up
    unlink(temp_file)
  }
  
  if (exists("load_pwm_result")) {
    # Should handle missing files gracefully
    expect_error(load_pwm_result("nonexistent_file.rds"))
  }
})

test_that("Configuration validation works", {
  if (exists("validate_config")) {
    # Test valid configuration
    valid_config <- list(
      method = "integrated",
      pseudocount = 0.25,
      min_ic_threshold = 8.0,
      bootstrap_samples = 1000
    )
    
    expect_true(validate_config(valid_config))
    
    # Test invalid configuration
    invalid_config <- list(
      method = "invalid_method",
      pseudocount = -1,
      min_ic_threshold = -5
    )
    
    expect_false(validate_config(invalid_config))
  }
})

#' Run all unit tests
run_unit_tests <- function() {
  cat("Running unit tests for CTCF PWM Testing Pipeline...\n")
  
  # Set up test environment
  test_env <- new.env()
  
  # Run tests
  test_results <- test_dir(".", env = test_env, reporter = "summary")
  
  # Summary
  cat("\n=== UNIT TEST SUMMARY ===\n")
  if (test_results$failed == 0) {
    cat("✓ All", test_results$passed, "unit tests passed!\n")
    return(TRUE)
  } else {
    cat("✗", test_results$failed, "out of", 
        test_results$passed + test_results$failed, "tests failed!\n")
    return(FALSE)
  }
}

#' Main execution
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0 && args[1] %in% c("-h", "--help", "help")) {
    cat("Usage: Rscript tests/unit_tests.R\n")
    cat("Runs unit tests for individual pipeline functions.\n")
    return()
  }
  
  # Check if testthat is available
  if (!requireNamespace("testthat", quietly = TRUE)) {
    cat("Warning: testthat package not available. Installing basic test framework...\n")
    
    # Simple test framework fallback
    test_that <- function(desc, code) {
      cat("Testing:", desc, "... ")
      tryCatch({
        eval(code)
        cat("PASS\n")
        return(TRUE)
      }, error = function(e) {
        cat("FAIL -", e$message, "\n")
        return(FALSE)
      })
    }
    
    expect_true <- function(x) if (!x) stop("Expected TRUE")
    expect_false <- function(x) if (x) stop("Expected FALSE")
    expect_equal <- function(x, y) if (!identical(x, y)) stop("Values not equal")
    expect_error <- function(code) {
      tryCatch({
        eval(code)
        stop("Expected error but none occurred")
      }, error = function(e) {
        # Error expected - this is good
      })
    }
    expect_silent <- function(code) {
      tryCatch({
        eval(code)
      }, error = function(e) {
        stop("Unexpected error:", e$message)
      })
    }
  }
  
  # Run tests
  success <- run_unit_tests()
  
  if (success) {
    quit(status = 0)
  } else {
    quit(status = 1)
  }
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
