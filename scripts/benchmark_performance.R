#!/usr/bin/env Rscript

#' Performance Benchmarking Script for CTCF PWM Testing Pipeline
#' 
#' This script performs comprehensive performance benchmarking of the PWM construction
#' pipeline across different dataset sizes, methods, and system configurations.
#' 
#' Usage:
#'   Rscript scripts/benchmark_performance.R \
#'     --dataset-sizes "100,500,1000,5000,10000" \
#'     --methods "center,consensus,integrated" \
#'     --repetitions 10 \
#'     --output results/performance_benchmark.html

suppressPackageStartupMessages({
  library(optparse)
  library(microbenchmark)
  library(Biostrings)
  library(pryr)
  library(ggplot2)
  library(dplyr)
  library(knitr)
  library(rmarkdown)
})

# Source required functions
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

#' Generate synthetic sequences for benchmarking
#' @param n_sequences Number of sequences to generate
#' @param seq_length Length of each sequence
#' @param ctcf_motif Include CTCF-like motif if TRUE
#' @return DNAStringSet of synthetic sequences
generate_benchmark_sequences <- function(n_sequences, seq_length = 150, ctcf_motif = TRUE) {
  sequences <- character(n_sequences)
  
  # CTCF-like consensus for realistic benchmarking
  ctcf_consensus <- "CCGCGGGGGCGCT"
  
  for (i in 1:n_sequences) {
    # Generate random DNA sequence
    seq <- paste(sample(c("A", "C", "G", "T"), seq_length, replace = TRUE), collapse = "")
    
    if (ctcf_motif && runif(1) > 0.3) {  # 70% contain motif
      # Insert CTCF-like motif at random position
      pos <- sample(1:(seq_length - nchar(ctcf_consensus) + 1), 1)
      
      # Add some variation to the motif
      motif_var <- ctcf_consensus
      if (runif(1) > 0.8) {  # 20% have 1-2 mutations
        mut_pos <- sample(1:nchar(motif_var), sample(1:2, 1))
        for (mp in mut_pos) {
          substr(motif_var, mp, mp) <- sample(c("A", "C", "G", "T"), 1)
        }
      }
      
      substr(seq, pos, pos + nchar(motif_var) - 1) <- motif_var
    }
    
    sequences[i] <- seq
  }
  
  return(DNAStringSet(sequences))
}

#' Benchmark PWM construction method
#' @param sequences Input sequences
#' @param method PWM construction method
#' @param repetitions Number of repetitions for timing
#' @return List with timing and memory usage results
benchmark_pwm_method <- function(sequences, method = "center", repetitions = 5) {
  
  # Define the benchmark function based on method
  benchmark_func <- switch(method,
    "center" = function() {
      if (exists("build_center_aligned_pwm")) {
        build_center_aligned_pwm(sequences)
      } else {
        # Fallback basic implementation
        build_basic_pwm(sequences)
      }
    },
    "consensus" = function() {
      if (exists("build_consensus_aligned_pwm")) {
        build_consensus_aligned_pwm(sequences)
      } else {
        build_basic_pwm(sequences)
      }
    },
    "integrated" = function() {
      if (exists("build_pwm_robust")) {
        build_pwm_robust(sequences)
      } else {
        build_basic_pwm(sequences)
      }
    },
    function() build_basic_pwm(sequences)
  )
  
  # Memory profiling function
  profile_memory <- function() {
    gc()  # Force garbage collection
    mem_before <- object_size(ls(envir = .GlobalEnv))
    
    result <- benchmark_func()
    
    gc()
    mem_after <- object_size(ls(envir = .GlobalEnv))
    
    list(
      result = result,
      memory_delta = as.numeric(mem_after - mem_before)
    )
  }
  
  # Timing benchmark
  timing_result <- microbenchmark(
    benchmark_func(),
    times = repetitions,
    unit = "s"
  )
  
  # Memory benchmark (single run due to overhead)
  memory_result <- profile_memory()
  
  # Calculate quality metrics if possible
  quality_metrics <- tryCatch({
    pwm_result <- memory_result$result
    if (is.list(pwm_result) && "information_content" %in% names(pwm_result)) {
      list(
        total_ic = sum(pwm_result$information_content),
        max_ic = max(pwm_result$information_content),
        mean_ic = mean(pwm_result$information_content)
      )
    } else {
      list(total_ic = NA, max_ic = NA, mean_ic = NA)
    }
  }, error = function(e) {
    list(total_ic = NA, max_ic = NA, mean_ic = NA)
  })
  
  return(list(
    method = method,
    n_sequences = length(sequences),
    timing = timing_result,
    memory_usage_mb = memory_result$memory_delta / (1024^2),
    quality_metrics = quality_metrics,
    median_time_s = median(timing_result$time) / 1e9,
    mean_time_s = mean(timing_result$time) / 1e9,
    sd_time_s = sd(timing_result$time) / 1e9
  ))
}

#' Fallback basic PWM builder for benchmarking
#' @param sequences Input sequences
#' @return Basic PWM result
build_basic_pwm <- function(sequences) {
  if (length(sequences) == 0) return(NULL)
  
  # Convert to character matrix
  seq_chars <- strsplit(as.character(sequences), "")
  seq_length <- max(lengths(seq_chars))
  
  # Count nucleotides at each position
  pwm_matrix <- matrix(0, nrow = 4, ncol = seq_length,
                       dimnames = list(c("A", "C", "G", "T"), 1:seq_length))
  
  for (i in 1:length(seq_chars)) {
    seq_char <- seq_chars[[i]]
    for (j in 1:length(seq_char)) {
      if (seq_char[j] %in% c("A", "C", "G", "T")) {
        pwm_matrix[seq_char[j], j] <- pwm_matrix[seq_char[j], j] + 1
      }
    }
  }
  
  # Convert to frequencies with pseudocounts
  pwm_matrix <- pwm_matrix + 0.25
  pwm_matrix <- t(t(pwm_matrix) / colSums(pwm_matrix))
  
  # Calculate information content
  ic <- apply(pwm_matrix, 2, function(col) {
    2 + sum(col * log2(col + 1e-10), na.rm = TRUE)
  })
  
  return(list(
    pwm = pwm_matrix,
    information_content = ic,
    total_ic = sum(ic),
    quality_grade = if(sum(ic) > 10) "GOOD" else "FAIR"
  ))
}

#' Run comprehensive performance benchmark
#' @param dataset_sizes Vector of dataset sizes to test
#' @param methods Vector of methods to test
#' @param repetitions Number of repetitions per test
#' @return Data frame with benchmark results
run_performance_benchmark <- function(dataset_sizes, methods, repetitions = 5) {
  results <- list()
  total_tests <- length(dataset_sizes) * length(methods)
  test_count <- 0
  
  cat("Running performance benchmark...\n")
  cat("Dataset sizes:", paste(dataset_sizes, collapse = ", "), "\n")
  cat("Methods:", paste(methods, collapse = ", "), "\n")
  cat("Repetitions per test:", repetitions, "\n\n")
  
  for (size in dataset_sizes) {
    cat("Generating", size, "synthetic sequences...\n")
    sequences <- generate_benchmark_sequences(size)
    
    for (method in methods) {
      test_count <- test_count + 1
      cat(sprintf("Test %d/%d: %s method with %d sequences...\n", 
                  test_count, total_tests, method, size))
      
      tryCatch({
        result <- benchmark_pwm_method(sequences, method, repetitions)
        results[[length(results) + 1]] <- result
      }, error = function(e) {
        cat("Error in benchmark:", e$message, "\n")
        # Add failed result
        results[[length(results) + 1]] <- list(
          method = method,
          n_sequences = size,
          median_time_s = NA,
          mean_time_s = NA,
          sd_time_s = NA,
          memory_usage_mb = NA,
          quality_metrics = list(total_ic = NA, max_ic = NA, mean_ic = NA),
          error = e$message
        )
      })
    }
  }
  
  # Convert to data frame
  df_results <- data.frame(
    method = sapply(results, function(x) x$method),
    n_sequences = sapply(results, function(x) x$n_sequences),
    median_time_s = sapply(results, function(x) x$median_time_s),
    mean_time_s = sapply(results, function(x) x$mean_time_s),
    sd_time_s = sapply(results, function(x) x$sd_time_s),
    memory_usage_mb = sapply(results, function(x) x$memory_usage_mb),
    total_ic = sapply(results, function(x) x$quality_metrics$total_ic),
    max_ic = sapply(results, function(x) x$quality_metrics$max_ic),
    mean_ic = sapply(results, function(x) x$quality_metrics$mean_ic),
    stringsAsFactors = FALSE
  )
  
  return(df_results)
}

#' Generate performance benchmark report
#' @param results Benchmark results data frame
#' @param output_file Output HTML file path
generate_benchmark_report <- function(results, output_file) {
  
  # Create temporary R Markdown file
  rmd_content <- '
---
title: "CTCF PWM Pipeline Performance Benchmark Report"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(dplyr)
library(knitr)
```

# Performance Benchmark Results

## Summary Statistics

```{r summary_table}
summary_stats <- results %>%
  group_by(method) %>%
  summarise(
    tests = n(),
    avg_time_s = mean(median_time_s, na.rm = TRUE),
    avg_memory_mb = mean(memory_usage_mb, na.rm = TRUE),
    avg_total_ic = mean(total_ic, na.rm = TRUE),
    .groups = "drop"
  )

kable(summary_stats, digits = 2, caption = "Summary Statistics by Method")
```

## Performance vs Dataset Size

```{r performance_plot, fig.width=12, fig.height=8}
if (nrow(results) > 0 && any(!is.na(results$median_time_s))) {
  p1 <- ggplot(results, aes(x = n_sequences, y = median_time_s, color = method)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Processing Time vs Dataset Size",
         x = "Number of Sequences",
         y = "Median Time (seconds)",
         color = "Method") +
    theme_minimal()
  
  print(p1)
}
```

## Memory Usage vs Dataset Size

```{r memory_plot, fig.width=12, fig.height=6}
if (nrow(results) > 0 && any(!is.na(results$memory_usage_mb))) {
  p2 <- ggplot(results, aes(x = n_sequences, y = memory_usage_mb, color = method)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    labs(title = "Memory Usage vs Dataset Size",
         x = "Number of Sequences",
         y = "Memory Usage (MB)",
         color = "Method") +
    theme_minimal()
  
  print(p2)
}
```

## Quality vs Performance Trade-off

```{r quality_performance, fig.width=12, fig.height=6}
if (nrow(results) > 0 && any(!is.na(results$total_ic)) && any(!is.na(results$median_time_s))) {
  p3 <- ggplot(results, aes(x = median_time_s, y = total_ic, 
                           color = method, size = n_sequences)) +
    geom_point(alpha = 0.7) +
    labs(title = "Quality vs Performance Trade-off",
         x = "Processing Time (seconds)",
         y = "Total Information Content (bits)",
         color = "Method",
         size = "Dataset Size") +
    theme_minimal()
  
  print(p3)
}
```

## Detailed Results Table

```{r detailed_table}
detailed_results <- results %>%
  select(method, n_sequences, median_time_s, memory_usage_mb, total_ic) %>%
  arrange(method, n_sequences)

kable(detailed_results, digits = 3, caption = "Detailed Benchmark Results")
```

## System Information

```{r system_info}
sys_info <- data.frame(
  Property = c("R Version", "Platform", "CPU Cores", "Memory"),
  Value = c(
    paste(R.version$major, R.version$minor, sep = "."),
    R.version$platform,
    parallel::detectCores(),
    paste(round(as.numeric(system("free -m | grep Mem: | awk \'{print $2}\'", intern = TRUE))/1024, 1), "GB")
  )
)

kable(sys_info, caption = "System Configuration")
```
  '
  
  # Write temporary RMD file
  temp_rmd <- tempfile(fileext = ".Rmd")
  writeLines(rmd_content, temp_rmd)
  
  # Render to HTML
  tryCatch({
    rmarkdown::render(temp_rmd, output_file = output_file, 
                     envir = list(results = results))
    cat("Benchmark report generated:", output_file, "\n")
  }, error = function(e) {
    cat("Error generating report:", e$message, "\n")
    # Save results as CSV fallback
    csv_file <- sub("\\.html$", ".csv", output_file)
    write.csv(results, csv_file, row.names = FALSE)
    cat("Results saved as CSV:", csv_file, "\n")
  })
  
  # Clean up
  unlink(temp_rmd)
}

#' Main execution function
main <- function() {
  # Parse command line arguments
  option_list <- list(
    make_option(c("--dataset-sizes"), type = "character", default = "100,500,1000",
                help = "Comma-separated list of dataset sizes to test [default: %default]"),
    make_option(c("--methods"), type = "character", default = "center,consensus,integrated",
                help = "Comma-separated list of methods to test [default: %default]"),
    make_option(c("--repetitions"), type = "integer", default = 5,
                help = "Number of repetitions per test [default: %default]"),
    make_option(c("--output"), type = "character", default = "results/performance_benchmark.html",
                help = "Output HTML report file [default: %default]"),
    make_option(c("--help"), action = "store_true", default = FALSE,
                help = "Show this help message and exit")
  )
  
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  opt <- parse_args(parser)
  
  if (opt$help) {
    print_help(parser)
    return()
  }
  
  # Parse parameters
  dataset_sizes <- as.numeric(strsplit(opt$`dataset-sizes`, ",")[[1]])
  methods <- trimws(strsplit(opt$methods, ",")[[1]])
  repetitions <- opt$repetitions
  output_file <- opt$output
  
  # Validate parameters
  if (any(is.na(dataset_sizes)) || any(dataset_sizes <= 0)) {
    stop("Invalid dataset sizes. Must be positive integers.")
  }
  
  if (repetitions <= 0) {
    stop("Repetitions must be a positive integer.")
  }
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Run benchmark
  cat("Starting performance benchmark...\n")
  start_time <- Sys.time()
  
  results <- run_performance_benchmark(dataset_sizes, methods, repetitions)
  
  end_time <- Sys.time()
  cat("Benchmark completed in", format(end_time - start_time), "\n")
  
  # Generate report
  generate_benchmark_report(results, output_file)
  
  # Print summary
  cat("\n=== BENCHMARK SUMMARY ===\n")
  if (nrow(results) > 0) {
    for (method in unique(results$method)) {
      method_results <- results[results$method == method, ]
      avg_time <- mean(method_results$median_time_s, na.rm = TRUE)
      avg_memory <- mean(method_results$memory_usage_mb, na.rm = TRUE)
      avg_ic <- mean(method_results$total_ic, na.rm = TRUE)
      
      cat(sprintf("%s: %.2fs avg time, %.1f MB avg memory, %.1f bits avg IC\n",
                  method, avg_time, avg_memory, avg_ic))
    }
  } else {
    cat("No successful benchmark results.\n")
  }
  
  cat("Report saved to:", output_file, "\n")
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
