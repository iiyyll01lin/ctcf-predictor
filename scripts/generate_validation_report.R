#!/usr/bin/env Rscript

#' Validation Report Generator for CTCF PWM Testing Pipeline
#' 
#' This script generates comprehensive validation reports including quality
#' certificates and summary statistics for PWM validation results.
#' 
#' Usage:
#'   Rscript scripts/generate_validation_report.R \
#'     --validation-results results/validation/ \
#'     --output results/validation_report.html

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(knitr)
  library(rmarkdown)
  library(ggplot2)
  library(dplyr)
})

#' Load validation results from directory
#' @param results_dir Directory containing validation result files
#' @return List of validation results
load_validation_results <- function(results_dir) {
  if (!dir.exists(results_dir)) {
    stop("Validation results directory does not exist: ", results_dir)
  }
  
  # Find all validation result files
  result_files <- list.files(results_dir, 
                            pattern = "\\.(rds|json|csv)$", 
                            full.names = TRUE, 
                            recursive = TRUE)
  
  if (length(result_files) == 0) {
    warning("No validation result files found in: ", results_dir)
    return(list())
  }
  
  results <- list()
  
  for (file in result_files) {
    tryCatch({
      ext <- tools::file_ext(file)
      result_name <- tools::file_path_sans_ext(basename(file))
      
      result_data <- switch(ext,
        "rds" = readRDS(file),
        "json" = fromJSON(file),
        "csv" = read.csv(file, stringsAsFactors = FALSE),
        NULL
      )
      
      if (!is.null(result_data)) {
        results[[result_name]] <- result_data
      }
    }, error = function(e) {
      warning("Failed to load ", file, ": ", e$message)
    })
  }
  
  return(results)
}

#' Extract key metrics from validation results
#' @param validation_results List of validation results
#' @return Data frame with key metrics
extract_key_metrics <- function(validation_results) {
  metrics <- data.frame(
    test_name = character(),
    status = character(),
    total_ic = numeric(),
    max_ic = numeric(),
    mean_ic = numeric(),
    quality_grade = character(),
    p_value = numeric(),
    confidence_level = numeric(),
    biological_similarity = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (test_name in names(validation_results)) {
    result <- validation_results[[test_name]]
    
    # Extract metrics with safe defaults
    row <- data.frame(
      test_name = test_name,
      status = extract_safely(result, "status", "UNKNOWN"),
      total_ic = extract_safely(result, "total_ic", NA),
      max_ic = extract_safely(result, "max_ic", NA),
      mean_ic = extract_safely(result, "mean_ic", NA),
      quality_grade = extract_safely(result, "quality_grade", "UNKNOWN"),
      p_value = extract_safely(result, "p_value", NA),
      confidence_level = extract_safely(result, "confidence_level", NA),
      biological_similarity = extract_safely(result, "biological_similarity", NA),
      stringsAsFactors = FALSE
    )
    
    metrics <- rbind(metrics, row)
  }
  
  return(metrics)
}

#' Safely extract value from nested list
#' @param obj Source object
#' @param key Key to extract
#' @param default Default value if key not found
#' @return Extracted value or default
extract_safely <- function(obj, key, default = NA) {
  tryCatch({
    if (is.list(obj) && key %in% names(obj)) {
      return(obj[[key]])
    } else if (is.list(obj)) {
      # Try to find in nested structures
      for (subkey in names(obj)) {
        if (is.list(obj[[subkey]]) && key %in% names(obj[[subkey]])) {
          return(obj[[subkey]][[key]])
        }
      }
    }
    return(default)
  }, error = function(e) {
    return(default)
  })
}

#' Determine overall quality grade
#' @param metrics Key metrics data frame
#' @return Overall quality grade
determine_overall_quality <- function(metrics) {
  if (nrow(metrics) == 0) return("UNKNOWN")
  
  # Count status types
  status_counts <- table(metrics$status)
  
  # Check for failures
  if ("FAIL" %in% names(status_counts) && status_counts["FAIL"] > 0) {
    return("POOR")
  }
  
  # Check total IC values
  ic_values <- metrics$total_ic[!is.na(metrics$total_ic)]
  if (length(ic_values) > 0) {
    avg_ic <- mean(ic_values)
    if (avg_ic >= 15) return("EXCELLENT")
    if (avg_ic >= 12) return("GOOD") 
    if (avg_ic >= 8) return("FAIR")
    return("POOR")
  }
  
  # Check quality grades
  grades <- metrics$quality_grade[metrics$quality_grade != "UNKNOWN"]
  if (length(grades) > 0) {
    if (any(grades == "EXCELLENT")) return("GOOD")
    if (any(grades == "GOOD")) return("FAIR")
  }
  
  return("UNKNOWN")
}

#' Generate quality certificate
#' @param validation_results Validation results
#' @param metrics Key metrics
#' @return Quality certificate list
generate_quality_certificate <- function(validation_results, metrics) {
  overall_quality <- determine_overall_quality(metrics)
  
  # Count validation passes
  passes <- sum(metrics$status == "PASS", na.rm = TRUE)
  total <- nrow(metrics)
  
  # Calculate certification level
  pass_rate <- if (total > 0) passes / total else 0
  cert_level <- if (pass_rate >= 0.9) "GOLD" else 
                if (pass_rate >= 0.8) "SILVER" else 
                if (pass_rate >= 0.7) "BRONZE" else "BASIC"
  
  certificate <- list(
    certificate_id = generate_certificate_id(),
    validation_date = Sys.Date(),
    overall_quality_grade = overall_quality,
    certification_level = cert_level,
    validation_passes = passes,
    total_validations = total,
    pass_rate = pass_rate,
    key_metrics = list(
      avg_total_ic = mean(metrics$total_ic, na.rm = TRUE),
      avg_max_ic = mean(metrics$max_ic, na.rm = TRUE),
      significant_tests = sum(metrics$p_value < 0.05, na.rm = TRUE),
      high_confidence_tests = sum(metrics$confidence_level > 0.95, na.rm = TRUE)
    ),
    validation_signature = generate_validation_signature(validation_results)
  )
  
  return(certificate)
}

#' Generate unique certificate ID
#' @return Certificate ID string
generate_certificate_id <- function() {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  random_suffix <- sample(1000:9999, 1)
  return(paste0("CTCF_", timestamp, "_", random_suffix))
}

#' Generate validation signature (checksum)
#' @param validation_results Validation results
#' @return Signature string
generate_validation_signature <- function(validation_results) {
  # Create a simple signature based on results
  content <- paste(names(validation_results), collapse = "|")
  signature <- digest::digest(content, algo = "md5")
  return(substr(signature, 1, 12))  # First 12 characters
}

#' Generate validation report
#' @param validation_results Validation results
#' @param output_file Output HTML file path
generate_validation_report <- function(validation_results, output_file) {
  
  # Extract metrics and generate certificate
  metrics <- extract_key_metrics(validation_results)
  certificate <- generate_quality_certificate(validation_results, metrics)
  
  # Create R Markdown content
  rmd_content <- '
---
title: "CTCF PWM Pipeline Validation Report"
date: "`r Sys.Date()`"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: bootstrap
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(dplyr)
library(knitr)
library(DT)
```

# Executive Summary

## Quality Certificate

```{r certificate, results="asis"}
cert <- certificate

cat("**Certificate ID:**", cert$certificate_id, "\\n\\n")
cat("**Validation Date:**", as.character(cert$validation_date), "\\n\\n")
cat("**Overall Quality Grade:**", cert$overall_quality_grade, "\\n\\n")
cat("**Certification Level:**", cert$certification_level, "\\n\\n")
cat("**Pass Rate:**", sprintf("%.1f%% (%d/%d)", cert$pass_rate * 100, cert$validation_passes, cert$total_validations), "\\n\\n")
```

## Key Performance Indicators

```{r kpis}
kpi_data <- data.frame(
  Metric = c(
    "Average Total IC",
    "Average Max IC", 
    "Significant Tests",
    "High Confidence Tests"
  ),
  Value = c(
    sprintf("%.2f bits", certificate$key_metrics$avg_total_ic),
    sprintf("%.2f bits", certificate$key_metrics$avg_max_ic),
    paste(certificate$key_metrics$significant_tests, "tests"),
    paste(certificate$key_metrics$high_confidence_tests, "tests")
  ),
  stringsAsFactors = FALSE
)

kable(kpi_data, caption = "Key Performance Indicators")
```

# Detailed Validation Results

## Test Results Summary

```{r summary_table}
summary_data <- metrics %>%
  group_by(status) %>%
  summarise(
    count = n(),
    avg_total_ic = mean(total_ic, na.rm = TRUE),
    avg_p_value = mean(p_value, na.rm = TRUE),
    .groups = "drop"
  )

kable(summary_data, digits = 3, caption = "Summary by Test Status")
```

## Information Content Distribution

```{r ic_distribution, fig.width=12, fig.height=6}
if (nrow(metrics) > 0 && any(!is.na(metrics$total_ic))) {
  p1 <- ggplot(metrics, aes(x = reorder(test_name, total_ic), y = total_ic, fill = status)) +
    geom_col() +
    coord_flip() +
    labs(title = "Total Information Content by Test",
         x = "Test Name",
         y = "Total IC (bits)",
         fill = "Status") +
    theme_minimal() +
    scale_fill_manual(values = c("PASS" = "green3", "FAIL" = "red3", "WARNING" = "orange"))
  
  print(p1)
}
```

## Quality Grade Distribution

```{r quality_distribution, fig.width=10, fig.height=6}
if (nrow(metrics) > 0) {
  quality_counts <- table(metrics$quality_grade)
  
  if (length(quality_counts) > 0) {
    quality_df <- data.frame(
      grade = names(quality_counts),
      count = as.numeric(quality_counts),
      stringsAsFactors = FALSE
    )
    
    p2 <- ggplot(quality_df, aes(x = grade, y = count, fill = grade)) +
      geom_col() +
      labs(title = "Distribution of Quality Grades",
           x = "Quality Grade",
           y = "Number of Tests") +
      theme_minimal() +
      scale_fill_manual(values = c(
        "EXCELLENT" = "darkgreen",
        "GOOD" = "green3", 
        "FAIR" = "yellow3",
        "POOR" = "orange",
        "UNKNOWN" = "gray"
      ))
    
    print(p2)
  }
}
```

## Statistical Significance

```{r significance_plot, fig.width=10, fig.height=6}
if (nrow(metrics) > 0 && any(!is.na(metrics$p_value))) {
  significance_data <- metrics[!is.na(metrics$p_value), ]
  
  p3 <- ggplot(significance_data, aes(x = p_value)) +
    geom_histogram(bins = 20, fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
    labs(title = "Distribution of P-values",
         x = "P-value",
         y = "Count",
         caption = "Red line indicates significance threshold (p < 0.05)") +
    theme_minimal()
  
  print(p3)
}
```

## Detailed Results Table

```{r detailed_table}
if (nrow(metrics) > 0) {
  display_metrics <- metrics %>%
    select(test_name, status, quality_grade, total_ic, p_value) %>%
    mutate(
      total_ic = round(total_ic, 2),
      p_value = ifelse(is.na(p_value), "-", sprintf("%.4f", p_value))
    )
  
  datatable(display_metrics, 
            caption = "Detailed Validation Results",
            options = list(pageLength = 15, scrollX = TRUE),
            filter = "top")
}
```

# Validation Methodology

## Test Descriptions

```{r test_descriptions, results="asis"}
test_descriptions <- list(
  "data_validation" = "Input data quality and format validation",
  "pwm_construction" = "PWM construction quality and robustness testing",
  "statistical_validation" = "Statistical significance and null model comparison",
  "bootstrap_validation" = "Bootstrap confidence interval validation",
  "biological_validation" = "CTCF-specific biological pattern validation",
  "alignment_quality" = "Sequence alignment quality assessment"
)

for (test_type in names(test_descriptions)) {
  if (any(grepl(test_type, metrics$test_name, ignore.case = TRUE))) {
    cat("**", stringr::str_to_title(gsub("_", " ", test_type)), ":**", 
        test_descriptions[[test_type]], "\\n\\n")
  }
}
```

## Quality Thresholds

```{r thresholds}
threshold_data <- data.frame(
  Metric = c("Total Information Content", "P-value", "Biological Similarity", "Bootstrap Confidence"),
  Excellent = c(">= 15 bits", "< 0.001", ">= 0.8", ">= 0.95"),
  Good = c(">= 12 bits", "< 0.01", ">= 0.7", ">= 0.90"),
  Fair = c(">= 8 bits", "< 0.05", ">= 0.6", ">= 0.80"),
  Poor = c("< 8 bits", ">= 0.05", "< 0.6", "< 0.80"),
  stringsAsFactors = FALSE
)

kable(threshold_data, caption = "Quality Assessment Thresholds")
```

# Recommendations

```{r recommendations, results="asis"}
# Generate recommendations based on results
recommendations <- c()

# Check for failed tests
failed_tests <- metrics[metrics$status == "FAIL", ]
if (nrow(failed_tests) > 0) {
  recommendations <- c(recommendations, 
    paste("- Address", nrow(failed_tests), "failed validation tests"))
}

# Check IC levels
low_ic_tests <- metrics[!is.na(metrics$total_ic) & metrics$total_ic < 8, ]
if (nrow(low_ic_tests) > 0) {
  recommendations <- c(recommendations,
    paste("- Improve information content for", nrow(low_ic_tests), "tests (< 8 bits)"))
}

# Check significance
non_sig_tests <- metrics[!is.na(metrics$p_value) & metrics$p_value >= 0.05, ]
if (nrow(non_sig_tests) > 0) {
  recommendations <- c(recommendations,
    paste("- Review", nrow(non_sig_tests), "non-significant tests"))
}

if (length(recommendations) == 0) {
  recommendations <- "- All validation tests meet quality standards"
}

for (rec in recommendations) {
  cat(rec, "\\n")
}
```

---

*Report generated on `r Sys.time()` by CTCF PWM Testing Pipeline Validation System*
  '
  
  # Write temporary RMD file
  temp_rmd <- tempfile(fileext = ".Rmd")
  writeLines(rmd_content, temp_rmd)
  
  # Render to HTML
  tryCatch({
    rmarkdown::render(temp_rmd, output_file = output_file, 
                     envir = list(validation_results = validation_results,
                                 metrics = metrics,
                                 certificate = certificate))
    cat("Validation report generated:", output_file, "\n")
  }, error = function(e) {
    cat("Error generating HTML report:", e$message, "\n")
    # Generate text report as fallback
    generate_text_report(validation_results, metrics, certificate, output_file)
  })
  
  # Clean up
  unlink(temp_rmd)
  
  return(certificate)
}

#' Generate text report fallback
#' @param validation_results Validation results
#' @param metrics Key metrics
#' @param certificate Quality certificate
#' @param output_file Output file path
generate_text_report <- function(validation_results, metrics, certificate, output_file) {
  text_file <- sub("\\.html$", ".txt", output_file)
  
  report_lines <- c(
    "CTCF PWM Pipeline Validation Report",
    paste("Generated:", Sys.time()),
    "",
    "=== QUALITY CERTIFICATE ===",
    paste("Certificate ID:", certificate$certificate_id),
    paste("Validation Date:", certificate$validation_date),
    paste("Overall Quality Grade:", certificate$overall_quality_grade),
    paste("Certification Level:", certificate$certification_level),
    paste("Pass Rate:", sprintf("%.1f%% (%d/%d)", 
                               certificate$pass_rate * 100,
                               certificate$validation_passes,
                               certificate$total_validations)),
    "",
    "=== KEY METRICS ===",
    paste("Average Total IC:", sprintf("%.2f bits", certificate$key_metrics$avg_total_ic)),
    paste("Average Max IC:", sprintf("%.2f bits", certificate$key_metrics$avg_max_ic)),
    paste("Significant Tests:", certificate$key_metrics$significant_tests),
    paste("High Confidence Tests:", certificate$key_metrics$high_confidence_tests),
    "",
    "=== DETAILED RESULTS ==="
  )
  
  if (nrow(metrics) > 0) {
    for (i in 1:nrow(metrics)) {
      row <- metrics[i, ]
      report_lines <- c(report_lines,
        paste("Test:", row$test_name),
        paste("  Status:", row$status),
        paste("  Quality:", row$quality_grade),
        paste("  Total IC:", ifelse(is.na(row$total_ic), "N/A", sprintf("%.2f bits", row$total_ic))),
        paste("  P-value:", ifelse(is.na(row$p_value), "N/A", sprintf("%.4f", row$p_value))),
        ""
      )
    }
  }
  
  writeLines(report_lines, text_file)
  cat("Text report saved to:", text_file, "\n")
}

#' Main execution function
main <- function() {
  # Parse command line arguments
  option_list <- list(
    make_option(c("--validation-results"), type = "character", 
                default = "results/validation/",
                help = "Directory containing validation results [default: %default]"),
    make_option(c("--output"), type = "character", 
                default = "results/validation_report.html",
                help = "Output HTML report file [default: %default]"),
    make_option(c("--certificate-only"), action = "store_true", default = FALSE,
                help = "Generate only quality certificate (JSON format)"),
    make_option(c("--help"), action = "store_true", default = FALSE,
                help = "Show this help message and exit")
  )
  
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  opt <- parse_args(parser)
  
  if (opt$help) {
    print_help(parser)
    return()
  }
  
  # Create output directory if needed
  output_dir <- dirname(opt$output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load validation results
  cat("Loading validation results from:", opt$`validation-results`, "\n")
  validation_results <- load_validation_results(opt$`validation-results`)
  
  if (length(validation_results) == 0) {
    stop("No validation results found")
  }
  
  cat("Loaded", length(validation_results), "validation results\n")
  
  # Generate report or certificate
  if (opt$`certificate-only`) {
    # Generate only certificate
    metrics <- extract_key_metrics(validation_results)
    certificate <- generate_quality_certificate(validation_results, metrics)
    
    cert_file <- sub("\\.html$", "_certificate.json", opt$output)
    writeLines(toJSON(certificate, pretty = TRUE), cert_file)
    cat("Quality certificate saved to:", cert_file, "\n")
  } else {
    # Generate full report
    certificate <- generate_validation_report(validation_results, opt$output)
    
    # Also save certificate separately
    cert_file <- sub("\\.html$", "_certificate.json", opt$output)
    writeLines(toJSON(certificate, pretty = TRUE), cert_file)
    cat("Quality certificate saved to:", cert_file, "\n")
  }
  
  cat("Validation report generation complete!\n")
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
