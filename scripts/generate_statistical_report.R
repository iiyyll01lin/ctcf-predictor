# Statistical Reporting Template for PWM Validation
# Author: PWM Improvement Team

library(Biostrings)

# Generate comprehensive statistical validation report
generate_statistical_report <- function(pwm_results, null_results, output_file = NULL) {
  cat("=== CTCF PWM Statistical Validation Report ===\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("Observed Results:\n")
  cat("  Information Content:", pwm_results$observed_ic, "bits\n")
  cat("  Method:", pwm_results$method %||% "Unknown", "\n")
  cat("  Sequences Used:", pwm_results$n_sequences %||% "Unknown", "\n\n")
  
  cat("Null Model Comparison:\n")
  cat("  Null Models Generated:", length(null_results), "\n")
  cat("  Mean Null IC:", round(mean(null_results), 3), "bits\n")
  cat("  Null SD:", round(sd(null_results), 3), "bits\n\n")
  
  cat("Statistical Significance:\n")
  cat("  P-value:", format(pwm_results$p_value, digits = 3), "\n")
  cat("  Effect Size (Cohen's d):", round(pwm_results$effect_size, 2), "\n")
  cat("  Interpretation:", pwm_results$interpretation %||% "Unknown", "\n\n")
  
  cat("Conclusion:\n")
  if (pwm_results$p_value < 0.001) {
    cat("  ✅ HIGHLY SIGNIFICANT: Strong evidence for biological signal\n")
  } else if (pwm_results$p_value < 0.05) {
    cat("  ✅ SIGNIFICANT: Evidence for biological signal\n")
  } else {
    cat("  ❌ NOT SIGNIFICANT: No evidence above random expectation\n")
  }
  
  # Additional quality metrics if available
  if (!is.null(pwm_results$conserved_positions)) {
    cat("\nQuality Metrics:\n")
    cat("  Conserved Positions (>1 bit):", pwm_results$conserved_positions, "\n")
    if (!is.null(pwm_results$avg_info)) {
      cat("  Average Information per Position:", round(pwm_results$avg_info, 3), "bits\n")
    }
  }
  
  # Recommendations based on results
  cat("\nRecommendations:\n")
  if (pwm_results$p_value < 0.001 && abs(pwm_results$effect_size) > 0.8) {
    cat("  ✅ HIGH CONFIDENCE: PWM suitable for publication and applications\n")
  } else if (pwm_results$p_value < 0.01 && abs(pwm_results$effect_size) > 0.5) {
    cat("  ✅ GOOD CONFIDENCE: PWM suitable for most applications with validation\n")
  } else if (pwm_results$p_value < 0.05 && abs(pwm_results$effect_size) > 0.2) {
    cat("  ⚠️  CAUTION: PWM shows significance but requires additional validation\n")
  } else {
    cat("  ❌ NOT RECOMMENDED: Insufficient evidence for biological signal\n")
  }
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
}

# Generate HTML statistical report
generate_html_statistical_report <- function(pwm_results, null_results, output_file) {
  cat("Generating HTML statistical validation report...\n")
  
  # Create HTML content
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>CTCF PWM Statistical Validation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }
        .header { background-color: #f0f8ff; padding: 15px; border-radius: 5px; margin-bottom: 20px; }
        .section { margin: 20px 0; padding: 15px; border-left: 4px solid #007acc; background-color: #f9f9f9; }
        .metric { margin: 10px 0; }
        .significant { color: #008000; font-weight: bold; }
        .not-significant { color: #cc0000; font-weight: bold; }
        .caution { color: #ff8c00; font-weight: bold; }
        .table { border-collapse: collapse; width: 100%; margin: 15px 0; }
        .table th, .table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        .table th { background-color: #f2f2f2; }
    </style>
</head>
<body>")
  
  # Header
  html_content <- paste0(html_content, "
    <div class='header'>
        <h1>CTCF PWM Statistical Validation Report</h1>
        <p><strong>Generated:</strong> ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>
    </div>")
  
  # Observed Results Section
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Observed Results</h2>
        <div class='metric'><strong>Information Content:</strong> ", pwm_results$observed_ic, " bits</div>
        <div class='metric'><strong>Method:</strong> ", pwm_results$method %||% "Unknown", "</div>
        <div class='metric'><strong>Sequences Used:</strong> ", pwm_results$n_sequences %||% "Unknown", "</div>
    </div>")
  
  # Null Model Comparison Section
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Null Model Comparison</h2>
        <div class='metric'><strong>Null Models Generated:</strong> ", length(null_results), "</div>
        <div class='metric'><strong>Mean Null IC:</strong> ", round(mean(null_results), 3), " bits</div>
        <div class='metric'><strong>Null SD:</strong> ", round(sd(null_results), 3), " bits</div>
    </div>")
  
  # Statistical Significance Section
  significance_class <- if (pwm_results$p_value < 0.001) "significant" else if (pwm_results$p_value < 0.05) "significant" else "not-significant"
  
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Statistical Significance</h2>
        <div class='metric'><strong>P-value:</strong> <span class='", significance_class, "'>", format(pwm_results$p_value, digits = 3), "</span></div>
        <div class='metric'><strong>Effect Size (Cohen's d):</strong> ", round(pwm_results$effect_size, 2), "</div>
        <div class='metric'><strong>Interpretation:</strong> ", pwm_results$interpretation %||% "Unknown", "</div>
    </div>")
  
  # Conclusion Section
  conclusion_class <- if (pwm_results$p_value < 0.001) "significant" else if (pwm_results$p_value < 0.05) "significant" else "not-significant"
  conclusion_text <- if (pwm_results$p_value < 0.001) {
    "✅ HIGHLY SIGNIFICANT: Strong evidence for biological signal"
  } else if (pwm_results$p_value < 0.05) {
    "✅ SIGNIFICANT: Evidence for biological signal"  
  } else {
    "❌ NOT SIGNIFICANT: No evidence above random expectation"
  }
  
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Conclusion</h2>
        <div class='metric ", conclusion_class, "'>", conclusion_text, "</div>
    </div>")
  
  # Quality Metrics Section (if available)
  if (!is.null(pwm_results$conserved_positions)) {
    html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Quality Metrics</h2>
        <div class='metric'><strong>Conserved Positions (>1 bit):</strong> ", pwm_results$conserved_positions, "</div>")
    
    if (!is.null(pwm_results$avg_info)) {
      html_content <- paste0(html_content, "
        <div class='metric'><strong>Average Information per Position:</strong> ", round(pwm_results$avg_info, 3), " bits</div>")
    }
    
    html_content <- paste0(html_content, "
    </div>")
  }
  
  # Recommendations Section
  if (pwm_results$p_value < 0.001 && abs(pwm_results$effect_size) > 0.8) {
    rec_class <- "significant"
    rec_text <- "✅ HIGH CONFIDENCE: PWM suitable for publication and applications"
  } else if (pwm_results$p_value < 0.01 && abs(pwm_results$effect_size) > 0.5) {
    rec_class <- "significant"
    rec_text <- "✅ GOOD CONFIDENCE: PWM suitable for most applications with validation"
  } else if (pwm_results$p_value < 0.05 && abs(pwm_results$effect_size) > 0.2) {
    rec_class <- "caution"
    rec_text <- "⚠️ CAUTION: PWM shows significance but requires additional validation"
  } else {
    rec_class <- "not-significant"
    rec_text <- "❌ NOT RECOMMENDED: Insufficient evidence for biological signal"
  }
  
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Recommendations</h2>
        <div class='metric ", rec_class, "'>", rec_text, "</div>
    </div>")
  
  # Close HTML
  html_content <- paste0(html_content, "
</body>
</html>")
  
  # Write HTML file
  writeLines(html_content, output_file)
  cat("HTML statistical report saved to:", output_file, "\n")
}

# Generate statistical report for multiple PWMs
generate_multiple_pwm_report <- function(statistical_results, null_summary, output_file = NULL) {
  if (!is.null(output_file)) {
    sink(output_file)
  }
  
  cat("=== MULTIPLE PWM STATISTICAL VALIDATION REPORT ===\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  pwm_count <- 0
  significant_count <- 0
  
  for (pwm_name in names(statistical_results)) {
    if (pwm_name == "multiple_testing_correction") next
    
    pwm_count <- pwm_count + 1
    cat("PWM:", pwm_name, "\n")
    cat(paste(rep("-", 40), collapse=""), "\n")
    
    # Get first null model type results
    null_types <- names(statistical_results[[pwm_name]])
    first_null_type <- null_types[1]
    pwm_results <- statistical_results[[pwm_name]][[first_null_type]]
    
    if (!is.null(null_summary) && first_null_type %in% names(null_summary)) {
      null_results <- null_summary[[first_null_type]]$total_info
    } else {
      null_results <- c(0.04, 0.041, 0.042)  # Default null baseline
    }
    
    # Add method and sequences info if available
    pwm_results$method <- pwm_name
    pwm_results$n_sequences <- "Unknown"
    
    # Generate report for this PWM
    generate_statistical_report(pwm_results, null_results)
    
    if (pwm_results$p_value < 0.05) {
      significant_count <- significant_count + 1
    }
    
    cat("\n")
  }
  
  # Summary statistics
  cat("=== SUMMARY ===\n")
  cat("Total PWMs evaluated:", pwm_count, "\n")
  cat("Statistically significant PWMs:", significant_count, "\n")
  cat("Success rate:", round(significant_count / pwm_count * 100, 1), "%\n")
  
  # Multiple testing correction results
  if ("multiple_testing_correction" %in% names(statistical_results)) {
    cat("\n=== MULTIPLE TESTING CORRECTION ===\n")
    correction_results <- statistical_results$multiple_testing_correction
    significant_after_correction <- sum(correction_results$significant_adjusted)
    cat("Tests after correction:", nrow(correction_results), "\n")
    cat("Significant after correction:", significant_after_correction, "\n")
    cat("Adjusted success rate:", round(significant_after_correction / nrow(correction_results) * 100, 1), "%\n")
  }
  
  if (!is.null(output_file)) {
    sink()
    cat("Multiple PWM statistical report saved to:", output_file, "\n")
  }
}

# Helper function for null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript generate_statistical_report.R <statistical_results.rds> <null_summary.rds> [output.txt] [output.html]\n")
    cat("  statistical_results.rds: Results from statistical_significance_test.R\n")
    cat("  null_summary.rds: Null model summary statistics\n")
    cat("  output.txt: Text report output file (optional)\n")
    cat("  output.html: HTML report output file (optional)\n")
    quit(status = 1)
  }
  
  statistical_results_file <- args[1]
  null_summary_file <- args[2]
  output_txt <- if (length(args) >= 3) args[3] else "results/statistical_validation_report.txt"
  output_html <- if (length(args) >= 4) args[4] else "results/statistical_validation_report.html"
  
  cat("PWM Statistical Reporting Tool\n")
  cat("==============================\n")
  cat("Statistical results file:", statistical_results_file, "\n")
  cat("Null summary file:", null_summary_file, "\n")
  cat("Text output file:", output_txt, "\n")
  cat("HTML output file:", output_html, "\n\n")
  
  # Load data
  if (!file.exists(statistical_results_file)) {
    stop("Statistical results file not found: ", statistical_results_file)
  }
  
  if (!file.exists(null_summary_file)) {
    stop("Null summary file not found: ", null_summary_file)
  }
  
  statistical_results <- readRDS(statistical_results_file)
  null_summary <- readRDS(null_summary_file)
  
  # Generate reports
  cat("Generating text report...\n")
  generate_multiple_pwm_report(statistical_results, null_summary, output_txt)
  
  cat("Generating HTML report...\n")
  # For HTML, use first PWM as example
  pwm_names <- names(statistical_results)[names(statistical_results) != "multiple_testing_correction"]
  if (length(pwm_names) > 0) {
    first_pwm <- pwm_names[1]
    first_null_type <- names(statistical_results[[first_pwm]])[1]
    pwm_results <- statistical_results[[first_pwm]][[first_null_type]]
    pwm_results$method <- first_pwm
    pwm_results$n_sequences <- "Unknown"
    
    null_results <- null_summary[[first_null_type]]$total_info
    generate_html_statistical_report(pwm_results, null_results, output_html)
  }
  
  cat("Statistical reporting completed!\n")
}
