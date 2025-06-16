# Comprehensive quality assessment for PWM statistical validation
# Author: PWM Improvement Team

library(Biostrings)

# Comprehensive quality assessment for PWM statistical validation
assess_statistical_quality <- function(pwm_results, null_results) {
  cat("Performing comprehensive statistical quality assessment...\n")
  
  checks <- list()
  
  # 1. Significance check
  if (!is.null(pwm_results$p_value)) {
    checks$is_significant <- pwm_results$p_value < 0.05
    checks$is_highly_significant <- pwm_results$p_value < 0.001
  } else {
    checks$is_significant <- FALSE
    checks$is_highly_significant <- FALSE
  }
  
  # 2. Effect size check
  if (!is.null(pwm_results$effect_size)) {
    checks$large_effect <- abs(pwm_results$effect_size) > 0.8
    checks$moderate_effect <- abs(pwm_results$effect_size) > 0.5
  } else {
    checks$large_effect <- FALSE
    checks$moderate_effect <- FALSE
  }
  
  # 3. Information content check
  if (!is.null(pwm_results$observed_ic)) {
    checks$sufficient_ic <- pwm_results$observed_ic > 8.0
    checks$high_ic <- pwm_results$observed_ic > 15.0
  } else {
    checks$sufficient_ic <- FALSE
    checks$high_ic <- FALSE
  }
  
  # 4. Null model validation
  checks$valid_nulls <- validate_null_models(null_results)
  
  # 5. Consistency check
  checks$consistent_results <- check_method_consistency(pwm_results)
  
  # Overall assessment
  checks$overall_quality <- all(unlist(checks[c("is_significant", "moderate_effect", 
                                               "sufficient_ic", "valid_nulls")]))
  
  # Quality score (0-100)
  quality_components <- c(
    significance = ifelse(checks$is_highly_significant, 25, ifelse(checks$is_significant, 15, 0)),
    effect_size = ifelse(checks$large_effect, 25, ifelse(checks$moderate_effect, 15, 0)),
    information_content = ifelse(checks$high_ic, 25, ifelse(checks$sufficient_ic, 15, 0)),
    null_validation = ifelse(checks$valid_nulls, 25, 0)
  )
  
  checks$quality_score <- sum(quality_components)
  checks$quality_breakdown <- quality_components
  
  return(checks)
}

# Validate null models
validate_null_models <- function(null_results) {
  if (is.null(null_results) || length(null_results) == 0) {
    return(FALSE)
  }
  
  # Check if null models have appropriate properties
  if (length(null_results) < 50) {
    return(FALSE)
  }
  
  # Check if null distribution is reasonable
  null_mean <- mean(null_results)
  null_sd <- sd(null_results)
  
  # Null models should have low IC and reasonable variance
  reasonable_mean <- null_mean < 2.0
  reasonable_variance <- null_sd > 0.01 && null_sd < 0.5
  
  return(reasonable_mean && reasonable_variance)
}

# Check method consistency (placeholder)
check_method_consistency <- function(pwm_results) {
  # This would compare results across multiple methods
  # For now, return TRUE as placeholder
  return(TRUE)
}

# Generate quality assessment report
generate_quality_report <- function(quality_results, pwm_name = "PWM") {
  cat("\n=== Quality Assessment Report for", pwm_name, "===\n")
  cat("Quality Score:", quality_results$quality_score, "/100\n")
  
  cat("\nDetailed Assessment:\n")
  cat("✓ Statistical Significance:", ifelse(quality_results$is_significant, "PASS", "FAIL"))
  if (quality_results$is_highly_significant) cat(" (highly significant)")
  cat("\n")
  
  cat("✓ Effect Size:", ifelse(quality_results$moderate_effect, "PASS", "FAIL"))
  if (quality_results$large_effect) cat(" (large effect)")
  cat("\n")
  
  cat("✓ Information Content:", ifelse(quality_results$sufficient_ic, "PASS", "FAIL"))
  if (quality_results$high_ic) cat(" (high IC)")
  cat("\n")
  
  cat("✓ Null Model Validation:", ifelse(quality_results$valid_nulls, "PASS", "FAIL"), "\n")
  cat("✓ Overall Quality:", ifelse(quality_results$overall_quality, "PASS", "FAIL"), "\n")
  
  cat("\nQuality Breakdown:\n")
  for (component in names(quality_results$quality_breakdown)) {
    cat(" -", component, ":", quality_results$quality_breakdown[component], "/25\n")
  }
  
  # Interpretation
  if (quality_results$quality_score >= 80) {
    cat("\nInterpretation: EXCELLENT quality PWM\n")
  } else if (quality_results$quality_score >= 60) {
    cat("\nInterpretation: GOOD quality PWM\n")
  } else if (quality_results$quality_score >= 40) {
    cat("\nInterpretation: MODERATE quality PWM - consider improvements\n")
  } else {
    cat("\nInterpretation: POOR quality PWM - significant improvements needed\n")
  }
}

# Assess multiple PWMs
assess_multiple_pwms <- function(statistical_results, null_summary) {
  cat("Assessing quality for multiple PWMs...\n")
  
  quality_assessments <- list()
  
  for (pwm_name in names(statistical_results)) {
    if (pwm_name == "multiple_testing_correction") next
    
    cat("\nAssessing", pwm_name, "...\n")
    
    # Get the first null model type results (assuming they're similar)
    null_types <- names(statistical_results[[pwm_name]])
    first_null_type <- null_types[1]
    
    pwm_results <- statistical_results[[pwm_name]][[first_null_type]]
    
    # Extract null results for this PWM type
    null_results <- NULL
    if (!is.null(null_summary) && first_null_type %in% names(null_summary)) {
      null_results <- null_summary[[first_null_type]]$total_info
    }
    
    # Assess quality
    quality_assessment <- assess_statistical_quality(pwm_results$total_info, null_results)
    quality_assessments[[pwm_name]] <- quality_assessment
    
    # Generate report
    generate_quality_report(quality_assessment, pwm_name)
  }
  
  return(quality_assessments)
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  statistical_results_file <- if (length(args) >= 1) args[1] else "results/statistical_significance_results.rds"
  null_summary_file <- if (length(args) >= 2) args[2] else "results/null_models/null_summary_statistics.rds"
  output_file <- if (length(args) >= 3) args[3] else "results/quality_assessment_report.txt"
  
  cat("PWM Quality Assessment Tool\n")
  cat("===========================\n")
  cat("Statistical results file:", statistical_results_file, "\n")
  cat("Null summary file:", null_summary_file, "\n")
  cat("Output file:", output_file, "\n\n")
  
  # Load statistical results
  if (!file.exists(statistical_results_file)) {
    stop("Statistical results file not found: ", statistical_results_file)
  }
  
  statistical_results <- readRDS(statistical_results_file)
  
  # Load null summary (optional)
  null_summary <- NULL
  if (file.exists(null_summary_file)) {
    null_summary <- readRDS(null_summary_file)
    cat("Loaded null model summary\n")
  } else {
    cat("Warning: Null summary file not found, proceeding without null validation\n")
  }
  
  # Redirect output to file
  sink(output_file)
  
  # Assess quality
  quality_assessments <- assess_multiple_pwms(statistical_results, null_summary)
  
  # Stop redirecting output
  sink()
  
  # Save quality assessments
  output_rds <- sub("\\.txt$", ".rds", output_file)
  saveRDS(quality_assessments, output_rds)
  
  cat("Quality assessment completed!\n")
  cat("Report saved to:", output_file, "\n")
  cat("Data saved to:", output_rds, "\n")
}
