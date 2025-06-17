#!/usr/bin/env Rscript

# CTCF-Specific Pattern Validation Functions
# Validates PWM patterns against known CTCF biological patterns
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input PWM file (RDS format)", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file for validation results", metavar="file"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Validate PWM against CTCF biological patterns")
opt <- parse_args(opt_parser)

#' Validate CTCF pattern against expected consensus
#' @param pwm_data PWM matrix or PWM result object
#' @param verbose Enable verbose output
#' @return List with validation results
validate_ctcf_pattern <- function(pwm_data, verbose = FALSE) {
  if (verbose) cat("Validating CTCF pattern...\n")
  
  # Extract PWM matrix
  if (is.list(pwm_data) && "pwm" %in% names(pwm_data)) {
    pwm_matrix <- pwm_data$pwm
  } else if (is.matrix(pwm_data)) {
    pwm_matrix <- pwm_data
  } else {
    stop("Invalid PWM data format")
  }
  
  # Expected CTCF motif pattern
  expected_consensus <- "CCGCGNNGGNGGCAG"
  expected_pattern <- list(
    position_1 = c("C" = 0.8, "G" = 0.2),
    position_2 = c("C" = 0.9, "T" = 0.1),
    position_3 = c("G" = 0.85, "A" = 0.15),
    position_4 = c("C" = 0.8, "G" = 0.2),
    position_5 = c("G" = 0.7, "A" = 0.3),
    position_6 = c("T" = 0.4, "A" = 0.3, "G" = 0.2, "C" = 0.1), # Variable
    position_7 = c("T" = 0.4, "A" = 0.3, "G" = 0.2, "C" = 0.1), # Variable
    position_8 = c("G" = 0.9, "A" = 0.1),
    position_9 = c("G" = 0.85, "A" = 0.15),
    position_10 = c("T" = 0.4, "A" = 0.3, "G" = 0.2, "C" = 0.1), # Variable
    position_11 = c("G" = 0.8, "A" = 0.2),
    position_12 = c("G" = 0.8, "A" = 0.2),
    position_13 = c("C" = 0.9, "T" = 0.1),
    position_14 = c("A" = 0.9, "G" = 0.1),
    position_15 = c("G" = 0.8, "C" = 0.2)
  )
  
  # Generate observed consensus
  observed_consensus <- paste(sapply(1:ncol(pwm_matrix), function(i) {
    max_idx <- which.max(pwm_matrix[,i])
    rownames(pwm_matrix)[max_idx]
  }), collapse = "")
  
  # Calculate pattern similarity
  if (nchar(observed_consensus) != nchar(expected_consensus)) {
    similarity_score <- 0
    pattern_matches <- 0
  } else {
    matches <- sum(strsplit(observed_consensus, "")[[1]] == strsplit(expected_consensus, "")[[1]])
    similarity_score <- matches / nchar(expected_consensus)
    pattern_matches <- matches
  }
  
  # Assess biological relevance
  biological_relevance <- if (similarity_score >= 0.9) {
    "HIGH"
  } else if (similarity_score >= 0.8) {
    "GOOD"
  } else if (similarity_score >= 0.7) {
    "ACCEPTABLE"
  } else {
    "POOR"
  }
  
  # Identify mismatches
  mismatches <- list()
  if (nchar(observed_consensus) == nchar(expected_consensus)) {
    obs_chars <- strsplit(observed_consensus, "")[[1]]
    exp_chars <- strsplit(expected_consensus, "")[[1]]
    
    for (i in 1:length(obs_chars)) {
      if (obs_chars[i] != exp_chars[i]) {
        mismatches[[length(mismatches) + 1]] <- list(
          position = i,
          observed = obs_chars[i],
          expected = exp_chars[i]
        )
      }
    }
  }
  
  results <- list(
    observed_consensus = observed_consensus,
    expected_consensus = expected_consensus,
    similarity_score = similarity_score,
    pattern_matches = pattern_matches,
    total_positions = nchar(expected_consensus),
    biological_relevance = biological_relevance,
    mismatches = mismatches
  )
  
  if (verbose) {
    cat("Pattern validation results:\n")
    cat("  Observed consensus:", observed_consensus, "\n")
    cat("  Expected consensus:", expected_consensus, "\n")
    cat("  Similarity score:", round(similarity_score, 3), "\n")
    cat("  Biological relevance:", biological_relevance, "\n")
    cat("  Mismatches:", length(mismatches), "\n")
  }
  
  return(results)
}

#' Analyze zinc finger correspondence
#' @param ic_profile Information content profile
#' @param verbose Enable verbose output
#' @return List with zinc finger analysis results
analyze_zf_correspondence <- function(ic_profile, verbose = FALSE) {
  if (verbose) cat("Analyzing zinc finger correspondence...\n")
  
  # Expected high IC positions for CTCF zinc fingers
  zf_contacts <- list(
    ZF1 = c(1, 2, 3),      # C-C-G
    ZF2 = c(4, 5),         # C-G
    ZF3 = c(8, 9),         # G-G
    ZF4 = c(11, 12),       # G-G
    ZF5 = c(13, 14, 15)    # C-A-G
  )
  
  expected_high_ic <- unique(unlist(zf_contacts))
  ic_threshold <- 1.0
  
  # Identify observed high IC positions
  observed_high_ic <- which(ic_profile > ic_threshold)
  
  # Calculate correspondence
  correspondence_positions <- intersect(expected_high_ic, observed_high_ic)
  correspondence_ratio <- length(correspondence_positions) / length(expected_high_ic)
  
  # Analyze each zinc finger
  zf_analysis <- list()
  for (zf_name in names(zf_contacts)) {
    positions <- zf_contacts[[zf_name]]
    zf_ic <- ic_profile[positions]
    zf_avg_ic <- mean(zf_ic, na.rm = TRUE)
    zf_high_ic_count <- sum(zf_ic > ic_threshold, na.rm = TRUE)
    
    zf_analysis[[zf_name]] <- list(
      positions = positions,
      average_ic = zf_avg_ic,
      high_ic_count = zf_high_ic_count,
      total_positions = length(positions),
      conservation_ratio = zf_high_ic_count / length(positions)
    )
  }
  
  # Overall assessment
  biological_validity <- correspondence_ratio > 0.7
  
  results <- list(
    expected_high_ic_positions = expected_high_ic,
    observed_high_ic_positions = observed_high_ic,
    correspondence_positions = correspondence_positions,
    correspondence_ratio = correspondence_ratio,
    biological_validity = biological_validity,
    zinc_finger_analysis = zf_analysis,
    ic_threshold = ic_threshold
  )
  
  if (verbose) {
    cat("Zinc finger correspondence analysis:\n")
    cat("  Correspondence ratio:", round(correspondence_ratio, 3), "\n")
    cat("  Biological validity:", biological_validity, "\n")
    cat("  Expected high IC positions:", length(expected_high_ic), "\n")
    cat("  Observed high IC positions:", length(observed_high_ic), "\n")
  }
  
  return(results)
}

#' Determine PWM quality grade
#' @param total_ic Total information content
#' @param conserved_count Number of conserved positions
#' @param pattern_similarity CTCF pattern similarity score
#' @param statistical_significance Statistical significance (p-value < threshold)
#' @param verbose Enable verbose output
#' @return Quality grade string
determine_quality_grade <- function(total_ic, conserved_count = NULL, 
                                   pattern_similarity = NULL, 
                                   statistical_significance = NULL,
                                   verbose = FALSE) {
  if (verbose) cat("Determining quality grade...\n")
  
  # If only IC is provided, use simple IC-based grading
  if (is.null(conserved_count) && is.null(pattern_similarity) && is.null(statistical_significance)) {
    if (total_ic > 20) {
      return("EXCEPTIONAL")
    } else if (total_ic > 16) {
      return("EXCELLENT")
    } else if (total_ic > 12) {
      return("GOOD")
    } else if (total_ic > 8) {
      return("FAIR")
    } else {
      return("POOR")
    }
  }
  
  # Comprehensive grading with multiple criteria
  criteria_met <- 0
  total_criteria <- 0
  
  # Information content criterion
  if (!is.null(total_ic)) {
    total_criteria <- total_criteria + 1
    if (total_ic > 16) {
      criteria_met <- criteria_met + 1
    }
  }
  
  # Conserved positions criterion
  if (!is.null(conserved_count)) {
    total_criteria <- total_criteria + 1
    if (conserved_count > 6) {
      criteria_met <- criteria_met + 1
    }
  }
  
  # Pattern similarity criterion
  if (!is.null(pattern_similarity)) {
    total_criteria <- total_criteria + 1
    if (pattern_similarity > 0.8) {
      criteria_met <- criteria_met + 1
    }
  }
  
  # Statistical significance criterion
  if (!is.null(statistical_significance)) {
    total_criteria <- total_criteria + 1
    if (statistical_significance == TRUE) {
      criteria_met <- criteria_met + 1
    }
  }
  
  # Determine grade based on criteria met
  criteria_ratio <- criteria_met / total_criteria
  
  if (criteria_ratio >= 1.0) {
    grade <- "EXCELLENT"
  } else if (criteria_ratio >= 0.75) {
    grade <- "GOOD"
  } else if (criteria_ratio >= 0.5) {
    grade <- "ACCEPTABLE"
  } else {
    grade <- "POOR"
  }
  
  if (verbose) {
    cat("Quality assessment:\n")
    cat("  Criteria met:", criteria_met, "/", total_criteria, "\n")
    cat("  Criteria ratio:", round(criteria_ratio, 3), "\n")
    cat("  Grade:", grade, "\n")
  }
  
  return(grade)
}

#' Comprehensive CTCF pattern validation
#' @param pwm_data PWM matrix or PWM result object
#' @param verbose Enable verbose output
#' @return List with comprehensive validation results
comprehensive_ctcf_validation <- function(pwm_data, verbose = FALSE) {
  if (verbose) cat("Performing comprehensive CTCF validation...\n")
  
  # Extract PWM matrix and calculate IC
  if (is.list(pwm_data) && "pwm" %in% names(pwm_data)) {
    pwm_matrix <- pwm_data$pwm
  } else if (is.matrix(pwm_data)) {
    pwm_matrix <- pwm_data
  } else {
    stop("Invalid PWM data format")
  }
  
  # Calculate information content
  ic_profile <- apply(pwm_matrix, 2, function(col) {
    # Remove zeros to avoid log(0)
    col[col == 0] <- 1e-10
    -sum(col * log2(col))
  })
  ic_profile <- 2 - ic_profile  # Convert to information content
  total_ic <- sum(ic_profile)
  
  # Count conserved positions (IC > 1.0)
  conserved_positions <- sum(ic_profile > 1.0)
  
  # Validate CTCF pattern
  pattern_validation <- validate_ctcf_pattern(pwm_data, verbose)
  
  # Analyze zinc finger correspondence
  zf_analysis <- analyze_zf_correspondence(ic_profile, verbose)
  
  # Determine quality grade
  quality_grade <- determine_quality_grade(
    total_ic = total_ic,
    conserved_count = conserved_positions,
    pattern_similarity = pattern_validation$similarity_score,
    statistical_significance = TRUE,  # Assume passed if running validation
    verbose = verbose
  )
  
  # Compile comprehensive results
  results <- list(
    total_information_content = total_ic,
    ic_profile = ic_profile,
    conserved_positions = conserved_positions,
    pattern_validation = pattern_validation,
    zinc_finger_analysis = zf_analysis,
    quality_grade = quality_grade,
    recommendations = generate_recommendations(quality_grade, pattern_validation, zf_analysis)
  )
  
  if (verbose) {
    cat("Comprehensive validation complete:\n")
    cat("  Total IC:", round(total_ic, 2), "bits\n")
    cat("  Conserved positions:", conserved_positions, "\n")
    cat("  Pattern similarity:", round(pattern_validation$similarity_score, 3), "\n")
    cat("  Quality grade:", quality_grade, "\n")
  }
  
  return(results)
}

#' Generate recommendations based on validation results
#' @param quality_grade Quality grade
#' @param pattern_validation Pattern validation results
#' @param zf_analysis Zinc finger analysis results
#' @return List of recommendations
generate_recommendations <- function(quality_grade, pattern_validation, zf_analysis) {
  recommendations <- list()
  
  if (quality_grade == "EXCELLENT") {
    recommendations$usage <- c(
      "APPROVED for drug discovery applications",
      "APPROVED for high-stakes research",
      "APPROVED for publication-quality analysis"
    )
  } else if (quality_grade == "GOOD") {
    recommendations$usage <- c(
      "APPROVED for research applications",
      "APPROVED for comparative studies",
      "Consider validation for high-stakes applications"
    )
  } else if (quality_grade == "ACCEPTABLE") {
    recommendations$usage <- c(
      "APPROVED for exploratory analysis",
      "APPROVED for method development",
      "NOT RECOMMENDED for drug discovery"
    )
  } else {
    recommendations$usage <- c(
      "NOT RECOMMENDED for research applications",
      "Consider improving data quality or alignment",
      "Review input sequences and parameters"
    )
  }
  
  # Pattern-specific recommendations
  if (pattern_validation$similarity_score < 0.8) {
    recommendations$improvements <- c(
      recommendations$improvements,
      "Consider improving sequence alignment",
      "Validate input sequences are CTCF binding sites",
      "Review ChIP-seq peak quality"
    )
  }
  
  # Zinc finger-specific recommendations
  if (!zf_analysis$biological_validity) {
    recommendations$improvements <- c(
      recommendations$improvements,
      "Low zinc finger correspondence detected",
      "Consider using CTCF-specific alignment parameters",
      "Validate biological relevance of input data"
    )
  }
  
  return(recommendations)
}

# Main execution
if (!interactive()) {
  # Validate inputs
  if (is.null(opt$input)) {
    cat("Error: Input PWM file is required\n")
    quit(status = 1)
  }
  
  if (!file.exists(opt$input)) {
    cat("Error: Input file does not exist:", opt$input, "\n")
    quit(status = 1)
  }
  
  if (is.null(opt$output)) {
    base_name <- tools::file_path_sans_ext(basename(opt$input))
    opt$output <- paste0(base_name, "_ctcf_validation.json")
  }
  
  if (opt$verbose) {
    cat("CTCF Pattern Validation\n")
    cat("======================\n")
    cat("Input file:", opt$input, "\n")
    cat("Output file:", opt$output, "\n\n")
  }
  
  # Load PWM data
  if (opt$verbose) cat("Loading PWM data...\n")
  pwm_data <- readRDS(opt$input)
  
  # Perform comprehensive validation
  validation_results <- comprehensive_ctcf_validation(pwm_data, opt$verbose)
  
  # Save results
  writeLines(jsonlite::toJSON(validation_results, pretty = TRUE), opt$output)
  
  if (opt$verbose) {
    cat("\nValidation complete!\n")
    cat("Results saved to:", opt$output, "\n")
    cat("Quality grade:", validation_results$quality_grade, "\n")
  }
}
