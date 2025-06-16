# PWM Quality Validation Script

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
pwm_file <- if (length(args) >= 1) args[1] else "results/generated_pwm.rds"
sequences_file <- if (length(args) >= 2) args[2] else NULL
output_file <- if (length(args) >= 3) args[3] else NULL

# --- Validation Functions ---

# Test information content
test_information_content <- function(pwm) {
  ic_per_pos <- apply(pwm, 2, function(pos_probs) {
    valid_probs <- pmax(pos_probs, 1e-10)
    sum(valid_probs * log2(valid_probs / 0.25))
  })
  
  total_ic <- sum(ic_per_pos)
  
  return(list(
    total_ic = total_ic,
    ic_per_position = ic_per_pos,
    mean_ic = mean(ic_per_pos),
    max_ic = max(ic_per_pos),
    conserved_positions = sum(ic_per_pos > 1.0),
    passes_minimum = total_ic >= 8.0,
    quality_level = if (total_ic >= 16.0) "excellent"
                   else if (total_ic >= 12.0) "good"
                   else if (total_ic >= 8.0) "acceptable"
                   else "poor"
  ))
}

# Test conservation pattern
test_conservation_pattern <- function(pwm) {
  # Find positions with high conservation (>1.5 bits)
  ic_per_pos <- apply(pwm, 2, function(pos_probs) {
    valid_probs <- pmax(pos_probs, 1e-10)
    sum(valid_probs * log2(valid_probs / 0.25))
  })
  
  high_conservation <- ic_per_pos > 1.5
  conserved_positions <- which(high_conservation)
  
  # Check for conserved core
  has_conserved_core <- any(high_conservation)
  core_length <- if (has_conserved_core) {
    # Find longest consecutive stretch of conserved positions
    max_consecutive <- 0
    current_consecutive <- 0
    for (i in seq_along(high_conservation)) {
      if (high_conservation[i]) {
        current_consecutive <- current_consecutive + 1
        max_consecutive <- max(max_consecutive, current_consecutive)
      } else {
        current_consecutive <- 0
      }
    }
    max_consecutive
  } else {
    0
  }
  
  return(list(
    conserved_positions = conserved_positions,
    n_conserved = length(conserved_positions),
    has_conserved_core = has_conserved_core,
    core_length = core_length,
    conservation_pattern = if (core_length >= 3) "strong_core"
                          else if (length(conserved_positions) >= 3) "distributed"
                          else "weak"
  ))
}

# Test biological relevance
test_biological_pattern <- function(pwm) {
  # Check for known CTCF-like patterns
  # Look for conserved G/C bases (common in CTCF motifs)
  gc_positions <- which(apply(pwm, 2, function(col) {
    (col["G"] + col["C"]) > 0.6
  }))
  
  # Check for palindromic patterns (CTCF often binds palindromic sequences)
  n_pos <- ncol(pwm)
  palindrome_score <- 0
  if (n_pos >= 4) {
    for (i in 1:(n_pos %/% 2)) {
      j <- n_pos - i + 1
      # Compare forward and reverse-complement positions
      forward <- pwm[, i]
      reverse_comp <- pwm[c("T", "G", "C", "A"), j]  # A<->T, C<->G
      similarity <- sum(abs(forward - reverse_comp))
      palindrome_score <- palindrome_score + (1 - similarity/2)
    }
    palindrome_score <- palindrome_score / (n_pos %/% 2)
  }
  
  return(list(
    gc_rich_positions = gc_positions,
    n_gc_rich = length(gc_positions),
    palindrome_score = palindrome_score,
    is_palindromic = palindrome_score > 0.7,
    biological_relevance = if (length(gc_positions) >= 3 && palindrome_score > 0.5) "high"
                          else if (length(gc_positions) >= 2 || palindrome_score > 0.3) "medium"
                          else "low"
  ))
}

# Test statistical significance
test_statistical_significance <- function(pwm, sequences = NULL) {
  if (is.null(sequences)) {
    return(list(
      significance_test = "skipped",
      reason = "no_sequences_provided"
    ))
  }
  
  # Calculate observed PWM score distribution
  observed_scores <- numeric(length(sequences))
  for (i in seq_along(sequences)) {
    seq_str <- as.character(sequences[i])
    if (nchar(seq_str) >= ncol(pwm)) {
      # Score the sequence at best position
      best_score <- -Inf
      for (pos in 1:(nchar(seq_str) - ncol(pwm) + 1)) {
        subseq <- substr(seq_str, pos, pos + ncol(pwm) - 1)
        score <- 0
        for (j in 1:ncol(pwm)) {
          base <- substr(subseq, j, j)
          if (base %in% c("A", "C", "G", "T")) {
            score <- score + log2(pwm[base, j] / 0.25)
          }
        }
        best_score <- max(best_score, score)
      }
      observed_scores[i] <- best_score
    }
  }
  
  # Simple significance test - compare to expected random score
  random_expected <- 0  # Expected score for random sequence
  random_sd <- 2  # Approximate standard deviation
  
  mean_observed <- mean(observed_scores, na.rm = TRUE)
  z_score <- (mean_observed - random_expected) / random_sd
  p_value <- 2 * pnorm(-abs(z_score))  # Two-tailed test
  
  return(list(
    mean_score = mean_observed,
    z_score = z_score,
    p_value = p_value,
    is_significant = p_value < 0.05,
    significance_level = if (p_value < 0.001) "highly_significant"
                        else if (p_value < 0.01) "significant"
                        else if (p_value < 0.05) "marginally_significant"
                        else "not_significant"
  ))
}

# Compile validation report
compile_validation_report <- function(tests) {
  # Calculate overall quality score
  ic_score <- if (tests$information_content$passes_minimum) 1 else 0
  conservation_score <- if (tests$conservation$has_conserved_core) 1 else 0
  biological_score <- if (tests$biological$biological_relevance != "low") 1 else 0
  statistical_score <- if (!is.null(tests$statistical$is_significant) && tests$statistical$is_significant) 1 else 0
  
  overall_score <- (ic_score + conservation_score + biological_score + statistical_score) / 4
  
  overall_quality <- if (overall_score >= 0.75) "excellent"
                    else if (overall_score >= 0.5) "good"
                    else if (overall_score >= 0.25) "fair"
                    else "poor"
  
  return(list(
    overall_quality = overall_quality,
    overall_score = overall_score,
    individual_tests = tests,
    recommendations = generate_recommendations(tests)
  ))
}

# Generate recommendations
generate_recommendations <- function(tests) {
  recommendations <- character(0)
  
  if (!tests$information_content$passes_minimum) {
    recommendations <- c(recommendations, "Increase sequence alignment quality or dataset size to improve information content")
  }
  
  if (!tests$conservation$has_conserved_core) {
    recommendations <- c(recommendations, "Consider better sequence alignment methods to identify conserved motif core")
  }
  
  if (tests$biological$biological_relevance == "low") {
    recommendations <- c(recommendations, "Verify that input sequences are genuine CTCF binding sites")
  }
  
  if (!is.null(tests$statistical$is_significant) && !tests$statistical$is_significant) {
    recommendations <- c(recommendations, "Increase dataset size or improve sequence quality for better statistical significance")
  }
  
  if (length(recommendations) == 0) {
    recommendations <- "PWM quality is good - no major improvements needed"
  }
  
  return(recommendations)
}

# Main validation function
validate_quality <- function(pwm, sequences = NULL) {
  tests <- list(
    information_content = test_information_content(pwm),
    conservation = test_conservation_pattern(pwm),
    biological_relevance = test_biological_pattern(pwm),
    statistical_significance = test_statistical_significance(pwm, sequences)
  )
  
  return(compile_validation_report(tests))
}

# --- Main Analysis ---
cat("=== PWM Quality Validation Report ===\n")
cat("PWM file:", pwm_file, "\n")
if (!is.null(sequences_file)) cat("Sequences file:", sequences_file, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Validate input files
if (!file.exists(pwm_file)) {
  stop("PWM file not found: ", pwm_file)
}

# Load PWM
cat("Loading PWM...\n")
pwm_data <- readRDS(pwm_file)

# Handle different PWM formats
if (is.list(pwm_data) && "pwm" %in% names(pwm_data)) {
  pwm <- pwm_data$pwm
} else if (is.matrix(pwm_data)) {
  pwm <- pwm_data
} else {
  stop("Unrecognized PWM format")
}

# Validate PWM format
if (!is.matrix(pwm)) {
  stop("PWM object is not a matrix")
}

if (nrow(pwm) != 4) {
  stop("PWM should have 4 rows (A, C, G, T)")
}

# Ensure row names are set
if (is.null(rownames(pwm))) {
  rownames(pwm) <- c("A", "C", "G", "T")
}

cat("PWM dimensions:", nrow(pwm), "x", ncol(pwm), "(bases x positions)\n")
cat("PWM length:", ncol(pwm), "bp\n\n")

# Load sequences if provided
sequences <- NULL
if (!is.null(sequences_file) && file.exists(sequences_file)) {
  cat("Loading sequences for statistical validation...\n")
  sequences <- readDNAStringSet(sequences_file)
  cat("Loaded", length(sequences), "sequences\n\n")
}

# Perform validation
cat("Performing comprehensive PWM validation...\n")
validation_results <- validate_quality(pwm, sequences)

# Display results
cat("\n=== VALIDATION RESULTS ===\n")
cat("Overall Quality:", toupper(validation_results$overall_quality), "\n")
cat("Overall Score:", round(validation_results$overall_score * 100, 1), "%\n\n")

# Information Content Results
ic_results <- validation_results$individual_tests$information_content
cat("1. INFORMATION CONTENT ANALYSIS\n")
cat("   Total IC:", round(ic_results$total_ic, 3), "bits\n")
cat("   Mean IC per position:", round(ic_results$mean_ic, 3), "bits\n")
cat("   Max IC:", round(ic_results$max_ic, 3), "bits\n")
cat("   Conserved positions (>1 bit):", ic_results$conserved_positions, "\n")
cat("   Quality level:", toupper(ic_results$quality_level), "\n")
cat("   Passes minimum threshold:", if (ic_results$passes_minimum) "✅ YES" else "❌ NO", "\n\n")

# Conservation Pattern Results
cons_results <- validation_results$individual_tests$conservation
cat("2. CONSERVATION PATTERN ANALYSIS\n")
cat("   Conserved positions:", length(cons_results$conserved_positions), "\n")
cat("   Has conserved core:", if (cons_results$has_conserved_core) "✅ YES" else "❌ NO", "\n")
cat("   Core length:", cons_results$core_length, "positions\n")
cat("   Pattern type:", toupper(cons_results$conservation_pattern), "\n\n")

# Biological Relevance Results
bio_results <- validation_results$individual_tests$biological_relevance
cat("3. BIOLOGICAL RELEVANCE ANALYSIS\n")
cat("   GC-rich positions:", bio_results$n_gc_rich, "\n")
cat("   Palindrome score:", round(bio_results$palindrome_score, 3), "\n")
cat("   Is palindromic:", if (bio_results$is_palindromic) "✅ YES" else "❌ NO", "\n")
cat("   Biological relevance:", toupper(bio_results$biological_relevance), "\n\n")

# Statistical Significance Results
stat_results <- validation_results$individual_tests$statistical_significance
cat("4. STATISTICAL SIGNIFICANCE\n")
if (stat_results$significance_test == "skipped") {
  cat("   Test: SKIPPED (no sequences provided)\n")
} else {
  cat("   Mean score:", round(stat_results$mean_score, 3), "\n")
  cat("   Z-score:", round(stat_results$z_score, 3), "\n")
  cat("   P-value:", format(stat_results$p_value, scientific = TRUE), "\n")
  cat("   Significance:", toupper(stat_results$significance_level), "\n")
  cat("   Is significant:", if (stat_results$is_significant) "✅ YES" else "❌ NO", "\n")
}

# Recommendations
cat("\n=== RECOMMENDATIONS ===\n")
for (i in seq_along(validation_results$recommendations)) {
  cat(i, ".", validation_results$recommendations[i], "\n")
}

# Save results if output file specified
if (!is.null(output_file)) {
  cat("\nSaving validation results to:", output_file, "\n")
  
  # Create a comprehensive report
  report_data <- list(
    pwm_file = pwm_file,
    sequences_file = sequences_file,
    analysis_time = as.character(Sys.time()),
    pwm_dimensions = paste(nrow(pwm), "x", ncol(pwm)),
    validation_results = validation_results
  )
  
  saveRDS(report_data, output_file)
}

cat("\n=== PWM Quality Validation Complete ===\n")
