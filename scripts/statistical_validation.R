# Statistical Validation Script - Statistical significance testing
# This script implements comprehensive statistical validation for PWM models

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Please install Biostrings: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
pwm_file <- if (length(args) >= 1) args[1] else "results/generated_pwm.rds"
sequences_file <- if (length(args) >= 2) args[2] else "data/training_sequences.fasta"
output_file <- if (length(args) >= 3) args[3] else "results/statistical_validation_report.txt"
n_permutations <- if (length(args) >= 4) as.numeric(args[4]) else 1000

cat("Statistical Validation of PWM Models\n")
cat("====================================\n")
cat("PWM file:", pwm_file, "\n")
cat("Sequences file:", sequences_file, "\n")
cat("Output file:", output_file, "\n")
cat("Permutations:", n_permutations, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Functions ---

# Load PWM model
load_pwm_model <- function(pwm_file) {
  cat("Loading PWM model from:", pwm_file, "\n")
  
  if (!file.exists(pwm_file)) {
    stop("PWM file not found: ", pwm_file)
  }
  
  pwm <- readRDS(pwm_file)
  
  # Validate PWM structure
  if (!is.matrix(pwm)) {
    stop("PWM must be a matrix")
  }
  
  if (nrow(pwm) != 4) {
    stop("PWM must have 4 rows (A, C, G, T)")
  }
  
  if (is.null(rownames(pwm))) {
    rownames(pwm) <- c("A", "C", "G", "T")
  }
  
  cat("PWM loaded successfully: ", nrow(pwm), "x", ncol(pwm), "\n")
  
  return(pwm)
}

# Calculate information content of PWM
calculate_pwm_information_content <- function(pwm) {
  total_ic <- 0
  position_ic <- numeric(ncol(pwm))
  
  for (pos in 1:ncol(pwm)) {
    # Get probabilities for this position
    probs <- pwm[, pos]
    
    # Remove zero probabilities for log calculation
    probs <- probs[probs > 0]
    
    if (length(probs) > 0) {
      # Calculate information content: IC = sum(p * log2(p/0.25))
      ic <- sum(probs * log2(probs / 0.25))
      position_ic[pos] <- ic
      total_ic <- total_ic + ic
    }
  }
  
  return(list(
    total_ic = total_ic,
    position_ic = position_ic,
    mean_ic = total_ic / ncol(pwm)
  ))
}

# Generate random sequences with same base composition
generate_random_sequences <- function(sequences, n_random = 1000) {
  cat("Generating", n_random, "random sequences...\n")
  
  # Calculate overall base composition
  all_bases <- unlist(strsplit(as.character(sequences), ""))
  base_composition <- table(all_bases[all_bases %in% c("A", "C", "G", "T")])
  base_probs <- base_composition / sum(base_composition)
  
  # Get sequence lengths
  seq_lengths <- width(sequences)
  
  # Generate random sequences
  random_seqs <- character(n_random)
  
  for (i in 1:n_random) {
    # Choose random length from original distribution
    target_length <- sample(seq_lengths, 1)
    
    # Generate random sequence with same base composition
    random_bases <- sample(names(base_probs), target_length, replace = TRUE, prob = base_probs)
    random_seqs[i] <- paste(random_bases, collapse = "")
  }
  
  return(DNAStringSet(random_seqs))
}

# Generate shuffled sequences (preserve individual compositions)
generate_shuffled_sequences <- function(sequences, n_shuffled = 1000) {
  cat("Generating", n_shuffled, "shuffled sequences...\n")
  
  shuffled_seqs <- character(n_shuffled)
  
  for (i in 1:n_shuffled) {
    # Pick a random sequence to shuffle
    original_seq <- sample(sequences, 1)
    seq_bases <- strsplit(as.character(original_seq), "")[[1]]
    
    # Shuffle the bases
    shuffled_bases <- sample(seq_bases)
    shuffled_seqs[i] <- paste(shuffled_bases, collapse = "")
  }
  
  return(DNAStringSet(shuffled_seqs))
}

# Score sequences against PWM
score_sequences_against_pwm <- function(sequences, pwm) {
  pwm_length <- ncol(pwm)
  scores <- numeric(length(sequences))
  
  for (i in 1:length(sequences)) {
    seq_char <- as.character(sequences[i])
    seq_length <- nchar(seq_char)
    
    if (seq_length < pwm_length) {
      scores[i] <- -Inf
      next
    }
    
    # Find best scoring position
    best_score <- -Inf
    
    for (pos in 1:(seq_length - pwm_length + 1)) {
      subseq <- substr(seq_char, pos, pos + pwm_length - 1)
      bases <- strsplit(subseq, "")[[1]]
      
      # Calculate PWM score
      position_score <- 0
      valid_bases <- TRUE
      
      for (j in 1:length(bases)) {
        if (bases[j] %in% c("A", "C", "G", "T")) {
          prob <- pwm[bases[j], j]
          if (prob > 0) {
            position_score <- position_score + log2(prob / 0.25)
          } else {
            valid_bases <- FALSE
            break
          }
        } else {
          valid_bases <- FALSE
          break
        }
      }
      
      if (valid_bases && position_score > best_score) {
        best_score <- position_score
      }
    }
    
    scores[i] <- best_score
  }
  
  return(scores)
}

# Perform permutation test
perform_permutation_test <- function(real_scores, null_scores) {
  cat("Performing permutation test...\n")
  
  # Calculate test statistics
  real_mean <- mean(real_scores[is.finite(real_scores)])
  null_mean <- mean(null_scores[is.finite(null_scores)])
  
  observed_difference <- real_mean - null_mean
  
  # Permutation test
  all_scores <- c(real_scores[is.finite(real_scores)], null_scores[is.finite(null_scores)])
  n_real <- length(real_scores[is.finite(real_scores)])
  n_total <- length(all_scores)
  
  if (n_total == 0) {
    return(list(p_value = 1.0, effect_size = 0, observed_difference = 0))
  }
  
  permutation_differences <- numeric(n_permutations)
  
  for (i in 1:n_permutations) {
    # Randomly assign labels
    shuffled_indices <- sample(n_total, n_real)
    perm_real <- all_scores[shuffled_indices]
    perm_null <- all_scores[-shuffled_indices]
    
    perm_difference <- mean(perm_real) - mean(perm_null)
    permutation_differences[i] <- perm_difference
  }
  
  # Calculate p-value (two-tailed)
  p_value <- sum(abs(permutation_differences) >= abs(observed_difference)) / n_permutations
  
  # Calculate effect size (Cohen's d)
  pooled_sd <- sqrt(((length(real_scores[is.finite(real_scores)]) - 1) * var(real_scores[is.finite(real_scores)]) + 
                     (length(null_scores[is.finite(null_scores)]) - 1) * var(null_scores[is.finite(null_scores)])) / 
                    (length(real_scores[is.finite(real_scores)]) + length(null_scores[is.finite(null_scores)]) - 2))
  
  effect_size <- observed_difference / pooled_sd
  
  return(list(
    p_value = p_value,
    effect_size = effect_size,
    observed_difference = observed_difference,
    real_mean = real_mean,
    null_mean = null_mean,
    permutation_differences = permutation_differences
  ))
}

# Perform Kolmogorov-Smirnov test
perform_ks_test <- function(real_scores, null_scores) {
  cat("Performing Kolmogorov-Smirnov test...\n")
  
  # Remove infinite values
  real_finite <- real_scores[is.finite(real_scores)]
  null_finite <- null_scores[is.finite(null_scores)]
  
  if (length(real_finite) == 0 || length(null_finite) == 0) {
    return(list(p_value = 1.0, statistic = 0))
  }
  
  # Perform KS test
  ks_result <- ks.test(real_finite, null_finite)
  
  return(list(
    p_value = ks_result$p.value,
    statistic = ks_result$statistic,
    method = ks_result$method
  ))
}

# Calculate bootstrap confidence intervals
calculate_bootstrap_ci <- function(scores, n_bootstrap = 1000, alpha = 0.05) {
  cat("Calculating bootstrap confidence intervals...\n")
  
  finite_scores <- scores[is.finite(scores)]
  
  if (length(finite_scores) == 0) {
    return(list(lower = NA, upper = NA, mean = NA))
  }
  
  bootstrap_means <- replicate(n_bootstrap, {
    sample_scores <- sample(finite_scores, length(finite_scores), replace = TRUE)
    mean(sample_scores)
  })
  
  ci_lower <- quantile(bootstrap_means, alpha/2)
  ci_upper <- quantile(bootstrap_means, 1 - alpha/2)
  
  return(list(
    lower = ci_lower,
    upper = ci_upper,
    mean = mean(finite_scores),
    bootstrap_means = bootstrap_means
  ))
}

# Generate comprehensive statistics
generate_comprehensive_statistics <- function(pwm, sequences, real_scores, random_scores, shuffled_scores) {
  cat("Generating comprehensive statistics...\n")
  
  # PWM information content
  ic_stats <- calculate_pwm_information_content(pwm)
  
  # Score distributions
  real_finite <- real_scores[is.finite(real_scores)]
  random_finite <- random_scores[is.finite(random_scores)]
  shuffled_finite <- shuffled_scores[is.finite(shuffled_scores)]
  
  # Permutation tests
  random_perm_test <- perform_permutation_test(real_scores, random_scores)
  shuffled_perm_test <- perform_permutation_test(real_scores, shuffled_scores)
  
  # KS tests
  random_ks_test <- perform_ks_test(real_scores, random_scores)
  shuffled_ks_test <- perform_ks_test(real_scores, shuffled_scores)
  
  # Bootstrap confidence intervals
  real_ci <- calculate_bootstrap_ci(real_scores)
  random_ci <- calculate_bootstrap_ci(random_scores)
  shuffled_ci <- calculate_bootstrap_ci(shuffled_scores)
  
  return(list(
    pwm_stats = list(
      dimensions = paste(nrow(pwm), "x", ncol(pwm)),
      total_ic = ic_stats$total_ic,
      mean_ic = ic_stats$mean_ic,
      position_ic = ic_stats$position_ic
    ),
    score_stats = list(
      real = list(
        n = length(real_finite),
        mean = mean(real_finite),
        median = median(real_finite),
        sd = sd(real_finite),
        min = min(real_finite),
        max = max(real_finite)
      ),
      random = list(
        n = length(random_finite),
        mean = mean(random_finite),
        median = median(random_finite),
        sd = sd(random_finite),
        min = min(random_finite),
        max = max(random_finite)
      ),
      shuffled = list(
        n = length(shuffled_finite),
        mean = mean(shuffled_finite),
        median = median(shuffled_finite),
        sd = sd(shuffled_finite),
        min = min(shuffled_finite),
        max = max(shuffled_finite)
      )
    ),
    statistical_tests = list(
      random_permutation = random_perm_test,
      shuffled_permutation = shuffled_perm_test,
      random_ks = random_ks_test,
      shuffled_ks = shuffled_ks_test
    ),
    confidence_intervals = list(
      real = real_ci,
      random = random_ci,
      shuffled = shuffled_ci
    )
  ))
}

# --- Main Processing ---

# Load PWM and sequences
pwm <- load_pwm_model(pwm_file)

if (!file.exists(sequences_file)) {
  stop("Sequences file not found: ", sequences_file)
}

cat("Loading sequences from:", sequences_file, "\n")
sequences <- readDNAStringSet(sequences_file)
cat("Sequences loaded:", length(sequences), "\n\n")

if (length(sequences) == 0) {
  stop("No sequences found in sequences file")
}

# Generate null models
cat("Generating null model sequences...\n")
random_sequences <- generate_random_sequences(sequences, n_permutations)
shuffled_sequences <- generate_shuffled_sequences(sequences, n_permutations)

# Score all sequence sets against PWM
cat("Scoring sequences against PWM...\n")
real_scores <- score_sequences_against_pwm(sequences, pwm)
cat("Real sequences scored\n")

random_scores <- score_sequences_against_pwm(random_sequences, pwm)
cat("Random sequences scored\n")

shuffled_scores <- score_sequences_against_pwm(shuffled_sequences, pwm)
cat("Shuffled sequences scored\n")

# Generate comprehensive statistics
stats <- generate_comprehensive_statistics(pwm, sequences, real_scores, random_scores, shuffled_scores)

# Create output directory
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Generate report
cat("\nGenerating statistical validation report...\n")

sink(output_file)
cat("Statistical Validation Report\n")
cat("=============================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("PWM file:", pwm_file, "\n")
cat("Sequences file:", sequences_file, "\n")
cat("Permutations:", n_permutations, "\n\n")

cat("PWM CHARACTERISTICS\n")
cat("===================\n")
cat("Dimensions:", stats$pwm_stats$dimensions, "\n")
cat("Total Information Content:", round(stats$pwm_stats$total_ic, 3), "bits\n")
cat("Mean Information Content:", round(stats$pwm_stats$mean_ic, 3), "bits/position\n\n")

cat("SCORE DISTRIBUTIONS\n")
cat("===================\n")
cat("Real Sequences (n =", stats$score_stats$real$n, "):\n")
cat("  Mean:", round(stats$score_stats$real$mean, 3), "\n")
cat("  Median:", round(stats$score_stats$real$median, 3), "\n")
cat("  SD:", round(stats$score_stats$real$sd, 3), "\n")
cat("  Range:", round(stats$score_stats$real$min, 3), "to", round(stats$score_stats$real$max, 3), "\n\n")

cat("Random Sequences (n =", stats$score_stats$random$n, "):\n")
cat("  Mean:", round(stats$score_stats$random$mean, 3), "\n")
cat("  Median:", round(stats$score_stats$random$median, 3), "\n")
cat("  SD:", round(stats$score_stats$random$sd, 3), "\n")
cat("  Range:", round(stats$score_stats$random$min, 3), "to", round(stats$score_stats$random$max, 3), "\n\n")

cat("Shuffled Sequences (n =", stats$score_stats$shuffled$n, "):\n")
cat("  Mean:", round(stats$score_stats$shuffled$mean, 3), "\n")
cat("  Median:", round(stats$score_stats$shuffled$median, 3), "\n")
cat("  SD:", round(stats$score_stats$shuffled$sd, 3), "\n")
cat("  Range:", round(stats$score_stats$shuffled$min, 3), "to", round(stats$score_stats$shuffled$max, 3), "\n\n")

cat("STATISTICAL TESTS\n")
cat("=================\n")
cat("Real vs Random Sequences:\n")
cat("  Permutation test p-value:", sprintf("%.2e", stats$statistical_tests$random_permutation$p_value), "\n")
cat("  Effect size (Cohen's d):", round(stats$statistical_tests$random_permutation$effect_size, 3), "\n")
cat("  Mean difference:", round(stats$statistical_tests$random_permutation$observed_difference, 3), "\n")
cat("  KS test p-value:", sprintf("%.2e", stats$statistical_tests$random_ks$p_value), "\n\n")

cat("Real vs Shuffled Sequences:\n")
cat("  Permutation test p-value:", sprintf("%.2e", stats$statistical_tests$shuffled_permutation$p_value), "\n")
cat("  Effect size (Cohen's d):", round(stats$statistical_tests$shuffled_permutation$effect_size, 3), "\n")
cat("  Mean difference:", round(stats$statistical_tests$shuffled_permutation$observed_difference, 3), "\n")
cat("  KS test p-value:", sprintf("%.2e", stats$statistical_tests$shuffled_ks$p_value), "\n\n")

cat("CONFIDENCE INTERVALS (95%)\n")
cat("==========================\n")
cat("Real sequences:", round(stats$confidence_intervals$real$lower, 3), "to", 
    round(stats$confidence_intervals$real$upper, 3), "\n")
cat("Random sequences:", round(stats$confidence_intervals$random$lower, 3), "to", 
    round(stats$confidence_intervals$random$upper, 3), "\n")
cat("Shuffled sequences:", round(stats$confidence_intervals$shuffled$lower, 3), "to", 
    round(stats$confidence_intervals$shuffled$upper, 3), "\n\n")

cat("INTERPRETATION\n")
cat("==============\n")
significance_level <- 0.05

if (stats$statistical_tests$random_permutation$p_value < significance_level) {
  cat("✓ PWM shows SIGNIFICANT discrimination against random sequences\n")
} else {
  cat("✗ PWM does NOT show significant discrimination against random sequences\n")
}

if (stats$statistical_tests$shuffled_permutation$p_value < significance_level) {
  cat("✓ PWM shows SIGNIFICANT discrimination against shuffled sequences\n")
} else {
  cat("✗ PWM does NOT show significant discrimination against shuffled sequences\n")
}

# Effect size interpretation
if (abs(stats$statistical_tests$random_permutation$effect_size) > 0.8) {
  cat("✓ LARGE effect size against random sequences\n")
} else if (abs(stats$statistical_tests$random_permutation$effect_size) > 0.5) {
  cat("✓ MODERATE effect size against random sequences\n")
} else if (abs(stats$statistical_tests$random_permutation$effect_size) > 0.2) {
  cat("✓ SMALL effect size against random sequences\n")
} else {
  cat("✗ NEGLIGIBLE effect size against random sequences\n")
}

if (stats$pwm_stats$total_ic > 10) {
  cat("✓ HIGH information content PWM (", round(stats$pwm_stats$total_ic, 1), "bits)\n")
} else if (stats$pwm_stats$total_ic > 5) {
  cat("✓ MODERATE information content PWM (", round(stats$pwm_stats$total_ic, 1), "bits)\n")
} else {
  cat("✗ LOW information content PWM (", round(stats$pwm_stats$total_ic, 1), "bits)\n")
}

sink()

# Save detailed results as RDS
results_file <- gsub("\\.txt$", "_results.rds", output_file)
cat("Saving detailed results to:", results_file, "\n")

validation_results <- list(
  pwm = pwm,
  sequences_file = sequences_file,
  n_sequences = length(sequences),
  n_permutations = n_permutations,
  scores = list(
    real = real_scores,
    random = random_scores,
    shuffled = shuffled_scores
  ),
  statistics = stats,
  timestamp = Sys.time()
)

saveRDS(validation_results, results_file)

cat("\nStatistical validation completed successfully!\n")
cat("Report saved to:", output_file, "\n")
cat("Detailed results saved to:", results_file, "\n")
