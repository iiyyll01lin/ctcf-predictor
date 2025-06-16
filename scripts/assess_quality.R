#!/usr/bin/env Rscript

# CTCF PWM Testing Pipeline - Quality Assessment Engine
# Comprehensive quality metrics for PWM evaluation
# Author: CTCF Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(Biostrings)
  library(seqLogo)
})

# Logging function
log_message <- function(level, message, script = "assess_quality.R") {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(sprintf("%s %s [%s] %s\n", timestamp, level, script, message))
}

# Main quality assessment function
assess_pwm_quality <- function(pwm_file, sequences_file = NULL, 
                              output_dir = "results/quality_assessment",
                              generate_plots = TRUE) {
  
  log_message("INFO", "Starting comprehensive PWM quality assessment")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    log_message("INFO", paste("Created output directory:", output_dir))
  }
  
  # Load PWM
  log_message("INFO", paste("Loading PWM from:", pwm_file))
  pwm_data <- load_pwm_data(pwm_file)
  
  # Load sequences if provided
  sequences <- NULL
  if (!is.null(sequences_file) && file.exists(sequences_file)) {
    log_message("INFO", paste("Loading sequences from:", sequences_file))
    sequences <- readDNAStringSet(sequences_file)
    log_message("INFO", paste("Loaded", length(sequences), "sequences"))
  }
  
  # Perform comprehensive quality assessment
  quality_results <- list()
  
  # 1. Information Content Analysis
  log_message("INFO", "Analyzing information content")
  quality_results$information_content <- analyze_information_content(pwm_data)
  
  # 2. Conservation Pattern Analysis
  log_message("INFO", "Analyzing conservation patterns")
  quality_results$conservation <- analyze_conservation_patterns(pwm_data)
  
  # 3. Motif Structure Analysis
  log_message("INFO", "Analyzing motif structure")
  quality_results$structure <- analyze_motif_structure(pwm_data)
  
  # 4. Statistical Properties
  log_message("INFO", "Analyzing statistical properties")
  quality_results$statistics <- analyze_statistical_properties(pwm_data)
  
  # 5. Biological Relevance (if sequences provided)
  if (!is.null(sequences)) {
    log_message("INFO", "Analyzing biological relevance")
    quality_results$biological <- analyze_biological_relevance(pwm_data, sequences)
  }
  
  # 6. Overall Quality Scoring
  log_message("INFO", "Computing overall quality score")
  quality_results$overall <- compute_overall_quality_score(quality_results)
  
  # Generate visualizations
  if (generate_plots) {
    log_message("INFO", "Generating quality visualizations")
    generate_quality_plots(pwm_data, quality_results, output_dir)
  }
  
  # Generate comprehensive report
  generate_quality_report(pwm_data, quality_results, output_dir)
  
  # Save results
  results_file <- file.path(output_dir, "quality_assessment_results.rds")
  saveRDS(quality_results, results_file)
  
  log_message("INFO", paste("Quality assessment completed. Results saved to:", results_file))
  return(quality_results)
}

# Load PWM data
load_pwm_data <- function(pwm_file) {
  if (grepl("\\.rds$", pwm_file, ignore.case = TRUE)) {
    pwm_data <- readRDS(pwm_file)
    if (is.list(pwm_data)) {
      return(pwm_data)
    } else if (is.matrix(pwm_data)) {
      return(list(pwm = pwm_data))
    }
  } else {
    pwm_matrix <- as.matrix(read.table(pwm_file, header = TRUE, row.names = 1))
    return(list(pwm = pwm_matrix))
  }
}

# Analyze information content
analyze_information_content <- function(pwm_data) {
  pwm <- pwm_data$pwm
  
  # Position-wise information content
  position_ic <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col))
  })
  
  # Total information content
  total_ic <- sum(position_ic)
  
  # Information content statistics
  ic_stats <- list(
    position_ic = position_ic,
    total_ic = total_ic,
    mean_ic = mean(position_ic),
    max_ic = max(position_ic),
    min_ic = min(position_ic),
    ic_variance = var(position_ic),
    high_ic_positions = which(position_ic > 1.0),
    moderate_ic_positions = which(position_ic > 0.5 & position_ic <= 1.0),
    low_ic_positions = which(position_ic <= 0.5)
  )
  
  # Quality scoring
  ic_stats$quality_score <- calculate_ic_quality_score(ic_stats)
  
  return(ic_stats)
}

# Calculate IC quality score
calculate_ic_quality_score <- function(ic_stats) {
  # Base score from total IC
  base_score <- min(100, ic_stats$total_ic * 5)
  
  # Bonus for high-IC positions
  high_ic_bonus <- length(ic_stats$high_ic_positions) * 5
  
  # Penalty for too many low-IC positions
  low_ic_penalty <- length(ic_stats$low_ic_positions) * 2
  
  # Uniformity bonus (prefer some variation but not too much)
  uniformity_score <- 1 / (1 + ic_stats$ic_variance)
  uniformity_bonus <- uniformity_score * 10
  
  final_score <- max(0, min(100, base_score + high_ic_bonus - low_ic_penalty + uniformity_bonus))
  
  return(final_score)
}

# Analyze conservation patterns
analyze_conservation_patterns <- function(pwm_data) {
  pwm <- pwm_data$pwm
  
  # Position-wise conservation (entropy-based)
  conservation_scores <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    entropy <- -sum(col * log2(col))
    conservation <- 2 - entropy  # Max entropy is 2 for 4 bases
    return(conservation)
  })
  
  # Identify conserved regions
  conserved_positions <- which(conservation_scores > 1.0)
  highly_conserved <- which(conservation_scores > 1.5)
  
  # Find contiguous conserved regions
  conserved_regions <- find_conserved_regions(conserved_positions)
  
  # Dominant nucleotides per position
  dominant_bases <- apply(pwm, 2, function(col) {
    rownames(pwm)[which.max(col)]
  })
  
  dominant_strengths <- apply(pwm, 2, max)
  
  conservation_analysis <- list(
    conservation_scores = conservation_scores,
    mean_conservation = mean(conservation_scores),
    max_conservation = max(conservation_scores),
    conserved_positions = conserved_positions,
    highly_conserved_positions = highly_conserved,
    conserved_regions = conserved_regions,
    dominant_bases = dominant_bases,
    dominant_strengths = dominant_strengths,
    conservation_pattern = classify_conservation_pattern(conservation_scores)
  )
  
  # Quality scoring
  conservation_analysis$quality_score <- calculate_conservation_quality_score(conservation_analysis)
  
  return(conservation_analysis)
}

# Find conserved regions
find_conserved_regions <- function(positions) {
  if (length(positions) == 0) return(list())
  
  regions <- list()
  current_start <- positions[1]
  current_end <- positions[1]
  
  for (i in 2:length(positions)) {
    if (positions[i] == current_end + 1) {
      current_end <- positions[i]
    } else {
      regions[[length(regions) + 1]] <- list(
        start = current_start,
        end = current_end,
        length = current_end - current_start + 1
      )
      current_start <- positions[i]
      current_end <- positions[i]
    }
  }
  
  # Add final region
  regions[[length(regions) + 1]] <- list(
    start = current_start,
    end = current_end,
    length = current_end - current_start + 1
  )
  
  return(regions)
}

# Classify conservation pattern
classify_conservation_pattern <- function(conservation_scores) {
  high_cons <- sum(conservation_scores > 1.0)
  mod_cons <- sum(conservation_scores > 0.5 & conservation_scores <= 1.0)
  low_cons <- sum(conservation_scores <= 0.5)
  
  total_positions <- length(conservation_scores)
  
  if (high_cons / total_positions > 0.5) {
    return("HIGHLY_CONSERVED")
  } else if (high_cons / total_positions > 0.3) {
    return("MODERATELY_CONSERVED")
  } else if (mod_cons / total_positions > 0.5) {
    return("WEAKLY_CONSERVED")
  } else {
    return("POORLY_CONSERVED")
  }
}

# Calculate conservation quality score
calculate_conservation_quality_score <- function(conservation_analysis) {
  # Base score from mean conservation
  base_score <- min(50, conservation_analysis$mean_conservation * 25)
  
  # Bonus for conserved positions
  conserved_bonus <- length(conservation_analysis$conserved_positions) * 3
  
  # Bonus for conserved regions
  region_bonus <- sum(sapply(conservation_analysis$conserved_regions, function(r) r$length)) * 2
  
  # Pattern bonus
  pattern_bonus <- switch(conservation_analysis$conservation_pattern,
                         "HIGHLY_CONSERVED" = 20,
                         "MODERATELY_CONSERVED" = 15,
                         "WEAKLY_CONSERVED" = 10,
                         "POORLY_CONSERVED" = 0)
  
  final_score <- max(0, min(100, base_score + conserved_bonus + region_bonus + pattern_bonus))
  
  return(final_score)
}

# Analyze motif structure
analyze_motif_structure <- function(pwm_data) {
  pwm <- pwm_data$pwm
  
  # Core motif identification
  ic_values <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col))
  })
  
  # Identify core positions (high IC)
  core_positions <- which(ic_values > quantile(ic_values, 0.7))
  flanking_positions <- setdiff(1:ncol(pwm), core_positions)
  
  # GC content analysis
  gc_content_per_position <- apply(pwm, 2, function(col) {
    sum(col[c("C", "G")])
  })
  
  # Palindrome analysis
  palindrome_score <- calculate_palindrome_score(pwm)
  
  # Repeat analysis
  repeat_score <- calculate_repeat_score(pwm)
  
  structure_analysis <- list(
    motif_length = ncol(pwm),
    core_positions = core_positions,
    flanking_positions = flanking_positions,
    core_length = length(core_positions),
    gc_content_positions = gc_content_per_position,
    mean_gc_content = mean(gc_content_per_position),
    palindrome_score = palindrome_score,
    repeat_score = repeat_score,
    structure_type = classify_motif_structure(core_positions, ncol(pwm))
  )
  
  # Quality scoring
  structure_analysis$quality_score <- calculate_structure_quality_score(structure_analysis)
  
  return(structure_analysis)
}

# Calculate palindrome score
calculate_palindrome_score <- function(pwm) {
  n_pos <- ncol(pwm)
  palindrome_score <- 0
  
  for (i in 1:floor(n_pos/2)) {
    left_pos <- i
    right_pos <- n_pos - i + 1
    
    # Compare A-T and C-G complementarity
    at_score <- pwm["A", left_pos] * pwm["T", right_pos] + pwm["T", left_pos] * pwm["A", right_pos]
    gc_score <- pwm["G", left_pos] * pwm["C", right_pos] + pwm["C", left_pos] * pwm["G", right_pos]
    
    palindrome_score <- palindrome_score + at_score + gc_score
  }
  
  return(palindrome_score / floor(n_pos/2))
}

# Calculate repeat score
calculate_repeat_score <- function(pwm) {
  # Simple repeat detection based on column similarity
  n_pos <- ncol(pwm)
  if (n_pos < 4) return(0)
  
  repeat_scores <- numeric(n_pos - 1)
  
  for (i in 1:(n_pos - 1)) {
    # Calculate correlation between adjacent positions
    repeat_scores[i] <- cor(pwm[, i], pwm[, i + 1])
  }
  
  return(mean(repeat_scores, na.rm = TRUE))
}

# Classify motif structure
classify_motif_structure <- function(core_positions, total_length) {
  core_length <- length(core_positions)
  core_ratio <- core_length / total_length
  
  if (core_ratio > 0.8) {
    return("UNIFORM_HIGH_IC")
  } else if (core_ratio > 0.5) {
    return("CORE_DOMINANT")
  } else if (core_length >= 3 && max(diff(core_positions)) == 1) {
    return("CONTIGUOUS_CORE")
  } else if (core_length >= 2) {
    return("DISPERSED_CORE")
  } else {
    return("WEAK_STRUCTURE")
  }
}

# Calculate structure quality score
calculate_structure_quality_score <- function(structure_analysis) {
  # Base score from structure type
  structure_score <- switch(structure_analysis$structure_type,
                           "UNIFORM_HIGH_IC" = 40,
                           "CORE_DOMINANT" = 35,
                           "CONTIGUOUS_CORE" = 30,
                           "DISPERSED_CORE" = 20,
                           "WEAK_STRUCTURE" = 10)
  
  # Core length bonus
  core_bonus <- min(20, structure_analysis$core_length * 3)
  
  # GC content bonus (balanced GC content is good)
  gc_balance_score <- 1 - abs(structure_analysis$mean_gc_content - 0.5)
  gc_bonus <- gc_balance_score * 20
  
  # Palindrome bonus
  palindrome_bonus <- structure_analysis$palindrome_score * 20
  
  final_score <- max(0, min(100, structure_score + core_bonus + gc_bonus + palindrome_bonus))
  
  return(final_score)
}

# Analyze statistical properties
analyze_statistical_properties <- function(pwm_data) {
  pwm <- pwm_data$pwm
  
  # Matrix properties
  row_sums <- rowSums(pwm)
  col_sums <- colSums(pwm)
  
  # Uniformity tests
  row_uniformity <- 1 - var(row_sums) / mean(row_sums)^2
  col_uniformity <- 1 - var(col_sums) / mean(col_sums)^2
  
  # Entropy measures
  total_entropy <- -sum(pwm * log2(pwm + 1e-10))
  position_entropies <- apply(pwm, 2, function(col) {
    -sum(col * log2(col + 1e-10))
  })
  
  statistical_props <- list(
    row_sums = row_sums,
    col_sums = col_sums,
    row_uniformity = row_uniformity,
    col_uniformity = col_uniformity,
    total_entropy = total_entropy,
    position_entropies = position_entropies,
    mean_position_entropy = mean(position_entropies),
    matrix_rank = qr(pwm)$rank,
    condition_number = kappa(pwm),
    determinant = det(pwm %*% t(pwm))
  )
  
  # Quality scoring
  statistical_props$quality_score <- calculate_statistical_quality_score(statistical_props)
  
  return(statistical_props)
}

# Calculate statistical quality score
calculate_statistical_quality_score <- function(statistical_props) {
  # Column sum uniformity (should be close to 1)
  col_uniformity_score <- statistical_props$col_uniformity * 30
  
  # Entropy balance
  entropy_balance <- 1 - var(statistical_props$position_entropies) / mean(statistical_props$position_entropies)^2
  entropy_score <- entropy_balance * 30
  
  # Matrix properties
  rank_score <- min(20, statistical_props$matrix_rank * 5)
  
  # Condition number (lower is better, but not too low)
  condition_score <- max(0, 20 - log10(statistical_props$condition_number + 1))
  
  final_score <- max(0, min(100, col_uniformity_score + entropy_score + rank_score + condition_score))
  
  return(final_score)
}

# Analyze biological relevance
analyze_biological_relevance <- function(pwm_data, sequences) {
  pwm <- pwm_data$pwm
  motif_length <- ncol(pwm)
  
  # Score sequences against PWM
  sequence_scores <- score_sequences_against_pwm(sequences, pwm)
  
  # Score distribution analysis
  score_stats <- list(
    mean_score = mean(sequence_scores, na.rm = TRUE),
    median_score = median(sequence_scores, na.rm = TRUE),
    max_score = max(sequence_scores, na.rm = TRUE),
    min_score = min(sequence_scores, na.rm = TRUE),
    score_variance = var(sequence_scores, na.rm = TRUE)
  )
  
  # High-scoring sequence analysis
  high_score_threshold <- quantile(sequence_scores, 0.9, na.rm = TRUE)
  high_scoring_sequences <- which(sequence_scores >= high_score_threshold)
  
  biological_analysis <- list(
    sequence_scores = sequence_scores,
    score_statistics = score_stats,
    high_scoring_sequences = high_scoring_sequences,
    n_high_scoring = length(high_scoring_sequences),
    score_distribution_type = classify_score_distribution(sequence_scores)
  )
  
  # Quality scoring
  biological_analysis$quality_score <- calculate_biological_quality_score(biological_analysis)
  
  return(biological_analysis)
}

# Score sequences against PWM
score_sequences_against_pwm <- function(sequences, pwm) {
  motif_length <- ncol(pwm)
  scores <- numeric(length(sequences))
  
  for (i in 1:length(sequences)) {
    seq_str <- as.character(sequences[i])
    seq_len <- nchar(seq_str)
    
    if (seq_len >= motif_length) {
      best_score <- -Inf
      
      for (pos in 1:(seq_len - motif_length + 1)) {
        subseq <- substr(seq_str, pos, pos + motif_length - 1)
        score <- 0
        
        for (j in 1:nchar(subseq)) {
          base <- substr(subseq, j, j)
          if (base %in% rownames(pwm)) {
            score <- score + log2(pwm[base, j] / 0.25)
          }
        }
        
        if (score > best_score) {
          best_score <- score
        }
      }
      
      scores[i] <- best_score
    } else {
      scores[i] <- NA
    }
  }
  
  return(scores)
}

# Classify score distribution
classify_score_distribution <- function(scores) {
  scores <- scores[!is.na(scores)]
  
  if (length(scores) == 0) return("UNDEFINED")
  
  # Normality test
  if (length(scores) > 3) {
    shapiro_test <- shapiro.test(sample(scores, min(5000, length(scores))))
    is_normal <- shapiro_test$p.value > 0.05
  } else {
    is_normal <- FALSE
  }
  
  # Skewness
  mean_score <- mean(scores)
  median_score <- median(scores)
  skewness <- (mean_score - median_score) / sd(scores)
  
  if (is_normal) {
    return("NORMAL")
  } else if (skewness > 0.5) {
    return("RIGHT_SKEWED")
  } else if (skewness < -0.5) {
    return("LEFT_SKEWED")
  } else {
    return("NON_NORMAL")
  }
}

# Calculate biological quality score
calculate_biological_quality_score <- function(biological_analysis) {
  # Mean score bonus
  mean_score_bonus <- min(30, biological_analysis$score_statistics$mean_score * 3)
  
  # Score variance (moderate variance is good)
  variance_score <- 1 / (1 + biological_analysis$score_statistics$score_variance)
  variance_bonus <- variance_score * 20
  
  # High-scoring sequences
  high_score_ratio <- biological_analysis$n_high_scoring / length(biological_analysis$sequence_scores)
  high_score_bonus <- min(25, high_score_ratio * 100)
  
  # Distribution bonus
  distribution_bonus <- switch(biological_analysis$score_distribution_type,
                              "NORMAL" = 15,
                              "RIGHT_SKEWED" = 10,
                              "LEFT_SKEWED" = 5,
                              "NON_NORMAL" = 5,
                              "UNDEFINED" = 0)
  
  final_score <- max(0, min(100, mean_score_bonus + variance_bonus + high_score_bonus + distribution_bonus))
  
  return(final_score)
}

# Compute overall quality score
compute_overall_quality_score <- function(quality_results) {
  # Weight different aspects
  weights <- list(
    information_content = 0.3,
    conservation = 0.25,
    structure = 0.2,
    statistics = 0.15,
    biological = 0.1
  )
  
  # Adjust weights if biological analysis is not available
  if (!"biological" %in% names(quality_results)) {
    weights$biological <- 0
    # Redistribute weight
    weights$information_content <- 0.35
    weights$conservation <- 0.3
    weights$structure <- 0.2
    weights$statistics <- 0.15
  }
  
  # Calculate weighted score
  weighted_score <- 0
  component_scores <- list()
  
  for (component in names(weights)) {
    if (component %in% names(quality_results) && weights[[component]] > 0) {
      component_score <- quality_results[[component]]$quality_score
      weighted_score <- weighted_score + weights[[component]] * component_score
      component_scores[[component]] <- component_score
    }
  }
  
  # Overall quality grade
  grade <- if (weighted_score >= 90) {
    "EXCELLENT"
  } else if (weighted_score >= 80) {
    "VERY_GOOD"
  } else if (weighted_score >= 70) {
    "GOOD"
  } else if (weighted_score >= 60) {
    "FAIR"
  } else if (weighted_score >= 50) {
    "POOR"
  } else {
    "VERY_POOR"
  }
  
  return(list(
    overall_score = weighted_score,
    grade = grade,
    component_scores = component_scores,
    weights = weights
  ))
}

# Generate quality plots
generate_quality_plots <- function(pwm_data, quality_results, output_dir) {
  # This would generate various plots - simplified for now
  log_message("INFO", "Generating sequence logo")
  
  # Create sequence logo (if seqLogo package is available)
  tryCatch({
    png(file.path(output_dir, "sequence_logo.png"), width = 800, height = 400)
    seqLogo(pwm_data$pwm)
    dev.off()
  }, error = function(e) {
    log_message("WARN", "Could not generate sequence logo")
  })
  
  # Information content plot
  png(file.path(output_dir, "information_content.png"), width = 800, height = 400)
  barplot(quality_results$information_content$position_ic,
          main = "Position-wise Information Content",
          xlab = "Position", ylab = "Information Content (bits)",
          col = "steelblue")
  dev.off()
  
  log_message("INFO", "Quality plots generated")
}

# Generate quality report
generate_quality_report <- function(pwm_data, quality_results, output_dir) {
  report_file <- file.path(output_dir, "quality_assessment_report.txt")
  
  sink(report_file)
  cat("CTCF PWM Quality Assessment Report\n")
  cat("=================================\n\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # PWM Information
  cat("PWM Information:\n")
  cat("  Dimensions:", nrow(pwm_data$pwm), "x", ncol(pwm_data$pwm), "\n\n")
  
  # Overall Quality
  overall <- quality_results$overall
  cat("OVERALL QUALITY ASSESSMENT:\n")
  cat("  Overall Score:", round(overall$overall_score, 1), "/100\n")
  cat("  Quality Grade:", overall$grade, "\n\n")
  
  # Component Scores
  cat("Component Scores:\n")
  for (component in names(overall$component_scores)) {
    score <- overall$component_scores[[component]]
    weight <- overall$weights[[component]]
    cat("  ", gsub("_", " ", toupper(component)), ":", round(score, 1), 
        "/100 (weight:", round(weight, 2), ")\n")
  }
  cat("\n")
  
  # Detailed Analysis
  cat("DETAILED ANALYSIS:\n")
  cat("==================\n\n")
  
  # Information Content
  ic <- quality_results$information_content
  cat("Information Content Analysis:\n")
  cat("  Total IC:", round(ic$total_ic, 3), "bits\n")
  cat("  Mean IC per position:", round(ic$mean_ic, 3), "bits\n")
  cat("  High-IC positions:", length(ic$high_ic_positions), "\n")
  cat("  Low-IC positions:", length(ic$low_ic_positions), "\n\n")
  
  # Conservation
  cons <- quality_results$conservation
  cat("Conservation Analysis:\n")
  cat("  Conservation pattern:", cons$conservation_pattern, "\n")
  cat("  Mean conservation:", round(cons$mean_conservation, 3), "\n")
  cat("  Conserved positions:", length(cons$conserved_positions), "\n")
  cat("  Conserved regions:", length(cons$conserved_regions), "\n\n")
  
  # Structure
  struct <- quality_results$structure
  cat("Structure Analysis:\n")
  cat("  Motif length:", struct$motif_length, "bp\n")
  cat("  Structure type:", struct$structure_type, "\n")
  cat("  Core positions:", length(struct$core_positions), "\n")
  cat("  Mean GC content:", round(struct$mean_gc_content, 3), "\n")
  cat("  Palindrome score:", round(struct$palindrome_score, 3), "\n\n")
  
  # Statistical Properties
  stats <- quality_results$statistics
  cat("Statistical Properties:\n")
  cat("  Column uniformity:", round(stats$col_uniformity, 3), "\n")
  cat("  Mean position entropy:", round(stats$mean_position_entropy, 3), "\n")
  cat("  Matrix rank:", stats$matrix_rank, "\n")
  cat("  Condition number:", round(stats$condition_number, 3), "\n\n")
  
  # Biological Relevance (if available)
  if ("biological" %in% names(quality_results)) {
    bio <- quality_results$biological
    cat("Biological Relevance:\n")
    cat("  Mean sequence score:", round(bio$score_statistics$mean_score, 3), "\n")
    cat("  High-scoring sequences:", bio$n_high_scoring, "\n")
    cat("  Score distribution:", bio$score_distribution_type, "\n\n")
  }
  
  # Recommendations
  cat("RECOMMENDATIONS:\n")
  cat("================\n")
  
  if (overall$overall_score >= 80) {
    cat("✅ EXCELLENT PWM quality. Ready for publication and biological applications.\n")
  } else if (overall$overall_score >= 70) {
    cat("✅ GOOD PWM quality. Suitable for most applications.\n")
  } else if (overall$overall_score >= 60) {
    cat("⚠️  FAIR PWM quality. Consider refinement for critical applications.\n")
  } else {
    cat("❌ POOR PWM quality. Significant improvement needed.\n")
  }
  
  if (ic$total_ic < 5) {
    cat("- Consider increasing information content through better alignment or filtering.\n")
  }
  
  if (length(cons$conserved_positions) < 3) {
    cat("- Low conservation suggests weak motif signal. Review input sequences.\n")
  }
  
  if (struct$structure_type == "WEAK_STRUCTURE") {
    cat("- Weak structural organization. Consider motif refinement.\n")
  }
  
  sink()
  
  log_message("INFO", paste("Quality assessment report written to:", report_file))
}

# Command line interface
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript assess_quality.R <pwm_file> [sequences_file] [output_dir]\n")
    cat("Example: Rscript assess_quality.R results/pwm_robust.rds data/filtered_sequences.fasta results/quality\n")
    quit(status = 1)
  }
  
  pwm_file <- args[1]
  sequences_file <- if (length(args) >= 2 && args[2] != "NULL") args[2] else NULL
  output_dir <- ifelse(length(args) >= 3, args[3], "results/quality_assessment")
  
  tryCatch({
    result <- assess_pwm_quality(pwm_file, sequences_file, output_dir)
    cat("\nQuality Assessment Results:\n")
    cat("Overall Score:", round(result$overall$overall_score, 1), "/100\n")
    cat("Quality Grade:", result$overall$grade, "\n")
  }, error = function(e) {
    log_message("ERROR", paste("Quality assessment failed:", e$message))
    quit(status = 1)
  })
}

# Run if called directly
if (!interactive()) {
  main()
}
