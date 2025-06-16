#!/usr/bin/env Rscript

# Biological Validation Script
# CTCF-specific biological validation of PWM patterns
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

# Command line options
option_list <- list(
  make_option(c("--sequences", "-s"), type="character", default=NULL,
              help="Input aligned sequences file (FASTA format) [required]", metavar="character"),
  make_option(c("--pwm-file"), type="character", default=NULL,
              help="PWM file (if available) [optional]", metavar="character"),
  make_option(c("--output", "-o"), type="character", default="results/biological_validation.html",
              help="Output validation report file [default: %default]", metavar="character"),
  make_option(c("--reference-motif"), type="character", default="CCGCGNGGNGGCAG",
              help="Reference CTCF motif pattern [default: %default]", metavar="character"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Validate PWM against known CTCF biological patterns")
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$sequences)) {
  cat("Error: Input sequences file is required\n")
  print_help(opt_parser)
  quit(status=1)
}

if (!file.exists(opt$sequences)) {
  cat("Error: Input sequences file does not exist:", opt$sequences, "\n")
  quit(status=1)
}

# Create output directory
output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (opt$verbose) {
  cat("Starting biological validation...\n")
  cat("Input sequences:", opt$sequences, "\n")
  cat("Reference motif:", opt$`reference-motif`, "\n")
  cat("Output file:", opt$output, "\n")
}

#' Perform CTCF-specific biological validation
#' @param sequences_file Path to aligned sequences file
#' @param pwm_file Path to PWM file (optional)
#' @param reference_motif Reference CTCF motif pattern
#' @param verbose Enable verbose output
ctcf_biological_validation <- function(sequences_file, pwm_file = NULL, 
                                     reference_motif = "CCGCGNGGNGGCAG", verbose = FALSE) {
  
  # Load sequences
  if (verbose) cat("Loading aligned sequences from:", sequences_file, "\n")
  sequences <- readDNAStringSet(sequences_file)
  
  if (verbose) cat("Loaded", length(sequences), "aligned sequences\n")
  
  # Build PWM if not provided
  if (!is.null(pwm_file) && file.exists(pwm_file)) {
    if (verbose) cat("Loading PWM from file:", pwm_file, "\n")
    pwm_data <- readRDS(pwm_file)
    if (is.matrix(pwm_data)) {
      pwm <- pwm_data
    } else if (is.list(pwm_data) && "pwm" %in% names(pwm_data)) {
      pwm <- pwm_data$pwm
    } else {
      stop("Invalid PWM file format")
    }
  } else {
    if (verbose) cat("Building PWM from sequences...\n")
    pwm <- build_pwm_from_sequences(sequences, verbose)
  }
  
  # Validate CTCF pattern
  if (verbose) cat("Validating CTCF motif pattern...\n")
  pattern_validation <- validate_ctcf_pattern(pwm, reference_motif, verbose)
  
  # Validate information content profile
  if (verbose) cat("Validating IC profile...\n")
  ic_validation <- validate_ic_profile(pwm, verbose)
  
  # Validate zinc finger binding positions
  if (verbose) cat("Validating zinc finger positions...\n")
  zf_validation <- validate_zinc_finger_positions(pwm, verbose)
  
  # Validate conserved positions
  if (verbose) cat("Validating conserved positions...\n")
  conservation_validation <- validate_conserved_positions(pwm, sequences, verbose)
  
  # Validate GC content bias
  if (verbose) cat("Validating GC content...\n")
  gc_validation <- validate_gc_content(sequences, verbose)
  
  # Validate against known CTCF sites
  if (verbose) cat("Validating against known CTCF characteristics...\n")
  known_sites_validation <- validate_known_ctcf_characteristics(pwm, sequences, verbose)
  
  # Calculate overall biological validity score
  overall_score <- calculate_biological_validity_score(
    pattern_validation, ic_validation, zf_validation, 
    conservation_validation, gc_validation, known_sites_validation, verbose
  )
  
  return(list(
    pattern_validation = pattern_validation,
    ic_validation = ic_validation,
    zf_validation = zf_validation,
    conservation_validation = conservation_validation,
    gc_validation = gc_validation,
    known_sites_validation = known_sites_validation,
    overall_score = overall_score,
    biological_grade = classify_biological_grade(overall_score),
    pwm = pwm,
    reference_motif = reference_motif
  ))
}

#' Build PWM from sequences
build_pwm_from_sequences <- function(sequences, verbose = FALSE) {
  # Convert sequences to matrix
  seq_matrix <- do.call(rbind, strsplit(as.character(sequences), ""))
  
  # Remove positions with too many Ns
  valid_positions <- apply(seq_matrix, 2, function(col) {
    sum(col %in% c("A", "T", "G", "C")) >= length(col) * 0.5
  })
  
  seq_matrix <- seq_matrix[, valid_positions, drop = FALSE]
  
  # Build PWM
  bases <- c("A", "C", "G", "T")
  n_positions <- ncol(seq_matrix)
  pwm <- matrix(0, nrow = 4, ncol = n_positions)
  rownames(pwm) <- bases
  
  pseudocount <- 0.1
  
  for (pos in 1:n_positions) {
    col <- seq_matrix[, pos]
    valid_bases <- col[col %in% bases]
    
    if (length(valid_bases) > 0) {
      counts <- table(factor(valid_bases, levels = bases))
      freqs <- (counts + pseudocount) / (sum(counts) + 4 * pseudocount)
      pwm[, pos] <- freqs
    } else {
      pwm[, pos] <- rep(0.25, 4)
    }
  }
  
  return(pwm)
}

#' Validate CTCF motif pattern
validate_ctcf_pattern <- function(pwm, reference_motif, verbose = FALSE) {
  
  # Generate consensus from PWM
  observed_consensus <- generate_consensus_from_pwm(pwm)
  
  if (verbose) cat("  Observed consensus:", observed_consensus, "\n")
  if (verbose) cat("  Reference motif:", reference_motif, "\n")
  
  # Calculate pattern similarity
  similarity <- calculate_pattern_similarity(observed_consensus, reference_motif)
  
  # Check for key CTCF positions
  key_positions <- identify_key_ctcf_positions(pwm, reference_motif)
  
  # Validate core binding sequence
  core_validation <- validate_core_binding_sequence(pwm, observed_consensus)
  
  return(list(
    observed_consensus = observed_consensus,
    reference_motif = reference_motif,
    similarity_score = similarity,
    key_positions = key_positions,
    core_validation = core_validation,
    pattern_match = similarity > 0.7,
    core_intact = core_validation$core_intact
  ))
}

#' Generate consensus sequence from PWM
generate_consensus_from_pwm <- function(pwm) {
  bases <- c("A", "C", "G", "T")
  consensus_chars <- apply(pwm, 2, function(col) {
    max_idx <- which.max(col)
    if (max(col) > 0.4) {  # Strong preference
      bases[max_idx]
    } else if (max(col) > 0.3) {  # Moderate preference
      bases[max_idx]
    } else {
      "N"  # Ambiguous
    }
  })
  
  return(paste(consensus_chars, collapse = ""))
}

#' Calculate pattern similarity between two sequences
calculate_pattern_similarity <- function(seq1, seq2) {
  # Pad sequences to same length
  max_len <- max(nchar(seq1), nchar(seq2))
  seq1_padded <- sprintf("%-*s", max_len, seq1)
  seq2_padded <- sprintf("%-*s", max_len, seq2)
  
  # Convert to character vectors
  chars1 <- strsplit(seq1_padded, "")[[1]]
  chars2 <- strsplit(seq2_padded, "")[[1]]
  
  # Calculate matches (allowing N to match anything)
  matches <- sum(chars1 == chars2 | chars1 == "N" | chars2 == "N")
  
  return(matches / max_len)
}

#' Identify key CTCF positions
identify_key_ctcf_positions <- function(pwm, reference_motif) {
  
  # Known key positions in CTCF motif (1-indexed)
  # These are highly conserved positions based on literature
  expected_key_positions <- c(1, 2, 3, 4, 8, 9, 11, 12, 13, 14)  # CC[GC][GC]---[GC][GC]-[GC][GC][GC][GC]
  
  # Calculate information content for each position
  ic_values <- apply(pwm, 2, function(col) {
    max(0, 2 + sum(col * log2(col + 1e-10)))
  })
  
  # Identify high IC positions (>1.0 bits)
  high_ic_positions <- which(ic_values > 1.0)
  
  # Check overlap with expected key positions
  expected_in_range <- expected_key_positions[expected_key_positions <= length(ic_values)]
  overlap <- intersect(high_ic_positions, expected_in_range)
  
  return(list(
    high_ic_positions = high_ic_positions,
    expected_key_positions = expected_in_range,
    overlap = overlap,
    overlap_fraction = length(overlap) / length(expected_in_range),
    key_positions_identified = length(overlap) >= length(expected_in_range) * 0.6
  ))
}

#' Validate core binding sequence
validate_core_binding_sequence <- function(pwm, consensus) {
  
  # CTCF core binding sequence is typically the first 4-6 positions
  core_length <- min(6, ncol(pwm))
  core_consensus <- substr(consensus, 1, core_length)
  
  # Expected core patterns for CTCF
  expected_cores <- c("CCGCG", "CCGC", "GCGC", "CCGG")
  
  # Check if observed core matches any expected pattern
  core_matches <- sapply(expected_cores, function(pattern) {
    calculate_pattern_similarity(core_consensus, pattern)
  })
  
  best_match <- max(core_matches)
  best_pattern <- expected_cores[which.max(core_matches)]
  
  return(list(
    core_consensus = core_consensus,
    best_matching_pattern = best_pattern,
    best_match_score = best_match,
    core_intact = best_match > 0.6
  ))
}

#' Validate information content profile
validate_ic_profile <- function(pwm, verbose = FALSE) {
  
  # Calculate IC for each position
  ic_values <- apply(pwm, 2, function(col) {
    max(0, 2 + sum(col * log2(col + 1e-10)))
  })
  
  # Expected CTCF IC profile characteristics
  total_ic <- sum(ic_values)
  mean_ic <- mean(ic_values)
  max_ic <- max(ic_values)
  
  # High IC positions (conserved)
  high_ic_positions <- which(ic_values > 1.0)
  n_conserved <- length(high_ic_positions)
  
  # IC profile shape analysis
  ic_profile_shape <- analyze_ic_profile_shape(ic_values)
  
  # Validate against expected CTCF characteristics
  ic_validation <- list(
    total_ic = total_ic,
    mean_ic = mean_ic,
    max_ic = max_ic,
    n_conserved_positions = n_conserved,
    conserved_positions = high_ic_positions,
    ic_profile_shape = ic_profile_shape,
    sufficient_ic = total_ic > 8.0,  # Minimum acceptable IC for CTCF
    adequate_conservation = n_conserved >= 3,  # At least 3 conserved positions
    good_profile = ic_profile_shape$has_peaks
  )
  
  return(ic_validation)
}

#' Analyze IC profile shape
analyze_ic_profile_shape <- function(ic_values) {
  
  # Find peaks (local maxima)
  peaks <- c()
  if (length(ic_values) > 2) {
    for (i in 2:(length(ic_values) - 1)) {
      if (ic_values[i] > ic_values[i-1] && ic_values[i] > ic_values[i+1] && ic_values[i] > 1.0) {
        peaks <- c(peaks, i)
      }
    }
  }
  
  # Check for characteristic CTCF profile
  # CTCF typically has high IC at beginning and end, lower in middle
  beginning_ic <- mean(ic_values[1:min(4, length(ic_values))])
  end_ic <- mean(ic_values[max(1, length(ic_values)-3):length(ic_values)])
  middle_ic <- if (length(ic_values) > 8) {
    mean(ic_values[5:(length(ic_values)-4)])
  } else {
    mean(ic_values)
  }
  
  return(list(
    peaks = peaks,
    n_peaks = length(peaks),
    has_peaks = length(peaks) >= 2,
    beginning_ic = beginning_ic,
    end_ic = end_ic,
    middle_ic = middle_ic,
    u_shaped = beginning_ic > middle_ic && end_ic > middle_ic
  ))
}

#' Validate zinc finger binding positions
validate_zinc_finger_positions <- function(pwm, verbose = FALSE) {
  
  # CTCF has 11 zinc fingers, with specific preferences
  # Zinc fingers 4-7 are most important for DNA binding
  
  # Calculate base preferences at each position
  base_preferences <- apply(pwm, 2, function(col) {
    bases <- c("A", "C", "G", "T")
    max_base <- bases[which.max(col)]
    max_prob <- max(col)
    list(base = max_base, probability = max_prob)
  })
  
  # Check for zinc finger-specific patterns
  # ZF typically prefer G/C in major groove contacts
  gc_positions <- sapply(base_preferences, function(x) x$base %in% c("G", "C"))
  strong_gc_positions <- which(gc_positions & sapply(base_preferences, function(x) x$probability > 0.5))
  
  # Validate spacing patterns (zinc fingers are typically spaced 3-4 bp apart)
  spacing_validation <- validate_zf_spacing(strong_gc_positions)
  
  return(list(
    base_preferences = base_preferences,
    gc_positions = which(gc_positions),
    strong_gc_positions = strong_gc_positions,
    n_strong_gc = length(strong_gc_positions),
    spacing_validation = spacing_validation,
    zf_pattern_likely = length(strong_gc_positions) >= 4 && spacing_validation$regular_spacing
  ))
}

#' Validate zinc finger spacing
validate_zf_spacing <- function(positions) {
  if (length(positions) < 2) {
    return(list(regular_spacing = FALSE, mean_spacing = NA))
  }
  
  spacings <- diff(positions)
  mean_spacing <- mean(spacings)
  
  # Zinc finger spacing is typically 3-4 bp
  regular_spacing <- all(spacings >= 2 & spacings <= 6) && mean_spacing >= 2.5 && mean_spacing <= 5.5
  
  return(list(
    spacings = spacings,
    mean_spacing = mean_spacing,
    regular_spacing = regular_spacing
  ))
}

#' Validate conserved positions
validate_conserved_positions <- function(pwm, sequences, verbose = FALSE) {
  
  # Calculate conservation scores from sequences
  seq_matrix <- do.call(rbind, strsplit(as.character(sequences), ""))
  
  # Remove positions with too many Ns
  valid_positions <- apply(seq_matrix, 2, function(col) {
    sum(col %in% c("A", "T", "G", "C")) >= length(col) * 0.5
  })
  
  seq_matrix <- seq_matrix[, valid_positions, drop = FALSE]
  
  # Calculate conservation for each position
  conservation_scores <- apply(seq_matrix, 2, function(col) {
    valid_bases <- col[col %in% c("A", "T", "G", "C")]
    if (length(valid_bases) < 5) return(0)
    
    freqs <- table(valid_bases)
    max_freq <- max(freqs)
    return(max_freq / length(valid_bases))
  })
  
  # Identify highly conserved positions
  highly_conserved <- which(conservation_scores > 0.8)
  moderately_conserved <- which(conservation_scores > 0.6 & conservation_scores <= 0.8)
  
  return(list(
    conservation_scores = conservation_scores,
    mean_conservation = mean(conservation_scores),
    highly_conserved_positions = highly_conserved,
    moderately_conserved_positions = moderately_conserved,
    n_highly_conserved = length(highly_conserved),
    n_moderately_conserved = length(moderately_conserved),
    adequate_conservation = length(highly_conserved) >= 3
  ))
}

#' Validate GC content
validate_gc_content <- function(sequences, verbose = FALSE) {
  
  # Calculate GC content for each sequence
  gc_contents <- letterFrequency(sequences, letters = "GC", as.prob = TRUE)[,1]
  
  # CTCF sites are typically GC-rich
  mean_gc <- mean(gc_contents)
  median_gc <- median(gc_contents)
  
  # Check if GC content is in expected range for CTCF (typically 50-80%)
  expected_gc_range <- c(0.5, 0.8)
  gc_in_range <- sum(gc_contents >= expected_gc_range[1] & gc_contents <= expected_gc_range[2])
  gc_fraction_in_range <- gc_in_range / length(gc_contents)
  
  return(list(
    gc_contents = gc_contents,
    mean_gc = mean_gc,
    median_gc = median_gc,
    expected_range = expected_gc_range,
    fraction_in_range = gc_fraction_in_range,
    gc_rich = mean_gc > 0.5,
    appropriate_gc_distribution = gc_fraction_in_range > 0.7
  ))
}

#' Validate against known CTCF characteristics
validate_known_ctcf_characteristics <- function(pwm, sequences, verbose = FALSE) {
  
  # Length validation - CTCF motifs are typically 15-20 bp
  motif_length <- ncol(pwm)
  appropriate_length <- motif_length >= 12 && motif_length <= 25
  
  # Palindromic tendency - CTCF has some palindromic characteristics
  palindromic_score <- assess_palindromic_tendency(pwm)
  
  # Binding strength distribution
  binding_strength <- assess_binding_strength_distribution(pwm, sequences)
  
  # Degeneracy assessment - CTCF tolerates some degeneracy
  degeneracy <- assess_motif_degeneracy(pwm)
  
  return(list(
    motif_length = motif_length,
    appropriate_length = appropriate_length,
    palindromic_score = palindromic_score,
    binding_strength = binding_strength,
    degeneracy = degeneracy,
    ctcf_like = appropriate_length && palindromic_score$somewhat_palindromic && binding_strength$reasonable_distribution
  ))
}

#' Assess palindromic tendency
assess_palindromic_tendency <- function(pwm) {
  
  n_pos <- ncol(pwm)
  if (n_pos < 4) {
    return(list(palindromic_score = 0, somewhat_palindromic = FALSE))
  }
  
  # Compare first half with reverse complement of second half
  mid_point <- floor(n_pos / 2)
  first_half <- pwm[, 1:mid_point]
  second_half <- pwm[, (n_pos - mid_point + 1):n_pos]
  
  # Reverse complement mapping: A<->T, C<->G
  rc_mapping <- c("A" = 4, "C" = 3, "G" = 2, "T" = 1)  # T, G, C, A indices
  
  # Reverse the second half and apply complement
  second_half_rc <- second_half[c(4, 3, 2, 1), ncol(second_half):1, drop = FALSE]
  
  # Calculate correlation
  if (ncol(first_half) == ncol(second_half_rc)) {
    correlation <- cor(as.vector(first_half), as.vector(second_half_rc))
    palindromic_score <- max(0, correlation)  # Only positive correlations count
  } else {
    palindromic_score <- 0
  }
  
  return(list(
    palindromic_score = palindromic_score,
    somewhat_palindromic = palindromic_score > 0.3
  ))
}

#' Assess binding strength distribution
assess_binding_strength_distribution <- function(pwm, sequences) {
  
  # Calculate binding scores for all sequences
  binding_scores <- calculate_binding_scores(pwm, sequences)
  
  # Analyze distribution
  mean_score <- mean(binding_scores)
  median_score <- median(binding_scores)
  score_range <- range(binding_scores)
  score_sd <- sd(binding_scores)
  
  # Check for reasonable distribution
  # Should have a range of scores, not all identical
  reasonable_distribution <- score_sd > 0.5 && diff(score_range) > 2.0
  
  return(list(
    binding_scores = binding_scores,
    mean_score = mean_score,
    median_score = median_score,
    score_range = score_range,
    score_sd = score_sd,
    reasonable_distribution = reasonable_distribution
  ))
}

#' Calculate binding scores for sequences
calculate_binding_scores <- function(pwm, sequences) {
  
  scores <- sapply(sequences, function(seq) {
    seq_chars <- strsplit(as.character(seq), "")[[1]]
    
    if (length(seq_chars) != ncol(pwm)) {
      return(NA)  # Length mismatch
    }
    
    score <- 0
    for (i in 1:length(seq_chars)) {
      base <- seq_chars[i]
      if (base %in% c("A", "C", "G", "T")) {
        base_idx <- match(base, c("A", "C", "G", "T"))
        score <- score + log2(pwm[base_idx, i] + 1e-10)
      }
    }
    
    return(score)
  })
  
  return(scores[!is.na(scores)])
}

#' Assess motif degeneracy
assess_motif_degeneracy <- function(pwm) {
  
  # Calculate entropy for each position
  entropies <- apply(pwm, 2, function(col) {
    -sum(col * log2(col + 1e-10))
  })
  
  # Classify positions by degeneracy
  highly_specific <- sum(entropies < 0.5)  # Low entropy = high specificity
  moderately_specific <- sum(entropies >= 0.5 & entropies < 1.0)
  degenerate <- sum(entropies >= 1.0)
  
  mean_entropy <- mean(entropies)
  
  # CTCF should have a mix of specific and degenerate positions
  balanced_degeneracy <- highly_specific >= 2 && degenerate >= 1
  
  return(list(
    entropies = entropies,
    mean_entropy = mean_entropy,
    highly_specific = highly_specific,
    moderately_specific = moderately_specific,
    degenerate = degenerate,
    balanced_degeneracy = balanced_degeneracy
  ))
}

#' Calculate overall biological validity score
calculate_biological_validity_score <- function(pattern_validation, ic_validation, zf_validation,
                                              conservation_validation, gc_validation, known_sites_validation, verbose = FALSE) {
  
  # Weight different aspects of biological validation
  weights <- list(
    pattern_match = 0.25,
    ic_profile = 0.20,
    zinc_finger = 0.15,
    conservation = 0.15,
    gc_content = 0.10,
    known_characteristics = 0.15
  )
  
  # Score each component (0-1 scale)
  scores <- list(
    pattern_match = ifelse(pattern_validation$pattern_match, 1.0, pattern_validation$similarity_score),
    ic_profile = calculate_ic_score(ic_validation),
    zinc_finger = ifelse(zf_validation$zf_pattern_likely, 1.0, 0.5),
    conservation = calculate_conservation_score(conservation_validation),
    gc_content = calculate_gc_score(gc_validation),
    known_characteristics = ifelse(known_sites_validation$ctcf_like, 1.0, 0.5)
  )
  
  # Calculate weighted average
  overall_score <- sum(mapply(function(w, s) w * s, weights, scores))
  
  return(list(
    component_scores = scores,
    weights = weights,
    overall_score = overall_score
  ))
}

#' Calculate IC validation score
calculate_ic_score <- function(ic_validation) {
  score <- 0
  
  if (ic_validation$sufficient_ic) score <- score + 0.4
  if (ic_validation$adequate_conservation) score <- score + 0.3
  if (ic_validation$good_profile) score <- score + 0.3
  
  return(score)
}

#' Calculate conservation validation score
calculate_conservation_score <- function(conservation_validation) {
  score <- 0
  
  if (conservation_validation$adequate_conservation) score <- score + 0.5
  if (conservation_validation$mean_conservation > 0.6) score <- score + 0.3
  if (conservation_validation$n_moderately_conserved >= 2) score <- score + 0.2
  
  return(score)
}

#' Calculate GC content validation score
calculate_gc_score <- function(gc_validation) {
  score <- 0
  
  if (gc_validation$gc_rich) score <- score + 0.5
  if (gc_validation$appropriate_gc_distribution) score <- score + 0.5
  
  return(score)
}

#' Classify biological validation grade
classify_biological_grade <- function(overall_score) {
  if (overall_score >= 0.8) {
    return("EXCELLENT")
  } else if (overall_score >= 0.7) {
    return("GOOD")
  } else if (overall_score >= 0.6) {
    return("ACCEPTABLE")
  } else {
    return("POOR")
  }
}

#' Generate biological validation report
generate_biological_report <- function(validation_results, output_file) {
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>CTCF Biological Validation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #f0f8ff; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }
        .excellent { background-color: #e8f5e8; }
        .good { background-color: #e8f0ff; }
        .acceptable { background-color: #fff8e8; }
        .poor { background-color: #ffeaea; }
        table { border-collapse: collapse; width: 100%; margin: 10px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .metric { margin: 10px 0; }
        .pass { color: green; font-weight: bold; }
        .fail { color: red; font-weight: bold; }
        .consensus { font-family: monospace; font-size: 16px; background-color: #f0f0f0; padding: 5px; }
    </style>
</head>
<body>
    <div class='header'>
        <h1>CTCF Biological Validation Report</h1>
        <p>Generated: ", Sys.time(), "</p>
        <p>Overall Biological Grade: <strong>", validation_results$biological_grade, "</strong></p>
        <p>Overall Score: ", round(validation_results$overall_score$overall_score, 3), "</p>
    </div>")
  
  # Pattern validation
  pattern <- validation_results$pattern_validation
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>CTCF Pattern Validation</h2>
        <div class='metric'>Observed Consensus: <span class='consensus'>", pattern$observed_consensus, "</span></div>
        <div class='metric'>Reference Motif: <span class='consensus'>", pattern$reference_motif, "</span></div>
        <div class='metric'>Similarity Score: ", round(pattern$similarity_score, 3), "</div>
        <div class='metric'>Pattern Match: <span class='", if(pattern$pattern_match) "pass'>PASS" else "fail'>FAIL", "</span></div>
        <div class='metric'>Core Binding Intact: <span class='", if(pattern$core_validation$core_intact) "pass'>PASS" else "fail'>FAIL", "</span></div>
        <div class='metric'>Key Positions Identified: ", pattern$key_positions$overlap_fraction * 100, "% (", length(pattern$key_positions$overlap), "/", length(pattern$key_positions$expected_key_positions), ")</div>
    </div>")
  
  # IC validation
  ic <- validation_results$ic_validation
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Information Content Profile</h2>
        <div class='metric'>Total IC: ", round(ic$total_ic, 2), " bits</div>
        <div class='metric'>Mean IC: ", round(ic$mean_ic, 2), " bits</div>
        <div class='metric'>Conserved Positions: ", ic$n_conserved_positions, " (", paste(ic$conserved_positions, collapse=", "), ")</div>
        <div class='metric'>Sufficient IC: <span class='", if(ic$sufficient_ic) "pass'>PASS" else "fail'>FAIL", "</span> (≥8.0 bits)</div>
        <div class='metric'>Adequate Conservation: <span class='", if(ic$adequate_conservation) "pass'>PASS" else "fail'>FAIL", "</span> (≥3 positions)</div>
        <div class='metric'>Profile Shape: ", if(ic$ic_profile_shape$has_peaks) "Has characteristic peaks" else "Flat profile", "</div>
    </div>")
  
  # Zinc finger validation
  zf <- validation_results$zf_validation
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Zinc Finger Binding Pattern</h2>
        <div class='metric'>Strong GC Positions: ", zf$n_strong_gc, " (", paste(zf$strong_gc_positions, collapse=", "), ")</div>
        <div class='metric'>Regular Spacing: <span class='", if(zf$spacing_validation$regular_spacing) "pass'>PASS" else "fail'>FAIL", "</span></div>
        <div class='metric'>Mean Spacing: ", round(zf$spacing_validation$mean_spacing, 1), " bp</div>
        <div class='metric'>Zinc Finger Pattern Likely: <span class='", if(zf$zf_pattern_likely) "pass'>PASS" else "fail'>FAIL", "</span></div>
    </div>")
  
  # Conservation validation
  cons <- validation_results$conservation_validation
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Sequence Conservation</h2>
        <div class='metric'>Mean Conservation: ", round(cons$mean_conservation, 3), "</div>
        <div class='metric'>Highly Conserved Positions: ", cons$n_highly_conserved, " (>80% conservation)</div>
        <div class='metric'>Moderately Conserved Positions: ", cons$n_moderately_conserved, " (60-80% conservation)</div>
        <div class='metric'>Adequate Conservation: <span class='", if(cons$adequate_conservation) "pass'>PASS" else "fail'>FAIL", "</span></div>
    </div>")
  
  # GC content validation
  gc <- validation_results$gc_validation
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>GC Content Analysis</h2>
        <div class='metric'>Mean GC Content: ", round(gc$mean_gc * 100, 1), "%</div>
        <div class='metric'>Expected Range: ", gc$expected_range[1] * 100, "-", gc$expected_range[2] * 100, "%</div>
        <div class='metric'>Sequences in Range: ", round(gc$fraction_in_range * 100, 1), "%</div>
        <div class='metric'>GC-Rich: <span class='", if(gc$gc_rich) "pass'>PASS" else "fail'>FAIL", "</span></div>
        <div class='metric'>Appropriate Distribution: <span class='", if(gc$appropriate_gc_distribution) "pass'>PASS" else "fail'>FAIL", "</span></div>
    </div>")
  
  # Known characteristics validation
  known <- validation_results$known_sites_validation
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Known CTCF Characteristics</h2>
        <div class='metric'>Motif Length: ", known$motif_length, " bp</div>
        <div class='metric'>Appropriate Length: <span class='", if(known$appropriate_length) "pass'>PASS" else "fail'>FAIL", "</span> (12-25 bp range)</div>
        <div class='metric'>Palindromic Score: ", round(known$palindromic_score$palindromic_score, 3), "</div>
        <div class='metric'>Somewhat Palindromic: <span class='", if(known$palindromic_score$somewhat_palindromic) "pass'>PASS" else "fail'>FAIL", "</span></div>
        <div class='metric'>Reasonable Binding Distribution: <span class='", if(known$binding_strength$reasonable_distribution) "pass'>PASS" else "fail'>FAIL", "</span></div>
        <div class='metric'>CTCF-like: <span class='", if(known$ctcf_like) "pass'>PASS" else "fail'>FAIL", "</span></div>
    </div>")
  
  # Component scores
  scores <- validation_results$overall_score$component_scores
  weights <- validation_results$overall_score$weights
  
  html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Validation Component Scores</h2>
        <table>
            <tr><th>Component</th><th>Weight</th><th>Score</th><th>Weighted Score</th></tr>
            <tr><td>Pattern Match</td><td>", weights$pattern_match, "</td><td>", round(scores$pattern_match, 3), "</td><td>", round(weights$pattern_match * scores$pattern_match, 3), "</td></tr>
            <tr><td>IC Profile</td><td>", weights$ic_profile, "</td><td>", round(scores$ic_profile, 3), "</td><td>", round(weights$ic_profile * scores$ic_profile, 3), "</td></tr>
            <tr><td>Zinc Finger</td><td>", weights$zinc_finger, "</td><td>", round(scores$zinc_finger, 3), "</td><td>", round(weights$zinc_finger * scores$zinc_finger, 3), "</td></tr>
            <tr><td>Conservation</td><td>", weights$conservation, "</td><td>", round(scores$conservation, 3), "</td><td>", round(weights$conservation * scores$conservation, 3), "</td></tr>
            <tr><td>GC Content</td><td>", weights$gc_content, "</td><td>", round(scores$gc_content, 3), "</td><td>", round(weights$gc_content * scores$gc_content, 3), "</td></tr>
            <tr><td>Known Characteristics</td><td>", weights$known_characteristics, "</td><td>", round(scores$known_characteristics, 3), "</td><td>", round(weights$known_characteristics * scores$known_characteristics, 3), "</td></tr>
            <tr><th>Total</th><th>1.0</th><th>-</th><th>", round(validation_results$overall_score$overall_score, 3), "</th></tr>
        </table>
    </div>")
  
  html_content <- paste0(html_content, "
</body>
</html>")
  
  writeLines(html_content, output_file)
}

# Main execution
validation_results <- ctcf_biological_validation(
  sequences_file = opt$sequences,
  pwm_file = opt$`pwm-file`,
  reference_motif = opt$`reference-motif`,
  verbose = opt$verbose
)

if (opt$verbose) cat("Generating biological validation report...\n")

generate_biological_report(validation_results, opt$output)

if (opt$verbose) {
  cat("Biological validation complete!\n")
  cat("Biological grade:", validation_results$biological_grade, "\n")
  cat("Overall score:", round(validation_results$overall_score$overall_score, 3), "\n")
  cat("Report saved to:", opt$output, "\n")
}

# Exit with appropriate status
if (validation_results$biological_grade %in% c("EXCELLENT", "GOOD")) {
  quit(status = 0)
} else {
  quit(status = 1)
}
