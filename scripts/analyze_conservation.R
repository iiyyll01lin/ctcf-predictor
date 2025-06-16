#!/usr/bin/env Rscript

# CTCF PWM Testing Pipeline - Conservation Pattern Analysis
# Detailed analysis of conservation patterns in PWM matrices
# Author: CTCF Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(Biostrings)
  library(cluster)
})

# Logging function
log_message <- function(level, message, script = "analyze_conservation.R") {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(sprintf("%s %s [%s] %s\n", timestamp, level, script, message))
}

# Main conservation analysis function
analyze_conservation_patterns <- function(pwm_file, sequences_file = NULL,
                                        output_dir = "results/conservation_analysis",
                                        conservation_threshold = 1.0) {
  
  log_message("INFO", "Starting conservation pattern analysis")
  
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
  
  # Perform conservation analyses
  conservation_results <- list()
  
  # 1. Position-wise Conservation Analysis
  log_message("INFO", "Analyzing position-wise conservation")
  conservation_results$position_wise <- analyze_position_conservation(pwm_data, conservation_threshold)
  
  # 2. Regional Conservation Analysis
  log_message("INFO", "Analyzing regional conservation patterns")
  conservation_results$regional <- analyze_regional_conservation(pwm_data, conservation_threshold)
  
  # 3. Comparative Conservation Analysis
  log_message("INFO", "Analyzing comparative conservation")
  conservation_results$comparative <- analyze_comparative_conservation(pwm_data)
  
  # 4. Structural Conservation Analysis
  log_message("INFO", "Analyzing structural conservation")
  conservation_results$structural <- analyze_structural_conservation(pwm_data)
  
  # 5. Evolutionary Conservation Analysis (if sequences provided)
  if (!is.null(sequences)) {
    log_message("INFO", "Analyzing evolutionary conservation")
    conservation_results$evolutionary <- analyze_evolutionary_conservation(pwm_data, sequences)
  }
  
  # 6. Conservation Clustering Analysis
  log_message("INFO", "Performing conservation clustering")
  conservation_results$clustering <- analyze_conservation_clustering(pwm_data)
  
  # Generate comprehensive report
  generate_conservation_report(pwm_data, conservation_results, output_dir)
  
  # Generate visualizations
  generate_conservation_plots(pwm_data, conservation_results, output_dir)
  
  # Save results
  results_file <- file.path(output_dir, "conservation_analysis_results.rds")
  saveRDS(conservation_results, results_file)
  
  log_message("INFO", paste("Conservation analysis completed. Results saved to:", results_file))
  return(conservation_results)
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

# Analyze position-wise conservation
analyze_position_conservation <- function(pwm_data, threshold = 1.0) {
  pwm <- pwm_data$pwm
  
  # Calculate conservation scores for each position
  conservation_scores <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    entropy <- -sum(col * log2(col))
    conservation <- 2 - entropy  # Maximum entropy is 2 for 4 bases
    return(conservation)
  })
  
  # Information content scores
  information_content <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col))
  })
  
  # Dominant base analysis
  dominant_bases <- apply(pwm, 2, function(col) {
    rownames(pwm)[which.max(col)]
  })
  
  dominant_probabilities <- apply(pwm, 2, max)
  
  # Conservation classification
  conservation_classes <- ifelse(conservation_scores >= threshold, "CONSERVED", "VARIABLE")
  conservation_classes[conservation_scores >= 1.5] <- "HIGHLY_CONSERVED"
  conservation_classes[conservation_scores < 0.5] <- "POORLY_CONSERVED"
  
  # Position-wise statistics
  position_stats <- data.frame(
    position = 1:ncol(pwm),
    conservation_score = conservation_scores,
    information_content = information_content,
    dominant_base = dominant_bases,
    dominant_probability = dominant_probabilities,
    conservation_class = conservation_classes,
    stringsAsFactors = FALSE
  )
  
  # Summary statistics
  summary_stats <- list(
    mean_conservation = mean(conservation_scores),
    median_conservation = median(conservation_scores),
    max_conservation = max(conservation_scores),
    min_conservation = min(conservation_scores),
    conservation_variance = var(conservation_scores),
    n_conserved = sum(conservation_scores >= threshold),
    n_highly_conserved = sum(conservation_scores >= 1.5),
    n_poorly_conserved = sum(conservation_scores < 0.5),
    conservation_ratio = sum(conservation_scores >= threshold) / length(conservation_scores)
  )
  
  return(list(
    position_stats = position_stats,
    summary_stats = summary_stats,
    conservation_scores = conservation_scores,
    threshold = threshold
  ))
}

# Analyze regional conservation patterns
analyze_regional_conservation <- function(pwm_data, threshold = 1.0) {
  pwm <- pwm_data$pwm
  position_conservation <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 - (-sum(col * log2(col)))
  })
  
  # Find conserved regions
  conserved_positions <- which(position_conservation >= threshold)
  conserved_regions <- find_conserved_regions(conserved_positions)
  
  # Analyze flanking regions
  flanking_analysis <- analyze_flanking_regions(pwm, conserved_regions)
  
  # Core-flanking pattern analysis
  core_flanking_pattern <- classify_core_flanking_pattern(position_conservation, ncol(pwm))
  
  # Regional statistics
  regional_stats <- list(
    n_regions = length(conserved_regions),
    total_conserved_length = sum(sapply(conserved_regions, function(r) r$length)),
    mean_region_length = if (length(conserved_regions) > 0) {
      mean(sapply(conserved_regions, function(r) r$length))
    } else { 0 },
    max_region_length = if (length(conserved_regions) > 0) {
      max(sapply(conserved_regions, function(r) r$length))
    } else { 0 },
    core_flanking_pattern = core_flanking_pattern
  )
  
  return(list(
    conserved_regions = conserved_regions,
    flanking_analysis = flanking_analysis,
    regional_stats = regional_stats,
    position_conservation = position_conservation
  ))
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
        length = current_end - current_start + 1,
        positions = current_start:current_end
      )
      current_start <- positions[i]
      current_end <- positions[i]
    }
  }
  
  # Add final region
  regions[[length(regions) + 1]] <- list(
    start = current_start,
    end = current_end,
    length = current_end - current_start + 1,
    positions = current_start:current_end
  )
  
  return(regions)
}

# Analyze flanking regions
analyze_flanking_regions <- function(pwm, conserved_regions) {
  if (length(conserved_regions) == 0) {
    return(list(
      left_flanking_conservation = numeric(0),
      right_flanking_conservation = numeric(0),
      flanking_pattern = "NO_CONSERVED_REGIONS"
    ))
  }
  
  # Calculate conservation for all positions
  position_conservation <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 - (-sum(col * log2(col)))
  })
  
  # Analyze flanking regions around conserved core
  left_flanking <- numeric(0)
  right_flanking <- numeric(0)
  
  for (region in conserved_regions) {
    # Left flanking
    left_start <- max(1, region$start - 3)
    left_end <- region$start - 1
    if (left_end >= left_start) {
      left_flanking <- c(left_flanking, position_conservation[left_start:left_end])
    }
    
    # Right flanking
    right_start <- region$end + 1
    right_end <- min(ncol(pwm), region$end + 3)
    if (right_end >= right_start) {
      right_flanking <- c(right_flanking, position_conservation[right_start:right_end])
    }
  }
  
  # Classify flanking pattern
  flanking_pattern <- classify_flanking_pattern(left_flanking, right_flanking)
  
  return(list(
    left_flanking_conservation = left_flanking,
    right_flanking_conservation = right_flanking,
    flanking_pattern = flanking_pattern,
    mean_left_flanking = if (length(left_flanking) > 0) mean(left_flanking) else 0,
    mean_right_flanking = if (length(right_flanking) > 0) mean(right_flanking) else 0
  ))
}

# Classify core-flanking pattern
classify_core_flanking_pattern <- function(conservation_scores, motif_length) {
  n_pos <- length(conservation_scores)
  
  # Identify core (high conservation) positions
  high_cons_positions <- which(conservation_scores > quantile(conservation_scores, 0.7))
  
  if (length(high_cons_positions) == 0) {
    return("NO_CLEAR_CORE")
  }
  
  # Analyze core distribution
  core_start <- min(high_cons_positions)
  core_end <- max(high_cons_positions)
  core_span <- core_end - core_start + 1
  
  # Classify pattern
  if (core_span / n_pos > 0.8) {
    return("UNIFORM_CONSERVATION")
  } else if (core_start <= 0.2 * n_pos && core_end >= 0.8 * n_pos) {
    return("CENTRAL_CORE_WITH_FLANKING")
  } else if (core_start <= 0.3 * n_pos) {
    return("LEFT_BIASED_CORE")
  } else if (core_end >= 0.7 * n_pos) {
    return("RIGHT_BIASED_CORE")
  } else {
    return("CENTRAL_CORE")
  }
}

# Classify flanking pattern
classify_flanking_pattern <- function(left_flanking, right_flanking) {
  if (length(left_flanking) == 0 && length(right_flanking) == 0) {
    return("NO_FLANKING")
  }
  
  left_mean <- if (length(left_flanking) > 0) mean(left_flanking) else 0
  right_mean <- if (length(right_flanking) > 0) mean(right_flanking) else 0
  
  if (left_mean > 0.8 && right_mean > 0.8) {
    return("CONSERVED_FLANKING")
  } else if (left_mean > 0.8) {
    return("LEFT_CONSERVED_FLANKING")
  } else if (right_mean > 0.8) {
    return("RIGHT_CONSERVED_FLANKING")
  } else {
    return("VARIABLE_FLANKING")
  }
}

# Analyze comparative conservation
analyze_comparative_conservation <- function(pwm_data) {
  pwm <- pwm_data$pwm
  
  # Compare conservation across different nucleotides
  nucleotide_conservation <- list()
  for (nuc in rownames(pwm)) {
    nuc_positions <- which(apply(pwm, 2, which.max) == which(rownames(pwm) == nuc))
    if (length(nuc_positions) > 0) {
      nuc_conservation <- apply(pwm[, nuc_positions, drop = FALSE], 2, function(col) {
        col[col == 0] <- 1e-10
        2 - (-sum(col * log2(col)))
      })
      nucleotide_conservation[[nuc]] <- list(
        positions = nuc_positions,
        conservation_scores = nuc_conservation,
        mean_conservation = mean(nuc_conservation),
        n_positions = length(nuc_positions)
      )
    }
  }
  
  # Purine vs Pyrimidine analysis
  purine_positions <- which(apply(pwm, 2, which.max) %in% which(rownames(pwm) %in% c("A", "G")))
  pyrimidine_positions <- which(apply(pwm, 2, which.max) %in% which(rownames(pwm) %in% c("C", "T")))
  
  purine_conservation <- if (length(purine_positions) > 0) {
    mean(apply(pwm[, purine_positions, drop = FALSE], 2, function(col) {
      col[col == 0] <- 1e-10
      2 - (-sum(col * log2(col)))
    }))
  } else { 0 }
  
  pyrimidine_conservation <- if (length(pyrimidine_positions) > 0) {
    mean(apply(pwm[, pyrimidine_positions, drop = FALSE], 2, function(col) {
      col[col == 0] <- 1e-10
      2 - (-sum(col * log2(col)))
    }))
  } else { 0 }
  
  return(list(
    nucleotide_conservation = nucleotide_conservation,
    purine_conservation = purine_conservation,
    pyrimidine_conservation = pyrimidine_conservation,
    purine_positions = purine_positions,
    pyrimidine_positions = pyrimidine_positions
  ))
}

# Analyze structural conservation
analyze_structural_conservation <- function(pwm_data) {
  pwm <- pwm_data$pwm
  
  # Palindrome conservation analysis
  palindrome_conservation <- analyze_palindrome_conservation(pwm)
  
  # Repeat conservation analysis
  repeat_conservation <- analyze_repeat_conservation(pwm)
  
  # Symmetry analysis
  symmetry_analysis <- analyze_motif_symmetry(pwm)
  
  # Spacing conservation (if applicable)
  spacing_conservation <- analyze_spacing_conservation(pwm)
  
  return(list(
    palindrome = palindrome_conservation,
    repeat = repeat_conservation,
    symmetry = symmetry_analysis,
    spacing = spacing_conservation
  ))
}

# Analyze palindrome conservation
analyze_palindrome_conservation <- function(pwm) {
  n_pos <- ncol(pwm)
  palindrome_scores <- numeric(floor(n_pos/2))
  
  for (i in 1:floor(n_pos/2)) {
    left_pos <- i
    right_pos <- n_pos - i + 1
    
    # Calculate complement match scores
    at_complement <- pwm["A", left_pos] * pwm["T", right_pos] + pwm["T", left_pos] * pwm["A", right_pos]
    gc_complement <- pwm["G", left_pos] * pwm["C", right_pos] + pwm["C", left_pos] * pwm["G", right_pos]
    
    palindrome_scores[i] <- at_complement + gc_complement
  }
  
  return(list(
    palindrome_scores = palindrome_scores,
    mean_palindrome_score = mean(palindrome_scores),
    max_palindrome_score = max(palindrome_scores),
    palindrome_strength = classify_palindrome_strength(palindrome_scores)
  ))
}

# Classify palindrome strength
classify_palindrome_strength <- function(palindrome_scores) {
  mean_score <- mean(palindrome_scores)
  
  if (mean_score > 0.8) {
    return("STRONG_PALINDROME")
  } else if (mean_score > 0.6) {
    return("MODERATE_PALINDROME")
  } else if (mean_score > 0.4) {
    return("WEAK_PALINDROME")
  } else {
    return("NO_PALINDROME")
  }
}

# Analyze repeat conservation
analyze_repeat_conservation <- function(pwm) {
  n_pos <- ncol(pwm)
  
  if (n_pos < 4) {
    return(list(
      repeat_scores = numeric(0),
      repeat_strength = "TOO_SHORT",
      best_repeat_length = 0
    ))
  }
  
  # Test different repeat lengths
  repeat_lengths <- 2:min(6, floor(n_pos/2))
  repeat_scores <- numeric(length(repeat_lengths))
  
  for (i in 1:length(repeat_lengths)) {
    repeat_len <- repeat_lengths[i]
    n_repeats <- floor(n_pos / repeat_len)
    
    if (n_repeats >= 2) {
      # Calculate correlation between repeats
      correlations <- numeric(n_repeats - 1)
      
      for (j in 1:(n_repeats - 1)) {
        start1 <- (j - 1) * repeat_len + 1
        end1 <- j * repeat_len
        start2 <- j * repeat_len + 1
        end2 <- (j + 1) * repeat_len
        
        if (end2 <= n_pos) {
          correlations[j] <- cor(as.vector(pwm[, start1:end1]), as.vector(pwm[, start2:end2]))
        }
      }
      
      repeat_scores[i] <- mean(correlations, na.rm = TRUE)
    }
  }
  
  best_repeat_idx <- which.max(repeat_scores)
  best_repeat_length <- repeat_lengths[best_repeat_idx]
  best_repeat_score <- repeat_scores[best_repeat_idx]
  
  repeat_strength <- if (best_repeat_score > 0.8) {
    "STRONG_REPEAT"
  } else if (best_repeat_score > 0.6) {
    "MODERATE_REPEAT"
  } else if (best_repeat_score > 0.4) {
    "WEAK_REPEAT"
  } else {
    "NO_REPEAT"
  }
  
  return(list(
    repeat_scores = repeat_scores,
    repeat_lengths = repeat_lengths,
    best_repeat_length = best_repeat_length,
    best_repeat_score = best_repeat_score,
    repeat_strength = repeat_strength
  ))
}

# Analyze motif symmetry
analyze_motif_symmetry <- function(pwm) {
  n_pos <- ncol(pwm)
  
  # Test for reflection symmetry
  reflection_symmetry <- cor(as.vector(pwm), as.vector(pwm[, ncol(pwm):1]))
  
  # Test for rotational symmetry (if applicable)
  rotational_symmetry <- 0
  if (n_pos >= 4) {
    mid_point <- ceiling(n_pos / 2)
    left_half <- pwm[, 1:mid_point]
    right_half <- pwm[, (n_pos - ncol(left_half) + 1):n_pos]
    
    if (ncol(left_half) == ncol(right_half)) {
      rotational_symmetry <- cor(as.vector(left_half), as.vector(right_half))
    }
  }
  
  return(list(
    reflection_symmetry = reflection_symmetry,
    rotational_symmetry = rotational_symmetry,
    symmetry_type = classify_symmetry_type(reflection_symmetry, rotational_symmetry)
  ))
}

# Classify symmetry type
classify_symmetry_type <- function(reflection_symmetry, rotational_symmetry) {
  if (reflection_symmetry > 0.8 && rotational_symmetry > 0.8) {
    return("HIGHLY_SYMMETRIC")
  } else if (reflection_symmetry > 0.6) {
    return("REFLECTION_SYMMETRIC")
  } else if (rotational_symmetry > 0.6) {
    return("ROTATIONALLY_SYMMETRIC")
  } else if (reflection_symmetry > 0.4 || rotational_symmetry > 0.4) {
    return("MODERATELY_SYMMETRIC")
  } else {
    return("ASYMMETRIC")
  }
}

# Analyze spacing conservation
analyze_spacing_conservation <- function(pwm) {
  n_pos <- ncol(pwm)
  
  # Calculate information content for each position
  ic_values <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col))
  })
  
  # Find peaks in information content
  peaks <- find_ic_peaks(ic_values)
  
  # Analyze spacing between peaks
  if (length(peaks) >= 2) {
    spacings <- diff(peaks)
    spacing_stats <- list(
      spacings = spacings,
      mean_spacing = mean(spacings),
      spacing_variance = var(spacings),
      regular_spacing = var(spacings) < 1
    )
  } else {
    spacing_stats <- list(
      spacings = numeric(0),
      mean_spacing = 0,
      spacing_variance = 0,
      regular_spacing = FALSE
    )
  }
  
  return(list(
    peaks = peaks,
    spacing_stats = spacing_stats,
    spacing_pattern = classify_spacing_pattern(spacing_stats)
  ))
}

# Find information content peaks
find_ic_peaks <- function(ic_values) {
  if (length(ic_values) < 3) return(integer(0))
  
  peaks <- integer(0)
  threshold <- quantile(ic_values, 0.7)
  
  for (i in 2:(length(ic_values) - 1)) {
    if (ic_values[i] > ic_values[i-1] && ic_values[i] > ic_values[i+1] && ic_values[i] > threshold) {
      peaks <- c(peaks, i)
    }
  }
  
  return(peaks)
}

# Classify spacing pattern
classify_spacing_pattern <- function(spacing_stats) {
  if (length(spacing_stats$spacings) == 0) {
    return("NO_PEAKS")
  } else if (spacing_stats$regular_spacing) {
    return("REGULAR_SPACING")
  } else {
    return("IRREGULAR_SPACING")
  }
}

# Analyze evolutionary conservation
analyze_evolutionary_conservation <- function(pwm_data, sequences) {
  # This is a simplified evolutionary analysis
  # In practice, this would require multiple species data
  
  pwm <- pwm_data$pwm
  motif_length <- ncol(pwm)
  
  # Score sequences and analyze conservation
  sequence_scores <- score_sequences_against_pwm(sequences, pwm)
  
  # Analyze conservation across high-scoring sequences
  high_score_threshold <- quantile(sequence_scores, 0.9, na.rm = TRUE)
  high_scoring_indices <- which(sequence_scores >= high_score_threshold)
  
  if (length(high_scoring_indices) > 10) {
    high_scoring_sequences <- sequences[high_scoring_indices]
    
    # Extract best-matching subsequences
    best_matches <- extract_best_matches(high_scoring_sequences, pwm)
    
    # Build conservation profile
    conservation_profile <- build_conservation_profile(best_matches)
    
    return(list(
      n_high_scoring = length(high_scoring_sequences),
      conservation_profile = conservation_profile,
      evolutionary_score = mean(conservation_profile, na.rm = TRUE)
    ))
  } else {
    return(list(
      n_high_scoring = length(high_scoring_indices),
      conservation_profile = NULL,
      evolutionary_score = 0
    ))
  }
}

# Score sequences against PWM (simplified)
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

# Extract best matches (simplified)
extract_best_matches <- function(sequences, pwm) {
  motif_length <- ncol(pwm)
  best_matches <- character(0)
  
  for (seq in sequences) {
    seq_str <- as.character(seq)
    seq_len <- nchar(seq_str)
    
    if (seq_len >= motif_length) {
      # Find center region (simplified)
      center_start <- max(1, floor((seq_len - motif_length) / 2) + 1)
      center_end <- center_start + motif_length - 1
      
      if (center_end <= seq_len) {
        best_matches <- c(best_matches, substr(seq_str, center_start, center_end))
      }
    }
  }
  
  return(best_matches)
}

# Build conservation profile
build_conservation_profile <- function(matches) {
  if (length(matches) == 0) return(NULL)
  
  motif_length <- nchar(matches[1])
  conservation_profile <- numeric(motif_length)
  
  for (pos in 1:motif_length) {
    position_bases <- substr(matches, pos, pos)
    base_counts <- table(position_bases)
    
    # Calculate conservation as 1 - normalized entropy
    if (length(base_counts) > 0) {
      base_probs <- base_counts / sum(base_counts)
      entropy <- -sum(base_probs * log2(base_probs + 1e-10))
      conservation_profile[pos] <- 1 - (entropy / log2(4))  # Normalize by max entropy
    }
  }
  
  return(conservation_profile)
}

# Analyze conservation clustering
analyze_conservation_clustering <- function(pwm_data) {
  pwm <- pwm_data$pwm
  
  # Transpose PWM for position clustering
  pwm_t <- t(pwm)
  
  # Calculate distance matrix
  distance_matrix <- dist(pwm_t, method = "euclidean")
  
  # Perform hierarchical clustering
  hc <- hclust(distance_matrix, method = "ward.D2")
  
  # Cut tree to get clusters
  k <- min(4, nrow(pwm_t) - 1)  # Number of clusters
  if (k >= 2) {
    clusters <- cutree(hc, k = k)
    
    # Analyze cluster characteristics
    cluster_analysis <- list()
    for (i in 1:k) {
      cluster_positions <- which(clusters == i)
      if (length(cluster_positions) > 0) {
        cluster_conservation <- apply(pwm[, cluster_positions, drop = FALSE], 2, function(col) {
          col[col == 0] <- 1e-10
          2 - (-sum(col * log2(col)))
        })
        
        cluster_analysis[[paste0("cluster_", i)]] <- list(
          positions = cluster_positions,
          size = length(cluster_positions),
          mean_conservation = mean(cluster_conservation),
          conservation_range = range(cluster_conservation)
        )
      }
    }
    
    return(list(
      clusters = clusters,
      n_clusters = k,
      cluster_analysis = cluster_analysis,
      clustering_quality = calculate_clustering_quality(pwm_t, clusters)
    ))
  } else {
    return(list(
      clusters = rep(1, nrow(pwm_t)),
      n_clusters = 1,
      cluster_analysis = list(),
      clustering_quality = 0
    ))
  }
}

# Calculate clustering quality (simplified silhouette score)
calculate_clustering_quality <- function(data, clusters) {
  if (length(unique(clusters)) <= 1) return(0)
  
  tryCatch({
    sil <- silhouette(clusters, dist(data))
    return(mean(sil[, 3]))
  }, error = function(e) {
    return(0)
  })
}

# Generate conservation plots
generate_conservation_plots <- function(pwm_data, conservation_results, output_dir) {
  # Conservation score plot
  if ("position_wise" %in% names(conservation_results)) {
    png(file.path(output_dir, "conservation_scores.png"), width = 800, height = 400)
    barplot(conservation_results$position_wise$conservation_scores,
            main = "Position-wise Conservation Scores",
            xlab = "Position", ylab = "Conservation Score",
            col = "darkgreen")
    abline(h = conservation_results$position_wise$threshold, col = "red", lty = 2)
    dev.off()
  }
  
  # Information content vs conservation plot
  if ("position_wise" %in% names(conservation_results)) {
    png(file.path(output_dir, "ic_vs_conservation.png"), width = 800, height = 400)
    plot(conservation_results$position_wise$position_stats$conservation_score,
         conservation_results$position_wise$position_stats$information_content,
         main = "Information Content vs Conservation",
         xlab = "Conservation Score", ylab = "Information Content",
         pch = 16, col = "blue")
    dev.off()
  }
  
  log_message("INFO", "Conservation plots generated")
}

# Generate conservation report
generate_conservation_report <- function(pwm_data, conservation_results, output_dir) {
  report_file <- file.path(output_dir, "conservation_analysis_report.txt")
  
  sink(report_file)
  cat("CTCF PWM Conservation Analysis Report\n")
  cat("====================================\n\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # PWM Information
  cat("PWM Information:\n")
  cat("  Dimensions:", nrow(pwm_data$pwm), "x", ncol(pwm_data$pwm), "\n\n")
  
  # Position-wise Conservation
  if ("position_wise" %in% names(conservation_results)) {
    pos_wise <- conservation_results$position_wise
    cat("POSITION-WISE CONSERVATION:\n")
    cat("  Mean conservation:", round(pos_wise$summary_stats$mean_conservation, 3), "\n")
    cat("  Conservation threshold:", pos_wise$threshold, "\n")
    cat("  Conserved positions:", pos_wise$summary_stats$n_conserved, "\n")
    cat("  Highly conserved positions:", pos_wise$summary_stats$n_highly_conserved, "\n")
    cat("  Conservation ratio:", round(pos_wise$summary_stats$conservation_ratio, 3), "\n\n")
    
    # Most conserved positions
    top_positions <- order(pos_wise$conservation_scores, decreasing = TRUE)[1:min(5, length(pos_wise$conservation_scores))]
    cat("  Top conserved positions:\n")
    for (pos in top_positions) {
      cat("    Position", pos, ":", round(pos_wise$conservation_scores[pos], 3), 
          "(", pos_wise$position_stats$dominant_base[pos], ")\n")
    }
    cat("\n")
  }
  
  # Regional Conservation
  if ("regional" %in% names(conservation_results)) {
    regional <- conservation_results$regional
    cat("REGIONAL CONSERVATION:\n")
    cat("  Number of conserved regions:", regional$regional_stats$n_regions, "\n")
    cat("  Total conserved length:", regional$regional_stats$total_conserved_length, "\n")
    cat("  Core-flanking pattern:", regional$regional_stats$core_flanking_pattern, "\n")
    
    if (length(regional$conserved_regions) > 0) {
      cat("  Conserved regions:\n")
      for (i in 1:length(regional$conserved_regions)) {
        region <- regional$conserved_regions[[i]]
        cat("    Region", i, ": positions", region$start, "-", region$end, 
            "(length", region$length, ")\n")
      }
    }
    cat("\n")
  }
  
  # Comparative Conservation
  if ("comparative" %in% names(conservation_results)) {
    comp <- conservation_results$comparative
    cat("COMPARATIVE CONSERVATION:\n")
    cat("  Purine conservation:", round(comp$purine_conservation, 3), "\n")
    cat("  Pyrimidine conservation:", round(comp$pyrimidine_conservation, 3), "\n")
    
    cat("  Nucleotide-specific conservation:\n")
    for (nuc in names(comp$nucleotide_conservation)) {
      nuc_data <- comp$nucleotide_conservation[[nuc]]
      cat("    ", nuc, ": mean =", round(nuc_data$mean_conservation, 3), 
          ", positions =", nuc_data$n_positions, "\n")
    }
    cat("\n")
  }
  
  # Structural Conservation
  if ("structural" %in% names(conservation_results)) {
    struct <- conservation_results$structural
    cat("STRUCTURAL CONSERVATION:\n")
    
    if ("palindrome" %in% names(struct)) {
      cat("  Palindrome strength:", struct$palindrome$palindrome_strength, "\n")
      cat("  Mean palindrome score:", round(struct$palindrome$mean_palindrome_score, 3), "\n")
    }
    
    if ("repeat" %in% names(struct)) {
      cat("  Repeat strength:", struct$repeat$repeat_strength, "\n")
      if (struct$repeat$best_repeat_length > 0) {
        cat("  Best repeat length:", struct$repeat$best_repeat_length, "\n")
        cat("  Best repeat score:", round(struct$repeat$best_repeat_score, 3), "\n")
      }
    }
    
    if ("symmetry" %in% names(struct)) {
      cat("  Symmetry type:", struct$symmetry$symmetry_type, "\n")
      cat("  Reflection symmetry:", round(struct$symmetry$reflection_symmetry, 3), "\n")
    }
    cat("\n")
  }
  
  # Evolutionary Conservation
  if ("evolutionary" %in% names(conservation_results)) {
    evol <- conservation_results$evolutionary
    cat("EVOLUTIONARY CONSERVATION:\n")
    cat("  High-scoring sequences:", evol$n_high_scoring, "\n")
    cat("  Evolutionary score:", round(evol$evolutionary_score, 3), "\n\n")
  }
  
  # Clustering Analysis
  if ("clustering" %in% names(conservation_results)) {
    clust <- conservation_results$clustering
    cat("CLUSTERING ANALYSIS:\n")
    cat("  Number of clusters:", clust$n_clusters, "\n")
    cat("  Clustering quality:", round(clust$clustering_quality, 3), "\n")
    
    if (length(clust$cluster_analysis) > 0) {
      cat("  Cluster details:\n")
      for (cluster_name in names(clust$cluster_analysis)) {
        cluster_data <- clust$cluster_analysis[[cluster_name]]
        cat("    ", cluster_name, ": size =", cluster_data$size, 
            ", conservation =", round(cluster_data$mean_conservation, 3), "\n")
      }
    }
    cat("\n")
  }
  
  # Summary and Recommendations
  cat("SUMMARY:\n")
  cat("========\n")
  
  if ("position_wise" %in% names(conservation_results)) {
    conservation_ratio <- conservation_results$position_wise$summary_stats$conservation_ratio
    
    if (conservation_ratio >= 0.7) {
      cat("✅ EXCELLENT conservation pattern - highly suitable for biological applications\n")
    } else if (conservation_ratio >= 0.5) {
      cat("✅ GOOD conservation pattern - suitable for most applications\n")
    } else if (conservation_ratio >= 0.3) {
      cat("⚠️  MODERATE conservation pattern - consider refinement\n")
    } else {
      cat("❌ POOR conservation pattern - significant improvement needed\n")
    }
  }
  
  sink()
  
  log_message("INFO", paste("Conservation analysis report written to:", report_file))
}

# Command line interface
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript analyze_conservation.R <pwm_file> [sequences_file] [output_dir] [threshold]\n")
    cat("Example: Rscript analyze_conservation.R results/pwm_robust.rds data/filtered_sequences.fasta results/conservation 1.0\n")
    quit(status = 1)
  }
  
  pwm_file <- args[1]
  sequences_file <- if (length(args) >= 2 && args[2] != "NULL") args[2] else NULL
  output_dir <- ifelse(length(args) >= 3, args[3], "results/conservation_analysis")
  threshold <- ifelse(length(args) >= 4, as.numeric(args[4]), 1.0)
  
  tryCatch({
    result <- analyze_conservation_patterns(pwm_file, sequences_file, output_dir, threshold)
    cat("\nConservation Analysis Completed!\n")
    
    if ("position_wise" %in% names(result)) {
      conservation_ratio <- result$position_wise$summary_stats$conservation_ratio
      cat("Conservation ratio:", round(conservation_ratio, 3), "\n")
      cat("Conserved positions:", result$position_wise$summary_stats$n_conserved, "\n")
    }
  }, error = function(e) {
    log_message("ERROR", paste("Conservation analysis failed:", e$message))
    quit(status = 1)
  })
}

# Run if called directly
if (!interactive()) {
  main()
}
