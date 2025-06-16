# Integrated Alignment Script - Hybrid approach combining center and consensus methods
# This script implements the integrated alignment strategy as defined in the system architecture

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Please install Biostrings: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/filtered_sequences.fasta"
output_file <- if (length(args) >= 2) args[2] else "data/integrated_aligned_sequences.fasta"
target_length <- if (length(args) >= 3) as.numeric(args[3]) else 200
center_weight <- if (length(args) >= 4) as.numeric(args[4]) else 0.5
consensus_weight <- if (length(args) >= 5) as.numeric(args[5]) else 0.5

cat("Integrated Sequence Alignment (Hybrid Approach)\n")
cat("================================================\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Target length:", target_length, "bp\n")
cat("Center weight:", center_weight, "\n")
cat("Consensus weight:", consensus_weight, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Load alignment methods ---
source_dir <- dirname(input_file)
scripts_dir <- dirname(sys.frame(1)$ofile)
if (is.null(scripts_dir)) scripts_dir <- "scripts"

# Try to source the other alignment scripts
align_center_path <- file.path(scripts_dir, "align_center.R")
align_consensus_path <- file.path(scripts_dir, "align_consensus.R")

# --- Functions ---

# Extract center region
extract_center_region <- function(sequence, target_len) {
  seq_len <- width(sequence)
  
  if (seq_len <= target_len) {
    return(sequence)
  }
  
  center_pos <- seq_len %/% 2
  start_pos <- max(1, center_pos - target_len %/% 2)
  end_pos <- min(seq_len, start_pos + target_len - 1)
  
  if (end_pos - start_pos + 1 < target_len) {
    if (start_pos == 1) {
      end_pos <- min(seq_len, target_len)
    } else {
      start_pos <- max(1, seq_len - target_len + 1)
    }
  }
  
  return(subseq(sequence, start_pos, end_pos))
}

# Find optimal center position with information content
find_optimal_center_position <- function(sequences, window_size = 50) {
  cat("Finding optimal center position using information content...\n")
  
  # Sample sequences for analysis
  sample_size <- min(1000, length(sequences))
  sample_idx <- sample(length(sequences), sample_size)
  sample_seqs <- sequences[sample_idx]
  
  # Find minimum length for analysis
  min_length <- min(width(sample_seqs))
  if (min_length < window_size) window_size <- min_length
  
  # Calculate IC at different positions
  max_positions <- min_length - window_size + 1
  if (max_positions <= 0) return(1)
  
  ic_scores <- numeric(max_positions)
  
  for (pos in 1:max_positions) {
    # Extract windows at this position
    windows <- sapply(sample_seqs, function(seq) {
      as.character(subseq(seq, pos, pos + window_size - 1))
    })
    
    # Calculate IC for this window
    ic_scores[pos] <- calculate_window_information_content(windows)
  }
  
  optimal_pos <- which.max(ic_scores)
  cat("Optimal center position:", optimal_pos, "with IC =", round(ic_scores[optimal_pos], 3), "\n")
  
  return(optimal_pos)
}

# Calculate information content for sequence windows
calculate_window_information_content <- function(windows) {
  window_matrix <- do.call(rbind, strsplit(windows, ""))
  
  if (ncol(window_matrix) == 0) return(0)
  
  total_ic <- 0
  for (pos in 1:ncol(window_matrix)) {
    base_counts <- table(window_matrix[, pos])
    valid_bases <- base_counts[names(base_counts) %in% c("A", "C", "G", "T")]
    
    if (length(valid_bases) > 0 && sum(valid_bases) > 0) {
      base_probs <- valid_bases / sum(valid_bases)
      ic <- sum(base_probs * log2(base_probs / 0.25))
      total_ic <- total_ic + ic
    }
  }
  
  return(total_ic)
}

# Find consensus motif using k-mer analysis
find_consensus_motif_advanced <- function(sequences, motif_length = 15) {
  cat("Finding consensus motif using k-mer analysis...\n")
  
  # Sample sequences for motif discovery
  sample_size <- min(2000, length(sequences))
  sample_idx <- sample(length(sequences), sample_size)
  sample_seqs <- sequences[sample_idx]
  
  # Extract k-mers
  kmer_counts <- list()
  
  for (seq in sample_seqs) {
    seq_char <- as.character(seq)
    seq_len <- nchar(seq_char)
    
    if (seq_len >= motif_length) {
      for (start in 1:(seq_len - motif_length + 1)) {
        kmer <- substr(seq_char, start, start + motif_length - 1)
        if (!grepl("N", kmer)) {
          kmer_counts[[kmer]] <- (kmer_counts[[kmer]] %||% 0) + 1
        }
      }
    }
  }
  
  if (length(kmer_counts) == 0) return(NULL)
  
  # Find most frequent k-mer
  kmer_frequencies <- unlist(kmer_counts)
  most_frequent <- names(kmer_frequencies)[which.max(kmer_frequencies)]
  max_freq <- max(kmer_frequencies) / length(sample_seqs)
  
  cat("Most frequent motif:", most_frequent, "frequency:", round(max_freq * 100, 1), "%\n")
  
  return(list(
    motif = most_frequent,
    frequency = max_freq,
    length = motif_length
  ))
}

# Build PWM from motif
build_motif_pwm <- function(sequences, motif, motif_length) {
  cat("Building PWM from motif instances...\n")
  
  motif_instances <- character()
  
  # Find motif instances in sequences
  for (seq in sequences) {
    seq_char <- as.character(seq)
    seq_len <- nchar(seq_char)
    
    for (start in 1:(seq_len - motif_length + 1)) {
      subseq <- substr(seq_char, start, start + motif_length - 1)
      
      # Check similarity to motif
      similarity <- sum(strsplit(subseq, "")[[1]] == strsplit(motif, "")[[1]]) / motif_length
      
      if (similarity >= 0.8 && !grepl("N", subseq)) {
        motif_instances <- c(motif_instances, subseq)
      }
    }
  }
  
  if (length(motif_instances) == 0) return(NULL)
  
  # Build PWM
  motif_matrix <- do.call(rbind, strsplit(motif_instances, ""))
  pwm <- matrix(0, nrow = 4, ncol = motif_length)
  rownames(pwm) <- c("A", "C", "G", "T")
  
  for (pos in 1:motif_length) {
    base_counts <- table(factor(motif_matrix[, pos], levels = c("A", "C", "G", "T")))
    pwm[, pos] <- as.numeric(base_counts) / sum(base_counts)
  }
  
  return(pwm)
}

# Score sequence against PWM
score_sequence_against_pwm <- function(sequence, pwm, position) {
  seq_char <- as.character(sequence)
  motif_length <- ncol(pwm)
  
  if (nchar(seq_char) < position + motif_length - 1) return(-Inf)
  
  subseq <- substr(seq_char, position, position + motif_length - 1)
  bases <- strsplit(subseq, "")[[1]]
  
  score <- 0
  for (i in 1:length(bases)) {
    if (bases[i] %in% c("A", "C", "G", "T")) {
      prob <- pwm[bases[i], i]
      if (prob > 0) {
        score <- score + log2(prob / 0.25)
      }
    }
  }
  
  return(score)
}

# Perform integrated alignment
perform_integrated_alignment <- function(sequences, target_len, center_wt, consensus_wt) {
  cat("Performing integrated alignment...\n")
  cat("Center weight:", center_wt, ", Consensus weight:", consensus_wt, "\n")
  
  # Step 1: Find optimal center position
  optimal_center <- find_optimal_center_position(sequences, window_size = 50)
  
  # Step 2: Find consensus motif
  consensus_info <- find_consensus_motif_advanced(sequences, motif_length = 15)
  
  # Step 3: Build consensus PWM if motif found
  consensus_pwm <- NULL
  if (!is.null(consensus_info)) {
    consensus_pwm <- build_motif_pwm(sequences, consensus_info$motif, consensus_info$length)
  }
  
  # Step 4: Align sequences using integrated approach
  aligned_sequences <- DNAStringSet()
  alignment_methods <- character(length(sequences))
  
  for (i in 1:length(sequences)) {
    if (i %% 1000 == 0) {
      cat("Aligning sequence", i, "/", length(sequences), "\n")
    }
    
    seq <- sequences[i]
    seq_len <- width(seq)
    
    # Method 1: Center-based alignment score
    center_score <- 0
    if (seq_len >= target_len) {
      # Score based on information content around optimal center
      center_start <- max(1, optimal_center - target_len %/% 2)
      center_end <- min(seq_len, center_start + target_len - 1)
      
      if (center_end - center_start + 1 >= target_len * 0.8) {
        center_score <- center_wt
      }
    }
    
    # Method 2: Consensus-based alignment score
    consensus_score <- 0
    best_consensus_pos <- NA
    
    if (!is.null(consensus_pwm) && seq_len >= consensus_info$length) {
      best_pwm_score <- -Inf
      
      for (pos in 1:(seq_len - consensus_info$length + 1)) {
        pwm_score <- score_sequence_against_pwm(seq, consensus_pwm, pos)
        if (pwm_score > best_pwm_score) {
          best_pwm_score <- pwm_score
          best_consensus_pos <- pos
        }
      }
      
      if (best_pwm_score > 0) {
        consensus_score <- consensus_wt * (best_pwm_score / 10)  # Normalize score
      }
    }
    
    # Choose alignment method based on scores
    if (consensus_score > center_score && !is.na(best_consensus_pos)) {
      # Use consensus-based alignment
      motif_center <- best_consensus_pos + consensus_info$length %/% 2
      align_start <- max(1, motif_center - target_len %/% 2)
      align_end <- min(seq_len, align_start + target_len - 1)
      
      if (align_end - align_start + 1 < target_len) {
        if (align_start == 1) {
          align_end <- min(seq_len, target_len)
        } else {
          align_start <- max(1, seq_len - target_len + 1)
        }
      }
      
      aligned_seq <- subseq(seq, align_start, align_end)
      alignment_methods[i] <- "consensus"
      
    } else {
      # Use center-based alignment
      aligned_seq <- extract_center_region(seq, target_len)
      alignment_methods[i] <- "center"
    }
    
    # Pad if necessary
    if (width(aligned_seq) < target_len) {
      aligned_seq <- pad_sequence(aligned_seq, target_len)
    }
    
    aligned_sequences <- c(aligned_sequences, aligned_seq)
  }
  
  # Copy names
  names(aligned_sequences) <- names(sequences)
  
  # Calculate method usage statistics
  method_stats <- table(alignment_methods)
  cat("Alignment method usage:\n")
  for (method in names(method_stats)) {
    cat("  ", method, ":", method_stats[method], "(", 
        round(method_stats[method] / length(sequences) * 100, 1), "%)\n")
  }
  
  return(list(
    sequences = aligned_sequences,
    optimal_center = optimal_center,
    consensus_motif = if (!is.null(consensus_info)) consensus_info$motif else NULL,
    consensus_frequency = if (!is.null(consensus_info)) consensus_info$frequency else 0,
    method_stats = method_stats
  ))
}

# Pad sequence to target length
pad_sequence <- function(sequence, target_len, pad_char = "N") {
  seq_len <- width(sequence)
  
  if (seq_len >= target_len) {
    return(sequence)
  }
  
  pad_needed <- target_len - seq_len
  left_pad <- pad_needed %/% 2
  right_pad <- pad_needed - left_pad
  
  left_padding <- paste(rep(pad_char, left_pad), collapse = "")
  right_padding <- paste(rep(pad_char, right_pad), collapse = "")
  
  padded_seq <- paste0(left_padding, as.character(sequence), right_padding)
  
  return(DNAString(padded_seq))
}

# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Calculate alignment quality metrics
calculate_integrated_alignment_quality <- function(aligned_sequences) {
  cat("Calculating integrated alignment quality...\n")
  
  lengths <- width(aligned_sequences)
  
  # Sample sequences for detailed analysis
  sample_size <- min(5000, length(aligned_sequences))
  sample_idx <- sample(length(aligned_sequences), sample_size)
  sample_seqs <- aligned_sequences[sample_idx]
  
  # Convert to matrix for position-wise analysis
  seq_matrix <- do.call(rbind, lapply(sample_seqs, function(x) strsplit(as.character(x), "")[[1]]))
  
  # Calculate position-wise information content
  n_positions <- ncol(seq_matrix)
  position_ic <- numeric(n_positions)
  
  for (pos in 1:n_positions) {
    base_counts <- table(seq_matrix[, pos])
    valid_bases <- base_counts[names(base_counts) %in% c("A", "C", "G", "T")]
    
    if (length(valid_bases) > 0 && sum(valid_bases) > 0) {
      base_probs <- valid_bases / sum(valid_bases)
      position_ic[pos] <- sum(base_probs * log2(base_probs / 0.25))
    }
  }
  
  return(list(
    length_consistency = length(unique(lengths)) == 1,
    mean_length = mean(lengths),
    total_information_content = sum(position_ic),
    mean_position_ic = mean(position_ic),
    max_position_ic = max(position_ic),
    high_ic_positions = which(position_ic > 1.0),
    position_ic_vector = position_ic
  ))
}

# --- Main Processing ---

# Check input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Validate weights
if (center_weight + consensus_weight != 1.0) {
  cat("Warning: Weights don't sum to 1.0, normalizing...\n")
  total_weight <- center_weight + consensus_weight
  center_weight <- center_weight / total_weight
  consensus_weight <- consensus_weight / total_weight
  cat("Normalized weights - Center:", center_weight, ", Consensus:", consensus_weight, "\n")
}

# Load sequences
cat("Loading sequences from:", input_file, "\n")
input_sequences <- readDNAStringSet(input_file)
cat("Input sequences loaded:", length(input_sequences), "\n\n")

if (length(input_sequences) == 0) {
  stop("No sequences found in input file")
}

# Display initial statistics
initial_lengths <- width(input_sequences)
cat("Initial sequence statistics:\n")
cat("  Count:", length(input_sequences), "\n")
cat("  Length: mean =", round(mean(initial_lengths), 1), 
    "bp, range =", min(initial_lengths), "-", max(initial_lengths), "bp\n")
cat("  Length SD:", round(sd(initial_lengths), 1), "bp\n\n")

# Perform integrated alignment
alignment_result <- perform_integrated_alignment(input_sequences, target_length, 
                                                center_weight, consensus_weight)
aligned_sequences <- alignment_result$sequences

# Calculate alignment quality
quality_metrics <- calculate_integrated_alignment_quality(aligned_sequences)

cat("\nIntegrated alignment quality metrics:\n")
cat("  Length consistency:", quality_metrics$length_consistency, "\n")
cat("  Mean length:", round(quality_metrics$mean_length, 1), "bp\n")
cat("  Total information content:", round(quality_metrics$total_information_content, 2), "bits\n")
cat("  Mean position IC:", round(quality_metrics$mean_position_ic, 3), "bits\n")
cat("  Max position IC:", round(quality_metrics$max_position_ic, 3), "bits\n")
cat("  High-IC positions:", length(quality_metrics$high_ic_positions), "\n")

# Create output directory if needed
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save aligned sequences
cat("\nSaving aligned sequences to:", output_file, "\n")
writeXStringSet(aligned_sequences, output_file)

# Generate comprehensive alignment report
report_file <- gsub("\\.fasta$", "_integrated_alignment_report.txt", output_file)
cat("Generating alignment report:", report_file, "\n")

sink(report_file)
cat("Integrated Alignment Report\n")
cat("===========================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n\n")

cat("Alignment Parameters:\n")
cat("  Target length:", target_length, "bp\n")
cat("  Center weight:", center_weight, "\n")
cat("  Consensus weight:", consensus_weight, "\n\n")

cat("Alignment Results:\n")
cat("  Input sequences:", length(input_sequences), "\n")
cat("  Aligned sequences:", length(aligned_sequences), "\n")
cat("  Optimal center position:", alignment_result$optimal_center, "\n")

if (!is.null(alignment_result$consensus_motif)) {
  cat("  Consensus motif:", alignment_result$consensus_motif, "\n")
  cat("  Motif frequency:", round(alignment_result$consensus_frequency * 100, 1), "%\n")
}

cat("\nMethod Usage:\n")
for (method in names(alignment_result$method_stats)) {
  cat("  ", method, ":", alignment_result$method_stats[method], 
      "(", round(alignment_result$method_stats[method] / length(input_sequences) * 100, 1), "%)\n")
}

cat("\nQuality Metrics:\n")
cat("  Length consistency:", quality_metrics$length_consistency, "\n")
cat("  Total information content:", round(quality_metrics$total_information_content, 2), "bits\n")
cat("  Mean position IC:", round(quality_metrics$mean_position_ic, 3), "bits\n")
cat("  Max position IC:", round(quality_metrics$max_position_ic, 3), "bits\n")
cat("  High-IC positions:", length(quality_metrics$high_ic_positions), "\n")

# Top information content positions
if (length(quality_metrics$position_ic_vector) > 0) {
  cat("\nTop 10 Information Content Positions:\n")
  top_positions <- order(quality_metrics$position_ic_vector, decreasing = TRUE)[1:min(10, length(quality_metrics$position_ic_vector))]
  for (pos in top_positions) {
    cat("  Position", pos, ":", round(quality_metrics$position_ic_vector[pos], 3), "bits\n")
  }
}

sink()

cat("\nIntegrated alignment completed successfully!\n")
cat("Aligned", length(aligned_sequences), "sequences saved to", output_file, "\n")
cat("This hybrid approach combines the best of center-based and consensus-based methods.\n")
