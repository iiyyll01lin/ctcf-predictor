# Center-based Alignment Script
# This script implements center-based alignment strategy for sequence alignment

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Please install Biostrings: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/filtered_sequences.fasta"
output_file <- if (length(args) >= 2) args[2] else "data/center_aligned_sequences.fasta"
target_length <- if (length(args) >= 3) as.numeric(args[3]) else 200
center_window <- if (length(args) >= 4) as.numeric(args[4]) else 50

cat("Center-based Sequence Alignment\n")
cat("===============================\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Target length:", target_length, "bp\n")
cat("Center window:", center_window, "bp\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Functions ---

# Extract center region of specified length
extract_center_region <- function(sequence, target_len) {
  seq_len <- width(sequence)
  
  if (seq_len <= target_len) {
    # Sequence is shorter than or equal to target, return as is
    return(sequence)
  }
  
  # Calculate center position
  center_pos <- seq_len %/% 2
  start_pos <- max(1, center_pos - target_len %/% 2)
  end_pos <- min(seq_len, start_pos + target_len - 1)
  
  # Adjust if we're at the boundaries
  if (end_pos - start_pos + 1 < target_len) {
    if (start_pos == 1) {
      end_pos <- min(seq_len, target_len)
    } else {
      start_pos <- max(1, seq_len - target_len + 1)
    }
  }
  
  return(subseq(sequence, start_pos, end_pos))
}

# Pad sequence to target length (if shorter)
pad_sequence <- function(sequence, target_len, pad_char = "N") {
  seq_len <- width(sequence)
  
  if (seq_len >= target_len) {
    return(sequence)
  }
  
  # Calculate padding needed
  pad_needed <- target_len - seq_len
  left_pad <- pad_needed %/% 2
  right_pad <- pad_needed - left_pad
  
  # Create padding strings
  left_padding <- paste(rep(pad_char, left_pad), collapse = "")
  right_padding <- paste(rep(pad_char, right_pad), collapse = "")
  
  # Combine sequence with padding
  padded_seq <- paste0(left_padding, as.character(sequence), right_padding)
  
  return(DNAString(padded_seq))
}

# Find optimal center position based on information content
find_optimal_center <- function(sequences, window_size = 50) {
  cat("Finding optimal center position...\n")
  
  # Convert sequences to matrix for analysis
  seq_lengths <- width(sequences)
  min_length <- min(seq_lengths)
  
  if (min_length < window_size) {
    cat("Warning: Some sequences shorter than window size, using minimum length\n")
    window_size <- min_length
  }
  
  # Sample sequences for analysis (for speed)
  sample_size <- min(1000, length(sequences))
  sample_idx <- sample(length(sequences), sample_size)
  sample_seqs <- sequences[sample_idx]
  
  # Calculate information content at each position
  max_pos <- min(seq_lengths) - window_size + 1
  if (max_pos <= 0) {
    return(1)  # Default to position 1 if sequences are too short
  }
  
  ic_scores <- numeric(max_pos)
  
  for (pos in 1:max_pos) {
    if (pos %% 50 == 0) {
      cat("Analyzing position", pos, "/", max_pos, "\n")
    }
    
    # Extract windows from all sequences at this position
    windows <- sapply(sample_seqs, function(seq) {
      as.character(subseq(seq, pos, pos + window_size - 1))
    })
    
    # Calculate information content for this position
    ic_scores[pos] <- calculate_window_ic(windows)
  }
  
  # Find position with maximum information content
  optimal_pos <- which.max(ic_scores)
  cat("Optimal center position:", optimal_pos, "with IC =", round(ic_scores[optimal_pos], 3), "\n")
  
  return(optimal_pos)
}

# Calculate information content for a set of sequence windows
calculate_window_ic <- function(windows) {
  # Convert to matrix
  window_matrix <- do.call(rbind, strsplit(windows, ""))
  
  if (ncol(window_matrix) == 0) {
    return(0)
  }
  
  total_ic <- 0
  
  for (pos in 1:ncol(window_matrix)) {
    base_counts <- table(window_matrix[, pos])
    
    # Only count ACGT bases
    valid_bases <- base_counts[names(base_counts) %in% c("A", "C", "G", "T")]
    
    if (length(valid_bases) == 0 || sum(valid_bases) == 0) {
      next
    }
    
    # Calculate probabilities
    base_probs <- valid_bases / sum(valid_bases)
    
    # Calculate information content (bits)
    # IC = sum(p * log2(p/0.25)) where 0.25 is background probability
    ic <- sum(base_probs * log2(base_probs / 0.25))
    total_ic <- total_ic + ic
  }
  
  return(total_ic)
}

# Perform center-based alignment
perform_center_alignment <- function(sequences, target_len = 200, method = "fixed") {
  cat("Performing center-based alignment...\n")
  cat("Method:", method, "\n")
  cat("Target length:", target_len, "bp\n")
  
  aligned_sequences <- DNAStringSet()
  
  if (method == "optimal") {
    # Find optimal center position first
    optimal_center <- find_optimal_center(sequences, center_window)
  }
  
  for (i in 1:length(sequences)) {
    if (i %% 1000 == 0) {
      cat("Aligning sequence", i, "/", length(sequences), "\n")
    }
    
    seq <- sequences[i]
    
    if (method == "fixed") {
      # Simple center extraction
      aligned_seq <- extract_center_region(seq, target_len)
    } else if (method == "optimal") {
      # Use optimal center position
      seq_len <- width(seq)
      
      if (seq_len <= target_len) {
        aligned_seq <- seq
      } else {
        # Center around optimal position
        center_start <- max(1, optimal_center - target_len %/% 2)
        center_end <- min(seq_len, center_start + target_len - 1)
        
        if (center_end - center_start + 1 < target_len) {
          # Adjust boundaries
          if (center_start == 1) {
            center_end <- min(seq_len, target_len)
          } else {
            center_start <- max(1, seq_len - target_len + 1)
            center_end <- seq_len
          }
        }
        
        aligned_seq <- subseq(seq, center_start, center_end)
      }
    }
    
    # Pad if necessary
    if (width(aligned_seq) < target_len) {
      aligned_seq <- pad_sequence(aligned_seq, target_len)
    }
    
    aligned_sequences <- c(aligned_sequences, aligned_seq)
  }
  
  # Copy names from original sequences
  names(aligned_sequences) <- names(sequences)
  
  cat("Alignment completed: ", length(aligned_sequences), "sequences aligned\n")
  
  return(aligned_sequences)
}

# Calculate alignment quality metrics
calculate_alignment_quality <- function(aligned_sequences) {
  cat("Calculating alignment quality metrics...\n")
  
  # Check length consistency
  lengths <- width(aligned_sequences)
  length_consistency <- length(unique(lengths)) == 1
  
  # Calculate position-wise information content
  n_positions <- max(lengths)
  position_ic <- numeric(n_positions)
  
  # Sample sequences for IC calculation (for speed)
  sample_size <- min(5000, length(aligned_sequences))
  sample_idx <- sample(length(aligned_sequences), sample_size)
  sample_seqs <- aligned_sequences[sample_idx]
  
  # Convert to matrix for analysis
  seq_matrix <- do.call(rbind, lapply(sample_seqs, function(x) strsplit(as.character(x), "")[[1]]))
  
  for (pos in 1:n_positions) {
    if (pos <= ncol(seq_matrix)) {
      base_counts <- table(seq_matrix[, pos])
      valid_bases <- base_counts[names(base_counts) %in% c("A", "C", "G", "T")]
      
      if (length(valid_bases) > 0 && sum(valid_bases) > 0) {
        base_probs <- valid_bases / sum(valid_bases)
        position_ic[pos] <- sum(base_probs * log2(base_probs / 0.25))
      }
    }
  }
  
  return(list(
    length_consistency = length_consistency,
    mean_length = mean(lengths),
    length_range = range(lengths),
    total_information_content = sum(position_ic),
    mean_position_ic = mean(position_ic),
    position_ic = position_ic,
    high_ic_positions = which(position_ic > 1.0)
  ))
}

# --- Main Processing ---

# Check input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
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

# Determine alignment method
alignment_method <- if (center_window > 0) "optimal" else "fixed"

# Perform center-based alignment
aligned_sequences <- perform_center_alignment(input_sequences, target_length, alignment_method)

# Calculate alignment quality
quality_metrics <- calculate_alignment_quality(aligned_sequences)

cat("\nAlignment quality metrics:\n")
cat("  Length consistency:", quality_metrics$length_consistency, "\n")
cat("  Mean length:", round(quality_metrics$mean_length, 1), "bp\n")
cat("  Total information content:", round(quality_metrics$total_information_content, 2), "bits\n")
cat("  Mean position IC:", round(quality_metrics$mean_position_ic, 3), "bits\n")
cat("  High-IC positions:", length(quality_metrics$high_ic_positions), "\n")

# Create output directory if needed
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save aligned sequences
cat("\nSaving aligned sequences to:", output_file, "\n")
writeXStringSet(aligned_sequences, output_file)

# Generate alignment report
report_file <- gsub("\\.fasta$", "_alignment_report.txt", output_file)
cat("Generating alignment report:", report_file, "\n")

sink(report_file)
cat("Center-based Alignment Report\n")
cat("=============================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n\n")

cat("Alignment Parameters:\n")
cat("  Method:", alignment_method, "\n")
cat("  Target length:", target_length, "bp\n")
cat("  Center window:", center_window, "bp\n\n")

cat("Alignment Results:\n")
cat("  Input sequences:", length(input_sequences), "\n")
cat("  Aligned sequences:", length(aligned_sequences), "\n")
cat("  Length consistency:", quality_metrics$length_consistency, "\n")
cat("  Mean aligned length:", round(quality_metrics$mean_length, 1), "bp\n\n")

cat("Quality Metrics:\n")
cat("  Total information content:", round(quality_metrics$total_information_content, 2), "bits\n")
cat("  Mean position IC:", round(quality_metrics$mean_position_ic, 3), "bits\n")
cat("  High-IC positions:", length(quality_metrics$high_ic_positions), "\n")
cat("  Information content per position (top 10):\n")

if (length(quality_metrics$position_ic) > 0) {
  top_positions <- order(quality_metrics$position_ic, decreasing = TRUE)[1:min(10, length(quality_metrics$position_ic))]
  for (pos in top_positions) {
    cat("    Position", pos, ":", round(quality_metrics$position_ic[pos], 3), "bits\n")
  }
}

sink()

cat("\nCenter-based alignment completed successfully!\n")
cat("Aligned", length(aligned_sequences), "sequences saved to", output_file, "\n")
