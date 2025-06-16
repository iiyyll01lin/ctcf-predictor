# Consensus-based Alignment Script
# This script implements consensus-driven alignment strategy for sequence alignment

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Please install Biostrings: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/filtered_sequences.fasta"
output_file <- if (length(args) >= 2) args[2] else "data/consensus_aligned_sequences.fasta"
target_length <- if (length(args) >= 3) as.numeric(args[3]) else 200
consensus_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 0.8

cat("Consensus-based Sequence Alignment\n")
cat("==================================\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Target length:", target_length, "bp\n")
cat("Consensus threshold:", consensus_threshold, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Functions ---

# Find consensus motif pattern in sequences
find_consensus_motif <- function(sequences, motif_length = 15, min_support = 0.3) {
  cat("Finding consensus motif pattern...\n")
  cat("Motif length:", motif_length, "bp\n")
  cat("Minimum support:", min_support * 100, "%\n")
  
  # Sample sequences for motif discovery (for computational efficiency)
  sample_size <- min(2000, length(sequences))
  sample_idx <- sample(length(sequences), sample_size)
  sample_seqs <- sequences[sample_idx]
  
  # Extract all possible k-mers of motif_length
  all_kmers <- character()
  
  for (i in 1:length(sample_seqs)) {
    seq_char <- as.character(sample_seqs[i])
    seq_len <- nchar(seq_char)
    
    if (seq_len >= motif_length) {
      for (start in 1:(seq_len - motif_length + 1)) {
        kmer <- substr(seq_char, start, start + motif_length - 1)
        # Only include k-mers without N bases
        if (!grepl("N", kmer)) {
          all_kmers <- c(all_kmers, kmer)
        }
      }
    }
  }
  
  if (length(all_kmers) == 0) {
    cat("Warning: No valid k-mers found\n")
    return(NULL)
  }
  
  # Count k-mer frequencies
  kmer_counts <- table(all_kmers)
  kmer_freqs <- kmer_counts / length(sample_seqs)
  
  # Find k-mers above threshold
  frequent_kmers <- kmer_freqs[kmer_freqs >= min_support]
  
  if (length(frequent_kmers) == 0) {
    cat("Warning: No k-mers found above support threshold\n")
    return(NULL)
  }
  
  # Select most frequent k-mer as consensus
  consensus_motif <- names(frequent_kmers)[which.max(frequent_kmers)]
  consensus_freq <- max(frequent_kmers)
  
  cat("Consensus motif found:", consensus_motif, "\n")
  cat("Motif frequency:", round(consensus_freq * 100, 1), "%\n")
  
  return(list(
    motif = consensus_motif,
    frequency = consensus_freq,
    length = motif_length
  ))
}

# Find motif position in a sequence
find_motif_position <- function(sequence, motif_pattern, allow_mismatches = 1) {
  seq_char <- as.character(sequence)
  pattern_len <- nchar(motif_pattern)
  seq_len <- nchar(seq_char)
  
  if (seq_len < pattern_len) {
    return(NA)
  }
  
  best_pos <- NA
  best_score <- -1
  
  # Scan sequence for best match
  for (start in 1:(seq_len - pattern_len + 1)) {
    subseq <- substr(seq_char, start, start + pattern_len - 1)
    
    # Calculate similarity score
    matches <- sum(strsplit(subseq, "")[[1]] == strsplit(motif_pattern, "")[[1]])
    score <- matches / pattern_len
    
    if (score > best_score && matches >= (pattern_len - allow_mismatches)) {
      best_score <- score
      best_pos <- start
    }
  }
  
  return(best_pos)
}

# Build position weight matrix from sequences
build_consensus_pwm <- function(sequences, motif_length = 15) {
  cat("Building consensus PWM...\n")
  
  # Sample sequences for PWM building
  sample_size <- min(5000, length(sequences))
  sample_idx <- sample(length(sequences), sample_size)
  sample_seqs <- sequences[sample_idx]
  
  # Extract motif regions and build PWM
  motif_regions <- character()
  
  for (seq in sample_seqs) {
    seq_char <- as.character(seq)
    seq_len <- nchar(seq_char)
    
    if (seq_len >= motif_length) {
      # Use center region as motif
      center_start <- max(1, (seq_len - motif_length) %/% 2 + 1)
      center_end <- center_start + motif_length - 1
      motif_region <- substr(seq_char, center_start, center_end)
      
      if (!grepl("N", motif_region)) {
        motif_regions <- c(motif_regions, motif_region)
      }
    }
  }
  
  if (length(motif_regions) == 0) {
    return(NULL)
  }
  
  # Convert to matrix
  motif_matrix <- do.call(rbind, strsplit(motif_regions, ""))
  
  # Build PWM
  pwm <- matrix(0, nrow = 4, ncol = motif_length)
  rownames(pwm) <- c("A", "C", "G", "T")
  
  for (pos in 1:motif_length) {
    base_counts <- table(factor(motif_matrix[, pos], levels = c("A", "C", "G", "T")))
    base_freqs <- base_counts / sum(base_counts)
    pwm[, pos] <- as.numeric(base_freqs)
  }
  
  return(pwm)
}

# Score sequence against PWM
score_sequence_pwm <- function(sequence, pwm, start_pos = 1) {
  seq_char <- as.character(sequence)
  motif_length <- ncol(pwm)
  
  if (nchar(seq_char) < start_pos + motif_length - 1) {
    return(-Inf)
  }
  
  motif_seq <- substr(seq_char, start_pos, start_pos + motif_length - 1)
  bases <- strsplit(motif_seq, "")[[1]]
  
  score <- 0
  for (i in 1:length(bases)) {
    if (bases[i] %in% c("A", "C", "G", "T")) {
      # Log-likelihood score
      prob <- pwm[bases[i], i]
      if (prob > 0) {
        score <- score + log2(prob / 0.25)  # Background probability = 0.25
      }
    }
  }
  
  return(score)
}

# Perform consensus-based alignment
perform_consensus_alignment <- function(sequences, target_len = 200, threshold = 0.8) {
  cat("Performing consensus-based alignment...\n")
  
  # Step 1: Find consensus motif
  consensus_info <- find_consensus_motif(sequences, motif_length = 15, min_support = 0.1)
  
  if (is.null(consensus_info)) {
    cat("Warning: Could not find consensus motif, falling back to center alignment\n")
    return(perform_center_fallback(sequences, target_len))
  }
  
  # Step 2: Build PWM from consensus
  consensus_pwm <- build_consensus_pwm(sequences, consensus_info$length)
  
  if (is.null(consensus_pwm)) {
    cat("Warning: Could not build consensus PWM, falling back to center alignment\n")
    return(perform_center_fallback(sequences, target_len))
  }
  
  # Step 3: Align sequences based on motif positions
  aligned_sequences <- DNAStringSet()
  alignment_stats <- list(
    motif_found = 0,
    center_fallback = 0
  )
  
  for (i in 1:length(sequences)) {
    if (i %% 1000 == 0) {
      cat("Aligning sequence", i, "/", length(sequences), "\n")
    }
    
    seq <- sequences[i]
    seq_len <- width(seq)
    
    # Find best motif position
    best_pos <- NA
    best_score <- -Inf
    
    if (seq_len >= consensus_info$length) {
      for (start in 1:(seq_len - consensus_info$length + 1)) {
        score <- score_sequence_pwm(seq, consensus_pwm, start)
        if (score > best_score) {
          best_score <- score
          best_pos <- start
        }
      }
    }
    
    # Align based on motif position
    if (!is.na(best_pos) && best_score > log2(threshold)) {
      # Calculate alignment window around motif
      motif_center <- best_pos + consensus_info$length %/% 2
      align_start <- max(1, motif_center - target_len %/% 2)
      align_end <- min(seq_len, align_start + target_len - 1)
      
      # Adjust if we don't have enough sequence
      if (align_end - align_start + 1 < target_len) {
        if (align_start == 1) {
          align_end <- min(seq_len, target_len)
        } else {
          align_start <- max(1, seq_len - target_len + 1)
          align_end <- seq_len
        }
      }
      
      aligned_seq <- subseq(seq, align_start, align_end)
      alignment_stats$motif_found <- alignment_stats$motif_found + 1
    } else {
      # Fallback to center alignment
      aligned_seq <- extract_center_region(seq, target_len)
      alignment_stats$center_fallback <- alignment_stats$center_fallback + 1
    }
    
    # Pad if necessary
    if (width(aligned_seq) < target_len) {
      aligned_seq <- pad_sequence(aligned_seq, target_len)
    }
    
    aligned_sequences <- c(aligned_sequences, aligned_seq)
  }
  
  # Copy names
  names(aligned_sequences) <- names(sequences)
  
  cat("Consensus alignment completed:\n")
  cat("  Motif-based alignments:", alignment_stats$motif_found, "\n")
  cat("  Center fallback alignments:", alignment_stats$center_fallback, "\n")
  
  return(list(
    sequences = aligned_sequences,
    consensus_motif = consensus_info$motif,
    motif_frequency = consensus_info$frequency,
    alignment_stats = alignment_stats
  ))
}

# Center alignment fallback
perform_center_fallback <- function(sequences, target_len) {
  aligned_sequences <- DNAStringSet()
  
  for (seq in sequences) {
    aligned_seq <- extract_center_region(seq, target_len)
    if (width(aligned_seq) < target_len) {
      aligned_seq <- pad_sequence(aligned_seq, target_len)
    }
    aligned_sequences <- c(aligned_sequences, aligned_seq)
  }
  
  names(aligned_sequences) <- names(sequences)
  
  return(list(
    sequences = aligned_sequences,
    consensus_motif = NULL,
    motif_frequency = 0,
    alignment_stats = list(motif_found = 0, center_fallback = length(sequences))
  ))
}

# Extract center region (from align_center.R)
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

# Pad sequence (from align_center.R)
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
    "bp, range =", min(initial_lengths), "-", max(initial_lengths), "bp\n\n")

# Perform consensus-based alignment
alignment_result <- perform_consensus_alignment(input_sequences, target_length, consensus_threshold)
aligned_sequences <- alignment_result$sequences

# Create output directory if needed
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save aligned sequences
cat("\nSaving aligned sequences to:", output_file, "\n")
writeXStringSet(aligned_sequences, output_file)

# Generate alignment report
report_file <- gsub("\\.fasta$", "_consensus_alignment_report.txt", output_file)
cat("Generating alignment report:", report_file, "\n")

sink(report_file)
cat("Consensus-based Alignment Report\n")
cat("================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n\n")

cat("Alignment Parameters:\n")
cat("  Target length:", target_length, "bp\n")
cat("  Consensus threshold:", consensus_threshold, "\n\n")

cat("Alignment Results:\n")
cat("  Input sequences:", length(input_sequences), "\n")
cat("  Aligned sequences:", length(aligned_sequences), "\n")

if (!is.null(alignment_result$consensus_motif)) {
  cat("  Consensus motif:", alignment_result$consensus_motif, "\n")
  cat("  Motif frequency:", round(alignment_result$motif_frequency * 100, 1), "%\n")
}

if (!is.null(alignment_result$alignment_stats)) {
  cat("  Motif-based alignments:", alignment_result$alignment_stats$motif_found, "\n")
  cat("  Center fallback alignments:", alignment_result$alignment_stats$center_fallback, "\n")
}

# Final sequence statistics
final_lengths <- width(aligned_sequences)
cat("\nFinal sequence statistics:\n")
cat("  Mean length:", round(mean(final_lengths), 1), "bp\n")
cat("  Length consistency:", length(unique(final_lengths)) == 1, "\n")

sink()

cat("\nConsensus-based alignment completed successfully!\n")
cat("Aligned", length(aligned_sequences), "sequences saved to", output_file, "\n")
