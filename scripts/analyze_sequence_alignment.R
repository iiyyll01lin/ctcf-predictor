# Sequence Alignment Analysis and Improvement Script

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/training_sequences.fasta"
output_file <- if (length(args) >= 2) args[2] else "data/aligned_sequences.fasta"
alignment_method <- if (length(args) >= 3) args[3] else "center"  # "center", "left", "right", "consensus"

# --- Utility Functions ---

# Extract central region of specified length
extract_center <- function(seq, target_length) {
  seq_len <- nchar(seq)
  if (seq_len <= target_length) {
    return(seq)
  }
  
  start_pos <- floor((seq_len - target_length) / 2) + 1
  end_pos <- start_pos + target_length - 1
  return(substr(seq, start_pos, end_pos))
}

# Find the most common length
find_consensus_length <- function(lengths) {
  length_table <- table(lengths)
  most_common <- as.numeric(names(length_table)[which.max(length_table)])
  return(most_common)
}

# Simple scoring function for motif finding
score_position <- function(sequences, pos, window = 11) {
  if (pos + window - 1 > min(nchar(sequences))) {
    return(0)
  }
  
  subseqs <- substr(sequences, pos, pos + window - 1)
  
  # Calculate information content of this window
  position_matrix <- matrix(0, nrow = 4, ncol = window)
  rownames(position_matrix) <- c("A", "C", "G", "T")
  
  for (i in 1:length(subseqs)) {
    seq_chars <- strsplit(subseqs[i], "")[[1]]
    for (j in 1:length(seq_chars)) {
      if (seq_chars[j] %in% rownames(position_matrix)) {
        position_matrix[seq_chars[j], j] <- position_matrix[seq_chars[j], j] + 1
      }
    }
  }
  
  # Convert to probabilities and calculate information content
  total_seqs <- length(sequences)
  info_content <- 0
  
  for (j in 1:window) {
    col_sum <- sum(position_matrix[, j])
    if (col_sum > 0) {
      probs <- position_matrix[, j] / col_sum
      probs[probs == 0] <- 1e-10
      info_content <- info_content + sum(probs * log2(probs / 0.25))
    }
  }
  
  return(info_content)
}

# Find optimal motif position in sequences
find_motif_positions <- function(sequences, motif_length = 11) {
  sequence_strings <- as.character(sequences)
  min_seq_len <- min(nchar(sequence_strings))
  
  if (min_seq_len < motif_length) {
    stop("Sequences too short for motif length ", motif_length)
  }
  
  best_positions <- numeric(length(sequences))
  
  # For each sequence, find the position that best matches a consensus
  for (i in 1:length(sequence_strings)) {
    seq_len <- nchar(sequence_strings[i])
    best_score <- -Inf
    best_pos <- 1
    
    # Try each possible position
    for (pos in 1:(seq_len - motif_length + 1)) {
      # Score this position by comparing with all other sequences
      total_score <- 0
      current_subseq <- substr(sequence_strings[i], pos, pos + motif_length - 1)
      
      for (j in 1:length(sequence_strings)) {
        if (i != j) {
          other_seq_len <- nchar(sequence_strings[j])
          # Find best match in other sequence
          for (other_pos in 1:(other_seq_len - motif_length + 1)) {
            other_subseq <- substr(sequence_strings[j], other_pos, other_pos + motif_length - 1)
            # Simple similarity score
            matches <- sum(strsplit(current_subseq, "")[[1]] == strsplit(other_subseq, "")[[1]])
            total_score <- total_score + matches
          }
        }
      }
      
      if (total_score > best_score) {
        best_score <- total_score
        best_pos <- pos
      }
    }
    
    best_positions[i] <- best_pos
  }
  
  return(best_positions)
}

# --- Main Analysis ---
cat("=== Sequence Alignment Analysis ===\n")
cat("Input file:", input_file, "\n")
cat("Alignment method:", alignment_method, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Validate input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Read sequences
cat("Reading sequences...\n")
seqs <- readDNAStringSet(input_file)
cat("Total sequences:", length(seqs), "\n")

if (length(seqs) == 0) {
  stop("No sequences found in the input file.")
}

sequence_strings <- as.character(seqs)
lengths <- nchar(sequence_strings)

# 1. Length Analysis
cat("\n1. LENGTH DISTRIBUTION ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

cat("Length statistics:\n")
cat("Min:", min(lengths), "bp\n")
cat("Max:", max(lengths), "bp\n")
cat("Mean:", round(mean(lengths), 1), "bp\n")
cat("Median:", median(lengths), "bp\n")
cat("Std Dev:", round(sd(lengths), 1), "bp\n")
cat("Coefficient of Variation:", round(sd(lengths)/mean(lengths), 3), "\n\n")

# Length distribution
unique_lengths <- unique(lengths)
cat("Number of unique lengths:", length(unique_lengths), "\n")

if (length(unique_lengths) <= 10) {
  cat("Length distribution:\n")
  length_table <- table(lengths)
  for (i in 1:length(length_table)) {
    len <- names(length_table)[i]
    count <- length_table[i]
    percent <- round(count / length(seqs) * 100, 1)
    cat("  ", len, "bp:", count, "sequences (", percent, "%)\n")
  }
} else {
  cat("Too many unique lengths to display. Top 5:\n")
  length_table <- sort(table(lengths), decreasing = TRUE)
  for (i in 1:min(5, length(length_table))) {
    len <- names(length_table)[i]
    count <- length_table[i]
    percent <- round(count / length(seqs) * 100, 1)
    cat("  ", len, "bp:", count, "sequences (", percent, "%)\n")
  }
}

# 2. Alignment Strategy Analysis
cat("\n2. ALIGNMENT STRATEGY ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

consensus_length <- find_consensus_length(lengths)
cat("Most common length:", consensus_length, "bp\n")

# Assess alignment need
cv <- sd(lengths) / mean(lengths)
if (cv < 0.1) {
  cat("✅ Length variation is minimal (CV < 0.1) - alignment may not be necessary\n")
  alignment_needed <- FALSE
} else if (cv < 0.3) {
  cat("⚠️  Moderate length variation (CV 0.1-0.3) - alignment recommended\n")
  alignment_needed <- TRUE
} else {
  cat("❌ High length variation (CV > 0.3) - alignment strongly recommended\n")
  alignment_needed <- TRUE
}

# 3. Perform Alignment Based on Method
cat("\n3. PERFORMING ALIGNMENT\n")
cat(paste(rep("=", 50), collapse=""), "\n")

aligned_sequences <- character(length(seqs))
target_length <- NULL

if (!alignment_needed && alignment_method == "center") {
  cat("Alignment not needed - using original sequences\n")
  aligned_sequences <- sequence_strings
  aligned_seqs <- seqs
} else {
  
  if (alignment_method == "center") {
    cat("Method: Extract central regions\n")
    target_length <- consensus_length
    cat("Target length:", target_length, "bp\n")
    
    for (i in 1:length(sequence_strings)) {
      aligned_sequences[i] <- extract_center(sequence_strings[i], target_length)
    }
    
  } else if (alignment_method == "left") {
    cat("Method: Left-align (take first N bases)\n")
    target_length <- min(lengths)
    cat("Target length:", target_length, "bp\n")
    
    for (i in 1:length(sequence_strings)) {
      aligned_sequences[i] <- substr(sequence_strings[i], 1, target_length)
    }
    
  } else if (alignment_method == "right") {
    cat("Method: Right-align (take last N bases)\n")
    target_length <- min(lengths)
    cat("Target length:", target_length, "bp\n")
    
    for (i in 1:length(sequence_strings)) {
      seq_len <- nchar(sequence_strings[i])
      start_pos <- seq_len - target_length + 1
      aligned_sequences[i] <- substr(sequence_strings[i], start_pos, seq_len)
    }
    
  } else if (alignment_method == "consensus") {
    cat("Method: Consensus motif-based alignment\n")
    target_length <- min(consensus_length, min(lengths))
    cat("Target length:", target_length, "bp\n")
    
    # Find best motif positions
    motif_positions <- find_motif_positions(seqs, target_length)
    
    for (i in 1:length(sequence_strings)) {
      start_pos <- motif_positions[i]
      end_pos <- start_pos + target_length - 1
      aligned_sequences[i] <- substr(sequence_strings[i], start_pos, end_pos)
    }
    
  } else {
    stop("Unknown alignment method: ", alignment_method)
  }
  
  # Create aligned DNAStringSet
  aligned_seqs <- DNAStringSet(aligned_sequences)
  names(aligned_seqs) <- paste0(names(seqs), "_aligned")
}

# 4. Validate Alignment Results
cat("\n4. ALIGNMENT VALIDATION\n")
cat(paste(rep("=", 50), collapse=""), "\n")

aligned_lengths <- nchar(aligned_sequences)
cat("Aligned sequence statistics:\n")
cat("Number of sequences:", length(aligned_sequences), "\n")
cat("Aligned length range:", min(aligned_lengths), "-", max(aligned_lengths), "bp\n")
cat("Unique lengths after alignment:", length(unique(aligned_lengths)), "\n")

if (length(unique(aligned_lengths)) == 1) {
  cat("✅ Perfect alignment achieved - all sequences same length\n")
} else {
  cat("⚠️  Imperfect alignment - still have length variation\n")
}

# Calculate information gain from alignment
if (alignment_needed) {
  cat("\nAlignment improvement:\n")
  original_cv <- sd(lengths) / mean(lengths)
  aligned_cv <- sd(aligned_lengths) / mean(aligned_lengths)
  cat("Original CV:", round(original_cv, 3), "\n")
  cat("Aligned CV:", round(aligned_cv, 3), "\n")
  
  if (aligned_cv < original_cv) {
    improvement <- ((original_cv - aligned_cv) / original_cv) * 100
    cat("✅ Alignment improved uniformity by", round(improvement, 1), "%\n")
  } else {
    cat("❌ Alignment did not improve uniformity\n")
  }
}

# 5. Save Aligned Sequences
cat("\n5. SAVING RESULTS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

if (alignment_needed || !file.exists(output_file)) {
  cat("Saving aligned sequences to:", output_file, "\n")
  writeXStringSet(aligned_seqs, filepath = output_file)
  cat("✅ Aligned sequences saved successfully\n")
} else {
  cat("No alignment needed - original sequences are already uniform\n")
}

# 6. Recommendations
cat("\n6. RECOMMENDATIONS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

cat("Based on alignment analysis:\n\n")

if (length(unique(aligned_lengths)) == 1) {
  cat("✅ ALIGNMENT: Perfect - proceed with PWM building\n")
  cat("   Use file:", output_file, "\n")
} else {
  cat("⚠️  ALIGNMENT: Imperfect but improved\n")
  cat("   Consider trying different alignment methods:\n")
  cat("   - 'consensus' for motif-based alignment\n")
  cat("   - 'center' for central region extraction\n")
  cat("   - Manual curation may be needed\n")
}

# Suggest next steps
cat("\nNext steps:\n")
cat("1. Validate alignment quality with: validate_pwm_quality.R\n")
cat("2. Build new PWM with aligned sequences\n")
cat("3. Compare PWM quality before and after alignment\n")

cat("\n=== Alignment Analysis Complete ===\n")
