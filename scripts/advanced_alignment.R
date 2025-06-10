# Advanced Sequence Alignment and PWM Building
# This script implements multiple alignment strategies and PWM building approaches

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/training_sequences.fasta"
output_prefix <- if (length(args) >= 2) args[2] else "results/advanced_alignment"
alignment_method <- if (length(args) >= 3) args[3] else "consensus"
min_coverage <- if (length(args) >= 4) as.numeric(args[4]) else 0.5

cat("Advanced Sequence Alignment and PWM Builder\n")
cat("===========================================\n")
cat("Input file:", input_file, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Alignment method:", alignment_method, "\n")
cat("Minimum coverage:", min_coverage, "\n\n")

# --- Functions ---

# Calculate sequence similarity matrix
calculate_similarity_matrix <- function(sequences, sample_size = 1000) {
  cat("Calculating sequence similarity matrix...\n")
  
  # Sample sequences if too many
  if (length(sequences) > sample_size) {
    sample_indices <- sample(length(sequences), sample_size)
    sequences <- sequences[sample_indices]
    cat("Using sample of", length(sequences), "sequences for similarity calculation\n")
  }
  
  n_seqs <- length(sequences)
  similarity_matrix <- matrix(0, nrow = n_seqs, ncol = n_seqs)
  
  for (i in 1:(n_seqs-1)) {
    for (j in (i+1):n_seqs) {
      # Calculate pairwise alignment score
      seq1 <- sequences[[i]]
      seq2 <- sequences[[j]]
      
      # Simple similarity: count matches at each position
      min_length <- min(length(seq1), length(seq2))
      matches <- sum(seq1[1:min_length] == seq2[1:min_length])
      similarity <- matches / min_length
      
      similarity_matrix[i, j] <- similarity_matrix[j, i] <- similarity
    }
    
    if (i %% 100 == 0) {
      cat("Processed", i, "of", n_seqs-1, "sequences\n")
    }
  }
  
  diag(similarity_matrix) <- 1
  return(similarity_matrix)
}

# Consensus-based alignment
align_by_consensus <- function(sequences) {
  cat("Performing consensus-based alignment...\n")
  
  # Find most common length
  lengths <- width(sequences)
  common_length <- as.numeric(names(sort(table(lengths), decreasing = TRUE))[1])
  cat("Most common length:", common_length, "\n")
  
  # Filter sequences by common length ± tolerance
  length_tolerance <- round(common_length * 0.1)  # 10% tolerance
  target_sequences <- sequences[abs(lengths - common_length) <= length_tolerance]
  
  cat("Sequences within tolerance:", length(target_sequences), "\n")
  
  if (length(target_sequences) < 100) {
    cat("Too few sequences, using length-based alignment instead\n")
    return(align_by_length(sequences))
  }
  
  # Trim or pad sequences to common length
  aligned_sequences <- DNAStringSet(lapply(target_sequences, function(seq) {
    current_length <- length(seq)
    if (current_length > common_length) {
      # Trim from center
      start_pos <- ceiling((current_length - common_length) / 2) + 1
      end_pos <- start_pos + common_length - 1
      return(seq[start_pos:end_pos])
    } else if (current_length < common_length) {
      # Pad with N's at both ends
      pad_total <- common_length - current_length
      pad_left <- floor(pad_total / 2)
      pad_right <- ceiling(pad_total / 2)
      padded <- c(rep("N", pad_left), as.character(seq), rep("N", pad_right))
      return(DNAString(paste(padded, collapse = "")))
    } else {
      return(seq)
    }
  }))
  
  return(aligned_sequences)
}

# Length-based alignment
align_by_length <- function(sequences) {
  cat("Performing length-based alignment...\n")
  
  lengths <- width(sequences)
  median_length <- median(lengths)
  cat("Median length:", median_length, "\n")
  
  # Align all sequences to median length
  aligned_sequences <- DNAStringSet(lapply(sequences, function(seq) {
    current_length <- length(seq)
    
    if (current_length > median_length) {
      # Trim from both ends equally
      trim_total <- current_length - median_length
      trim_left <- floor(trim_total / 2)
      start_pos <- trim_left + 1
      end_pos <- start_pos + median_length - 1
      return(seq[start_pos:end_pos])
    } else if (current_length < median_length) {
      # Pad with N's
      pad_total <- median_length - current_length
      pad_left <- floor(pad_total / 2)
      pad_right <- ceiling(pad_total / 2)
      padded <- c(rep("N", pad_left), as.character(seq), rep("N", pad_right))
      return(DNAString(paste(padded, collapse = "")))
    } else {
      return(seq)
    }
  }))
  
  return(aligned_sequences)
}

# Progressive alignment using similarity
align_progressive <- function(sequences, max_seqs = 5000) {
  cat("Performing progressive alignment...\n")
  
  # Limit sequences for computational efficiency
  if (length(sequences) > max_seqs) {
    sample_indices <- sample(length(sequences), max_seqs)
    sequences <- sequences[sample_indices]
    cat("Using sample of", length(sequences), "sequences\n")
  }
  
  # Calculate similarity matrix
  similarity_matrix <- calculate_similarity_matrix(sequences)
  
  # Find most representative sequence (highest average similarity)
  avg_similarities <- rowMeans(similarity_matrix)
  reference_idx <- which.max(avg_similarities)
  reference_seq <- sequences[[reference_idx]]
  reference_length <- length(reference_seq)
  
  cat("Reference sequence index:", reference_idx, "length:", reference_length, "\n")
  
  # Align all sequences to reference
  aligned_sequences <- DNAStringSet(lapply(sequences, function(seq) {
    current_length <- length(seq)
    
    if (current_length != reference_length) {
      if (current_length > reference_length) {
        # Find best alignment position
        best_score <- -Inf
        best_start <- 1
        
        for (start_pos in 1:(current_length - reference_length + 1)) {
          end_pos <- start_pos + reference_length - 1
          subseq <- seq[start_pos:end_pos]
          
          # Calculate match score
          matches <- sum(subseq == reference_seq)
          score <- matches / reference_length
          
          if (score > best_score) {
            best_score <- score
            best_start <- start_pos
          }
        }
        
        # Extract best matching subsequence
        best_end <- best_start + reference_length - 1
        return(seq[best_start:best_end])
      } else {
        # Pad shorter sequence
        pad_total <- reference_length - current_length
        pad_left <- floor(pad_total / 2)
        pad_right <- ceiling(pad_total / 2)
        padded <- c(rep("N", pad_left), as.character(seq), rep("N", pad_right))
        return(DNAString(paste(padded, collapse = "")))
      }
    } else {
      return(seq)
    }
  }))
  
  return(aligned_sequences)
}

# Filter positions by coverage
filter_low_coverage_positions <- function(aligned_sequences, min_coverage = 0.5) {
  cat("Filtering positions by coverage (min:", min_coverage, ")...\n")
  
  # Calculate coverage per position
  consensus_matrix <- consensusMatrix(aligned_sequences, as.prob = FALSE)
  total_sequences <- length(aligned_sequences)
  
  # Count non-N bases per position
  n_counts <- if ("N" %in% rownames(consensus_matrix)) consensus_matrix["N", ] else rep(0, ncol(consensus_matrix))
  coverage <- (total_sequences - n_counts) / total_sequences
  
  # Find positions with sufficient coverage
  good_positions <- which(coverage >= min_coverage)
  cat("Positions with coverage >=", min_coverage, ":", length(good_positions), "out of", length(coverage), "\n")
  
  if (length(good_positions) == 0) {
    warning("No positions meet coverage criteria. Using all positions.")
    return(aligned_sequences)
  }
  
  # Extract high-coverage positions
  filtered_sequences <- DNAStringSet(lapply(aligned_sequences, function(seq) {
    seq[good_positions]
  }))
  
  return(filtered_sequences)
}

# Calculate information content
calculate_info_content <- function(pwm) {
  apply(pwm, 2, function(x) {
    x[x == 0] <- 1e-10
    2 + sum(x * log2(x))
  })
}

# Build PWM with advanced features
build_advanced_pwm <- function(aligned_sequences, pseudocount = 0.1, filter_coverage = TRUE, min_coverage = 0.5) {
  cat("Building advanced PWM...\n")
  
  original_length <- width(aligned_sequences)[1]
  
  # Filter by coverage if requested
  if (filter_coverage) {
    aligned_sequences <- filter_low_coverage_positions(aligned_sequences, min_coverage)
  }
  
  final_length <- width(aligned_sequences)[1]
  cat("PWM length after filtering:", final_length, "(was", original_length, ")\n")
  
  # Build consensus matrix
  consensus_matrix <- consensusMatrix(aligned_sequences, as.prob = FALSE)
  nucleotides <- c("A", "C", "G", "T")
  count_matrix <- consensus_matrix[nucleotides, , drop = FALSE]
  
  # Add pseudocounts and normalize
  count_matrix <- count_matrix + pseudocount
  pwm <- apply(count_matrix, 2, function(x) x / sum(x))
  
  # Calculate information content
  info_content <- calculate_info_content(pwm)
  total_info <- sum(info_content)
  
  cat("Total information content:", round(total_info, 3), "bits\n")
  cat("Average per position:", round(total_info / ncol(pwm), 3), "bits\n")
  
  # Analyze motif structure
  conserved_positions <- which(info_content > 1.0)
  highly_conserved <- which(info_content > 1.5)
  
  cat("Conserved positions (>1 bit):", length(conserved_positions), "\n")
  cat("Highly conserved positions (>1.5 bit):", length(highly_conserved), "\n")
  
  # Find contiguous conserved regions
  if (length(conserved_positions) > 0) {
    conserved_regions <- find_conserved_regions(conserved_positions)
    cat("Conserved regions found:", length(conserved_regions), "\n")
  } else {
    conserved_regions <- list()
  }
  
  return(list(
    pwm = pwm,
    count_matrix = count_matrix,
    info_content = info_content,
    total_info = total_info,
    conserved_positions = conserved_positions,
    highly_conserved_positions = highly_conserved,
    conserved_regions = conserved_regions,
    num_sequences = length(aligned_sequences),
    original_length = original_length,
    final_length = final_length,
    pseudocount = pseudocount,
    min_coverage = min_coverage,
    method = "advanced_alignment"
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
      # Extend current region
      current_end <- positions[i]
    } else {
      # Save current region and start new one
      regions[[length(regions) + 1]] <- list(start = current_start, end = current_end, length = current_end - current_start + 1)
      current_start <- positions[i]
      current_end <- positions[i]
    }
  }
  
  # Save the last region
  regions[[length(regions) + 1]] <- list(start = current_start, end = current_end, length = current_end - current_start + 1)
  
  return(regions)
}

# --- Main Execution ---

# Check input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Read sequences
cat("Reading sequences from:", input_file, "\n")
sequences <- readDNAStringSet(input_file)

if (length(sequences) == 0) {
  stop("No sequences found in input file")
}

cat("Loaded", length(sequences), "sequences\n")

# Perform alignment based on selected method
aligned_sequences <- switch(alignment_method,
  "consensus" = align_by_consensus(sequences),
  "length" = align_by_length(sequences),
  "progressive" = align_progressive(sequences),
  {
    cat("Unknown alignment method. Using consensus alignment.\n")
    align_by_consensus(sequences)
  }
)

cat("Alignment completed. Aligned", length(aligned_sequences), "sequences\n")

# Build PWM with multiple approaches
approaches <- list(
  "basic" = list(filter_coverage = FALSE, min_coverage = 0.5),
  "filtered" = list(filter_coverage = TRUE, min_coverage = 0.5),
  "strict" = list(filter_coverage = TRUE, min_coverage = 0.8)
)

results <- list()

for (approach_name in names(approaches)) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("Building PWM with approach:", approach_name, "\n")
  
  approach_params <- approaches[[approach_name]]
  
  pwm_result <- build_advanced_pwm(
    aligned_sequences,
    filter_coverage = approach_params$filter_coverage,
    min_coverage = approach_params$min_coverage
  )
  
  pwm_result$approach = approach_name
  pwm_result$alignment_method = alignment_method
  pwm_result$creation_time = Sys.time()
  
  # Save individual result
  output_file <- paste0(output_prefix, "_", alignment_method, "_", approach_name, ".rds")
  saveRDS(pwm_result, output_file)
  cat("PWM saved to:", output_file, "\n")
  
  results[[approach_name]] <- pwm_result
}

# Save aligned sequences
aligned_output <- paste0(output_prefix, "_", alignment_method, "_aligned.fasta")
writeXStringSet(aligned_sequences, aligned_output)
cat("Aligned sequences saved to:", aligned_output, "\n")

# Save combined results
combined_output <- paste0(output_prefix, "_", alignment_method, "_all.rds")
saveRDS(results, combined_output)
cat("All results saved to:", combined_output, "\n")

# Generate comparison report
report_file <- paste0(output_prefix, "_", alignment_method, "_comparison.txt")
sink(report_file)

cat("=== Advanced Alignment PWM Comparison Report ===\n")
cat("Created:", as.character(Sys.time()), "\n")
cat("Input file:", input_file, "\n")
cat("Alignment method:", alignment_method, "\n")
cat("Number of sequences:", length(sequences), "\n")
cat("Aligned sequences:", length(aligned_sequences), "\n\n")

cat("Approach Comparison:\n")
cat(sprintf("%-10s | %-12s | %-10s | %-12s | %-8s\n", 
           "Approach", "Total Info", "Avg Info", "Conserved", "Length"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (approach_name in names(results)) {
  result <- results[[approach_name]]
  cat(sprintf("%-10s | %12.3f | %10.3f | %12d | %8d\n",
             approach_name,
             result$total_info,
             result$total_info / ncol(result$pwm),
             length(result$conserved_positions),
             result$final_length))
}

cat("\nDetailed Results:\n")
for (approach_name in names(results)) {
  result <- results[[approach_name]]
  cat("\n--- ", approach_name, " ---\n")
  cat("PWM dimensions:", nrow(result$pwm), "x", ncol(result$pwm), "\n")
  cat("Original length:", result$original_length, "\n")
  cat("Final length:", result$final_length, "\n")
  cat("Total information:", round(result$total_info, 3), "bits\n")
  cat("Average per position:", round(result$total_info / ncol(result$pwm), 3), "bits\n")
  cat("Conserved positions:", length(result$conserved_positions), "\n")
  cat("Highly conserved positions:", length(result$highly_conserved_positions), "\n")
  cat("Conserved regions:", length(result$conserved_regions), "\n")
  
  if (length(result$conserved_regions) > 0) {
    cat("Region details:\n")
    for (i in seq_along(result$conserved_regions)) {
      region <- result$conserved_regions[[i]]
      cat("  Region", i, ": positions", region$start, "-", region$end, 
          "(length", region$length, ")\n")
    }
  }
}

sink()

cat("\nComparison report saved to:", report_file, "\n")

# Print final summary
cat("\n=== Final Summary ===\n")
cat("Alignment method:", alignment_method, "\n")
cat("Input sequences:", length(sequences), "\n")
cat("Aligned sequences:", length(aligned_sequences), "\n")

# Find best approach
best_approach <- names(results)[which.max(sapply(results, function(x) x$total_info))]
best_result <- results[[best_approach]]

cat("Best approach:", best_approach, "\n")
cat("Best total information:", round(best_result$total_info, 3), "bits\n")
cat("Best average per position:", round(best_result$total_info / ncol(best_result$pwm), 3), "bits\n")

cat("\nScript completed successfully!\n")