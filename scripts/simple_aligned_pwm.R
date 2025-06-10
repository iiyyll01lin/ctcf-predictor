# Simple Aligned PWM Builder - Memory Efficient Implementation
# This script builds a PWM from aligned sequences with minimal memory usage

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
aligned_file <- if (length(args) >= 1) args[1] else "data/aligned_sequences.fasta"
output_file <- if (length(args) >= 2) args[2] else "results/simple_aligned_pwm.rds"
pseudocount <- if (length(args) >= 3) as.numeric(args[3]) else 0.1

cat("Simple Aligned PWM Builder\n")
cat("==========================\n")
cat("Input file:", aligned_file, "\n")
cat("Output file:", output_file, "\n")
cat("Pseudocount:", pseudocount, "\n\n")

# --- Functions ---

# Calculate information content for a position
calculate_position_info <- function(freqs) {
  # Remove zeros to avoid log(0)
  freqs[freqs == 0] <- 1e-10
  # Information content = 2 + sum(freq * log2(freq))
  return(2 + sum(freqs * log2(freqs)))
}

# Build PWM from aligned sequences
build_simple_aligned_pwm <- function(sequences, pseudocount = 0.1) {
  cat("Building PWM from", length(sequences), "aligned sequences...\n")
  
  # Get sequence lengths
  seq_lengths <- width(sequences)
  cat("Sequence length range:", min(seq_lengths), "-", max(seq_lengths), "\n")
  
  # Use consensus matrix for efficiency
  consensus_matrix <- consensusMatrix(sequences, as.prob = FALSE)
  
  # Keep only standard nucleotides (A, C, G, T)
  nucleotides <- c("A", "C", "G", "T")
  count_matrix <- consensus_matrix[nucleotides, , drop = FALSE]
  
  cat("PWM dimensions:", nrow(count_matrix), "x", ncol(count_matrix), "\n")
  
  # Add pseudocounts
  count_matrix <- count_matrix + pseudocount
  
  # Convert to frequencies
  pwm <- apply(count_matrix, 2, function(x) x / sum(x))
  
  # Calculate information content per position
  info_content <- apply(pwm, 2, calculate_position_info)
  total_info <- sum(info_content)
  
  cat("Total information content:", round(total_info, 3), "bits\n")
  cat("Average per position:", round(total_info / ncol(pwm), 3), "bits\n")
  
  # Find conserved positions (>1 bit)
  conserved_positions <- which(info_content > 1.0)
  cat("Conserved positions (>1 bit):", length(conserved_positions), "\n")
  if (length(conserved_positions) > 0) {
    cat("Conserved position indices:", paste(conserved_positions, collapse = ", "), "\n")
  }
  
  # Create result list
  result <- list(
    pwm = pwm,
    info_content = info_content,
    total_info = total_info,
    conserved_positions = conserved_positions,
    num_sequences = length(sequences),
    pseudocount = pseudocount,
    method = "simple_aligned",
    creation_time = Sys.time()
  )
  
  return(result)
}

# --- Main Execution ---

# Check if input file exists
if (!file.exists(aligned_file)) {
  stop("Input file not found: ", aligned_file)
}

# Read aligned sequences
cat("Reading aligned sequences from:", aligned_file, "\n")
sequences <- readDNAStringSet(aligned_file)

# Check if sequences were loaded
if (length(sequences) == 0) {
  stop("No sequences found in input file")
}

cat("Loaded", length(sequences), "sequences\n")

# Build PWM
pwm_result <- build_simple_aligned_pwm(sequences, pseudocount)

# Create output directory if needed
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Save PWM
saveRDS(pwm_result, output_file)
cat("PWM saved to:", output_file, "\n")

# Print summary
cat("\n=== PWM Summary ===\n")
cat("Dimensions:", nrow(pwm_result$pwm), "x", ncol(pwm_result$pwm), "\n")
cat("Total information content:", round(pwm_result$total_info, 3), "bits\n")
cat("Average per position:", round(pwm_result$total_info / ncol(pwm_result$pwm), 3), "bits\n")
cat("Conserved positions:", length(pwm_result$conserved_positions), "\n")
cat("Method: Simple aligned PWM\n")
cat("Pseudocount:", pwm_result$pseudocount, "\n")

cat("\nScript completed successfully!\n")