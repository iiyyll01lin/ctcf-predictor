# Robust PWM Building Script with Error Handling and Multiple Data Sources
# This script implements improved PWM building with better data handling

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/training_sequences.fasta"
output_prefix <- if (length(args) >= 2) args[2] else "results/robust_pwm"
min_sequences <- if (length(args) >= 3) as.numeric(args[3]) else 100

# --- Utility Functions ---

# Calculate information content for a PWM
calculate_info_content <- function(pwm) {
  apply(pwm, 2, function(x) {
    x[x == 0] <- 1e-10
    sum(x * log2(x/0.25))
  })
}

# Filter sequences by quality metrics
filter_high_quality_sequences <- function(sequences) {
  cat("Filtering sequences for quality...\n")
  
  # Remove sequences with too many N's
  n_counts <- sapply(sequences, function(seq) {
    sum(unlist(strsplit(as.character(seq), "")) == "N")
  })
  n_percent <- n_counts / width(sequences) * 100
  quality_seqs <- sequences[n_percent <= 5]  # Max 5% N's
  
  cat("Removed", length(sequences) - length(quality_seqs), "sequences with >5% N bases\n")
  
  # Filter by length consistency (keep most common length ± 10bp)
  lengths <- width(quality_seqs)
  most_common <- as.numeric(names(sort(table(lengths), decreasing = TRUE))[1])
  length_filtered <- quality_seqs[abs(width(quality_seqs) - most_common) <= 10]
  
  cat("Kept", length(length_filtered), "sequences with consistent length (", 
      most_common, "±10bp)\n")
  
  return(length_filtered)
}

# Build PWM with pseudocounts
build_robust_pwm <- function(sequences, pseudocount = 0.1) {
  cat("Building consensus matrix...\n")
  consensus_mat <- consensusMatrix(sequences, as.prob = FALSE)
  
  # Keep only standard bases
  valid_bases <- c("A", "C", "G", "T")
  consensus_mat <- consensus_mat[rownames(consensus_mat) %in% valid_bases, , drop = FALSE]
  
  # Add pseudocounts
  consensus_mat_pseudo <- consensus_mat + pseudocount
  
  # Convert to probabilities
  pwm <- prop.table(consensus_mat_pseudo, margin = 2)
  
  return(pwm)
}

# --- Main Script ---
cat("=== Robust PWM Building ===\n")
cat("Input file:", input_file, "\n")
cat("Minimum sequences required:", min_sequences, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Check input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Load sequences
cat("Loading sequences...\n")
all_sequences <- readDNAStringSet(input_file)
cat("Total sequences loaded:", length(all_sequences), "\n")

if (length(all_sequences) < min_sequences) {
  stop("Insufficient sequences: need at least ", min_sequences, 
       ", but only ", length(all_sequences), " available")
}

# Filter for high quality sequences
filtered_sequences <- filter_high_quality_sequences(all_sequences)

if (length(filtered_sequences) < min_sequences) {
  cat("WARNING: After quality filtering, only", length(filtered_sequences), 
      "sequences remain\n")
  cat("Using relaxed filtering to retain more sequences...\n")
  
  # Relaxed filtering: only remove sequences with >15% N's
  n_counts <- sapply(all_sequences, function(seq) {
    sum(unlist(strsplit(as.character(seq), "")) == "N")
  })
  n_percent <- n_counts / width(all_sequences) * 100
  filtered_sequences <- all_sequences[n_percent <= 15]
}

cat("Final sequence count for PWM building:", length(filtered_sequences), "\n")

# Build PWM with different pseudocount values
pseudocount_values <- c(0.01, 0.1, 0.5, 1.0)
best_pwm <- NULL
best_info <- 0
best_pseudocount <- 0.1

cat("\nTesting different pseudocount values:\n")
for (pc in pseudocount_values) {
  pwm <- build_robust_pwm(filtered_sequences, pseudocount = pc)
  info_content <- calculate_info_content(pwm)
  total_info <- sum(info_content)
  
  cat("Pseudocount", pc, ": Total info =", round(total_info, 3), "bits\n")
  
  if (total_info > best_info) {
    best_info <- total_info
    best_pwm <- pwm
    best_pseudocount <- pc
  }
}

cat("\nBest pseudocount:", best_pseudocount, "with", round(best_info, 3), "total bits\n")

# Calculate detailed statistics
info_content <- calculate_info_content(best_pwm)
high_info_positions <- sum(info_content > 1.0)

cat("\nPWM Quality Metrics:\n")
cat("- PWM dimensions:", nrow(best_pwm), "x", ncol(best_pwm), "\n")
cat("- Total information content:", round(sum(info_content), 3), "bits\n")
cat("- Average per position:", round(mean(info_content), 3), "bits\n")
cat("- High-information positions (>1 bit):", high_info_positions, "\n")
cat("- Maximum position info:", round(max(info_content), 3), "bits\n")

# Quality assessment
if (sum(info_content) > 15) {
  cat("✅ EXCELLENT PWM quality\n")
} else if (sum(info_content) > 10) {
  cat("✅ GOOD PWM quality\n")
} else if (sum(info_content) > 5) {
  cat("⚠️  FAIR PWM quality - consider more data\n")
} else {
  cat("❌ POOR PWM quality - need better alignment or more data\n")
}

# Save results
cat("\nSaving results...\n")
saveRDS(best_pwm, paste0(output_prefix, ".rds"))
write.table(best_pwm, paste0(output_prefix, ".txt"), 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Save metadata
metadata <- list(
  creation_time = as.character(Sys.time()),
  input_file = input_file,
  total_sequences = length(all_sequences),
  filtered_sequences = length(filtered_sequences),
  best_pseudocount = best_pseudocount,
  total_information = round(sum(info_content), 3),
  average_information = round(mean(info_content), 3),
  high_info_positions = high_info_positions,
  pwm_dimensions = paste(nrow(best_pwm), "x", ncol(best_pwm))
)

# Create simple JSON manually instead of using jsonlite
json_content <- paste0('{\n',
  '  "creation_time": "', metadata$creation_time, '",\n',
  '  "input_file": "', metadata$input_file, '",\n',
  '  "total_sequences": ', metadata$total_sequences, ',\n',
  '  "filtered_sequences": ', metadata$filtered_sequences, ',\n',
  '  "best_pseudocount": ', metadata$best_pseudocount, ',\n',
  '  "total_information": ', metadata$total_information, ',\n',
  '  "average_information": ', metadata$average_information, ',\n',
  '  "high_info_positions": ', metadata$high_info_positions, ',\n',
  '  "pwm_dimensions": "', metadata$pwm_dimensions, '"\n',
  '}'
)

writeLines(json_content, paste0(output_prefix, "_metadata.json"))

cat("Robust PWM saved to:", paste0(output_prefix, ".rds"), "\n")
cat("Metadata saved to:", paste0(output_prefix, "_metadata.json"), "\n")
cat("=== Robust PWM Building Complete ===\n")
