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
output_file <- if (length(args) >= 2) args[2] else "results/ctcf_pwm.meme"
output_format <- if (length(args) >= 3) args[3] else "meme"
pseudocount <- if (length(args) >= 4) as.numeric(args[4]) else 0.1
min_ic <- if (length(args) >= 5) as.numeric(args[5]) else 8.0
min_sequences <- if (length(args) >= 6) as.numeric(args[6]) else 100

# --- Utility Functions ---

# Calculate information content for a PWM
calculate_info_content <- function(pwm) {
  apply(pwm, 2, function(x) {
    x[x == 0] <- 1e-10
    sum(x * log2(x/0.25))
  })
}

# Calculate information content with custom background
calculate_information_content <- function(prob_matrix, background = NULL) {
  if (is.null(background)) {
    background <- c(0.25, 0.25, 0.25, 0.25)  # uniform background
  }
  
  ic_per_pos <- apply(prob_matrix, 2, function(pos_probs) {
    valid_probs <- pmax(pos_probs, 1e-10)  # avoid log(0)
    sum(valid_probs * log2(valid_probs / background))
  })
  
  return(ic_per_pos)
}

# Add pseudocounts to frequency matrix
add_pseudocounts <- function(freq_matrix, pseudocount) {
  return(freq_matrix + pseudocount)
}

# Normalize frequencies to probabilities
normalize_probabilities <- function(pseudo_matrix) {
  return(prop.table(pseudo_matrix, margin = 2))
}

# Calculate frequencies from sequences
calculate_frequencies <- function(sequences) {
  consensus_mat <- consensusMatrix(sequences, as.prob = FALSE)
  valid_bases <- c("A", "C", "G", "T")
  consensus_mat <- consensus_mat[rownames(consensus_mat) %in% valid_bases, , drop = FALSE]
  return(consensus_mat)
}

# Assess PWM quality based on multiple criteria
assess_pwm_quality <- function(ic_vector, config) {
  total_ic <- sum(ic_vector)
  conserved_positions <- sum(ic_vector > 1.0)
  max_ic <- max(ic_vector)
  mean_ic <- mean(ic_vector)
  
  quality_metrics <- list(
    total_ic = total_ic,
    conserved_positions = conserved_positions,
    max_ic = max_ic,
    mean_ic = mean_ic,
    passes_min_ic = total_ic >= config$min_ic,
    quality_level = if (total_ic >= 16.0) "excellent" 
                   else if (total_ic >= 12.0) "good"
                   else if (total_ic >= config$min_ic) "acceptable"
                   else "poor"
  )
  
  return(quality_metrics)
}

# Main PWM construction function
build_robust_pwm <- function(sequences, config) {
  cat("Building robust PWM...\n")
  
  # Calculate frequency matrix
  freq_matrix <- calculate_frequencies(sequences)
  
  # Add pseudocounts
  pseudo_matrix <- add_pseudocounts(freq_matrix, config$pseudocount)
  
  # Convert to probabilities
  prob_matrix <- normalize_probabilities(pseudo_matrix)
  
  # Calculate information content
  ic_vector <- calculate_information_content(prob_matrix, config$background)
  
  # Quality assessment
  quality <- assess_pwm_quality(ic_vector, config)
  
  return(list(
    pwm = prob_matrix,
    information_content = ic_vector,
    quality_metrics = quality,
    sequences_used = length(sequences),
    pseudocount_used = config$pseudocount
  ))
}

# Export PWM to different formats
export_pwm <- function(pwm_result, output_file, format = "meme") {
  pwm <- pwm_result$pwm
  
  if (format == "meme") {
    export_meme_format(pwm, output_file, pwm_result)
  } else if (format == "jaspar") {
    export_jaspar_format(pwm, output_file)
  } else if (format == "transfac") {
    export_transfac_format(pwm, output_file)
  } else {
    # Default: RDS format
    saveRDS(pwm_result, output_file)
  }
}

# Export to MEME format
export_meme_format <- function(pwm, output_file, pwm_result) {
  meme_content <- paste0(
    "MEME version 4\n\n",
    "ALPHABET= ACGT\n\n",
    "strands: + -\n\n",
    "Background letter frequencies (from uniform background):\n",
    "A 0.25000 C 0.25000 G 0.25000 T 0.25000\n\n",
    "MOTIF CTCF_PWM CTCF\n\n",
    "letter-probability matrix: alength= 4 w= ", ncol(pwm), 
    " nsites= ", pwm_result$sequences_used, 
    " E= 0\n"
  )
  
  # Add probability matrix
  for (i in 1:ncol(pwm)) {
    meme_content <- paste0(meme_content, 
                          sprintf("%.6f\t%.6f\t%.6f\t%.6f\n", 
                                  pwm["A", i], pwm["C", i], 
                                  pwm["G", i], pwm["T", i]))
  }
  
  writeLines(meme_content, output_file)
}

# Export to JASPAR format
export_jaspar_format <- function(pwm, output_file) {
  jaspar_content <- paste0(
    ">CTCF_PWM\tCTCF\n",
    "A  [", paste(sprintf("%.6f", pwm["A", ]), collapse=" "), "]\n",
    "C  [", paste(sprintf("%.6f", pwm["C", ]), collapse=" "), "]\n",
    "G  [", paste(sprintf("%.6f", pwm["G", ]), collapse=" "), "]\n",
    "T  [", paste(sprintf("%.6f", pwm["T", ]), collapse=" "), "]\n"
  )
  
  writeLines(jaspar_content, output_file)
}

# Export to TRANSFAC format
export_transfac_format <- function(pwm, output_file) {
  transfac_content <- paste0(
    "ID CTCF_PWM\n",
    "BF T00000\n",
    "P0\tA\tC\tG\tT\n"
  )
  
  for (i in 1:ncol(pwm)) {
    transfac_content <- paste0(transfac_content,
                              sprintf("%02d\t%.3f\t%.3f\t%.3f\t%.3f\n",
                                      i, pwm["A", i], pwm["C", i], 
                                      pwm["G", i], pwm["T", i]))
  }
  
  transfac_content <- paste0(transfac_content, "//\n")
  writeLines(transfac_content, output_file)
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

# --- Main Script ---
cat("=== Robust PWM Building ===\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")  
cat("Output format:", output_format, "\n")
cat("Pseudocount:", pseudocount, "\n")
cat("Min IC threshold:", min_ic, "\n")
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

# Prepare configuration
config <- list(
  pseudocount = pseudocount,
  min_ic = min_ic,
  background = c(0.25, 0.25, 0.25, 0.25)  # uniform background
)

# Build robust PWM
pwm_result <- build_robust_pwm(filtered_sequences, config)

# Display results
cat("\nPWM Quality Metrics:\n")
cat("- PWM dimensions:", nrow(pwm_result$pwm), "x", ncol(pwm_result$pwm), "\n")
cat("- Total information content:", round(pwm_result$quality_metrics$total_ic, 3), "bits\n")
cat("- Average per position:", round(pwm_result$quality_metrics$mean_ic, 3), "bits\n")
cat("- Conserved positions (>1 bit):", pwm_result$quality_metrics$conserved_positions, "\n")
cat("- Maximum position info:", round(pwm_result$quality_metrics$max_ic, 3), "bits\n")
cat("- Quality level:", toupper(pwm_result$quality_metrics$quality_level), "\n")

# Quality check
if (pwm_result$quality_metrics$passes_min_ic) {
  cat("✅ PWM meets minimum IC threshold\n")
} else {
  cat("⚠️  PWM below minimum IC threshold - consider improving alignment\n")
}

# Export PWM in requested format
cat("\nExporting PWM...\n")
export_pwm(pwm_result, output_file, output_format)

# Save metadata
output_dir <- dirname(output_file)
base_name <- tools::file_path_sans_ext(basename(output_file))
metadata_file <- file.path(output_dir, paste0(base_name, "_metadata.json"))

metadata <- list(
  creation_time = as.character(Sys.time()),
  input_file = input_file,
  output_file = output_file,
  output_format = output_format,
  total_sequences = length(all_sequences),
  filtered_sequences = length(filtered_sequences),
  pseudocount_used = pseudocount,
  total_information = round(pwm_result$quality_metrics$total_ic, 3),
  average_information = round(pwm_result$quality_metrics$mean_ic, 3),
  conserved_positions = pwm_result$quality_metrics$conserved_positions,
  quality_level = pwm_result$quality_metrics$quality_level,
  pwm_dimensions = paste(nrow(pwm_result$pwm), "x", ncol(pwm_result$pwm)),
  passes_threshold = pwm_result$quality_metrics$passes_min_ic
)

# Create JSON manually
json_content <- paste0('{\n',
  '  "creation_time": "', metadata$creation_time, '",\n',
  '  "input_file": "', metadata$input_file, '",\n',
  '  "output_file": "', metadata$output_file, '",\n',
  '  "output_format": "', metadata$output_format, '",\n',
  '  "total_sequences": ', metadata$total_sequences, ',\n',
  '  "filtered_sequences": ', metadata$filtered_sequences, ',\n',
  '  "pseudocount_used": ', metadata$pseudocount_used, ',\n',
  '  "total_information": ', metadata$total_information, ',\n',
  '  "average_information": ', metadata$average_information, ',\n',
  '  "conserved_positions": ', metadata$conserved_positions, ',\n',
  '  "quality_level": "', metadata$quality_level, '",\n',
  '  "pwm_dimensions": "', metadata$pwm_dimensions, '",\n',
  '  "passes_threshold": ', tolower(as.character(metadata$passes_threshold)), '\n',
  '}'
)

writeLines(json_content, metadata_file)

cat("PWM saved to:", output_file, "\n")
cat("Metadata saved to:", metadata_file, "\n")
cat("=== Robust PWM Building Complete ===\n")
