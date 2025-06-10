# High-Quality Subset PWM Builder
# This script builds PWMs from carefully selected high-quality sequence subsets

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/training_sequences.fasta"
output_prefix <- if (length(args) >= 2) args[2] else "results/subset_pwm"
subset_size <- if (length(args) >= 3) as.numeric(args[3]) else 5000
quality_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 0.01

cat("High-Quality Subset PWM Builder\n")
cat("===============================\n")
cat("Input file:", input_file, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Subset size:", subset_size, "\n")
cat("Quality threshold (max N ratio):", quality_threshold, "\n\n")

# --- Functions ---

# Calculate sequence quality metrics
calculate_sequence_quality <- function(seq_string) {
  seq_length <- nchar(seq_string)
  n_count <- length(gregexpr("N", seq_string, ignore.case = TRUE)[[1]])
  if (n_count == 1 && gregexpr("N", seq_string, ignore.case = TRUE)[[1]][1] == -1) {
    n_count <- 0
  }
  
  # Calculate nucleotide frequencies
  a_count <- length(gregexpr("A", seq_string, ignore.case = TRUE)[[1]])
  c_count <- length(gregexpr("C", seq_string, ignore.case = TRUE)[[1]])
  g_count <- length(gregexpr("G", seq_string, ignore.case = TRUE)[[1]])
  t_count <- length(gregexpr("T", seq_string, ignore.case = TRUE)[[1]])
  
  # Handle cases where pattern is not found
  if (a_count == 1 && gregexpr("A", seq_string, ignore.case = TRUE)[[1]][1] == -1) a_count <- 0
  if (c_count == 1 && gregexpr("C", seq_string, ignore.case = TRUE)[[1]][1] == -1) c_count <- 0
  if (g_count == 1 && gregexpr("G", seq_string, ignore.case = TRUE)[[1]][1] == -1) g_count <- 0
  if (t_count == 1 && gregexpr("T", seq_string, ignore.case = TRUE)[[1]][1] == -1) t_count <- 0
  
  gc_content <- (g_count + c_count) / seq_length
  n_ratio <- n_count / seq_length
  
  # Calculate sequence complexity (entropy)
  freqs <- c(a_count, c_count, g_count, t_count) / seq_length
  freqs[freqs == 0] <- 1e-10
  complexity <- -sum(freqs * log2(freqs))
  
  return(list(
    length = seq_length,
    n_ratio = n_ratio,
    gc_content = gc_content,
    complexity = complexity
  ))
}

# Select high-quality sequences
select_high_quality_sequences <- function(sequences, max_n_ratio = 0.01, 
                                        target_length = NULL, length_tolerance = 0.1) {
  cat("Analyzing sequence quality...\n")
  
  seq_strings <- as.character(sequences)
  quality_metrics <- lapply(seq_strings, calculate_sequence_quality)
  
  # Extract metrics
  lengths <- sapply(quality_metrics, function(x) x$length)
  n_ratios <- sapply(quality_metrics, function(x) x$n_ratio)
  gc_contents <- sapply(quality_metrics, function(x) x$gc_content)
  complexities <- sapply(quality_metrics, function(x) x$complexity)
  
  cat("Original sequences:", length(sequences), "\n")
  cat("Length range:", min(lengths), "-", max(lengths), "\n")
  cat("Mean GC content:", round(mean(gc_contents), 3), "\n")
  cat("Mean complexity:", round(mean(complexities), 3), "\n")
  
  # Filter by N content
  good_n_filter <- n_ratios <= max_n_ratio
  cat("Sequences with N ratio <=", max_n_ratio, ":", sum(good_n_filter), "\n")
  
  # Filter by length consistency
  if (is.null(target_length)) {
    target_length <- median(lengths)
  }
  length_filter <- abs(lengths - target_length) <= (target_length * length_tolerance)
  cat("Sequences within", length_tolerance*100, "% of target length", target_length, ":", sum(length_filter), "\n")
  
  # Filter by complexity (remove low complexity sequences)
  complexity_threshold <- quantile(complexities, 0.25)  # Bottom quartile
  complexity_filter <- complexities > complexity_threshold
  cat("Sequences above complexity threshold:", sum(complexity_filter), "\n")
  
  # Combine filters
  high_quality_filter <- good_n_filter & length_filter & complexity_filter
  cat("High-quality sequences (all filters):", sum(high_quality_filter), "\n")
  
  if (sum(high_quality_filter) == 0) {
    warning("No sequences passed all quality filters. Relaxing criteria...")
    high_quality_filter <- good_n_filter & length_filter
    cat("High-quality sequences (relaxed):", sum(high_quality_filter), "\n")
  }
  
  return(list(
    sequences = sequences[high_quality_filter],
    indices = which(high_quality_filter),
    metrics = list(
      original_count = length(sequences),
      filtered_count = sum(high_quality_filter),
      target_length = target_length,
      mean_gc = mean(gc_contents[high_quality_filter]),
      mean_complexity = mean(complexities[high_quality_filter])
    )
  ))
}

# Calculate information content for PWM
calculate_info_content <- function(pwm) {
  apply(pwm, 2, function(x) {
    x[x == 0] <- 1e-10
    2 + sum(x * log2(x))
  })
}

# Build PWM from high-quality subset
build_subset_pwm <- function(sequences, subset_size = 5000, pseudocount = 0.1) {
  cat("\nBuilding PWM from high-quality subset...\n")
  
  # Select subset if needed
  if (length(sequences) > subset_size) {
    # Random sampling from high-quality sequences
    sample_indices <- sample(length(sequences), subset_size)
    sequences <- sequences[sample_indices]
    cat("Random subset selected:", length(sequences), "sequences\n")
  } else {
    cat("Using all", length(sequences), "high-quality sequences\n")
  }
  
  # Get consensus matrix
  consensus_matrix <- consensusMatrix(sequences, as.prob = FALSE)
  
  # Keep only standard nucleotides
  nucleotides <- c("A", "C", "G", "T")
  count_matrix <- consensus_matrix[nucleotides, , drop = FALSE]
  
  cat("PWM dimensions:", nrow(count_matrix), "x", ncol(count_matrix), "\n")
  
  # Add pseudocounts and normalize
  count_matrix <- count_matrix + pseudocount
  pwm <- apply(count_matrix, 2, function(x) x / sum(x))
  
  # Calculate information content
  info_content <- calculate_info_content(pwm)
  total_info <- sum(info_content)
  
  cat("Total information content:", round(total_info, 3), "bits\n")
  cat("Average per position:", round(total_info / ncol(pwm), 3), "bits\n")
  
  # Find conserved positions
  conserved_positions <- which(info_content > 1.0)
  cat("Conserved positions (>1 bit):", length(conserved_positions), "\n")
  
  return(list(
    pwm = pwm,
    info_content = info_content,
    total_info = total_info,
    conserved_positions = conserved_positions,
    num_sequences = length(sequences),
    pseudocount = pseudocount,
    method = "high_quality_subset"
  ))
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

# Select high-quality sequences
quality_result <- select_high_quality_sequences(sequences, 
                                               max_n_ratio = quality_threshold)

if (length(quality_result$sequences) == 0) {
  stop("No high-quality sequences found")
}

# Build multiple PWMs with different subset sizes
subset_sizes <- c(1000, 2000, 5000, min(10000, length(quality_result$sequences)))
subset_sizes <- subset_sizes[subset_sizes <= length(quality_result$sequences)]

pwm_results <- list()

for (size in subset_sizes) {
  cat("\n", paste(rep("=", 50), collapse=""), "\n")
  cat("Building PWM with subset size:", size, "\n")
  
  pwm_result <- build_subset_pwm(quality_result$sequences, subset_size = size)
  pwm_result$subset_size <- size
  pwm_result$quality_metrics <- quality_result$metrics
  pwm_result$creation_time <- Sys.time()
  
  # Save individual PWM
  output_file <- paste0(output_prefix, "_size", size, ".rds")
  saveRDS(pwm_result, output_file)
  cat("PWM saved to:", output_file, "\n")
  
  pwm_results[[paste0("size_", size)]] <- pwm_result
}

# Save combined results
combined_file <- paste0(output_prefix, "_all_sizes.rds")
saveRDS(pwm_results, combined_file)
cat("\nAll PWM results saved to:", combined_file, "\n")

# Print comparison summary
cat("\n=== PWM Comparison Summary ===\n")
cat("Subset Size | Total Info | Avg Info | Conserved Pos\n")
cat(paste(rep("-", 50), collapse=""), "\n")
for (size in subset_sizes) {
  result <- pwm_results[[paste0("size_", size)]]
  cat(sprintf("%10d | %9.3f | %8.3f | %12d\n", 
              size, result$total_info, 
              result$total_info / ncol(result$pwm),
              length(result$conserved_positions)))
}

cat("\nScript completed successfully!\n")