# Data Preparation Script - Initial data preprocessing
# This script implements the data preparation layer as defined in the system architecture

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Please install Biostrings: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/extracted_sequences.fasta"
output_file <- if (length(args) >= 2) args[2] else "data/prepared_sequences.fasta"
min_length <- if (length(args) >= 3) as.numeric(args[3]) else 50
max_length <- if (length(args) >= 4) as.numeric(args[4]) else 500

cat("Data Preparation Pipeline\n")
cat("========================\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Length range:", min_length, "-", max_length, "bp\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Functions ---

# Validate sequence quality
validate_sequence_quality <- function(sequences) {
  cat("Validating sequence quality...\n")
  
  # Check for empty sequences
  empty_seqs <- width(sequences) == 0
  if (any(empty_seqs)) {
    cat("Warning: Found", sum(empty_seqs), "empty sequences\n")
  }
  
  # Check for N-base content
  n_content <- letterFrequency(sequences, "N", as.prob = TRUE)
  high_n <- n_content > 0.1  # More than 10% N bases
  if (any(high_n)) {
    cat("Warning: Found", sum(high_n), "sequences with >10% N bases\n")
  }
  
  return(list(
    empty_sequences = empty_seqs,
    high_n_content = high_n,
    n_content_stats = summary(as.vector(n_content))
  ))
}

# Filter sequences by length
filter_by_length <- function(sequences, min_len = 50, max_len = 500) {
  cat("Filtering sequences by length (", min_len, "-", max_len, "bp)...\n")
  
  lengths <- width(sequences)
  valid_length <- lengths >= min_len & lengths <= max_len
  
  cat("Sequences passing length filter:", sum(valid_length), "/", length(sequences), "\n")
  
  return(sequences[valid_length])
}

# Remove duplicate sequences
remove_duplicates <- function(sequences) {
  cat("Removing duplicate sequences...\n")
  
  # Convert to character for duplicate detection
  seq_chars <- as.character(sequences)
  duplicated_seqs <- duplicated(seq_chars)
  
  cat("Duplicate sequences removed:", sum(duplicated_seqs), "\n")
  
  return(sequences[!duplicated_seqs])
}

# Basic sequence statistics
calculate_basic_stats <- function(sequences) {
  lengths <- width(sequences)
  gc_content <- letterFrequency(sequences, "GC", as.prob = TRUE)
  n_content <- letterFrequency(sequences, "N", as.prob = TRUE)
  
  return(list(
    count = length(sequences),
    length_stats = list(
      mean = mean(lengths),
      median = median(lengths),
      min = min(lengths),
      max = max(lengths),
      sd = sd(lengths)
    ),
    gc_content_stats = list(
      mean = mean(gc_content),
      median = median(gc_content),
      min = min(gc_content),
      max = max(gc_content)
    ),
    n_content_stats = list(
      mean = mean(n_content),
      max = max(n_content),
      sequences_with_n = sum(n_content > 0)
    )
  ))
}

# --- Main Processing ---

# Check input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Load sequences
cat("Loading sequences from:", input_file, "\n")
raw_sequences <- readDNAStringSet(input_file)
cat("Raw sequences loaded:", length(raw_sequences), "\n\n")

if (length(raw_sequences) == 0) {
  stop("No sequences found in input file")
}

# Calculate initial statistics
cat("Initial sequence statistics:\n")
initial_stats <- calculate_basic_stats(raw_sequences)
cat("  Count:", initial_stats$count, "\n")
cat("  Length: mean =", round(initial_stats$length_stats$mean, 1), 
    "bp, range =", initial_stats$length_stats$min, "-", initial_stats$length_stats$max, "bp\n")
cat("  GC content: mean =", round(initial_stats$gc_content_stats$mean * 100, 1), "%\n")
cat("  N content: mean =", round(initial_stats$n_content_stats$mean * 100, 2), 
    "%, sequences with N =", initial_stats$n_content_stats$sequences_with_n, "\n\n")

# Validate sequence quality
quality_check <- validate_sequence_quality(raw_sequences)

# Step 1: Length filtering
filtered_sequences <- filter_by_length(raw_sequences, min_length, max_length)

# Step 2: Remove duplicates
unique_sequences <- remove_duplicates(filtered_sequences)

# Step 3: Final quality check
if (length(unique_sequences) == 0) {
  stop("No sequences remain after filtering")
}

# Calculate final statistics
cat("\nFinal sequence statistics:\n")
final_stats <- calculate_basic_stats(unique_sequences)
cat("  Count:", final_stats$count, "\n")
cat("  Length: mean =", round(final_stats$length_stats$mean, 1), 
    "bp, range =", final_stats$length_stats$min, "-", final_stats$length_stats$max, "bp\n")
cat("  GC content: mean =", round(final_stats$gc_content_stats$mean * 100, 1), "%\n")
cat("  N content: mean =", round(final_stats$n_content_stats$mean * 100, 2), "%\n")

# Create output directory if needed
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save prepared sequences
cat("\nSaving prepared sequences to:", output_file, "\n")
writeXStringSet(unique_sequences, output_file)

# Generate summary report
summary_file <- gsub("\\.fasta$", "_summary.txt", output_file)
cat("\nGenerating summary report:", summary_file, "\n")

sink(summary_file)
cat("Data Preparation Summary Report\n")
cat("===============================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n\n")

cat("Processing Parameters:\n")
cat("  Min length:", min_length, "bp\n")
cat("  Max length:", max_length, "bp\n\n")

cat("Initial Dataset:\n")
cat("  Sequences:", initial_stats$count, "\n")
cat("  Mean length:", round(initial_stats$length_stats$mean, 1), "bp\n")
cat("  Length range:", initial_stats$length_stats$min, "-", initial_stats$length_stats$max, "bp\n")
cat("  Mean GC content:", round(initial_stats$gc_content_stats$mean * 100, 1), "%\n\n")

cat("Final Dataset:\n")
cat("  Sequences:", final_stats$count, "\n")
cat("  Mean length:", round(final_stats$length_stats$mean, 1), "bp\n")
cat("  Length range:", final_stats$length_stats$min, "-", final_stats$length_stats$max, "bp\n")
cat("  Mean GC content:", round(final_stats$gc_content_stats$mean * 100, 1), "%\n\n")

cat("Filtering Results:\n")
cat("  Sequences retained:", final_stats$count, "/", initial_stats$count, 
    "(", round(final_stats$count / initial_stats$count * 100, 1), "%)\n")
cat("  Sequences removed by length filter:", initial_stats$count - length(filtered_sequences), "\n")
cat("  Duplicate sequences removed:", length(filtered_sequences) - final_stats$count, "\n")

sink()

cat("\nData preparation completed successfully!\n")
cat("Prepared", final_stats$count, "sequences saved to", output_file, "\n")
