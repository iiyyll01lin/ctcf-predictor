# R script to prepare training and testing datasets from extracted sequences

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
library(Biostrings)

# --- Parameters ---
# input_fasta <- "../data/extracted_sequences.fasta"
input_fasta <- "data/preprocessed_sequences_optimized.fasta" 
output_train_fasta <- "data/training_sequences.fasta"
output_test_fasta <- "data/test_sequences.fasta"

# Desired sequence length for PWM building and evaluation
# Set to NULL to skip length filtering
target_length <- NULL
# target_length <- 11 # Example: Match the original PWM length

# Proportion of data to use for the training set (e.g., 0.8 = 80%)
train_proportion <- 0.8

# Negative example generation parameters
generate_negatives <- TRUE # Set to FALSE to disable automatic negative generation
negative_generation_method <- "shuffle" # Options: "shuffle", "dinucleotide_shuffle", "random"
negative_to_positive_ratio <- 1.0 # Ratio of negative to positive examples in test set

# Random seed for reproducibility
set.seed(123)

# --- Negative Example Generation Functions ---

# Simple sequence shuffling (maintains nucleotide composition)
shuffle_sequence <- function(seq) {
  chars <- unlist(strsplit(as.character(seq), ""))
  paste0(sample(chars), collapse="")
}

# Generate a completely random DNA sequence of specified length
generate_random_sequence <- function(length) {
  nucleotides <- c("A", "C", "G", "T")
  paste0(sample(nucleotides, size=length, replace=TRUE), collapse="")
}

# Dinucleotide shuffling (better preserves sequence complexity)
# This is a simplified implementation of the Altschul-Erikson algorithm
dinucleotide_shuffle <- function(seq) {
  seq_str <- as.character(seq)
  seq_len <- nchar(seq_str)
  
  if (seq_len <= 2) return(seq_str) # Cannot shuffle sequences of length 1-2
  
  # Extract dinucleotides
  dinucs <- character(seq_len - 1)
  for (i in 1:(seq_len - 1)) {
    dinucs[i] <- substr(seq_str, i, i+1)
  }
  
  # Keep first nucleotide fixed (starting point)
  result <- substr(seq_str, 1, 1)
  available_dinucs <- dinucs
  current_nuc <- substr(seq_str, 1, 1)
  
  # Build the shuffled sequence
  while (length(available_dinucs) > 0) {
    # Find dinucleotides that start with current_nuc
    candidates <- grep(paste0("^", current_nuc), available_dinucs, value=TRUE)
    
    if (length(candidates) == 0) {
      # If no candidates, just pick a random remaining dinucleotide
      chosen <- sample(available_dinucs, 1)
    } else {
      # Select a random matching dinucleotide
      chosen <- sample(candidates, 1)
    }
    
    # Add the second nucleotide of chosen dinucleotide
    next_nuc <- substr(chosen, 2, 2)
    result <- paste0(result, next_nuc)
    current_nuc <- next_nuc
    
    # Remove the used dinucleotide
    available_dinucs <- available_dinucs[available_dinucs != chosen]
  }
  
  return(result)
}

# Generate negative examples based on positive examples
generate_negative_examples <- function(positive_seqs, method="shuffle", ratio=1.0) {
  n_positives <- length(positive_seqs)
  n_negatives <- ceiling(n_positives * ratio)
  
  cat("Generating", n_negatives, "negative examples using method:", method, "\n")
  
  negative_seqs <- DNAStringSet()
  
  for (i in 1:n_negatives) {
    # Choose a positive sequence to base the negative on
    pos_idx <- ((i - 1) %% n_positives) + 1
    pos_seq <- positive_seqs[[pos_idx]]
    seq_len <- nchar(as.character(pos_seq))
    
    # Apply the chosen method to generate a negative example
    if (method == "shuffle") {
      new_seq <- shuffle_sequence(pos_seq)
    } else if (method == "dinucleotide_shuffle") {
      new_seq <- dinucleotide_shuffle(pos_seq)
    } else if (method == "random") {
      new_seq <- generate_random_sequence(seq_len)
    } else {
      stop("Unknown negative generation method:", method)
    }
    
    # Create a named negative sequence
    neg_seq <- DNAStringSet(new_seq)
    names(neg_seq) <- paste0("negative_", i, " | class=0")
    
    # Add to our collection
    negative_seqs <- c(negative_seqs, neg_seq)
  }
  
  return(negative_seqs)
}

# --- Genomic Split Helper Functions ---

# Validate and show examples of sequence names for debugging
validate_sequence_names <- function(seq_names, show_examples = 5) {
  cat("Validating sequence name formats for chromosome extraction...\n")
  
  # Show first few examples
  n_examples <- min(show_examples, length(seq_names))
  cat("Example sequence names (first", n_examples, "):\n")
  for (i in 1:n_examples) {
    cat("  ", i, ":", seq_names[i], "\n")
  }
  
  # Test chromosome extraction on examples
  chr_patterns <- c(
    "chr[0-9XY]+:[0-9]+-[0-9]+",    # chr1:12345-67890
    "chr[0-9XY]+_[0-9]+_[0-9]+",    # chr1_12345_67890
    "chr[0-9XY]+\\.[0-9]+\\.[0-9]+", # chr1.12345.67890
    "chr[0-9XY]+"                   # just chr1 (minimal)
  )
  
  cat("Testing common chromosome patterns:\n")
  for (pattern in chr_patterns) {
    matches <- sum(grepl(pattern, seq_names[1:n_examples]))
    cat("  Pattern '", pattern, "': ", matches, "/", n_examples, " matches\n")
  }
  
  return(invisible(NULL))
}

# Extract chromosome information from sequence names
extract_chromosome <- function(seq_name) {
  # Try multiple common formats
  patterns <- list(
    chr_colon = "chr[0-9XY]+(?=:)",           # chr1:start-end
    chr_underscore = "chr[0-9XY]+(?=_)",      # chr1_start_end  
    chr_dot = "chr[0-9XY]+(?=\\.)",           # chr1.start.end
    chr_simple = "chr[0-9XY]+"                # chr1 (anywhere)
  )
  
  for (pattern in patterns) {
    chr_match <- regexpr(pattern, seq_name, perl = TRUE)
    if (chr_match > 0) {
      return(substr(seq_name, chr_match, chr_match + attr(chr_match, "match.length") - 1))
    }
  }
  
  return("unknown")
}

# --- Main Script ---

cat("Starting dataset preparation...\n")

# 1. Read Input Sequences
cat("Reading sequences from:", input_fasta, "\n")
if (!file.exists(input_fasta)) {
  stop("Input FASTA file not found: ", input_fasta, ". Run download_data.sh first.")
}
all_sequences <- readDNAStringSet(input_fasta)
cat("Read", length(all_sequences), "sequences.\n")

# 2. Filter by Length (Optional)
if (!is.null(target_length)) {
  cat("Filtering sequences to length:", target_length, "\n")
  original_count <- length(all_sequences)
  all_sequences <- all_sequences[width(all_sequences) == target_length]
  filtered_count <- length(all_sequences)
  cat("Kept", filtered_count, "out of", original_count, "sequences matching length", target_length, "\n")
  if (filtered_count == 0) {
    stop("No sequences remaining after length filtering. Check target_length or input data.")
  }
} else {
  cat("Skipping length filtering.\n")
  filtered_count <- length(all_sequences)
}

# 3. Chromosome-based Split (prevents genomic data leakage)
cat("Performing chromosome-based split (", train_proportion * 100, "% training)...\n", sep="")

# Validate sequence names first
validate_sequence_names(names(all_sequences))

# Get chromosome for each sequence
chromosomes <- sapply(names(all_sequences), extract_chromosome)
unique_chrs <- unique(chromosomes)
cat("Found chromosomes:", paste(sort(unique_chrs), collapse=", "), "\n")

# Handle case where chromosome extraction fails
if (length(unique_chrs) == 1 && unique_chrs[1] == "unknown") {
  cat("WARNING: Could not extract chromosome information from sequence names.\n")
  cat("Falling back to random split. Consider checking sequence name format.\n")
  cat("Expected formats: 'chr1:start-end', 'chr1_start_end', etc.\n")
  
  # Fallback to original random split
  shuffled_indices <- sample(filtered_count)
  num_train <- floor(train_proportion * filtered_count)
  train_indices <- shuffled_indices[1:num_train]
  test_indices <- shuffled_indices[(num_train + 1):filtered_count]
  
  train_sequences <- all_sequences[train_indices]
  test_sequences <- all_sequences[test_indices]
  
} else {
  # Split chromosomes between train/test (not individual sequences)
  set.seed(123)  # For reproducible chromosome assignment
  shuffled_chrs <- sample(unique_chrs)
  n_train_chrs <- max(1, ceiling(length(unique_chrs) * train_proportion))
  
  train_chrs <- shuffled_chrs[1:n_train_chrs]
  test_chrs <- shuffled_chrs[(n_train_chrs + 1):length(unique_chrs)]
  
  cat("Training chromosomes (", length(train_chrs), "):", paste(sort(train_chrs), collapse=", "), "\n")
  cat("Testing chromosomes (", length(test_chrs), "):", paste(sort(test_chrs), collapse=", "), "\n")
  
  # Split sequences based on chromosome assignment
  train_mask <- chromosomes %in% train_chrs
  test_mask <- chromosomes %in% test_chrs
  
  train_sequences <- all_sequences[train_mask]
  test_sequences <- all_sequences[test_mask]
  
  # Verify no overlap
  if (length(intersect(train_chrs, test_chrs)) > 0) {
    warning("Chromosome overlap detected between train/test sets!")
  }
  
  # Report split statistics
  cat("Split quality check:\n")
  cat("- Unique training chromosomes:", length(unique(chromosomes[train_mask])), "\n")
  cat("- Unique testing chromosomes:", length(unique(chromosomes[test_mask])), "\n")
}

cat("Training set size:", length(train_sequences), "\n")
cat("Test set size:", length(test_sequences), "\n")

# 4. Write Training Set
cat("Writing training set to:", output_train_fasta, "\n")
writeXStringSet(train_sequences, filepath = output_train_fasta)

# 5. Prepare Test Set (Adding Positive Labels)
cat("Preparing test set with positive/negative labels...\n")
# Add "| class=1" to headers, indicating these are positive examples
names(test_sequences) <- paste0(names(test_sequences), " | class=1")

# 6. Generate Negative Examples (if enabled)
if (generate_negatives) {
  cat("Automatic negative example generation is enabled.\n")
  negative_sequences <- generate_negative_examples(
    test_sequences, 
    method = negative_generation_method,
    ratio = negative_to_positive_ratio
  )
  
  # Add negative examples to test set
  test_sequences <- c(test_sequences, negative_sequences)
  cat("Final test set size (positives + negatives):", length(test_sequences), "\n")
} else {
  cat("Automatic negative example generation is disabled.\n")
  cat("IMPORTANT: The test set currently only contains positive examples (class=1).\n")
  cat("For meaningful evaluation, generate and include appropriate negative examples (class=0).\n")
}

# 7. Write Test Set
cat("Writing test set to:", output_test_fasta, "\n")
writeXStringSet(test_sequences, filepath = output_test_fasta)

cat("\nDataset preparation complete.\n")
