# R script to prepare training and testing datasets from extracted sequences

# --- Dependencies ---


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", force = TRUE)
}




library(Biostrings)

# --- Parameters ---
input_fasta <- "data/extracted_sequences.fasta"
output_train_fasta <- "data/training_sequences.fasta"
output_test_fasta <- "data/test_sequences.fasta"

# Desired sequence length for PWM building and evaluation
# Set to NULL to skip length filtering
target_length <- 82 # Example: Match the original PWM length

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

# --- Main Script ---

cat("Starting dataset preparation...\n")

# 1. Read Input Sequences
cat("Reading sequences from:", input_fasta, "\n")
if (!file.exists(input_fasta)) {
  stop("Input FASTA file not found: ", input_fasta, ". Run download_data.sh first.")
}
all_sequences <- readDNAStringSet(input_fasta)
if (!inherits(all_sequences, "DNAStringSet")) {
  all_sequences <- DNAStringSet(all_sequences)
}
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

# 3. Shuffle and Split Data
cat("Shuffling and splitting data (", train_proportion * 100, "% training)...\n", sep="")
shuffled_indices <- sample(filtered_count)
num_train <- floor(train_proportion * filtered_count)
num_test <- filtered_count - num_train

train_indices <- shuffled_indices[1:num_train]
test_indices <- shuffled_indices[(num_train + 1):filtered_count]

train_sequences <- all_sequences[train_indices]
test_sequences <- all_sequences[test_indices]

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


