# Test Script for Chromosome-based Split Implementation
# This script validates the new chromosome-based splitting functionality

library(Biostrings)

cat("=== Chromosome-based Split Test ===\n")
cat("Testing time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Test Data Preparation ---

# Check if input data exists
input_file <- "data/preprocessed_sequences_optimized.fasta"
if (!file.exists(input_file)) {
  cat("❌ Input file not found:", input_file, "\n")
  cat("Please run the preprocessing pipeline first.\n")
  quit(status = 1)
}

# --- Test 1: Load and Inspect Data ---
cat("1. LOADING AND INSPECTING DATA\n")
cat(paste(rep("=", 40), collapse=""), "\n")

sequences <- readDNAStringSet(input_file)
cat("Total sequences:", length(sequences), "\n")

# Show first few sequence names
cat("\nFirst 10 sequence names:\n")
for (i in 1:min(10, length(sequences))) {
  cat("  ", i, ":", names(sequences)[i], "\n")
}

# --- Test 2: Chromosome Extraction ---
cat("\n2. CHROMOSOME EXTRACTION TEST\n")
cat(paste(rep("=", 40), collapse=""), "\n")

# Source the new extraction function
source("scripts/prepare_datasets.R", local = TRUE)

# Test chromosome extraction on sample names
test_names <- names(sequences)[1:min(20, length(sequences))]
chromosomes <- sapply(test_names, extract_chromosome)

cat("Chromosome extraction results:\n")
for (i in 1:length(test_names)) {
  cat("  ", chromosomes[i], " <- ", test_names[i], "\n")
}

# Count chromosomes
chr_counts <- table(chromosomes)
cat("\nChromosome distribution:\n")
print(chr_counts)

# --- Test 3: Split Quality Analysis ---
cat("\n3. SPLIT QUALITY ANALYSIS\n")
cat(paste(rep("=", 40), collapse=""), "\n")

# Check if we have enough chromosomes for a good split
unique_chrs <- unique(chromosomes)
cat("Unique chromosomes found:", length(unique_chrs), "\n")
cat("Chromosomes:", paste(sort(unique_chrs), collapse=", "), "\n")

if (length(unique_chrs) < 3) {
  cat("⚠️  WARNING: Only", length(unique_chrs), "chromosomes found.\n")
  cat("   This may result in imbalanced train/test splits.\n")
} else {
  cat("✅ Good chromosome diversity for splitting.\n")
}

# --- Test 4: Run the Actual Split ---
cat("\n4. RUNNING CHROMOSOME-BASED SPLIT\n")
cat(paste(rep("=", 40), collapse=""), "\n")

# Run prepare_datasets.R and capture output
cat("Running prepare_datasets.R...\n")
system("Rscript scripts/prepare_datasets.R", intern = FALSE)

# --- Test 5: Validate Output Files ---
cat("\n5. VALIDATING OUTPUT FILES\n")
cat(paste(rep("=", 40), collapse=""), "\n")

train_file <- "data/training_sequences.fasta"
test_file <- "data/test_sequences.fasta"

if (file.exists(train_file)) {
  train_seqs <- readDNAStringSet(train_file)
  cat("✅ Training file created:", length(train_seqs), "sequences\n")
  
  # Check chromosomes in training set
  train_chrs <- sapply(names(train_seqs), extract_chromosome)
  train_unique <- unique(train_chrs)
  cat("   Training chromosomes:", paste(sort(train_unique), collapse=", "), "\n")
} else {
  cat("❌ Training file not found\n")
}

if (file.exists(test_file)) {
  test_seqs <- readDNAStringSet(test_file)
  cat("✅ Test file created:", length(test_seqs), "sequences\n")
  
  # Check chromosomes in test set (need to remove class labels first)
  test_names_clean <- gsub(" \\| class=[01]", "", names(test_seqs))
  test_chrs <- sapply(test_names_clean, extract_chromosome)
  test_unique <- unique(test_chrs)
  cat("   Test chromosomes:", paste(sort(test_unique), collapse=", "), "\n")
  
  # Check for overlap
  if (exists("train_unique") && length(intersect(train_unique, test_unique)) > 0) {
    cat("❌ CHROMOSOME OVERLAP DETECTED between train/test!\n")
    cat("   Overlapping chromosomes:", paste(intersect(train_unique, test_unique), collapse=", "), "\n")
  } else if (exists("train_unique")) {
    cat("✅ No chromosome overlap between train/test sets\n")
  }
  
  # Check positive/negative balance in test set
  pos_count <- sum(grepl("class=1", names(test_seqs)))
  neg_count <- sum(grepl("class=0", names(test_seqs)))
  cat("   Test set composition: ", pos_count, " positives, ", neg_count, " negatives\n")
  
} else {
  cat("❌ Test file not found\n")
}

# --- Test 6: Data Quality Checks ---
cat("\n6. DATA QUALITY CHECKS\n")
cat(paste(rep("=", 40), collapse=""), "\n")

if (exists("train_seqs") && exists("test_seqs")) {
  # Check sequence length distribution
  train_lengths <- width(train_seqs)
  test_lengths <- width(test_seqs[grepl("class=1", names(test_seqs))])  # Only positive test examples
  
  cat("Training sequence lengths - Min:", min(train_lengths), "Max:", max(train_lengths), "Mean:", round(mean(train_lengths), 1), "\n")
  cat("Test sequence lengths - Min:", min(test_lengths), "Max:", max(test_lengths), "Mean:", round(mean(test_lengths), 1), "\n")
  
  # Check for sufficient data
  if (length(train_seqs) < 1000) {
    cat("⚠️  WARNING: Training set may be too small (", length(train_seqs), " sequences)\n")
  } else {
    cat("✅ Training set size adequate\n")
  }
  
  if (pos_count < 100) {
    cat("⚠️  WARNING: Test set may be too small (", pos_count, " positive examples)\n")
  } else {
    cat("✅ Test set size adequate\n")
  }
}

cat("\n=== Test Complete ===\n")
cat("If all checks pass (✅), the chromosome-based split is working correctly.\n")
cat("If there are warnings (⚠️) or errors (❌), review the implementation.\n")
