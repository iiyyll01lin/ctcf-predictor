# Sequence Characteristics Analysis Script for CTCF Preprocessing Optimization
# This script analyzes DNA sequences to help optimize preprocessing parameters

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/extracted_sequences.fasta"
output_report <- if (length(args) >= 2) args[2] else "results/sequence_analysis_report.txt"

# --- Utility Functions ---

# Calculate sequence entropy (measure of complexity)
calculate_entropy <- function(sequence, window_size = 10) {
  seq_string <- as.character(sequence)
  seq_length <- nchar(seq_string)
  
  if (seq_length < window_size) {
    window_size <- seq_length
  }
  
  entropies <- numeric(seq_length - window_size + 1)
  
  for (i in 1:(seq_length - window_size + 1)) {
    window <- substr(seq_string, i, i + window_size - 1)
    base_counts <- table(strsplit(window, "")[[1]])
    base_probs <- base_counts / window_size
    entropy <- -sum(base_probs * log2(base_probs))
    entropies[i] <- entropy
  }
  
  return(min(entropies))
}

# Detect simple repeats in sequence
detect_repeats <- function(sequence, min_length = 5) {
  seq_string <- as.character(sequence)
  
  # Check for homopolymer repeats (AAAAA, TTTTT, etc.)
  homopolymer_count <- 0
  for (base in c("A", "C", "G", "T")) {
    pattern <- paste0(rep(base, min_length), collapse = "")
    if (grepl(pattern, seq_string)) {
      homopolymer_count <- homopolymer_count + 1
    }
  }
  
  # Check for dinucleotide repeats (ATATATATAT, etc.)
  dinuc_count <- 0
  dinucs <- c("AT", "TA", "GC", "CG", "AC", "CA", "GT", "TG", "AG", "GA", "CT", "TC")
  for (dinuc in dinucs) {
    repeat_count <- floor(min_length / 2)
    if (repeat_count >= 2) {
      pattern <- paste0(rep(dinuc, repeat_count), collapse = "")
      if (grepl(pattern, seq_string)) {
        dinuc_count <- dinuc_count + 1
      }
    }
  }
  
  return(c(homopolymer_count, dinuc_count))
}

# Simulate filtering effects with different parameter combinations
simulate_filtering <- function(sequences, scenarios) {
  lengths <- width(sequences)
  
  # Calculate N base percentages
  n_counts <- sapply(sequences, function(seq) {
    seq_str <- as.character(seq)
    lengths(gregexpr("N", seq_str))
  })
  n_percentages <- (n_counts / lengths) * 100
  
  # Calculate entropies for a sample (to save time)
  sample_size <- min(1000, length(sequences))
  sample_indices <- sample(length(sequences), sample_size)
  entropies <- sapply(sequences[sample_indices], calculate_entropy)
  
  # Extrapolate entropy filtering to full dataset (rough estimate)
  low_entropy_rate <- sum(entropies < 1.5) / length(entropies)
  
  results <- data.frame()
  
  for (i in 1:length(scenarios)) {
    scenario <- scenarios[[i]]
    name <- names(scenarios)[i]
    
    # Apply length filter
    length_pass <- lengths >= scenario$min_len & lengths <= scenario$max_len
    
    # Apply N base filter
    n_pass <- n_percentages <= scenario$max_n
    
    # Apply entropy filter (estimated)
    entropy_pass <- if (scenario$entropy_filter) {
      rep(c(TRUE, FALSE), times = c(round(length(sequences) * (1 - low_entropy_rate)), 
                                   round(length(sequences) * low_entropy_rate)))
    } else {
      rep(TRUE, length(sequences))
    }
    entropy_pass <- entropy_pass[1:length(sequences)]  # Ensure correct length
    
    # Combined filtering
    combined_pass <- length_pass & n_pass & entropy_pass
    
    results <- rbind(results, data.frame(
      Scenario = name,
      Length_Filter = sum(length_pass),
      N_Filter = sum(n_pass),
      Entropy_Filter = sum(entropy_pass),
      Combined = sum(combined_pass),
      Percentage = round(sum(combined_pass) / length(sequences) * 100, 2),
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# --- Main Analysis Script ---

cat("=== CTCF Sequence Characteristics Analysis Report ===\n")
cat("Input file:", input_file, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Redirect output to file if specified
if (output_report != "") {
  # Create results directory if it doesn't exist
  dir.create("results", showWarnings = FALSE)
  sink(output_report)
  cat("=== CTCF Sequence Characteristics Analysis Report ===\n")
  cat("Input file:", input_file, "\n")
  cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
}

# Validate input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Read sequences
cat("Reading sequences...\n")
seqs <- readDNAStringSet(input_file)
cat("Total sequences loaded:", length(seqs), "\n\n")

if (length(seqs) == 0) {
  stop("No sequences found in input file")
}

# 1. Length Analysis
cat("1. SEQUENCE LENGTH ANALYSIS\n")
cat("=" , rep("=", 50), "\n", sep = "")
lengths <- width(seqs)
cat("Average length:", round(mean(lengths), 2), "bp\n")
cat("Median length:", median(lengths), "bp\n")
cat("Shortest sequence:", min(lengths), "bp\n")
cat("Longest sequence:", max(lengths), "bp\n")
cat("Standard deviation:", round(sd(lengths), 2), "bp\n\n")

# Length distribution
cat("Length Distribution:\n")
breaks <- c(0, 50, 100, 150, 200, 300, 500, 1000, max(lengths))
length_groups <- cut(lengths, breaks, right = FALSE)
length_table <- table(length_groups)
for (i in 1:length(length_table)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", 
              names(length_table)[i], 
              length_table[i], 
              length_table[i]/length(seqs)*100))
}

# Key length thresholds
cat("\nLength Filter Impact:\n")
cat("Sequences <= 100bp:", sum(lengths <= 100), "(", round(sum(lengths <= 100)/length(seqs)*100, 2), "%)\n")
cat("Sequences <= 150bp:", sum(lengths <= 150), "(", round(sum(lengths <= 150)/length(seqs)*100, 2), "%)\n")
cat("Sequences <= 200bp:", sum(lengths <= 200), "(", round(sum(lengths <= 200)/length(seqs)*100, 2), "%)\n")
cat("Sequences <= 500bp:", sum(lengths <= 500), "(", round(sum(lengths <= 500)/length(seqs)*100, 2), "%)\n\n")

# 2. N Base Analysis
cat("2. N BASE ANALYSIS\n")
cat("=" , rep("=", 50), "\n", sep = "")
n_counts <- sapply(seqs, function(seq) {
  seq_str <- as.character(seq)
  lengths(gregexpr("N", seq_str))
})
n_percentages <- (n_counts / lengths) * 100

cat("Sequences containing N bases:", sum(n_counts > 0), "(", round(sum(n_counts > 0)/length(seqs)*100, 2), "%)\n")
cat("Average N percentage:", round(mean(n_percentages), 2), "%\n")
cat("Median N percentage:", round(median(n_percentages), 2), "%\n")
cat("Maximum N percentage:", round(max(n_percentages), 2), "%\n\n")

cat("N Base Filter Impact:\n")
cat("N percentage <= 10%:", sum(n_percentages <= 10), "(", round(sum(n_percentages <= 10)/length(seqs)*100, 2), "%)\n")
cat("N percentage <= 15%:", sum(n_percentages <= 15), "(", round(sum(n_percentages <= 15)/length(seqs)*100, 2), "%)\n")
cat("N percentage <= 20%:", sum(n_percentages <= 20), "(", round(sum(n_percentages <= 20)/length(seqs)*100, 2), "%)\n")
cat("N percentage <= 30%:", sum(n_percentages <= 30), "(", round(sum(n_percentages <= 30)/length(seqs)*100, 2), "%)\n\n")

# 3. Sequence Complexity Analysis
cat("3. SEQUENCE COMPLEXITY ANALYSIS\n")
cat("=" , rep("=", 50), "\n", sep = "")
cat("Calculating entropy for sample sequences (this may take a moment)...\n")

# Use a sample for entropy calculation to save time
sample_size <- min(1000, length(seqs))
sample_indices <- sample(length(seqs), sample_size)
sample_seqs <- seqs[sample_indices]

entropies <- sapply(sample_seqs, calculate_entropy)

cat("Entropy statistics (sample size:", sample_size, "):\n")
cat("Average entropy:", round(mean(entropies), 3), "\n")
cat("Median entropy:", round(median(entropies), 3), "\n")
cat("Minimum entropy:", round(min(entropies), 3), "\n")
cat("Maximum entropy:", round(max(entropies), 3), "\n\n")

cat("Complexity Filter Impact (estimated):\n")
cat("Entropy < 1.0:", sum(entropies < 1.0), "/", length(entropies), "(", round(sum(entropies < 1.0)/length(entropies)*100, 2), "%)\n")
cat("Entropy < 1.5:", sum(entropies < 1.5), "/", length(entropies), "(", round(sum(entropies < 1.5)/length(entropies)*100, 2), "%)\n")
cat("Entropy < 2.0:", sum(entropies < 2.0), "/", length(entropies), "(", round(sum(entropies < 2.0)/length(entropies)*100, 2), "%)\n\n")

# 4. Repeat Analysis
cat("4. REPEAT SEQUENCE ANALYSIS\n")
cat("=" , rep("=", 50), "\n", sep = "")
cat("Detecting repeats in sample sequences...\n")

repeat_results <- t(sapply(sample_seqs, detect_repeats))
colnames(repeat_results) <- c("Homopolymer", "Dinucleotide")

cat("Repeat statistics (sample size:", sample_size, "):\n")
cat("Sequences with homopolymer repeats:", sum(repeat_results[,1] > 0), "/", nrow(repeat_results), 
    "(", round(sum(repeat_results[,1] > 0)/nrow(repeat_results)*100, 2), "%)\n")
cat("Sequences with dinucleotide repeats:", sum(repeat_results[,2] > 0), "/", nrow(repeat_results), 
    "(", round(sum(repeat_results[,2] > 0)/nrow(repeat_results)*100, 2), "%)\n")
cat("Sequences with any repeats:", sum(rowSums(repeat_results) > 0), "/", nrow(repeat_results), 
    "(", round(sum(rowSums(repeat_results) > 0)/nrow(repeat_results)*100, 2), "%)\n\n")

# 5. Filtering Scenario Simulation
cat("5. PREPROCESSING PARAMETER OPTIMIZATION\n")
cat("=" , rep("=", 50), "\n", sep = "")

# Define different filtering scenarios
scenarios <- list(
  "Current_Strict (11-100bp, N<=15%, entropy)" = list(
    min_len = 11, max_len = 100, max_n = 15, entropy_filter = TRUE
  ),
  "Relaxed_Length (11-200bp, N<=15%, entropy)" = list(
    min_len = 11, max_len = 200, max_n = 15, entropy_filter = TRUE
  ),
  "Relaxed_N (11-200bp, N<=25%, entropy)" = list(
    min_len = 11, max_len = 200, max_n = 25, entropy_filter = TRUE
  ),
  "No_Entropy (11-200bp, N<=25%, no entropy)" = list(
    min_len = 11, max_len = 200, max_n = 25, entropy_filter = FALSE
  ),
  "Very_Relaxed (11-500bp, N<=30%, no entropy)" = list(
    min_len = 11, max_len = 500, max_n = 30, entropy_filter = FALSE
  ),
  "Length_Only (11-200bp, no other filters)" = list(
    min_len = 11, max_len = 200, max_n = 100, entropy_filter = FALSE
  )
)

cat("Simulating different preprocessing scenarios:\n\n")
results <- simulate_filtering(seqs, scenarios)
print(results, row.names = FALSE)

# 6. Recommendations
cat("\n\n6. RECOMMENDATIONS\n")
cat("=" , rep("=", 50), "\n", sep = "")

# Find the best compromise scenario
best_scenario <- results[results$Combined >= 1000 & results$Combined == max(results$Combined[results$Combined >= 1000]), ]
if (nrow(best_scenario) == 0) {
  best_scenario <- results[which.max(results$Combined), ]
}

cat("Analysis Summary:\n")
cat("- Total input sequences:", length(seqs), "\n")
cat("- Current strict filtering retains:", results$Combined[1], "sequences (", results$Percentage[1], "%)\n")
cat("- Most sequences are longer than 100bp (", round(sum(lengths > 100)/length(seqs)*100, 1), "%)\n")
cat("- N base content varies widely (max ", round(max(n_percentages), 1), "%)\n\n")

cat("Recommended preprocessing parameters:\n")
if (best_scenario$Combined[1] >= 1000) {
  cat("✓ Recommended scenario:", best_scenario$Scenario[1], "\n")
  cat("  - This will retain", best_scenario$Combined[1], "sequences (", best_scenario$Percentage[1], "%)\n")
  cat("  - Sufficient for reliable PWM model training\n\n")
} else {
  cat("⚠ Warning: Even the most relaxed filtering retains <1000 sequences\n")
  cat("  - Consider using original sequences with minimal preprocessing\n")
  cat("  - Or generate synthetic negative examples to increase dataset size\n\n")
}

# Specific parameter recommendations
cat("Suggested config.json updates:\n")
if (sum(lengths <= 200) / length(seqs) > 0.8) {
  cat('  "max_length": 200,  // Increase from 100 to retain more sequences\n')
} else {
  cat('  "max_length": 500,  // Most sequences are long, increase significantly\n')
}

if (sum(n_percentages <= 20) / length(seqs) > 0.9) {
  cat('  "max_n_percent": 20,  // Increase from 15% to be more permissive\n')
} else {
  cat('  "max_n_percent": 30,  // High N content detected, be very permissive\n')
}

cat('  "low_complexity_filter": false,  // Disable initially to retain more sequences\n')
cat('  "mask_repeats": false,  // Disable initially, can re-enable later if needed\n\n')

cat("Model Training Considerations:\n")
cat("- Minimum recommended sequences for PWM: 200-500\n")
cat("- Ideal range: 1000-5000 sequences\n")
cat("- Current analysis suggests retaining", max(results$Combined), "sequences is achievable\n")

# Close output file if writing to file
if (output_report != "") {
  sink()
  cat("Analysis complete. Report saved to:", output_report, "\n")
} else {
  cat("\nAnalysis complete.\n")
}