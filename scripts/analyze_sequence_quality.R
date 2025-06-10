# Complete Sequence Quality Analysis Script for PWM Building

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/training_sequences.fasta"

# --- Main Analysis ---
cat("=== CTCF Sequence Quality Analysis Report ===\n")
cat("Input file:", input_file, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Validate input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Read sequences
cat("Reading sequences...\n")
seqs <- readDNAStringSet(input_file)
cat("Total sequences loaded:", length(seqs), "\n\n")

if (length(seqs) == 0) {
  stop("No sequences found in the input file.")
}

# 1. Length Analysis
cat("1. SEQUENCE LENGTH ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")
lengths <- width(seqs)
cat("Average length:", round(mean(lengths), 2), "bp\n")
cat("Median length:", median(lengths), "bp\n")
cat("Shortest sequence:", min(lengths), "bp\n")
cat("Longest sequence:", max(lengths), "bp\n")
cat("Standard deviation:", round(sd(lengths), 2), "bp\n")
cat("Length coefficient of variation:", round(sd(lengths)/mean(lengths), 3), "\n\n")

# Length distribution
cat("Length Distribution:\n")
length_table <- table(lengths)
if (length(length_table) <= 10) {
  print(length_table)
} else {
  # Show top 10 most common lengths
  cat("Top 10 most common lengths:\n")
  print(head(sort(length_table, decreasing = TRUE), 10))
}

cat("\nLength Uniformity Assessment:\n")
unique_lengths <- length(unique(lengths))
cat("Number of unique lengths:", unique_lengths, "\n")
if (unique_lengths == 1) {
  cat("✅ All sequences have identical length - EXCELLENT for PWM\n")
} else if (unique_lengths <= 5) {
  cat("⚠️  Few different lengths - GOOD, but may need alignment\n")
} else {
  cat("❌ Many different lengths - POOR, alignment strongly recommended\n")
}

# 2. N Base Analysis
cat("\n2. N BASE ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")
n_counts <- sapply(seqs, function(seq) {
  seq_str <- as.character(seq)
  length(gregexpr("N", seq_str, fixed = TRUE)[[1]])
})
# Fix negative counts (when no N found, gregexpr returns -1)
n_counts[n_counts < 0] <- 0
n_percentages <- (n_counts / lengths) * 100

cat("Sequences containing N bases:", sum(n_counts > 0), 
    "(", round(sum(n_counts > 0)/length(seqs)*100, 2), "%)\n")
if (sum(n_counts > 0) > 0) {
  cat("Average N percentage (in affected sequences):", 
      round(mean(n_percentages[n_counts > 0]), 2), "%\n")
  cat("Maximum N percentage:", round(max(n_percentages), 2), "%\n")
} else {
  cat("✅ No N bases found in any sequence - EXCELLENT\n")
}

# 3. Nucleotide Composition Analysis
cat("\n3. NUCLEOTIDE COMPOSITION ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")
all_seq_string <- paste(as.character(seqs), collapse = "")
total_bases <- nchar(all_seq_string)
a_count <- length(gregexpr("A", all_seq_string, fixed = TRUE)[[1]])
c_count <- length(gregexpr("C", all_seq_string, fixed = TRUE)[[1]])
g_count <- length(gregexpr("G", all_seq_string, fixed = TRUE)[[1]])
t_count <- length(gregexpr("T", all_seq_string, fixed = TRUE)[[1]])

# Fix negative counts
a_count <- max(0, a_count)
c_count <- max(0, c_count)
g_count <- max(0, g_count)
t_count <- max(0, t_count)

cat("Overall nucleotide composition:\n")
cat("A:", round(a_count/total_bases*100, 2), "%\n")
cat("C:", round(c_count/total_bases*100, 2), "%\n")
cat("G:", round(g_count/total_bases*100, 2), "%\n")
cat("T:", round(t_count/total_bases*100, 2), "%\n")

gc_content <- (g_count + c_count) / total_bases * 100
cat("GC content:", round(gc_content, 2), "%\n")

if (gc_content >= 40 && gc_content <= 60) {
  cat("✅ GC content is balanced - GOOD for PWM\n")
} else if (gc_content < 30 || gc_content > 70) {
  cat("⚠️  Extreme GC content - may affect PWM quality\n")
} else {
  cat("⚠️  Moderate GC bias - acceptable for PWM\n")
}

# 4. Sequence Complexity Analysis
cat("\n4. SEQUENCE COMPLEXITY ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# Simple complexity measure: count of unique 3-mers
calculate_complexity <- function(seq_str) {
  if (nchar(seq_str) < 3) return(0)
  trimers <- character()
  for (i in 1:(nchar(seq_str) - 2)) {
    trimers <- c(trimers, substr(seq_str, i, i + 2))
  }
  return(length(unique(trimers)) / length(trimers))
}

complexities <- sapply(as.character(seqs), calculate_complexity)
cat("Average sequence complexity (unique 3-mers ratio):", round(mean(complexities), 3), "\n")
cat("Complexity range:", round(min(complexities), 3), "-", round(max(complexities), 3), "\n")

low_complexity_count <- sum(complexities < 0.5)
if (low_complexity_count > 0) {
  cat("⚠️  Low complexity sequences:", low_complexity_count, 
      "(", round(low_complexity_count/length(seqs)*100, 1), "%)\n")
} else {
  cat("✅ All sequences have good complexity\n")
}

# 5. Recommendations
cat("\n5. RECOMMENDATIONS FOR PWM BUILDING\n")
cat(paste(rep("=", 50), collapse=""), "\n")

cat("Based on the analysis:\n\n")

# Length recommendations
if (unique_lengths == 1) {
  cat("✅ LENGTH: Sequences are uniform - proceed with PWM building\n")
} else if (sd(lengths)/mean(lengths) > 0.3) {
  cat("❌ LENGTH: High variability detected\n")
  cat("   Recommendation: Align sequences or extract core regions\n")
} else {
  cat("⚠️  LENGTH: Moderate variability\n")
  cat("   Recommendation: Consider alignment for better results\n")
}

# Quality recommendations
if (sum(n_counts > 0) == 0) {
  cat("✅ QUALITY: No N bases - excellent quality\n")
} else if (mean(n_percentages) > 10) {
  cat("❌ QUALITY: High N content\n")
  cat("   Recommendation: Filter sequences with >10% N bases\n")
} else {
  cat("⚠️  QUALITY: Some N bases present but acceptable\n")
}

# Complexity recommendations
if (low_complexity_count/length(seqs) > 0.2) {
  cat("❌ COMPLEXITY: Many low-complexity sequences\n")
  cat("   Recommendation: Apply complexity filtering\n")
} else {
  cat("✅ COMPLEXITY: Good sequence complexity\n")
}

# Sample size recommendations
if (length(seqs) < 100) {
  cat("❌ SAMPLE SIZE: Too few sequences (<100)\n")
  cat("   Recommendation: Collect more training data\n")
} else if (length(seqs) < 500) {
  cat("⚠️  SAMPLE SIZE: Moderate number of sequences\n")
  cat("   Recommendation: More sequences would improve PWM reliability\n")
} else {
  cat("✅ SAMPLE SIZE: Good number of sequences (>=500)\n")
}

cat("\n=== Analysis Complete ===\n")
