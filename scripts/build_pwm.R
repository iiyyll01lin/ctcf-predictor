# R script to build a Position Weight Matrix (PWM) from training sequences

# --- Dependencies ---
# Check if Biostrings is installed, if not, provide installation instructions
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed for this script. Please install it by running: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
library(Biostrings)

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}
library(Biostrings)
library(optparse)

# --- Parse CLI options ---
option_list <- list(
  make_option("--input",  type = "character", default = "data/training_sequences.fasta"),
  make_option("--output", type = "character", default = "results/generated_pwm.rds"),
  make_option("--only_positive", action = "store_true", default = FALSE,
              help = "Keep only sequences whose name contains 'class=1'")
)
opt <- parse_args(OptionParser(option_list = option_list))

fasta_file         <- opt$input
output_pwm_file    <- opt$output
output_pwm_txt_file <- sub("\\.rds$", ".txt", output_pwm_file)



# --- Main Script ---

cat("Reading sequences from:", fasta_file, "\n")

# Read sequences
sequences <- readDNAStringSet(fasta_file)

# 把 only-positive 過濾段插在這裡
if (opt$only_positive) {
  keep <- grepl("class=1", names(sequences))
  sequences <- sequences[keep]
  cat("Keeping", length(sequences), "positive sequences (class=1).\n")
  if (length(sequences) == 0)
    stop("No positive sequences found after filtering.")
}

# Basic validation: Check if sequences are non-empty and have consistent length
if (length(sequences) == 0) {
  stop("No sequences found in the input file.")
}
seq_lengths <- width(sequences)
if (length(unique(seq_lengths)) > 1) {
  warning("Sequences have varying lengths. Using the length of the first sequence: ", seq_lengths[1])
  # Optional: Add logic here to handle varying lengths (e.g., trimming, padding, or erroring)
  # For simplicity, we'll proceed assuming alignment or focusing on a core region
  # Let's filter for the most common length or the first sequence's length
  target_length <- seq_lengths[1]
  sequences <- sequences[width(sequences) == target_length]
  if (length(sequences) == 0) {
      stop("No sequences remaining after filtering for length ", target_length)
  }
  cat("Filtered sequences to length:", target_length, "\n")
} else {
    target_length <- seq_lengths[1]
    cat("All sequences have length:", target_length, "\n")
}


cat("Calculating consensus matrix...\n")
# Calculate the consensus matrix (counts of A, C, G, T at each position)
# Note: This handles IUPAC ambiguity codes by default, but our example doesn't have them beyond 'N'
# We might want to exclude sequences with too many Ns or handle them specifically.
# For simplicity, consensusMatrix ignores Ns by default.
consensus_mat <- consensusMatrix(sequences, as.prob = FALSE)

# Remove rows other than A, C, G, T if they exist (like '-', '+', '.')
valid_bases <- c("A", "C", "G", "T")
consensus_mat <- consensus_mat[rownames(consensus_mat) %in% valid_bases, , drop = FALSE]

cat("Converting counts to probabilities (PWM)...\n")
# Add pseudocounts to avoid zero probabilities (common practice)
pseudocount <- 1
consensus_mat_pseudo <- consensus_mat + pseudocount

# Calculate probabilities (PWM)
pwm <- prop.table(consensus_mat_pseudo, margin = 2) # Normalize columns

cat("Generated PWM:\n")
print(pwm)

# Save the PWM object to an RDS file for use in other R scripts
cat("\nSaving PWM object to:", output_pwm_file, "\n")
saveRDS(pwm, file = output_pwm_file)

# Save the PWM as a text file (optional)
cat("Saving PWM matrix to text file:", output_pwm_txt_file, "\n")
write.table(pwm, file = output_pwm_txt_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("\nPWM generation complete.\n")
