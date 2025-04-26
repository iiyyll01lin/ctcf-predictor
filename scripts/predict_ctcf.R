# CTCF Binding Site Prediction Script

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript predict_ctcf.R <input_fasta> <output_tsv> <pwm_rds_file> <score_threshold>", call. = FALSE)
}

input_fasta_file <- args[1]
output_tsv_file <- args[2]
pwm_file <- args[3]
prediction_threshold <- as.numeric(args[4])

if (is.na(prediction_threshold)) {
    stop("Score threshold must be a numeric value.", call. = FALSE)
}

# --- Validate Inputs ---
if (!file.exists(input_fasta_file)) {
  stop("Input FASTA file not found: ", input_fasta_file, call. = FALSE)
}
if (!file.exists(pwm_file)) {
  stop("PWM file not found: ", pwm_file, call. = FALSE)
}

# --- Load PWM ---
cat("Loading PWM from:", pwm_file, "\n")
pwm <- readRDS(pwm_file)
pwm_length <- ncol(pwm)
if (is.null(pwm_length)) {
    stop("Could not determine PWM length from the loaded object.")
}
cat("Using PWM of length:", pwm_length, "\n")

# --- Functions ---

# Function to calculate the score of a sequence against the PWM
# Uses log-likelihood scoring for simplicity, assuming background probability of 0.25 for each base
calculate_pwm_score <- function(sequence, pwm) {
  seq_len <- nchar(sequence)
  pwm_len <- ncol(pwm)
  
  # This function expects sequence length to match PWM length
  if (seq_len != pwm_len) {
    # This case should ideally not be hit if called from scan_sequence
    warning(paste("Internal Error: calculate_pwm_score called with mismatched length:", seq_len, "vs", pwm_len))
    return(NA_real_)
  }
  
  score <- 0
  for (i in 1:seq_len) {
    base <- substr(sequence, i, i)
    if (base %in% rownames(pwm)) {
      prob <- pwm[base, i]
      if (prob > 0) { 
          score <- score + log2(prob / 0.25)
      } else {
          score <- score - 10 # Penalize zero probability
      }
    } else {
      # Handle non-ACGT characters (including N)
      # Option 1: Penalize
      score <- score - 10 
      # Option 2: Treat as neutral (score += 0) - less common for PWMs
      # Option 3: Distribute probability (complex)
    }
  }
  return(score)
}

# Function to scan a sequence with a sliding window
scan_sequence <- function(long_sequence, pwm, threshold = 0) {
  pwm_len <- ncol(pwm)
  long_seq_len <- nchar(long_sequence)
  results <- data.frame(start = integer(), end = integer(), sequence = character(), score = numeric(), stringsAsFactors = FALSE)
  
  if (long_seq_len < pwm_len) {
      # Return empty frame, warning issued by caller if needed
      return(results)
  }

  for (i in 1:(long_seq_len - pwm_len + 1)) {
    sub_sequence <- substr(long_sequence, i, i + pwm_len - 1)
    # Skip scoring if subsequence contains non-ACGT characters? 
    # Or handle in calculate_pwm_score (current approach)
    score <- calculate_pwm_score(sub_sequence, pwm)
    
    if (!is.na(score) && score >= threshold) {
      results <- rbind(results, data.frame(start = i, end = i + pwm_len - 1, sequence = sub_sequence, score = score))
    }
  }
  
  return(results)
}

# --- Main Processing ---

cat("Reading input sequences from:", input_fasta_file, "\n")
input_sequences <- readDNAStringSet(input_fasta_file)
cat("Read", length(input_sequences), "sequences.\n")

all_results <- list()

cat("Scanning sequences with threshold:", prediction_threshold, "\n")
for (i in 1:length(input_sequences)) {
  seq_name <- names(input_sequences)[i]
  sequence <- as.character(input_sequences[[i]])
  
  if (nchar(sequence) < pwm_length) {
      cat("  Skipping sequence", seq_name, "(length", nchar(sequence), ") - shorter than PWM length (", pwm_length, ").\n")
      next
  }
  
  cat("  Scanning sequence:", seq_name, "(length", nchar(sequence), ")...\n")
  predicted_sites <- scan_sequence(sequence, pwm, threshold = prediction_threshold)
  
  if (nrow(predicted_sites) > 0) {
    predicted_sites$sequence_name <- seq_name
    # Reorder columns
    predicted_sites <- predicted_sites[, c("sequence_name", "start", "end", "sequence", "score")]
    all_results[[seq_name]] <- predicted_sites
  }
}

# Combine results from all sequences
if (length(all_results) > 0) {
  final_results_df <- do.call(rbind, all_results)
  rownames(final_results_df) <- NULL # Reset row names
  cat("Found", nrow(final_results_df), "potential binding sites in total.\n")
  
  # Write results to output file
  cat("Writing results to:", output_tsv_file, "\n")
  write.table(final_results_df, file = output_tsv_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
} else {
  cat("No potential binding sites found above the threshold in any input sequence.\n")
  # Create an empty file with headers
  cat("Writing empty results file with headers to:", output_tsv_file, "\n")
  write.table(data.frame(sequence_name=character(), start=integer(), end=integer(), sequence=character(), score=numeric()), 
              file = output_tsv_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

cat("Prediction script finished.\n")

# --- Old Example Usage Removed ---
# (The code previously here using hardcoded sequences is removed)

