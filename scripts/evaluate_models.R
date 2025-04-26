# R script to evaluate PWM models using a labeled test set and ROC analysis

# --- Dependencies ---
# Check if Biostrings and pROC are installed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Package 'pROC' is needed. Please install it: install.packages('pROC')", call. = FALSE)
}
library(Biostrings)
library(pROC)

# --- Parameters ---
# Input file containing labeled test sequences in FASTA format
test_fasta_file <- "../data/test_sequences.fasta"
# Directory containing the trained PWM models (.rds files)
models_dir <- "../results/"
# Pattern to identify PWM model files
model_file_pattern <- "\.rds$" # Matches files ending in .rds

# --- Functions ---

# Function to extract class labels from FASTA headers
# Assumes header format like ">name | class=1" or ">name | class=0"
get_labels_from_fasta <- function(fasta_file) {
  headers <- names(readDNAStringSet(fasta_file))
  labels <- sub(".*\| *class=([01]).*", "\\1", headers)
  # Check if parsing worked, handle cases where pattern doesn't match
  if (any(!labels %in% c("0", "1"))) {
      warning("Could not parse class labels (0 or 1) from all FASTA headers. Check format.")
      # Attempt to return numeric, NAs for failures
      labels <- suppressWarnings(as.numeric(labels))
  } else {
      labels <- as.numeric(labels)
  }
  if (any(is.na(labels))) {
      stop("Failed to extract numeric labels from some FASTA headers. Ensure format is '>name | class=[0 or 1]'")
  }
  return(labels)
}

# Function to calculate the score of a single sequence against the PWM
# (Adapted from predict_ctcf.R - consider moving to a shared utils script later)
calculate_pwm_score <- function(sequence, pwm) {
  seq_len <- nchar(sequence)
  pwm_len <- ncol(pwm)
  
  # Ensure sequence length matches PWM length for direct scoring
  if (seq_len != pwm_len) {
    # Handle mismatch: return NA, or error, or implement sliding window max score?
    # For now, return NA as test sequences should match the model's expected length.
    warning(paste("Sequence length (", seq_len, ") does not match PWM length (", pwm_len, "). Returning NA score."))
    return(NA_real_)
  }
  
  score <- 0
  for (i in 1:seq_len) {
    base <- substr(sequence, i, i)
    if (base %in% rownames(pwm)) {
      prob <- pwm[base, i]
      if (prob > 0) { 
          score <- score + log2(prob / 0.25) # Log-likelihood ratio vs background
      } else {
          score <- score - 10 # Penalize zero probability
      }
    } else {
      score <- score - 10 # Penalize non-ACGT bases
    }
  }
  return(score)
}

# --- Main Evaluation Script ---

cat("Starting PWM Model Evaluation...\n")

# 1. Read Test Data and Labels
cat("Reading test sequences from:", test_fasta_file, "\n")
if (!file.exists(test_fasta_file)) {
  stop("Test FASTA file not found: ", test_fasta_file)
}
test_sequences <- readDNAStringSet(test_fasta_file)
true_labels <- get_labels_from_fasta(test_fasta_file)

if (length(test_sequences) != length(true_labels)) {
    stop("Mismatch between number of sequences and extracted labels.")
}
cat("Found", length(test_sequences), "test sequences with labels.\n")

# 2. Find PWM Model Files
cat("Looking for PWM models (", model_file_pattern, ") in:", models_dir, "\n")
model_files <- list.files(path = models_dir, pattern = model_file_pattern, full.names = TRUE)

if (length(model_files) == 0) {
  stop("No PWM model files (.rds) found in ", models_dir)
}
cat("Found", length(model_files), "model(s) to evaluate.\n")

# 3. Evaluate Each Model
results <- data.frame(model_file = character(), auc = numeric(), num_sequences_scored = integer(), pwm_length = integer(), stringsAsFactors = FALSE)

for (model_file in model_files) {
  cat("\nEvaluating model:", basename(model_file), "\n")
  
  # Load PWM
  pwm <- readRDS(model_file)
  pwm_len <- ncol(pwm)
  cat("  PWM length:", pwm_len, "\n")
  
  # Score test sequences
  scores <- sapply(as.character(test_sequences), calculate_pwm_score, pwm = pwm)
  
  # Filter out sequences that couldn't be scored (e.g., length mismatch)
  valid_indices <- !is.na(scores)
  scores_valid <- scores[valid_indices]
  labels_valid <- true_labels[valid_indices]
  num_scored <- length(scores_valid)
  cat("  Scored", num_scored, "out of", length(test_sequences), "sequences (matching PWM length).\n")

  if (num_scored < 2 || length(unique(labels_valid)) < 2) {
      cat("  Skipping ROC/AUC calculation: Need at least two scored sequences with both positive and negative labels.\n")
      auc_value <- NA
  } else {
      # Calculate ROC and AUC
      roc_obj <- roc(labels_valid, scores_valid, quiet = TRUE)
      auc_value <- auc(roc_obj)
      cat("  AUC:", sprintf("%.4f", auc_value), "\n")
      
      # Optional: Plot ROC curve
      # plot_file <- file.path(models_dir, paste0(gsub(".rds", "", basename(model_file)), "_roc.png"))
      # png(plot_file)
      # plot(roc_obj, main = paste("ROC Curve - ", basename(model_file)), print.auc = TRUE)
      # dev.off()
      # cat("  ROC plot saved to:", plot_file, "\n")
  }
  
  # Store results
  results <- rbind(results, data.frame(model_file = basename(model_file), 
                                        auc = auc_value, 
                                        num_sequences_scored = num_scored,
                                        pwm_length = pwm_len))
}

# 4. Display Summary
cat("\n--- Evaluation Summary ---\n")
print(results[order(-results$auc), ], row.names = FALSE)
cat("-------------------------\n")

cat("Evaluation complete.\n")
