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

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}
library(optparse)

# ---- Parse CLI arguments ----
option_list <- list(
  make_option("--pwm",      type = "character", help = "PWM model path"),
  make_option("--test",     type = "character", help = "FASTA file of test sequences"),
  make_option("--out",      type = "character", default = NULL, help = "Output metrics CSV path"),
  make_option("--strategy", type = "character", default = "unknown", help = "Split strategy name")
)
opt <- parse_args(OptionParser(option_list = option_list))


# --- Parameters ---
# Input file containing labeled test sequences in FASTA format
test_fasta_file <- opt$test
# Directory containing the trained PWM models (.rds files)
models_dir <- dirname(opt$pwm)
# Pattern to identify PWM model files
model_file_pattern <- paste0("^", basename(opt$pwm), "$") # Matches files ending in .rds

# --- Functions ---

# Function to extract class labels from FASTA headers
# Assumes header format like ">name | class=1" or ">name | class=0"
get_labels_from_fasta <- function(fasta_file) {
  headers <- names(readDNAStringSet(fasta_file))
  labels <- sub(".*\\| *class=([01]).*", "\\1", headers)
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
      fraction_sig <- NA
  } else {
      # Calculate ROC and AUC
      roc_obj <- roc(labels_valid, scores_valid, quiet = TRUE)
      auc_obj <- auc(roc_obj)
      auc_value <- as.numeric(auc_obj)
      ## ====== 再計 GC-matched p/q =====================================
      pos_idx <- labels_valid == 1
      neg_idx <- labels_valid == 0
      if (sum(neg_idx) == 0) {
          stop("No class=0 sequences ⇒ 無 GC 背景可用")
      }
      S_pos   <- scores_valid[pos_idx]
      S_neg   <- scores_valid[neg_idx]
      pvals   <- sapply(S_pos, function(s) mean(S_neg >= s))
      qvals   <- p.adjust(pvals, "BH")
      fraction_sig <- mean(qvals < 0.05)
  }
  
  # Store results
  results <- rbind(results, data.frame(
       model_file           = basename(model_file),
       auc                  = auc_value,
       num_sequences_scored = num_scored,
       pwm_length           = pwm_len,
       sig_frac_q05         = fraction_sig      # <-- 新欄位
   ))
}

# 4. Display Summary
cat("\n--- Evaluation Summary ---\n")
print(results[order(-results$auc), ], row.names = FALSE)
cat("-------------------------\n")

cat("Evaluation complete.\n")



# Write metrics if --out is specified
if (!is.null(opt$out)) {
  results <- dplyr::rename(results, AUC = auc)
  write.csv(results, opt$out, row.names = FALSE)   # 把整張 results 直接輸出
  cat("✔ All model metrics written to", opt$out, "\n")
}
