# filepath: /mnt/d/workspace/data-science/scripts/evaluate_models_with_cv.R
# R script to evaluate PWM models using stratified k-fold cross-validation

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Package 'pROC' is needed. Please install: install.packages('pROC')",
       call. = FALSE)
}
library(Biostrings)
library(pROC)

# --- Parameters ---
# File paths
input_fasta <- "../data/extracted_sequences.fasta"
pwm_dir <- "../results/"
pwm_file_pattern <- ".*\\.rds$"  # Regular expression matching PWM files

# Cross-validation parameters
k_folds <- 5    # Number of folds for cross-validation
num_repeats <- 3  # Number of repeated cross-validations for more reliable estimates
set.seed(123)   # For reproducibility

# --- Utility Functions ---

# Extract motif matches from sequences using PWM
scan_sequences_with_pwm <- function(sequences, pwm) {
  pwm_width <- ncol(pwm)
  scores <- numeric(length(sequences))
  
  for (i in 1:length(sequences)) {
    seq <- sequences[[i]]
    seq_length <- width(seq)
    
    # Skip if sequence is shorter than PWM
    if (seq_length < pwm_width) {
      scores[i] <- NA
      next
    }
    
    # Scan sequence with PWM
    max_score <- -Inf
    for (j in 1:(seq_length - pwm_width + 1)) {
      subseq <- subseq(seq, j, j + pwm_width - 1)
      score <- calculate_pwm_score(subseq, pwm)
      if (score > max_score) {
        max_score <- score
      }
    }
    scores[i] <- max_score
  }
  
  return(scores)
}

# Calculate PWM score for a sequence
calculate_pwm_score <- function(seq, pwm) {
  seq_string <- as.character(seq)
  score <- 0
  for (i in 1:nchar(seq_string)) {
    base <- substr(seq_string, i, i)
    if (base %in% c("A", "C", "G", "T")) {
      score <- score + pwm[base, i]
    }
  }
  return(score)
}

# Extract class labels from sequence names
extract_class_labels <- function(seq_names) {
  class_labels <- integer(length(seq_names))
  
  for (i in 1:length(seq_names)) {
    name <- seq_names[i]
    # Look for "class=1" or "class=0" in the header
    if (grepl("class=1", name)) {
      class_labels[i] <- 1
    } else if (grepl("class=0", name)) {
      class_labels[i] <- 0
    } else {
      warning("No class label found in sequence name: ", name)
      class_labels[i] <- NA
    }
  }
  
  return(class_labels)
}

# Create stratified folds - ensures each fold has same proportion of positive and negative examples
create_stratified_folds <- function(labels, k) {
  # Indices for positive and negative examples
  pos_idx <- which(labels == 1)
  neg_idx <- which(labels == 0)
  
  # Shuffle indices
  pos_idx <- sample(pos_idx)
  neg_idx <- sample(neg_idx)
  
  # Create folds
  pos_folds <- split(pos_idx, cut(seq_along(pos_idx), k, labels = FALSE))
  neg_folds <- split(neg_idx, cut(seq_along(neg_idx), k, labels = FALSE))
  
  # Combine positive and negative indices for each fold
  folds <- vector("list", k)
  for (i in 1:k) {
    folds[[i]] <- c(pos_folds[[i]], neg_folds[[i]])
  }
  
  return(folds)
}

# Perform cross-validation evaluation for a single PWM
evaluate_pwm_with_cv <- function(sequences, labels, pwm, k_folds, repeats) {
  n_samples <- length(sequences)
  
  # Store AUC values for each repeat and fold
  all_auc_values <- numeric(k_folds * repeats)
  all_thresholds <- list()
  
  fold_count <- 1
  
  # Repeat cross-validation multiple times for more reliable estimates
  for (repeat_idx in 1:repeats) {
    # Create stratified folds
    folds <- create_stratified_folds(labels, k_folds)
    
    # For each fold
    for (fold_idx in 1:k_folds) {
      # Indices for test set (current fold)
      test_indices <- folds[[fold_idx]]
      # Indices for training set (all other folds)
      train_indices <- setdiff(1:n_samples, test_indices)
      
      # Split data into training and test sets
      train_sequences <- sequences[train_indices]
      test_sequences <- sequences[test_indices]
      train_labels <- labels[train_indices]
      test_labels <- labels[test_indices]
      
      # Calculate scores for test set using PWM
      test_scores <- scan_sequences_with_pwm(test_sequences, pwm)
      
      # Calculate ROC and AUC
      if (length(unique(test_labels)) > 1 && !any(is.na(test_scores))) {
        roc_result <- pROC::roc(test_labels, test_scores, quiet = TRUE)
        all_auc_values[fold_count] <- as.numeric(pROC::auc(roc_result))
        
        # Find optimal threshold (Youden's J statistic - maximizes sensitivity + specificity - 1)
        coords <- pROC::coords(roc_result, "best", best.method = "youden")
        all_thresholds[[fold_count]] <- list(
          threshold = coords$threshold,
          sensitivity = coords$sensitivity,
          specificity = coords$specificity
        )
      } else {
        warning("Skipping fold ", fold_idx, " in repeat ", repeat_idx, " due to insufficient class diversity or NA scores.")
        all_auc_values[fold_count] <- NA
      }
      
      fold_count <- fold_count + 1
    }
  }
  
  # Remove any NA values
  all_auc_values <- all_auc_values[!is.na(all_auc_values)]
  
  # Calculate mean and standard deviation of AUC values
  mean_auc <- mean(all_auc_values)
  sd_auc <- sd(all_auc_values)
  
  # Calculate average optimal threshold
  thresholds <- sapply(all_thresholds, function(x) ifelse(is.null(x), NA, x$threshold))
  thresholds <- thresholds[!is.na(thresholds)]
  mean_threshold <- mean(thresholds)
  
  return(list(
    mean_auc = mean_auc,
    sd_auc = sd_auc,
    all_auc_values = all_auc_values,
    mean_threshold = mean_threshold,
    all_thresholds = all_thresholds
  ))
}

# --- Main Script ---

cat("Starting PWM model evaluation with cross-validation...\n")

# 1. Read Input Sequences
cat("Reading sequences from:", input_fasta, "\n")
if (!file.exists(input_fasta)) {
  stop("Input FASTA file not found: ", input_fasta, 
       ". Run download_data.sh and prepare_datasets.R first.")
}
all_sequences <- readDNAStringSet(input_fasta)
cat("Read", length(all_sequences), "sequences.\n")

# 2. Extract class labels
labels <- extract_class_labels(names(all_sequences))

# Check if we have both positive and negative examples
if (length(unique(labels[!is.na(labels)])) < 2) {
  stop("Both positive (class=1) and negative (class=0) examples are required for evaluation.",
       " Check that your input file has sequences with both '| class=1' and '| class=0' in headers.")
}

# Remove sequences without valid labels
valid_idx <- !is.na(labels)
all_sequences <- all_sequences[valid_idx]
labels <- labels[valid_idx]

cat("Using", sum(labels == 1), "positive sequences and", 
    sum(labels == 0), "negative sequences for evaluation.\n")

# 3. Get list of PWM files
pwm_files <- list.files(pwm_dir, pattern = pwm_file_pattern, full.names = TRUE)
if (length(pwm_files) == 0) {
  stop("No PWM files found in directory: ", pwm_dir, 
       ". Run build_pwm.R first to generate PWMs.")
}

# 4. Evaluate each PWM with cross-validation
cat("\nPerforming", k_folds, "fold cross-validation with", num_repeats, "repeats:\n\n")

results <- data.frame(
  PWM_File = character(),
  Mean_AUC = numeric(),
  SD_AUC = numeric(),
  Mean_Threshold = numeric(),
  stringsAsFactors = FALSE
)

for (pwm_file in pwm_files) {
  cat("Evaluating", basename(pwm_file), "...\n")
  
  # Load PWM
  pwm <- readRDS(pwm_file)
  
  # Run cross-validation
  cv_results <- evaluate_pwm_with_cv(all_sequences, labels, pwm, k_folds, num_repeats)
  
  # Print results
  cat("  Mean AUC:", round(cv_results$mean_auc, 4), 
      "±", round(cv_results$sd_auc, 4), "\n")
  cat("  Suggested score threshold:", round(cv_results$mean_threshold, 4), "\n")
  
  # Add to results table
  results <- rbind(results, data.frame(
    PWM_File = basename(pwm_file),
    Mean_AUC = cv_results$mean_auc,
    SD_AUC = cv_results$sd_auc,
    Mean_Threshold = cv_results$mean_threshold,
    stringsAsFactors = FALSE
  ))
}

# 5. Print summary table
cat("\nSummary of Evaluation Results:\n")
print(results)

# 6. Write results to file
results_file <- file.path(pwm_dir, "cv_evaluation_results.csv")
write.csv(results, file = results_file, row.names = FALSE)
cat("\nResults saved to:", results_file, "\n")

cat("\nCross-validation evaluation complete!\n")