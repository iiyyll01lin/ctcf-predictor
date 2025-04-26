# filepath: /mnt/d/workspace/data-science/scripts/optimize_threshold.R
# R script to find the optimal score threshold based on ROC analysis

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

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || length(args) > 4) {
  stop("Usage: Rscript optimize_threshold.R <test_fasta> <pwm_file> [optimization_method] [output_json]
       \n\t<test_fasta>: FASTA file with labeled sequences ('class=0' or 'class=1' in headers)
       \n\t<pwm_file>: RDS file containing the PWM model
       \n\t[optimization_method]: Optional. One of 'youden' (default), 'sensitivity_specificity', 'closest_topleft', or 'balanced'
       \n\t[output_json]: Optional. JSON file to save optimization results. Default: '../results/threshold_optimization.json'", 
       call. = FALSE)
}

test_fasta_file <- args[1]
pwm_file <- args[2]
optimization_method <- if (length(args) >= 3) args[3] else "youden"
output_json <- if (length(args) >= 4) args[4] else "../results/threshold_optimization.json"

# Validate method
valid_methods <- c("youden", "sensitivity_specificity", "closest_topleft", "balanced")
if (!(optimization_method %in% valid_methods)) {
  stop("Invalid optimization method. Choose one of: ", paste(valid_methods, collapse=", "))
}

# --- Utility Functions ---

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

# Find optimal threshold based on ROC analysis using different methods
find_optimal_threshold <- function(roc_obj, method = "youden") {
  if (method == "youden") {
    # Youden's J statistic (maximizes sensitivity + specificity - 1)
    coords <- coords(roc_obj, "best", best.method = "youden")
    threshold <- coords$threshold
    sensitivity <- coords$sensitivity
    specificity <- coords$specificity
    description <- "最大化 Youden's J 統計量 (敏感性 + 特異性 - 1)"
  } else if (method == "sensitivity_specificity") {
    # Equal sensitivity and specificity
    coords <- coords(roc_obj, "best", best.method = "closest.topleft")
    threshold <- coords$threshold
    sensitivity <- coords$sensitivity
    specificity <- coords$specificity
    description <- "敏感性和特異性最接近相等"
  } else if (method == "closest_topleft") {
    # Closest to top-left corner (0,1)
    coords <- coords(roc_obj, "best", best.method = "closest.topleft")
    threshold <- coords$threshold
    sensitivity <- coords$sensitivity
    specificity <- coords$specificity
    description <- "最接近 ROC 曲線左上角"
  } else if (method == "balanced") {
    # Find threshold where F1 score is maximized
    coords_list <- coords(roc_obj, "all")
    
    # Calculate F1 scores for all thresholds
    TP <- coords_list$tp
    TN <- coords_list$tn
    FP <- coords_list$fp
    FN <- coords_list$fn
    
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    
    # Handle division by zero
    precision[is.nan(precision)] <- 0
    recall[is.nan(recall)] <- 0
    
    F1 <- 2 * (precision * recall) / (precision + recall)
    F1[is.nan(F1)] <- 0
    
    # Find threshold with max F1
    best_idx <- which.max(F1)
    
    threshold <- coords_list$threshold[best_idx]
    sensitivity <- coords_list$sensitivity[best_idx]
    specificity <- coords_list$specificity[best_idx]
    description <- "平衡精準度和召回率 (最大化 F1 分數)"
  } else {
    stop("Unknown optimization method: ", method)
  }
  
  return(list(
    threshold = threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    method = method,
    description = description
  ))
}

# Create a JSON string from threshold results
create_json_output <- function(results, method) {
  # Calculate additional metrics
  ppv <- results$tp / (results$tp + results$fp)  # Positive Predictive Value (Precision)
  npv <- results$tn / (results$tn + results$fn)  # Negative Predictive Value
  accuracy <- (results$tp + results$tn) / (results$tp + results$tn + results$fp + results$fn)
  f1_score <- 2 * (ppv * results$sensitivity) / (ppv + results$sensitivity)
  
  # Handle NaN values
  ppv <- if(is.nan(ppv)) 0 else ppv
  npv <- if(is.nan(npv)) 0 else npv
  f1_score <- if(is.nan(f1_score)) 0 else f1_score
  
  json_str <- paste0(
    '{\n',
    '  "threshold": ', results$threshold, ',\n',
    '  "optimization_method": "', method, '",\n',
    '  "description": "', results$description, '",\n',
    '  "metrics": {\n',
    '    "sensitivity": ', results$sensitivity, ',\n',
    '    "specificity": ', results$specificity, ',\n',
    '    "precision": ', ppv, ',\n',
    '    "negative_predictive_value": ', npv, ',\n',
    '    "accuracy": ', accuracy, ',\n',
    '    "f1_score": ', f1_score, '\n',
    '  },\n',
    '  "confusion_matrix": {\n',
    '    "true_positive": ', results$tp, ',\n',
    '    "false_positive": ', results$fp, ',\n',
    '    "true_negative": ', results$tn, ',\n',
    '    "false_negative": ', results$fn, '\n',
    '  },\n',
    '  "roc_auc": ', results$auc, '\n',
    '}'
  )
  
  return(json_str)
}

# --- Main Script ---
cat("Starting threshold optimization...\n")

# 1. Validate inputs
if (!file.exists(test_fasta_file)) {
  stop("Test FASTA file not found: ", test_fasta_file)
}
if (!file.exists(pwm_file)) {
  stop("PWM file not found: ", pwm_file)
}

# 2. Read test sequences and labels
cat("Reading test sequences from:", test_fasta_file, "\n")
test_sequences <- readDNAStringSet(test_fasta_file)
labels <- extract_class_labels(names(test_sequences))

# Remove sequences without valid labels
valid_idx <- !is.na(labels)
test_sequences <- test_sequences[valid_idx]
labels <- labels[valid_idx]

# Check if we have both positive and negative examples
if (length(unique(labels)) < 2) {
  stop("Both positive (class=1) and negative (class=0) examples are required for threshold optimization.",
       " Check that your input file has sequences with both '| class=1' and '| class=0' in headers.")
}

cat("Using", sum(labels == 1), "positive sequences and", 
    sum(labels == 0), "negative sequences for optimization.\n")

# 3. Load PWM and calculate scores
cat("Loading PWM from:", pwm_file, "\n")
pwm <- readRDS(pwm_file)
pwm_width <- ncol(pwm)
cat("PWM width:", pwm_width, "\n")

cat("Calculating scores for all sequences...\n")
scores <- scan_sequences_with_pwm(test_sequences, pwm)

# Remove sequences with NA scores
valid_score_idx <- !is.na(scores)
if (sum(!valid_score_idx) > 0) {
  cat("Warning:", sum(!valid_score_idx), "sequences skipped due to length mismatch with PWM.\n")
}
scores <- scores[valid_score_idx]
labels <- labels[valid_score_idx]

# 4. Perform ROC analysis
cat("Performing ROC analysis...\n")
roc_result <- roc(labels, scores, quiet = TRUE)
auc_value <- auc(roc_result)
cat("AUC:", sprintf("%.4f", auc_value), "\n")

# 5. Find optimal threshold using the specified method
cat("Finding optimal threshold using method:", optimization_method, "\n")
threshold_result <- find_optimal_threshold(roc_result, method = optimization_method)

# 6. Calculate confusion matrix at optimal threshold
predicted_positive <- scores >= threshold_result$threshold
true_positive <- labels == 1

tp <- sum(predicted_positive & true_positive)
tn <- sum(!predicted_positive & !true_positive)
fp <- sum(predicted_positive & !true_positive)
fn <- sum(!predicted_positive & true_positive)

threshold_result$tp <- tp
threshold_result$tn <- tn
threshold_result$fp <- fp
threshold_result$fn <- fn
threshold_result$auc <- auc_value

# 7. Print results
cat("\n--- Threshold Optimization Results ---\n")
cat("Optimization Method:", optimization_method, "\n")
cat("Description:", threshold_result$description, "\n")
cat("Optimal Threshold:", sprintf("%.4f", threshold_result$threshold), "\n")
cat("Sensitivity:", sprintf("%.4f", threshold_result$sensitivity), "\n")
cat("Specificity:", sprintf("%.4f", threshold_result$specificity), "\n")
cat("Confusion Matrix:\n")
cat("  True Positive:", tp, "\n")
cat("  False Positive:", fp, "\n")
cat("  True Negative:", tn, "\n")
cat("  False Negative:", fn, "\n")
cat("AUC:", sprintf("%.4f", auc_value), "\n")
cat("--------------------------------------\n")

# 8. Save results to JSON
cat("Saving results to:", output_json, "\n")
json_output <- create_json_output(threshold_result, optimization_method)
write(json_output, file = output_json)

cat("\nUse this threshold in the prediction script as follows:\n")
cat("Rscript scripts/predict_ctcf.R <input_fasta> <output_tsv> <pwm_file>", threshold_result$threshold, "\n\n")

cat("Threshold optimization complete.\n")