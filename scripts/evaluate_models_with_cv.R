# filepath: scripts/evaluate_models_with_cv.R
# R script to evaluate PWM models using stratified k-fold cross-validation,
# 並依照傳入的 strategy_name（如 "cv" / "cv_uniform"）輸出 per-fold AUC

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("Biostrings")
}
if (!requireNamespace("pROC", quietly = TRUE)) {
  install.packages("pROC", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}

library(Biostrings)
library(pROC)
library(optparse)

# --- Parameters ---
option_list <- list(
  make_option("--fasta",      type="character", default="../data/extracted_sequences.fasta"),
  make_option("--pwm_dir",    type="character", default="../results/"),
  make_option("--k",          type="integer",   default=5),
  make_option("--repeats",    type="integer",   default=3, dest = "repeats"),
  make_option("--pwm_files",  type="character", default=NULL,
               help = "逗號分隔的 RDS 檔案清單（若指定，則只讀這些檔案，不掃整個 pwm_dir）")
)
opt <- parse_args(OptionParser(option_list = option_list))

input_fasta <- opt$fasta
pwm_dir     <- opt$pwm_dir
k_folds     <- opt$k
num_repeats <- opt$repeats
wanted      <- opt$pwm_files

pwm_file_pattern <- ".*\\.rds$"
set.seed(123)   # For reproducibility

# --- Utility Functions ---

# Extract motif matches from sequences using PWM
scan_sequences_with_pwm <- function(sequences, pwm) {
  pwm_width <- ncol(pwm)
  scores <- numeric(length(sequences))
  
  for (i in seq_along(sequences)) {
    seq <- sequences[[i]]
    seq_length <- nchar(seq)
    
    # Skip if sequence is shorter than PWM
    if (seq_length < pwm_width) {
      scores[i] <- NA
      next
    }
    
    # Scan sequence with PWM (maximum sub-window score)
    max_score <- -Inf
    for (j in 1:(seq_length - pwm_width + 1)) {
      subseq_window <- subseq(seq, j, j + pwm_width - 1)
      score <- calculate_pwm_score(subseq_window, pwm)
      if (score > max_score) {
        max_score <- score
      }
    }
    scores[i] <- max_score
  }
  
  return(scores)
}

# Calculate PWM score for a sequence (simple sum of probabilities)
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
  
  for (i in seq_along(seq_names)) {
    name <- seq_names[i]
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
  pos_idx <- which(labels == 1)
  neg_idx <- which(labels == 0)
  
  pos_idx <- sample(pos_idx)
  neg_idx <- sample(neg_idx)
  
  pos_folds <- split(pos_idx, cut(seq_along(pos_idx), k, labels = FALSE))
  neg_folds <- split(neg_idx, cut(seq_along(neg_idx), k, labels = FALSE))
  
  folds <- vector("list", k)
  for (i in 1:k) {
    folds[[i]] <- c(pos_folds[[i]], neg_folds[[i]])
  }
  
  return(folds)
}

# Perform cross-validation evaluation for a single PWM, 同時回傳 per-fold 結果
evaluate_pwm_with_cv <- function(sequences, labels, pwm, k_folds, repeats, strategy_name) {
  n_samples <- length(sequences)
  
  all_auc_values <- numeric(k_folds * repeats)
  all_thresholds <- list()
  
  # Collect per‐fold records in a local data.frame, 然後再回傳
  fold_records <- data.frame(
    strategy = character(),
    repeats  = integer(),
    fold     = integer(),
    AUC      = numeric(),
    stringsAsFactors = FALSE
  )
  
  fold_count <- 1
  
  for (repeat_idx in 1:repeats) {
    folds <- create_stratified_folds(labels, k_folds)
    
    for (fold_idx in 1:k_folds) {
      test_indices <- folds[[fold_idx]]
      train_indices <- setdiff(seq_len(n_samples), test_indices)
      
      train_sequences <- sequences[train_indices]
      test_sequences  <- sequences[test_indices]
      test_labels     <- labels[test_indices]
      
      test_scores <- scan_sequences_with_pwm(test_sequences, pwm)
      
      if (length(unique(test_labels)) > 1 && !any(is.na(test_scores))) {
        roc_result <- pROC::roc(test_labels, test_scores, quiet = TRUE)
        auc_val    <- as.numeric(pROC::auc(roc_result))
        all_auc_values[fold_count] <- auc_val
        
        # 收集單一 fold 的結果，把 strategy 設為 strategy_name
        fold_records <- rbind(
          fold_records,
          data.frame(
            strategy = strategy_name,
            repeats  = repeat_idx,
            fold     = fold_idx,
            AUC      = auc_val,
            stringsAsFactors = FALSE
          )
        )
        
        coords <- pROC::coords(roc_result, "best", best.method = "youden")
        all_thresholds[[fold_count]] <- list(
          threshold   = coords$threshold,
          sensitivity = coords$sensitivity,
          specificity = coords$specificity
        )
      } else {
        warning("Skipping fold ", fold_idx,
                " in repeat ", repeat_idx,
                " (strategy=", strategy_name, ") due to insufficient class diversity or NA scores.")
        all_auc_values[fold_count] <- NA
      }
      
      fold_count <- fold_count + 1
    }
  }
  
  # Remove NA values
  valid_auc <- all_auc_values[!is.na(all_auc_values)]
  
  mean_auc      <- mean(valid_auc)
  sd_auc        <- sd(valid_auc)
  
  thresholds    <- sapply(all_thresholds, function(x) ifelse(is.null(x), NA, x$threshold))
  valid_thresh  <- thresholds[!is.na(thresholds)]
  mean_threshold <- if (length(valid_thresh) > 0) mean(valid_thresh) else NA
  
  return(list(
    mean_auc      = mean_auc,
    sd_auc        = sd_auc,
    all_auc_values = valid_auc,
    mean_threshold = mean_threshold,
    all_thresholds = all_thresholds,
    fold_records  = fold_records
  ))
}

# --- Main Script ---

# 用來收集「每 repeat / 每 fold 的 AUC」
all_fold_results <- data.frame(
  strategy = character(),
  repeats  = integer(),
  fold     = integer(),
  AUC      = numeric(),
  stringsAsFactors = FALSE
)

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
if (length(unique(labels[!is.na(labels)])) < 2) {
  stop("Both positive (class=1) and negative (class=0) examples are required for evaluation.")
}
valid_idx <- !is.na(labels)
all_sequences <- all_sequences[valid_idx]
labels        <- labels[valid_idx]
cat("Using", sum(labels == 1), "positive and", sum(labels == 0), "negative sequences.\n")

# 3. Get list of PWM files (可選指定清單，否則掃整個 pwm_dir)
if (!is.null(wanted) && nzchar(wanted)) {
  pwm_files <- strsplit(wanted, ",")[[1]]
  pwm_files <- sapply(pwm_files, function(f) {
    if (!file.exists(f) && file.exists(file.path(pwm_dir, f))) {
      return(file.path(pwm_dir, f))
    } else {
      return(f)
    }
  })
} else {
  pwm_files <- list.files(pwm_dir, pattern = pwm_file_pattern, full.names = TRUE)
}
if (length(pwm_files) == 0) {
  stop("No PWM files found: ", ifelse(!is.null(wanted), wanted, pwm_dir))
}

# 4. Evaluate each PWM with cross-validation
cat("\nPerforming", k_folds, "fold cross-validation with", num_repeats, "repeats:\n\n")

results <- data.frame(
  PWM_File      = character(),
  Mean_AUC      = numeric(),
  SD_AUC        = numeric(),
  Mean_Threshold = numeric(),
  stringsAsFactors = FALSE
)

for (pwm_file in pwm_files) {
  file_basename <- basename(pwm_file)          # e.g. "pwm_cv.rds" 或 "null_cv_uniform.rds"
  cat("Evaluating", file_basename, "...\n")
  
  # Decide strategy_name: remove ".rds" 後的剩餘文字
  # 你可以自行定義要用 "cv" 或 "cv_uniform"；示例：
  if (grepl("^pwm_cv\\.rds$", file_basename)) {
    strategy_name <- "cv"
  } else if (grepl("^null_cv_uniform\\.rds$", file_basename)) {
    strategy_name <- "cv_uniform"
  } else {
    # 如果未命中上述兩者，就直接取檔名去除副檔名
    strategy_name <- sub("\\.rds$", "", file_basename)
  }
  
  # Load PWM
  pwm <- readRDS(pwm_file)
  
  # Run cross-validation, 並取得 per-fold 資料
  cv_out <- evaluate_pwm_with_cv(all_sequences, labels, pwm, k_folds, num_repeats, strategy_name)
  
  # 把 per-fold 的結果合併到全域 all_fold_results
  all_fold_results <- rbind(all_fold_results, cv_out$fold_records)
  
  # Print mean ± sd
  cat("  Mean AUC:", round(cv_out$mean_auc, 4),
      "±", round(cv_out$sd_auc, 4), "\n")
  cat("  Suggested threshold:", round(cv_out$mean_threshold, 4), "\n")
  
  # 加入 summary results
  results <- rbind(results, data.frame(
    PWM_File       = file_basename,
    Mean_AUC       = cv_out$mean_auc,
    SD_AUC         = cv_out$sd_auc,
    Mean_Threshold = cv_out$mean_threshold,
    stringsAsFactors = FALSE
  ))
}

# 5. Print summary table (mean ± sd for each model)
cat("\nSummary of Evaluation Results:\n")
print(results)

# 6. Write summary (Mean_AUC, SD_AUC) to file
results_file <- file.path(pwm_dir, "cv_evaluation_results.csv")
write.csv(results, file = results_file, row.names = FALSE)
cat("\nResults saved to:", results_file, "\n\n")

cat("Cross-validation evaluation complete!\n")

# 6b. 另外輸出 per-fold AUC
metrics_out <- file.path(pwm_dir, "metrics_cv.csv")
write.csv(all_fold_results, file = metrics_out, row.names = FALSE)
cat("Per-fold AUC (with strategy) saved to:", metrics_out, "\n")
