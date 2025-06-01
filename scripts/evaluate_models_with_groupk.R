# evaluate_models_with_groupk.R
# -------------------------------------------------------------
# 對 data/splits/groupk/fold*/ 中的每 fold：
#   - 用 train_sequences.fasta 建 PWM
#   - 用 test_sequences.fasta 評估 AUC
# 輸出每 fold 的 AUC 與平均 ± SD，格式對齊 CV 評估版本
# -------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(pROC)
})

option_list <- list(
  make_option("--fold_dir", type = "character", default = "data/splits/groupk",
              help = "根資料夾，包含 fold1 到 foldK"),
  make_option("--out_dir", type = "character", default = "results",
              help = "輸出結果的資料夾")
)
opt <- parse_args(OptionParser(option_list = option_list))

fold_root <- opt$fold_dir
out_dir   <- opt$out_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

calculate_pwm <- function(seqs) {
  cm <- consensusMatrix(seqs, as.prob = FALSE)
  cm <- cm[rownames(cm) %in% c("A", "C", "G", "T"), , drop = FALSE]
  pwm <- prop.table(cm + 1, margin = 2)  # 加 pseudocount
  return(pwm)
}

calculate_score <- function(seq, pwm) {
  if (nchar(seq) != ncol(pwm)) return(NA)  # ⚠️ 加入長度檢查
  score <- 0
  for (i in 1:nchar(seq)) {
    base <- substr(seq, i, i)
    if (base %in% rownames(pwm)) {
      prob <- pwm[base, i]
      score <- score + log2(prob / 0.25)  # log-likelihood ratio
    } else {
      score <- score - 10
    }
  }
  return(score)
}

get_label <- function(name) {
  m <- regexec("class=([01])", name)
  regmatches(name, m)[[1]][2]
}

# --- 收集每 fold 的 AUC 結果（統一欄位格式）---
all_results <- data.frame(
  strategy = character(),
  repeats = integer(),
  fold = integer(),
  AUC = numeric(),
  stringsAsFactors = FALSE
)

fold_dirs <- list.dirs(fold_root, full.names = TRUE, recursive = FALSE)
fold_index <- 1

for (fold_path in sort(fold_dirs)) {
  fold_name <- basename(fold_path)
  cat("\n▶ Evaluating", fold_name, "...\n")

  train_file <- file.path(fold_path, "train_sequences.fasta")
  test_file  <- file.path(fold_path, "test_sequences.fasta")

  if (!file.exists(train_file) || !file.exists(test_file)) {
    warning("Missing files in", fold_path)
    next
  }

  train_seqs <- readDNAStringSet(train_file)
  test_seqs  <- readDNAStringSet(test_file)

  pwm <- calculate_pwm(train_seqs)

  test_labels <- as.numeric(sapply(names(test_seqs), get_label))
  test_scores <- sapply(as.character(test_seqs), calculate_score, pwm = pwm)

  keep <- !is.na(test_scores) & !is.na(test_labels)
  test_scores <- test_scores[keep]
  test_labels <- test_labels[keep]

  if (length(unique(test_labels)) < 2) {
    warning("Test set in", fold_name, "沒有正負樣本，跳過")
    next
  }

  auc_val <- as.numeric(auc(roc(test_labels, test_scores, quiet = TRUE)))

  cat("  AUC =", round(auc_val, 4), "\n")

  all_results <- rbind(all_results, data.frame(
    strategy = "groupk",
    repeats = 1,
    fold = fold_index,
    AUC = auc_val
  ))

  fold_index <- fold_index + 1
}

# 輸出 metrics_groupk.csv（格式與 CV 一致）
metrics_out <- file.path(out_dir, "metrics_groupk.csv")
write.csv(all_results, metrics_out, row.names = FALSE)

# 額外輸出 summary（非 visualize_metrics 所用，但方便人工參考）
summary_df <- data.frame(
  mean_auc = mean(all_results$AUC, na.rm = TRUE),
  sd_auc   = sd(all_results$AUC, na.rm = TRUE)
)
write.csv(summary_df, file.path(out_dir, "summary_groupk.csv"), row.names = FALSE)

cat("\n✅ Done. AUC 平均值:", round(summary_df$mean_auc, 4), "+/-", round(summary_df$sd_auc, 4), "\n")
cat("✔ metrics written to:", metrics_out, "\n")
