#!/usr/bin/env Rscript
# evaluate_models_with_groupk_uniform.R
# -------------------------------------------------------------
# 用同一套 GroupKFold 拆分 (data/splits/groupk/fold*)：
#   • 針對每個 fold，不從 train_sequences 建 PWM，
#     而是直接產生 Uniform PWM (4×L，值全 0.25)。
#   • 在同一個 test_sequences 上計算 AUC。
# 輸出每個 fold 的 AUC、以及平均 ± 標準差，存成 metrics_groupk_uniform.csv
# -------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(pROC)
})

option_list <- list(
  make_option("--fold_dir", type = "character", default = "data/splits/groupk",
              help = "GroupKFold 拆分的根目錄 (e.g. data/splits/groupk)"),
  make_option("--out_dir", type = "character", default = "results",
              help = "輸出結果資料夾 (metrics 檔會放這裡)")
)
opt <- parse_args(OptionParser(option_list = option_list))

fold_root <- opt$fold_dir
out_dir   <- opt$out_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 產生 uniform PWM 的函式 ----
# 傳入 motif 長度 L (每條 seq 的長度假設一致)，回傳 4×L 機率全 0.25 的矩陣
make_uniform_pwm <- function(L) {
  bases <- c("A", "C", "G", "T")
  pwm <- matrix(0.25, nrow = 4, ncol = L,
                dimnames = list(bases, paste0("Pos", seq_len(L))))
  return(pwm)
}

# ---- 打分函式 (同原本) ----
calculate_score <- function(seq, pwm) {
  if (nchar(seq) != ncol(pwm)) return(NA)
  score <- 0
  for (i in 1:nchar(seq)) {
    base <- substr(seq, i, i)
    if (base %in% rownames(pwm)) {
      prob <- pwm[base, i]
      # log2(p_pwm / p_background)，背景用 0.25
      score <- score + log2(prob / 0.25)
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

# ---- 收集每個 fold 的結果 ----
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
  cat("\n▶ Evaluating (Uniform) ", fold_name, "...\n", sep = "")

  train_file <- file.path(fold_path, "train_sequences.fasta")
  test_file  <- file.path(fold_path, "test_sequences.fasta")

  if (!file.exists(train_file) || !file.exists(test_file)) {
    warning("Missing files in ", fold_path)
    next
  }

  # 只讀 train_sequences 主要是為了知道 motif 長度 L
  train_seqs <- readDNAStringSet(train_file)
  if (length(train_seqs) == 0) {
    warning("Train set in ", fold_name, " 沒有序列，跳過")
    next
  }
  L <- width(train_seqs)[1]
  # 如果 train 的序列長度不一致，就取第一條當作 L
  # （假設所有都是相同長度）
  cat("  motif length (L) = ", L, "\n", sep = "")

  # 產生 Uniform PWM
  pwm_uniform <- make_uniform_pwm(L)

  # 讀 test set
  test_seqs  <- readDNAStringSet(test_file)
  test_labels <- as.numeric(sapply(names(test_seqs), get_label))
  test_scores <- sapply(as.character(test_seqs), calculate_score, pwm = pwm_uniform)

  keep <- !is.na(test_scores) & !is.na(test_labels)
  test_scores <- test_scores[keep]
  test_labels <- test_labels[keep]

  if (length(unique(test_labels)) < 2) {
    warning("Test set in ", fold_name, " 沒有正負樣本，跳過")
    next
  }

  auc_val <- as.numeric(auc(roc(test_labels, test_scores, quiet = TRUE)))
  cat("  AUC =", round(auc_val, 4), "\n")

  all_results <- rbind(all_results, data.frame(
    strategy = "groupk",
    repeats = 1,
    fold = fold_index,
    AUC = auc_val,
    stringsAsFactors = FALSE
  ))
  fold_index <- fold_index + 1
}

# 輸出 metrics_groupk_uniform.csv（格式與真實 PWM 版相同）
metrics_out <- file.path(out_dir, "metrics_groupk_uniform.csv")
write.csv(all_results, metrics_out, row.names = FALSE)

# 額外輸出平均與 SD（方便自己看）
summary_df <- data.frame(
  mean_auc = mean(all_results$AUC, na.rm = TRUE),
  sd_auc   = sd(all_results$AUC, na.rm = TRUE)
)
write.csv(summary_df, file.path(out_dir, "summary_groupk_uniform.csv"), row.names = FALSE)

cat("\n✅ Done (Uniform). AUC 平均:", round(summary_df$mean_auc, 4),
    "+/-", round(summary_df$sd_auc, 4), "\n")
cat("✔ metrics written to:", metrics_out, "\n")
