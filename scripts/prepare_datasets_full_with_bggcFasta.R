#!/usr/bin/env Rscript
# prepare_datasets_full.R — 產生「完整」標記 FASTA（不做 train/test split）
#   另外固定輸出 GC-matched 背景檔 <basename>_bg_gc.fasta
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Package 'Biostrings' is required. Install via BiocManager::install('Biostrings')")
  library(Biostrings)
  library(tools)        # 取副檔名 / 去副檔名用
})

# --------------------------- helper functions -----------------------------------
shuffle_sequence <- function(seq) {
  paste0(sample(strsplit(as.character(seq), "")[[1]]), collapse = "")
}

dinuc_shuffle <- function(seq) {
  seq_str <- as.character(seq)
  if (nchar(seq_str) <= 2) return(seq_str)
  dinucs <- substring(seq_str, 1:(nchar(seq_str) - 1), 2:nchar(seq_str))
  first_nuc <- substring(seq_str, 1, 1)
  result <- first_nuc
  pool <- dinucs
  current <- first_nuc
  while (length(pool) > 0) {
    idx <- grep(paste0("^", current), pool)
    chosen <- if (length(idx)) sample(idx, 1) else sample(seq_along(pool), 1)
    dinuc <- pool[chosen]
    pool <- pool[-chosen]
    next_nuc <- substring(dinuc, 2, 2)
    result <- paste0(result, next_nuc)
    current <- next_nuc
  }
  result
}

generate_random_sequence <- function(len) {
  paste0(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
}

generate_negatives <- function(pos_seqs, method = "shuffle", ratio = 1) {
  n_pos <- length(pos_seqs)
  n_neg <- ceiling(n_pos * ratio)
  neg <- DNAStringSet()
  for (i in seq_len(n_neg)) {
    pos_seq <- pos_seqs[[ ((i - 1) %% n_pos) + 1 ]]
    new_seq_str <- switch(tolower(method),
                          shuffle = shuffle_sequence(pos_seq),
                          dinuc   = dinuc_shuffle(pos_seq),
                          random  = generate_random_sequence(width(pos_seq)),
                          stop("Unknown neg_method: ", method))
    neg_seq <- DNAStringSet(new_seq_str)
    names(neg_seq) <- paste0("negative_", i, " | class=0 | type=", tolower(method))
    neg <- c(neg, neg_seq)
  }
  neg
}

# ----- GC-matched 背景：每條正例對應一條長度相同、GC 比例相近 (±1bp 誤差) -----
generate_gc_matched_negatives <- function(pos_seqs, seed = 123) {
  set.seed(seed)
  pos_gc <- rowSums(letterFrequency(pos_seqs, c("G", "C"), as.prob = TRUE))  # GC 比例
  gc_neg <- DNAStringSet()
  for (i in seq_along(pos_seqs)) {
    len <- width(pos_seqs[i])
    n_gc <- round(len * pos_gc[i])
    n_at <- len - n_gc
    vec  <- c(rep("G", floor(n_gc / 2)), rep("C", ceiling(n_gc / 2)),
              rep("A", floor(n_at / 2)), rep("T", ceiling(n_at / 2)))
    new  <- DNAStringSet(paste(sample(vec), collapse = ""))
    names(new) <- paste0("gc_bg_", i, " | class=0 | type=gc")
    gc_neg <- c(gc_neg, new)
  }
  gc_neg
}

# --------------------------- argument parsing -----------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法：Rscript prepare_datasets_full.R <output_fasta> [options]\n",
       "  --input <path>       ：正例 FASTA (預設 data/extracted_sequences.fasta)\n",
       "  --target_len <n|NULL>：長度過濾；NULL 代表關閉 (預設 82)\n",
       "  --neg_method <m>     ：shuffle / dinuc / random (預設 shuffle)\n",
       "  --neg_ratio <r>      ：負例/正例 比例 (預設 1)\n",
       "  --no_neg             ：不產生 shuffle/dinuc/random 負例\n",
       "  --gc_bg_out <path>   ：覆寫 GC 背景檔輸出路徑\n")
}

input_fasta  <- "data/extracted_sequences.fasta"
output_fasta <- args[1]

# 根據 output_fasta 自動組成 <basename>_bg_gc.<ext>
out_dir       <- dirname(output_fasta)
base_no_ext   <- file_path_sans_ext(basename(output_fasta))
ext           <- file_ext(output_fasta)
ext           <- if (nzchar(ext)) paste0(".", ext) else ""
gc_bg_out     <- file.path(out_dir, paste0(base_no_ext, "_bg_gc", ext))

neg_method    <- "shuffle"
neg_ratio     <- 1.0
add_negatives <- TRUE
target_length <- 82  # NULL 表示不過濾

# 解析額外旗標
if (length(args) > 1) {
  idx <- 2
  while (idx <= length(args)) {
    flag <- args[idx]
    if (flag == "--input")           { input_fasta  <- args[idx+1]; idx <- idx + 2 }
    else if (flag == "--target_len") {
      val <- args[idx+1]
      if (tolower(val) %in% c("null", "none", "na", "0")) target_length <- NULL
      else {
        target_length <- as.numeric(val)
        if (is.na(target_length)) stop("--target_len expects a number or 'NULL'.")
      }
      idx <- idx + 2
    }
    else if (flag == "--neg_method") { neg_method   <- args[idx+1]; idx <- idx + 2 }
    else if (flag == "--neg_ratio")  { neg_ratio    <- as.numeric(args[idx+1]); idx <- idx + 2 }
    else if (flag == "--no_neg")     { add_negatives <- FALSE; idx <- idx + 1 }
    else if (flag == "--gc_bg_out")  { gc_bg_out    <- args[idx+1]; idx <- idx + 2 }
    else stop("Unknown option: ", flag)
  }
}

# --------------------------- 主流程 ---------------------------------------------
cat("[prepare_datasets_full] Reading positives from", input_fasta, "...\n")
if (!file.exists(input_fasta)) stop("Input FASTA not found: ", input_fasta)
pos <- readDNAStringSet(input_fasta)

# 標上 class=1
if (!all(grepl("class=", names(pos)))) {
  names(pos) <- paste0(names(pos), " | class=1")
}

# 長度過濾
if (!is.null(target_length)) {
  cat("  • Length filter:", target_length, "bp\n")
  kept <- width(pos) == target_length
  cat("    keeping", sum(kept), "of", length(pos), "sequences.\n")
  pos <- pos[kept]
  if (!length(pos)) stop("No sequences remain after length filtering.")
} else {
  cat("  • Length filter: **disabled** (keeping all", length(pos), "sequences)\n")
}

# ----- 產生 shuffle/dinuc/random 負例（若有） -----
neg <- DNAStringSet()
if (add_negatives) {
  cat("[prepare_datasets_full] Generating", neg_ratio, "× negatives using", neg_method, "...\n")
  neg <- generate_negatives(pos, method = neg_method, ratio = neg_ratio)
}

# ----- 永久產生 GC-matched 背景，分檔輸出 -----
cat("[prepare_datasets_full] Generating GC-matched negatives …\n")
gc_neg <- generate_gc_matched_negatives(pos)
cat("  • GC negatives:", length(gc_neg), "\n")
cat("  • Writing GC negatives to", gc_bg_out, "\n")
writeXStringSet(gc_neg, filepath = gc_bg_out)

# ----- 主 FASTA（正例 + 非 GC 背景） -----
all_seqs <- c(pos, neg)   # 不含 gc_neg
cat("[prepare_datasets_full] Writing", length(all_seqs), "sequences to", output_fasta, "...\n")
writeXStringSet(all_seqs, filepath = output_fasta)

cat("[prepare_datasets_full] Done.\n")
