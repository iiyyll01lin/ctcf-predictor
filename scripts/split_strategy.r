# split_strategy.R — create train/test (or CV folds) FASTA files under data/splits/
# --------------------------------------------------------------------------------
# 支援三種策略：
#   • random  —— 隨機抽取一定比例做 test（預設 80/20）
#   • chrom   —— 依染色體劃分（使用者指定哪些 chr 屬於 train/test）
#   • groupk  —— 以染色體為 group 做 K-fold CV，輸出 fold{i}/{train|test}.fasta
#
# 執行範例
# ----------
# 1) 隨機切分 80% train / 20% test
#    Rscript split_strategy.R --strategy random \
#            --input data/complete_labeled.fasta --out_dir data/splits/random \
#            --train_ratio 0.8 --seed 42
#
# 2) 染色體切分：chr1-6 為 train，chr7 為 val，chr8-9 為 test
#    Rscript split_strategy.R --strategy chrom \
#            --train_chr chr1,chr2,chr3,chr4,chr5,chr6 \
#            --test_chr  chr8,chr9  --val_chr chr7 \
#            --out_dir data/splits/chrom
#    （val_chr 可省略；若同時指定 test_chr/val_chr 以外的 chr，會自動歸入 train）
#
# 3) Group K‑Fold（以染色體為 group）
#    Rscript split_strategy.R --strategy groupk \
#            --k 5 --seed 123 \
#            --out_dir data/splits/groupk
# --------------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install via BiocManager::install('Biostrings')")
  }
  library(Biostrings)
})

# --------------------------- util functions --------------------------------------
msg <- function(...) cat("[split_strategy]", ..., "\n")

# 解析 chr 名稱；假設 header 形如 "chr1:123-456 | class=1"
extract_chr <- function(nm) {
  sub("^([a-zA-Z0-9_]+):.*$", "\\1", nm)   # 抓 ':' 前字串
}

write_set <- function(seqs, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeXStringSet(seqs, filepath = path)
  msg("  ↳", length(seqs), "seqs →", path)
}

# --------------------------- argument parsing ------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default settings
opt <- list(
  input        = "data/complete_labeled.fasta",
  out_dir      = "data/splits/random",
  strategy     = "random",      # random / chrom / groupk
  train_ratio  = 0.8,            # for random
  seed         = 42,
  train_chr    = NULL,           # for chrom
  test_chr     = NULL,
  val_chr      = NULL,
  k            = 5               # for groupk
)

# quick getopt‑like parse
if (length(args)) {
  i <- 1
  while (i <= length(args)) {
    flag <- args[i];
    val  <- if (i + 1 <= length(args)) args[i+1] else NA
    switch(flag,
      "--input"        = { opt$input   <- val; i <- i + 2 },
      "--out_dir"      = { opt$out_dir <- val; i <- i + 2 },
      "--strategy"     = { opt$strategy <- tolower(val); i <- i + 2 },
      "--train_ratio"  = { opt$train_ratio <- as.numeric(val); i <- i + 2 },
      "--seed"         = { opt$seed  <- as.integer(val); i <- i + 2 },
      "--train_chr"    = { opt$train_chr <- strsplit(val, ",")[[1]]; i <- i + 2 },
      "--test_chr"     = { opt$test_chr  <- strsplit(val, ",")[[1]]; i <- i + 2 },
      "--val_chr"      = { opt$val_chr   <- strsplit(val, ",")[[1]]; i <- i + 2 },
      "--k"            = { opt$k <- as.integer(val); i <- i + 2 },
      stop("Unknown flag: ", flag)
    )
  }
}

# --------------------------- load data -------------------------------------------
msg("Reading", opt$input, "...")
if (!file.exists(opt$input)) stop("Input FASTA not found.")
all_seqs <- readDNAStringSet(opt$input)
chr_vec  <- extract_chr(names(all_seqs))

# --------------------------- strategy handlers -----------------------------------
set.seed(opt$seed)
strategy <- opt$strategy

if (strategy == "random") {
  msg("Strategy: random (train_ratio", opt$train_ratio, ")")
  idx <- sample(seq_along(all_seqs))
  n_train <- floor(length(idx) * opt$train_ratio)
  train_seqs <- all_seqs[idx[1:n_train]]
  test_seqs  <- all_seqs[idx[(n_train + 1):length(idx)]]
  write_set(train_seqs, file.path(opt$out_dir, "train_sequences.fasta"))
  write_set(test_seqs,  file.path(opt$out_dir, "test_sequences.fasta"))

} else if (strategy == "chrom") {
  msg("Strategy: chrom")
  # 如果未指定 train_chr，預設把 test_chr / val_chr 以外的全部歸 train
  uniq_chr <- unique(chr_vec)
  if (is.null(opt$test_chr)) stop("--test_chr must be specified for chrom strategy")
  train_chr <- if (!is.null(opt$train_chr)) opt$train_chr else setdiff(uniq_chr, c(opt$test_chr, opt$val_chr))
  val_chr   <- opt$val_chr
  msg("  train_chr:", paste(train_chr, collapse=","))
  if (length(val_chr)) msg("  val_chr:",   paste(val_chr, collapse=","))
  msg("  test_chr:",  paste(opt$test_chr, collapse=","))

  train_seqs <- all_seqs[chr_vec %in% train_chr]
  test_seqs  <- all_seqs[chr_vec %in% opt$test_chr]
  if (length(val_chr)) {
    val_seqs <- all_seqs[chr_vec %in% val_chr]
    write_set(val_seqs,   file.path(opt$out_dir, "val_sequences.fasta"))
  }
  write_set(train_seqs, file.path(opt$out_dir, "train_sequences.fasta"))
  write_set(test_seqs,  file.path(opt$out_dir, "test_sequences.fasta"))

} else if (strategy == "groupk") {
  msg("Strategy: Group K‑Fold (k =", opt$k, ")")
  uniq_chr <- sample(unique(chr_vec))  # shuffle chr order for reproducibility
  folds <- split(uniq_chr, rep_len(seq_len(opt$k), length(uniq_chr)))
  for (fold_id in seq_len(opt$k)) {
    test_chr <- folds[[fold_id]]
    train_chr <- setdiff(uniq_chr, test_chr)
    msg(" Fold", fold_id, "→ test_chr:", paste(test_chr, collapse=","))
    train_seqs <- all_seqs[chr_vec %in% train_chr]
    test_seqs  <- all_seqs[chr_vec %in% test_chr]
    subdir <- file.path(opt$out_dir, paste0("fold", fold_id))
    write_set(train_seqs, file.path(subdir, "train_sequences.fasta"))
    write_set(test_seqs,  file.path(subdir, "test_sequences.fasta"))
  }

} else {
  stop("Unknown strategy: ", strategy)
}

msg("Done.")
