# prepare_datasets_full.R — generate a *complete* labeled dataset (no train/test split)
# ----------------------------------------------------------------------------------
# (…其餘註解同前，省略…)

suppressPackageStartupMessages({
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install via BiocManager::install('Biostrings')")
  }
  library(Biostrings)
})

# --------------------------- helper functions ------------------------------------
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
    names(neg_seq) <- paste0("negative_", i, " | class=0")
    neg <- c(neg, neg_seq)
  }
  neg
}

# --------------------------- argument parsing ------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("See header for usage. At minimum supply an output FASTA file.")
}

input_fasta      <- "data/extracted_sequences.fasta"
output_fasta     <- args[1]
neg_method       <- "shuffle"
neg_ratio        <- 1.0
add_negatives    <- TRUE
target_length    <- NULL   # default; set NULL to disable

if (length(args) > 1) {
  idx <- 2
  while (idx <= length(args)) {
    flag <- args[idx]
    if (flag == "--input")           { input_fasta   <- args[idx+1]; idx <- idx + 2 }
    else if (flag == "--target_len")  {
      val <- args[idx+1];
      if (tolower(val) %in% c("null", "none", "na", "0"))   target_length <- NULL
      else {
        target_length <- as.numeric(val)
        if (is.na(target_length)) stop("--target_len expects a number or 'NULL'.")
      }
      idx <- idx + 2
    }
    else if (flag == "--neg_method")  { neg_method    <- args[idx+1]; idx <- idx + 2 }
    else if (flag == "--neg_ratio")   { neg_ratio     <- as.numeric(args[idx+1]); idx <- idx + 2 }
    else if (flag == "--no_neg")      { add_negatives <- FALSE; idx <- idx + 1 }
    else { stop("Unknown option: ", flag) }
  }
}

cat("[prepare_datasets_full] Reading positives from", input_fasta, "...\n")
if (!file.exists(input_fasta)) stop("Input FASTA not found: ", input_fasta)
pos <- readDNAStringSet(input_fasta)

if (!all(grepl("class=", names(pos)))) {
  names(pos) <- paste0(names(pos), " | class=1")
}

if (!is.null(target_length)) {
  cat("  • Length filter:", target_length, "bp\n")
  kept <- width(pos) == target_length
  cat("    keeping", sum(kept), "of", length(pos), "sequences.\n")
  pos <- pos[kept]
  if (!length(pos)) stop("No sequences remain after length filtering.")
} else {
  cat("  • Length filter: **disabled** (keeping all", length(pos), "sequences)\n")
}

neg <- DNAStringSet()
if (add_negatives) {
  cat("[prepare_datasets_full] Generating", neg_ratio, "× negatives using", neg_method, "...\n")
  neg <- generate_negatives(pos, method = neg_method, ratio = neg_ratio)
}

all_seqs <- c(pos, neg)
cat("[prepare_datasets_full] Writing", length(all_seqs), "sequences to", output_fasta, "...\n")
writeXStringSet(all_seqs, filepath = output_fasta)
cat("[prepare_datasets_full] Done.\n")
