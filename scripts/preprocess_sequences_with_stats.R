# filepath: /mnt/d/workspace/data-science/scripts/preprocess_sequences.R
# R script for DNA sequence preprocessing: length filtering, N-base handling, and masking of low-complexity regions

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
library(Biostrings)
# s1 <- DNAString("ACGT")            # 單一 DNAString
# s2 <- DNAStringSet("ACGT")         # DNAStringSet 包裝

# names(s1) <- "wrong"               # ❌ 錯誤，會報錯
# names(s2) <- "correct"             # ✅ 正確

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript preprocess_sequences.R <input_fasta> <output_fasta> [config_file]
       \n\t<input_fasta>: Input FASTA file with sequences to preprocess
       \n\t<output_fasta>: Output file to save preprocessed sequences
       \n\t[config_file]: Optional JSON file with preprocessing parameters", 
       call. = FALSE)
}

input_fasta_file <- args[1]
output_fasta_file <- args[2]
config_file <- if (length(args) >= 3) args[3] else NULL

# --- Default Parameters ---
params <- list(
  # Length filtering
  length_filter = TRUE,
  min_length = 11,
  max_length = 100,
  target_length = NULL,  # If set, extract sequences of this exact length
  
  # N-base handling
  n_handling = "mask",   # Options: "mask", "remove", "keep"
  max_n_percent = 15,    # Maximum percentage of N bases allowed (sequences above this are removed)
  n_replacement = "A",   # Character to replace N with if n_handling is "mask"
  
  # Low-complexity filtering
  low_complexity_filter = TRUE,
  entropy_threshold = 1.5,  # Sequences with entropy below this are removed/masked
  complexity_handling = "remove",  # Options: "mask", "remove"
  complexity_window = 10,  # Window size for entropy calculation
  
  # Repeat masking
  mask_repeats = TRUE,
  repeat_min_length = 5,  # Minimum length of repeat to mask
  
  # Other options
  force_uppercase = TRUE,
  label_preprocessed = TRUE,  # Add preprocessing info to sequence headers
  verbose = TRUE
)

# --- Load Configuration from File (if provided) ---
if (!is.null(config_file)) {
  if (!file.exists(config_file)) {
    stop("Config file not found: ", config_file)
  }
  
  # Read JSON config file (can use jsonlite or rjson package)
  # For simplicity in this example, we'll just use a basic approach
  config_text <- readLines(config_file, warn = FALSE)
  config_text <- paste(config_text, collapse = "")
  
  # Very simple parsing - in a real application, use a proper JSON parser
  if (grepl("\"min_length\"[^0-9]*(\\d+)", config_text)) {
    params$min_length <- as.numeric(gsub(".*\"min_length\"[^0-9]*(\\d+).*", "\\1", config_text))
  }
  # Extract other parameters similarly (simplified for the example)
  
  if (params$verbose) {
    cat("Loaded configuration from:", config_file, "\n")
  }
}

# --- Utility Functions ---

# Calculate sequence entropy (measure of complexity)
calculate_entropy <- function(sequence, window_size) {
  seq_string <- as.character(sequence)
  seq_length <- nchar(seq_string)
  
  if (seq_length < window_size) {
    # For short sequences, calculate on the whole sequence
    window_size <- seq_length
  }
  
  # Prepare for sliding window entropy calculation
  entropies <- numeric(seq_length - window_size + 1)
  
  for (i in 1:(seq_length - window_size + 1)) {
    window <- substr(seq_string, i, i + window_size - 1)
    
    # Count frequency of each base in window
    base_counts <- table(strsplit(window, "")[[1]])
    base_probs <- base_counts / window_size
    
    # Calculate Shannon entropy: -Σ(p_i * log2(p_i))
    entropy <- -sum(base_probs * log2(base_probs))
    entropies[i] <- entropy
  }
  
  # Return minimum entropy found in any window (most low-complexity region)
  return(min(entropies))
}

# Detect and mask simple repeats
mask_repeats <- function(sequence, min_repeat_length) {
  seq_string <- as.character(sequence)
  seq_length <- nchar(seq_string)
  
  # Simple detection of homopolymer runs (e.g., AAAAA, TTTTT)
  for (base in c("A", "C", "G", "T")) {
    repeat_pattern <- paste0(rep(base, min_repeat_length), collapse = "")
    seq_string <- gsub(repeat_pattern, paste0(rep("N", min_repeat_length), collapse = ""), seq_string)
  }
  
  # Simple detection of dinucleotide repeats (e.g., ATATATATAT)
  dinucs <- c("AT", "TA", "GC", "CG", "AC", "CA", "GT", "TG", "AG", "GA", "CT", "TC")
  for (dinuc in dinucs) {
    repeat_count <- floor(min_repeat_length / 2)
    if (repeat_count >= 2) {  # At least 2 repetitions
      repeat_pattern <- paste0(rep(dinuc, repeat_count), collapse = "")
      replacement <- paste0(rep("N", nchar(repeat_pattern)), collapse = "")
      seq_string <- gsub(repeat_pattern, replacement, seq_string)
    }
  }
  
  return(DNAString(seq_string))
}


# --- Main Processing Function with reason tracking ---
preprocess_sequence <- function(sequence, seq_name, params) {
  seq_string <- as.character(sequence)
  original_length <- nchar(seq_string)
  modifications <- c()
  
  # Force uppercase
  if (params$force_uppercase) {
    seq_string <- toupper(seq_string)
  }

  # Length filter
  # Length filtering
  if (params$length_filter) {
  # 如果有指定 target_length，就裁切（但只允許長度足夠的）
    if (!is.null(params$target_length)) {
      if (original_length >= params$target_length) {
        start_pos <- floor((original_length - params$target_length) / 2) + 1
        seq_string <- substr(seq_string, start_pos, start_pos + params$target_length - 1)
        modifications <- c(modifications, paste0("trimmed_to_", params$target_length))
      } else {
        return(list(seq = NULL, reason = "length"))
      }
    }

  # 不論是否裁切，都要檢查長度是否符合 min_length 與 max_length 範圍
    final_length <- nchar(seq_string)
    if (final_length < params$min_length || final_length > params$max_length) {
      return(list(seq = NULL, reason = "length"))
    }
  }


  # N content filter
  n_count <- lengths(gregexpr("N", seq_string))
  n_percent <- (n_count / nchar(seq_string)) * 100

  if (n_percent > params$max_n_percent) {
    return(list(seq = NULL, reason = "n_content"))
  }

  if (n_count > 0 && params$n_handling == "mask") {
    seq_string <- gsub("N", params$n_replacement, seq_string)
    modifications <- c(modifications, paste0("masked_", n_count, "_Ns"))
  } else if (n_count > 0 && params$n_handling == "remove") {
    return(list(seq = NULL, reason = "n_content"))
  }

  # Entropy filter
  if (params$low_complexity_filter) {
    entropy <- calculate_entropy(seq_string, params$complexity_window)
    if (entropy < params$entropy_threshold) {
      if (params$complexity_handling == "remove") {
        return(list(seq = NULL, reason = "entropy"))
      } else if (params$complexity_handling == "mask") {
        modifications <- c(modifications, paste0("low_entropy_", round(entropy, 2)))
      }
    }
  }

  # Repeat masking
  if (params$mask_repeats) {
    masked_seq <- mask_repeats(DNAString(seq_string), params$repeat_min_length)
    seq_string <- as.character(masked_seq)
  }

  # Final output
  new_seq <- DNAStringSet(DNAString(seq_string))
  if (params$label_preprocessed && length(modifications) > 0) {
    names(new_seq) <- paste0(seq_name, " | preprocessed=", paste(modifications, collapse = ","))
  } else {
    names(new_seq) <- seq_name
  }
  return(list(seq = new_seq, reason = NULL))
}

# --- Main Execution ---
cat("Starting sequence preprocessing...\n")

if (!file.exists(input_fasta_file)) stop("Input FASTA file not found: ", input_fasta_file)

input_sequences <- readDNAStringSet(input_fasta_file)
# --- Debug: 顯示序列長度統計分布 ---
seq_lengths <- width(input_sequences)
cat("\n📊 原始序列長度摘要：\n")
print(summary(seq_lengths))
cat("⏳ 長度 ≤ 50bp：", sum(seq_lengths <= 50), "\n")
cat("⏳ 51 ~ 100bp：", sum(seq_lengths > 50 & seq_lengths < 101), "\n")
cat("✅ ≥ 101bp：", sum(seq_lengths >= 101), "\n\n")

cat("Read", length(input_sequences), "sequences.\n")
print(params)

processed_sequences <- DNAStringSet()
filter_stats <- list(length = 0, n_content = 0, entropy = 0, passed = 0)

for (i in seq_along(input_sequences)) {
  seq_name <- names(input_sequences)[i]
  sequence <- input_sequences[[i]]

  result <- preprocess_sequence(sequence, seq_name, params)

  if (!is.null(result$seq)) {
    processed_sequences <- c(processed_sequences, result$seq)
    filter_stats$passed <- filter_stats$passed + 1
  } else {
    filter_stats[[result$reason]] <- filter_stats[[result$reason]] + 1
    if (params$verbose) {
      cat("  Filtered:", seq_name, "| Reason:", result$reason, "\n")
    }
  }

  if (params$verbose && i %% 1000 == 0) {
    cat("  Processed", i, "/", length(input_sequences), "...\n")
  }
}

# Write results
writeXStringSet(processed_sequences, filepath = output_fasta_file)
cat("Preprocessing complete.\n")

# Output stats
cat("Filtering summary:\n")
cat("  Length filtered:  ", filter_stats$length, "\n")
cat("  N filtered:       ", filter_stats$n_content, "\n")
cat("  Entropy filtered: ", filter_stats$entropy, "\n")
cat("  Passed:           ", filter_stats$passed, "\n")
