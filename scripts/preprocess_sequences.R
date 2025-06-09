# filepath: /mnt/d/workspace/data-science/scripts/preprocess_sequences.R
# R script for DNA sequence preprocessing: length filtering, N-base handling, and masking of low-complexity regions

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
library(Biostrings)

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

# Handle case where output is a directory
if (dir.exists(output_fasta_file) || endsWith(output_fasta_file, "/") || endsWith(output_fasta_file, "\\")) {
  # If output is a directory, create a default filename
  output_fasta_file <- file.path(output_fasta_file, "preprocessed_sequences.fasta")
  cat("Output path is a directory. Using default filename:", output_fasta_file, "\n")
}

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
  
  # Read JSON config file - enhanced parsing
  config_text <- readLines(config_file, warn = FALSE)
  config_text <- paste(config_text, collapse = "")
  
  # Parse key parameters with improved regex
  if (grepl("\"min_length\"[^0-9]*(\\d+)", config_text)) {
    params$min_length <- as.numeric(gsub(".*\"min_length\"[^0-9]*(\\d+).*", "\\1", config_text))
  }
  if (grepl("\"max_length\"[^0-9]*(\\d+)", config_text)) {
    params$max_length <- as.numeric(gsub(".*\"max_length\"[^0-9]*(\\d+).*", "\\1", config_text))
  }
  if (grepl("\"max_n_percent\"[^0-9]*(\\d+)", config_text)) {
    params$max_n_percent <- as.numeric(gsub(".*\"max_n_percent\"[^0-9]*(\\d+).*", "\\1", config_text))
  }
  if (grepl("\"low_complexity_filter\"[^:]*:[^a-z]*(true|false)", config_text)) {
    filter_val <- gsub(".*\"low_complexity_filter\"[^:]*:[^a-z]*(true|false).*", "\\1", config_text)
    params$low_complexity_filter <- (filter_val == "true")
  }
  if (grepl("\"mask_repeats\"[^:]*:[^a-z]*(true|false)", config_text)) {
    mask_val <- gsub(".*\"mask_repeats\"[^:]*:[^a-z]*(true|false).*", "\\1", config_text)
    params$mask_repeats <- (mask_val == "true")
  }
  
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

# --- Main Processing Function ---
preprocess_sequence <- function(sequence, seq_name, params) {
  seq_string <- as.character(sequence)
  original_length <- nchar(seq_string)
  modifications <- c()
  
  # Force uppercase if requested
  if (params$force_uppercase) {
    seq_string <- toupper(seq_string)
    if (seq_string != as.character(sequence)) {
      modifications <- c(modifications, "uppercase")
    }
  }
  
  # Length filtering
  if (params$length_filter) {
    if (!is.null(params$target_length)) {
      # For exact length extraction (e.g., for PWM)
      if (original_length != params$target_length) {
        if (original_length > params$target_length) {
          # If sequence is longer, take the center portion
          start_pos <- floor((original_length - params$target_length) / 2) + 1
          seq_string <- substr(seq_string, start_pos, start_pos + params$target_length - 1)
          modifications <- c(modifications, paste0("trimmed_to_", params$target_length))
        } else {
          # Sequence too short, can't meet target length
          return(NULL)
        }
      }
    } else {
      # Apply min/max length filters
      if (original_length < params$min_length || original_length > params$max_length) {
        return(NULL)  # Sequence doesn't meet length criteria
      }
    }
  }
  
  # N-base handling
  n_count <- lengths(gregexpr("N", seq_string))
  n_percent <- (n_count / nchar(seq_string)) * 100
  
  if (n_percent > params$max_n_percent) {
    return(NULL)  # Too many Ns, reject sequence
  }
  
  if (n_count > 0 && params$n_handling == "mask") {
    old_seq <- seq_string
    seq_string <- gsub("N", params$n_replacement, seq_string)
    if (old_seq != seq_string) {
      modifications <- c(modifications, paste0("masked_", n_count, "_Ns"))
    }
  } else if (n_count > 0 && params$n_handling == "remove") {
    return(NULL)  # Any N presence causes rejection
  }
  
  # Low-complexity filtering
  if (params$low_complexity_filter) {
    entropy <- calculate_entropy(seq_string, params$complexity_window)
    
    if (entropy < params$entropy_threshold) {
      if (params$complexity_handling == "remove") {
        return(NULL)  # Low complexity, reject sequence
      } else if (params$complexity_handling == "mask") {
        # For simplicity, we'll just mark it in the header
        modifications <- c(modifications, paste0("low_complexity_", round(entropy, 2)))
      }
    }
  }
  
  # Repeat masking
  if (params$mask_repeats) {
    old_seq <- seq_string
    masked_seq <- mask_repeats(DNAString(seq_string), params$repeat_min_length)
    seq_string <- as.character(masked_seq)
    
    if (old_seq != seq_string) {
      n_masked <- lengths(gregexpr("N", seq_string)) - n_count
      if (n_masked > 0) {
        modifications <- c(modifications, paste0("masked_", n_masked, "_repeat_bases"))
      }
    }
  }
  
  # Create new sequence object
  new_seq <- DNAString(seq_string)
  
  # Update sequence name if modifications were made
  if (params$label_preprocessed && length(modifications) > 0) {
    new_name <- paste0(seq_name, " | preprocessed=", paste(modifications, collapse=","))
  } else {
    new_name <- seq_name
  }
  
  # Return preprocessed sequence with new name as DNAStringSet
  result <- DNAStringSet(new_seq)
  names(result) <- new_name
  return(result)
}

# --- Main Script ---
cat("Starting sequence preprocessing...\n")

# 1. Validate inputs
if (!file.exists(input_fasta_file)) {
  stop("Input FASTA file not found: ", input_fasta_file)
}

# 2. Read input sequences
cat("Reading sequences from:", input_fasta_file, "\n")
input_sequences <- readDNAStringSet(input_fasta_file)
cat("Read", length(input_sequences), "sequences.\n")

# 3. Process each sequence
cat("Applying preprocessing with parameters:\n")
print(params)

processed_sequences <- DNAStringSet()
filtered_count <- 0

for (i in 1:length(input_sequences)) {
  seq_name <- names(input_sequences)[i]
  sequence <- input_sequences[[i]]
  
  processed_seq <- preprocess_sequence(sequence, seq_name, params)
  
  if (!is.null(processed_seq)) {
    processed_sequences <- c(processed_sequences, processed_seq)
  } else {
    filtered_count <- filtered_count + 1
    if (params$verbose) {
      cat("  Filtered sequence:", seq_name, "\n")
    }
  }
  
  # Optional progress reporting for large datasets
  if (params$verbose && i %% 1000 == 0) {
    cat("  Processed", i, "of", length(input_sequences), "sequences...\n")
  }
}

# 4. Write results
cat("Writing", length(processed_sequences), "preprocessed sequences to:", output_fasta_file, "\n")
cat("Filtered out", filtered_count, "sequences.\n")
writeXStringSet(processed_sequences, filepath = output_fasta_file)

cat("Preprocessing complete.\n")