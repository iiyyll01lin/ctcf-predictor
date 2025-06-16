# R script to prepare training and testing datasets from extracted sequences

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager: \n",
       "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')",
       call. = FALSE)
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/preprocessed_sequences_optimized.fasta"
method <- if (length(args) >= 2) args[2] else "chromosome_split"
output_dir <- if (length(args) >= 3) args[3] else "data/datasets/"
train_proportion <- if (length(args) >= 4) as.numeric(args[4]) else 0.8
k_folds <- if (length(args) >= 5) as.numeric(args[5]) else 5

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Dataset Preparation ===\n")
cat("Input file:", input_file, "\n")
cat("Method:", method, "\n")
cat("Output directory:", output_dir, "\n")
cat("Train proportion:", train_proportion, "\n")
if (method == "cross_validation") cat("K-folds:", k_folds, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Random seed for reproducibility
set.seed(123)

# --- Dataset Preparation Functions ---

# Extract chromosome information from sequence names
extract_chromosome <- function(seq_names) {
  chr_patterns <- c(
    "chr([0-9]+|[XY])",  # chr1, chr2, ..., chrX, chrY
    "chromosome[_ ]([0-9]+|[XY])",  # chromosome 1, chromosome_X
    "([0-9]+|[XY]):"  # 1:, X:, etc.
  )
  
  chromosomes <- rep(NA, length(seq_names))
  
  for (i in seq_along(seq_names)) {
    name <- seq_names[i]
    for (pattern in chr_patterns) {
      match <- regexpr(pattern, name, ignore.case = TRUE)
      if (match > 0) {
        matched_text <- regmatches(name, match)
        chr_num <- gsub(".*([0-9]+|[XY]).*", "\\1", matched_text, ignore.case = TRUE)
        chromosomes[i] <- chr_num
        break
      }
    }
    
    # If no pattern matched, try to extract any number or X/Y
    if (is.na(chromosomes[i])) {
      simple_match <- regexpr("([0-9]+|[XY])", name, ignore.case = TRUE)
      if (simple_match > 0) {
        chromosomes[i] <- regmatches(name, simple_match)
      }
    }
  }
  
  return(chromosomes)
}

# Prepare chromosome split datasets
prepare_chromosome_split <- function(sequences, seq_names, output_dir) {
  cat("Preparing chromosome split datasets...\n")
  
  # Extract chromosome information
  chromosomes <- extract_chromosome(seq_names)
  
  # Remove sequences without chromosome information
  valid_indices <- !is.na(chromosomes)
  if (sum(valid_indices) == 0) {
    stop("No chromosome information found in sequence names")
  }
  
  sequences <- sequences[valid_indices]
  chromosomes <- chromosomes[valid_indices]
  seq_names <- seq_names[valid_indices]
  
  cat("Found sequences from chromosomes:", paste(unique(chromosomes), collapse = ", "), "\n")
  
  # Define training chromosomes (1-15) and testing chromosomes (16-22, X, Y)
  train_chrs <- as.character(1:15)
  test_chrs <- c(as.character(16:22), "X", "Y")
  
  # Split sequences
  train_indices <- chromosomes %in% train_chrs
  test_indices <- chromosomes %in% test_chrs
  
  if (sum(train_indices) == 0) {
    stop("No training sequences found (chromosomes 1-15)")
  }
  if (sum(test_indices) == 0) {
    stop("No testing sequences found (chromosomes 16-22, X, Y)")
  }
  
  train_sequences <- sequences[train_indices]
  test_sequences <- sequences[test_indices]
  
  cat("Training sequences:", length(train_sequences), "(chromosomes 1-15)\n")
  cat("Testing sequences:", length(test_sequences), "(chromosomes 16-22, X, Y)\n")
  
  # Save datasets
  train_file <- file.path(output_dir, "training_sequences.fasta")
  test_file <- file.path(output_dir, "test_sequences.fasta")
  
  writeXStringSet(train_sequences, train_file)
  writeXStringSet(test_sequences, test_file)
  
  # Save metadata
  metadata <- list(
    method = "chromosome_split",
    train_chromosomes = train_chrs,
    test_chromosomes = test_chrs,
    train_count = length(train_sequences),
    test_count = length(test_sequences),
    train_file = train_file,
    test_file = test_file
  )
  
  return(metadata)
}

# Prepare random split datasets
prepare_random_split <- function(sequences, seq_names, output_dir, train_prop = 0.8) {
  cat("Preparing random split datasets...\n")
  
  n_total <- length(sequences)
  n_train <- floor(n_total * train_prop)
  
  # Random sampling
  train_indices <- sample(n_total, n_train)
  test_indices <- setdiff(1:n_total, train_indices)
  
  train_sequences <- sequences[train_indices]
  test_sequences <- sequences[test_indices]
  
  cat("Training sequences:", length(train_sequences), "\n")
  cat("Testing sequences:", length(test_sequences), "\n")
  
  # Save datasets
  train_file <- file.path(output_dir, "training_sequences.fasta")
  test_file <- file.path(output_dir, "test_sequences.fasta")
  
  writeXStringSet(train_sequences, train_file)
  writeXStringSet(test_sequences, test_file)
  
  # Save metadata
  metadata <- list(
    method = "random_split",
    train_proportion = train_prop,
    train_count = length(train_sequences),
    test_count = length(test_sequences),
    train_file = train_file,
    test_file = test_file
  )
  
  return(metadata)
}

# Prepare cross-validation datasets
prepare_cross_validation <- function(sequences, seq_names, output_dir, k = 5) {
  cat("Preparing", k, "-fold cross-validation datasets...\n")
  
  n_total <- length(sequences)
  fold_size <- floor(n_total / k)
  
  # Create folds
  indices <- sample(n_total)  # randomize order
  folds <- list()
  
  for (i in 1:k) {
    if (i < k) {
      fold_indices <- indices[((i-1)*fold_size + 1):(i*fold_size)]
    } else {
      # Last fold gets remaining sequences
      fold_indices <- indices[((i-1)*fold_size + 1):n_total]
    }
    folds[[i]] <- fold_indices
  }
  
  # Create directories for each fold
  cv_dir <- file.path(output_dir, "cross_validation")
  dir.create(cv_dir, recursive = TRUE, showWarnings = FALSE)
  
  fold_metadata <- list()
  
  for (i in 1:k) {
    fold_dir <- file.path(cv_dir, paste0("fold_", i))
    dir.create(fold_dir, showWarnings = FALSE)
    
    # Test set is current fold, training set is all other folds
    test_indices <- folds[[i]]
    train_indices <- unlist(folds[-i])
    
    train_sequences <- sequences[train_indices]
    test_sequences <- sequences[test_indices]
    
    # Save fold datasets
    train_file <- file.path(fold_dir, "training_sequences.fasta")
    test_file <- file.path(fold_dir, "test_sequences.fasta")
    
    writeXStringSet(train_sequences, train_file)
    writeXStringSet(test_sequences, test_file)
    
    fold_metadata[[i]] <- list(
      fold = i,
      train_count = length(train_sequences),
      test_count = length(test_sequences),
      train_file = train_file,
      test_file = test_file
    )
    
    cat("Fold", i, ": Train =", length(train_sequences), ", Test =", length(test_sequences), "\n")
  }
  
  # Overall metadata
  metadata <- list(
    method = "cross_validation",
    k_folds = k,
    total_sequences = n_total,
    folds = fold_metadata,
    cv_directory = cv_dir
  )
  
  return(metadata)
}

# --- Main Execution ---

# Check input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Load sequences
cat("Loading sequences...\n")
sequences <- readDNAStringSet(input_file)
seq_names <- names(sequences)

if (length(sequences) == 0) {
  stop("No sequences found in input file")
}

cat("Total sequences loaded:", length(sequences), "\n")

# Prepare datasets based on selected method
metadata <- switch(method,
  "chromosome_split" = prepare_chromosome_split(sequences, seq_names, output_dir),
  "random_split" = prepare_random_split(sequences, seq_names, output_dir, train_proportion),
  "cross_validation" = prepare_cross_validation(sequences, seq_names, output_dir, k_folds),
  stop("Unknown method: ", method, ". Valid methods: chromosome_split, random_split, cross_validation")
)

# Save metadata
metadata_file <- file.path(output_dir, "dataset_metadata.json")

# Create JSON manually
json_lines <- c("{")
for (i in seq_along(metadata)) {
  name <- names(metadata)[i]
  value <- metadata[[i]]
  
  if (is.character(value)) {
    json_line <- paste0('  "', name, '": "', value, '"')
  } else if (is.numeric(value)) {
    json_line <- paste0('  "', name, '": ', value)
  } else if (is.logical(value)) {
    json_line <- paste0('  "', name, '": ', tolower(as.character(value)))
  } else if (is.list(value)) {
    # Skip complex nested structures for now
    json_line <- paste0('  "', name, '": "complex_object"')
  } else {
    json_line <- paste0('  "', name, '": "', as.character(value), '"')
  }
  
  if (i < length(metadata)) {
    json_line <- paste0(json_line, ",")
  }
  json_lines <- c(json_lines, json_line)
}
json_lines <- c(json_lines, "}")

writeLines(json_lines, metadata_file)

cat("\nDataset preparation completed!\n")
cat("Method:", method, "\n")
cat("Output directory:", output_dir, "\n")
cat("Metadata saved to:", metadata_file, "\n")
cat("=== Dataset Preparation Complete ===\n")

