#!/usr/bin/env Rscript

# Chromosome Split Validation for PWM Testing Pipeline
# Performs train/test validation by splitting data by chromosomes
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

#' Main function for chromosome split validation
#' @param input_file Path to input sequences (FASTA format)
#' @param output_dir Directory for output files
#' @param test_chromosomes Chromosomes to use for testing (comma-separated)
#' @param config_file Configuration file path
#' @param verbose Enable verbose output
chromosome_split_validation <- function(input_file, output_dir = "results", 
                                       test_chromosomes = "chr21,chr22", 
                                       config_file = NULL, verbose = FALSE) {
  
  if (verbose) cat("Starting chromosome split validation...\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load configuration
  config <- load_config(config_file)
  
  # Read sequences
  if (verbose) cat("Reading sequences from:", input_file, "\n")
  sequences <- readDNAStringSet(input_file)
  
  if (length(sequences) == 0) {
    stop("No sequences found in input file")
  }
  
  # Extract chromosome information from sequence names
  seq_names <- names(sequences)
  chromosomes <- extract_chromosomes(seq_names)
  
  if (verbose) {
    cat("Total sequences:", length(sequences), "\n")
    cat("Chromosomes found:", length(unique(chromosomes)), "\n")
  }
  
  # Parse test chromosomes
  test_chrs <- trimws(strsplit(test_chromosomes, ",")[[1]])
  
  # Split data into training and testing sets
  test_idx <- chromosomes %in% test_chrs
  train_idx <- !test_idx
  
  train_sequences <- sequences[train_idx]
  test_sequences <- sequences[test_idx]
  
  if (verbose) {
    cat("Training sequences:", length(train_sequences), "\n")
    cat("Test sequences:", length(test_sequences), "\n")
    cat("Test chromosomes:", paste(test_chrs, collapse = ", "), "\n")
  }
  
  # Validate split quality
  if (length(train_sequences) == 0) {
    stop("No training sequences found. Check chromosome names.")
  }
  
  if (length(test_sequences) == 0) {
    stop("No test sequences found. Check chromosome specification.")
  }
  
  # Build PWM on training data
  if (verbose) cat("Building PWM on training data...\n")
  train_pwm <- build_pwm_from_sequences(train_sequences, config)
  
  # Build PWM on test data for comparison
  if (verbose) cat("Building PWM on test data...\n")
  test_pwm <- build_pwm_from_sequences(test_sequences, config)
  
  # Calculate quality metrics
  if (verbose) cat("Calculating validation metrics...\n")
  validation_results <- calculate_validation_metrics(train_pwm, test_pwm, 
                                                    train_sequences, test_sequences, 
                                                    config)
  
  # Generate validation report
  report <- generate_validation_report(validation_results, test_chrs, config)
  
  # Save results
  save_validation_results(report, train_pwm, test_pwm, output_dir, verbose)
  
  # Return validation summary
  return(report$summary)
}

#' Load configuration file
load_config <- function(config_file) {
  default_config <- list(
    min_sequence_length = 50,
    max_sequence_length = 200,
    pseudocount = 0.01,
    background_probs = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
    min_information_content = 8.0,
    quality_threshold = 0.7
  )
  
  if (!is.null(config_file) && file.exists(config_file)) {
    if (grepl("\\.json$", config_file)) {
      user_config <- fromJSON(config_file)
    } else if (grepl("\\.yml$|\\.yaml$", config_file)) {
      user_config <- yaml::read_yaml(config_file)
    } else {
      warning("Unsupported config format. Using defaults.")
      return(default_config)
    }
    
    # Merge configs
    for (name in names(user_config)) {
      default_config[[name]] <- user_config[[name]]
    }
  }
  
  return(default_config)
}

#' Extract chromosome information from sequence names
extract_chromosomes <- function(seq_names) {
  # Try different patterns to extract chromosome names
  patterns <- c(
    "chr[0-9XY]+",           # chr1, chr2, chrX, chrY
    "chromosome_[0-9XY]+",   # chromosome_1, chromosome_2
    "(?<=chr)[0-9XY]+",      # Extract number/letter after chr
    "[0-9XY]+(?=:)"          # Extract before colon
  )
  
  chromosomes <- character(length(seq_names))
  
  for (i in seq_along(seq_names)) {
    name <- seq_names[i]
    found <- FALSE
    
    for (pattern in patterns) {
      matches <- regmatches(name, regexpr(pattern, name, perl = TRUE))
      if (length(matches) > 0 && matches != "") {
        chromosomes[i] <- matches[1]
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      # Fallback: use sequence index
      chromosomes[i] <- paste0("seq_", i)
    }
  }
  
  return(chromosomes)
}

#' Build PWM from sequences
build_pwm_from_sequences <- function(sequences, config) {
  # Filter sequences by length
  seq_lengths <- width(sequences)
  valid_idx <- seq_lengths >= config$min_sequence_length & 
               seq_lengths <= config$max_sequence_length
  
  if (sum(valid_idx) == 0) {
    stop("No sequences meet length criteria")
  }
  
  filtered_sequences <- sequences[valid_idx]
  
  # Align sequences to the same length (center-based alignment)
  target_length <- median(width(filtered_sequences))
  aligned_sequences <- align_sequences_to_center(filtered_sequences, target_length)
  
  # Build frequency matrix
  freq_matrix <- build_frequency_matrix(aligned_sequences)
  
  # Convert to PWM with pseudocounts
  pwm <- freq_matrix_to_pwm(freq_matrix, config$pseudocount, config$background_probs)
  
  return(pwm)
}

#' Align sequences to center
align_sequences_to_center <- function(sequences, target_length) {
  aligned <- lapply(sequences, function(seq) {
    seq_len <- width(seq)
    
    if (seq_len == target_length) {
      return(seq)
    } else if (seq_len > target_length) {
      # Trim from both ends
      trim_each_side <- (seq_len - target_length) %/% 2
      start_pos <- trim_each_side + 1
      end_pos <- seq_len - trim_each_side
      return(subseq(seq, start_pos, end_pos))
    } else {
      # Pad with N's (or skip short sequences)
      return(NULL)
    }
  })
  
  # Remove NULL sequences
  aligned <- aligned[!sapply(aligned, is.null)]
  return(DNAStringSet(aligned))
}

#' Build frequency matrix from aligned sequences
build_frequency_matrix <- function(sequences) {
  if (length(sequences) == 0) {
    stop("No sequences provided for frequency matrix")
  }
  
  seq_length <- unique(width(sequences))
  if (length(seq_length) > 1) {
    stop("Sequences must be of equal length")
  }
  
  # Initialize frequency matrix
  bases <- c("A", "C", "G", "T")
  freq_matrix <- matrix(0, nrow = 4, ncol = seq_length, 
                       dimnames = list(bases, paste0("pos", 1:seq_length)))
  
  # Count nucleotides at each position
  for (pos in 1:seq_length) {
    pos_chars <- as.character(subseq(sequences, pos, pos))
    for (base in bases) {
      freq_matrix[base, pos] <- sum(pos_chars == base)
    }
  }
  
  return(freq_matrix)
}

#' Convert frequency matrix to PWM
freq_matrix_to_pwm <- function(freq_matrix, pseudocount = 0.01, background_probs = NULL) {
  if (is.null(background_probs)) {
    background_probs <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
  }
  
  # Add pseudocounts
  pwm <- freq_matrix + pseudocount
  
  # Normalize to probabilities
  pwm <- sweep(pwm, 2, colSums(pwm), FUN = "/")
  
  # Calculate information content
  ic <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col / background_probs[names(col)]))
  })
  
  # Add information content as attribute
  attr(pwm, "information_content") <- ic
  attr(pwm, "total_ic") <- sum(ic)
  
  return(pwm)
}

#' Calculate validation metrics
calculate_validation_metrics <- function(train_pwm, test_pwm, train_sequences, 
                                       test_sequences, config) {
  
  train_ic <- attr(train_pwm, "total_ic")
  test_ic <- attr(test_pwm, "total_ic")
  
  # PWM similarity (correlation between position-wise information content)
  train_pos_ic <- attr(train_pwm, "information_content")
  test_pos_ic <- attr(test_pwm, "information_content")
  
  ic_correlation <- cor(train_pos_ic, test_pos_ic, use = "complete.obs")
  
  # PWM matrix correlation
  pwm_correlation <- cor(as.vector(train_pwm), as.vector(test_pwm), use = "complete.obs")
  
  # Conservation consistency
  train_conserved <- which(train_pos_ic > quantile(train_pos_ic, 0.75))
  test_conserved <- which(test_pos_ic > quantile(test_pos_ic, 0.75))
  conserved_overlap <- length(intersect(train_conserved, test_conserved)) / 
                       length(union(train_conserved, test_conserved))
  
  # Quality metrics
  train_quality <- calculate_pwm_quality(train_pwm, config)
  test_quality <- calculate_pwm_quality(test_pwm, config)
  
  return(list(
    train_ic = train_ic,
    test_ic = test_ic,
    ic_difference = abs(train_ic - test_ic),
    ic_correlation = ic_correlation,
    pwm_correlation = pwm_correlation,
    conserved_overlap = conserved_overlap,
    train_quality = train_quality,
    test_quality = test_quality,
    train_size = length(train_sequences),
    test_size = length(test_sequences)
  ))
}

#' Calculate PWM quality score
calculate_pwm_quality <- function(pwm, config) {
  ic <- attr(pwm, "total_ic")
  pos_ic <- attr(pwm, "information_content")
  
  # Quality components
  ic_score <- min(1.0, ic / 20.0)  # Normalize to max expected IC
  conservation_score <- sum(pos_ic > 1.0) / length(pos_ic)  # Fraction of well-conserved positions
  uniformity_score <- 1 - (sd(pos_ic) / mean(pos_ic))  # Lower variation = higher score
  
  # Combined quality score
  quality <- (ic_score + conservation_score + uniformity_score) / 3
  
  return(list(
    overall = quality,
    information_content = ic_score,
    conservation = conservation_score,
    uniformity = uniformity_score
  ))
}

#' Generate validation report
generate_validation_report <- function(metrics, test_chromosomes, config) {
  # Determine validation status
  ic_difference_threshold <- 2.0
  correlation_threshold <- 0.7
  quality_threshold <- config$quality_threshold
  
  ic_pass <- metrics$ic_difference < ic_difference_threshold
  correlation_pass <- metrics$ic_correlation > correlation_threshold
  quality_pass <- metrics$train_quality$overall > quality_threshold && 
                  metrics$test_quality$overall > quality_threshold
  
  overall_pass <- ic_pass && correlation_pass && quality_pass
  
  # Generate summary
  summary <- list(
    validation_status = if(overall_pass) "PASS" else "FAIL",
    test_chromosomes = paste(test_chromosomes, collapse = ", "),
    train_sequences = metrics$train_size,
    test_sequences = metrics$test_size,
    train_ic = round(metrics$train_ic, 3),
    test_ic = round(metrics$test_ic, 3),
    ic_difference = round(metrics$ic_difference, 3),
    ic_correlation = round(metrics$ic_correlation, 3),
    pwm_correlation = round(metrics$pwm_correlation, 3),
    conserved_overlap = round(metrics$conserved_overlap, 3),
    train_quality = round(metrics$train_quality$overall, 3),
    test_quality = round(metrics$test_quality$overall, 3)
  )
  
  # Generate detailed report
  detailed <- list(
    validation_criteria = list(
      ic_difference_threshold = ic_difference_threshold,
      correlation_threshold = correlation_threshold,
      quality_threshold = quality_threshold
    ),
    test_results = list(
      ic_difference_pass = ic_pass,
      correlation_pass = correlation_pass,
      quality_pass = quality_pass
    ),
    metrics = metrics
  )
  
  return(list(
    summary = summary,
    detailed = detailed,
    timestamp = Sys.time()
  ))
}

#' Save validation results
save_validation_results <- function(report, train_pwm, test_pwm, output_dir, verbose) {
  if (verbose) cat("Saving validation results to:", output_dir, "\n")
  
  # Save summary report
  summary_file <- file.path(output_dir, "chromosome_split_validation_summary.json")
  writeLines(toJSON(report$summary, pretty = TRUE), summary_file)
  
  # Save detailed report
  detailed_file <- file.path(output_dir, "chromosome_split_validation_detailed.json")
  writeLines(toJSON(report$detailed, pretty = TRUE), detailed_file)
  
  # Save PWMs
  train_pwm_file <- file.path(output_dir, "train_pwm.rds")
  test_pwm_file <- file.path(output_dir, "test_pwm.rds")
  
  saveRDS(train_pwm, train_pwm_file)
  saveRDS(test_pwm, test_pwm_file)
  
  # Save PWMs in MEME format
  save_pwm_meme_format(train_pwm, file.path(output_dir, "train_pwm.meme"))
  save_pwm_meme_format(test_pwm, file.path(output_dir, "test_pwm.meme"))
  
  if (verbose) {
    cat("Results saved:\n")
    cat("  Summary:", summary_file, "\n")
    cat("  Detailed:", detailed_file, "\n")
    cat("  Train PWM:", train_pwm_file, "\n")
    cat("  Test PWM:", test_pwm_file, "\n")
  }
}

#' Save PWM in MEME format
save_pwm_meme_format <- function(pwm, filename) {
  ic <- attr(pwm, "total_ic")
  
  lines <- c(
    "MEME version 4",
    "",
    "ALPHABET= ACGT",
    "",
    "Background letter frequencies",
    "A 0.25 C 0.25 G 0.25 T 0.25",
    "",
    paste0("MOTIF CTCF_PWM_IC_", round(ic, 2)),
    paste0("letter-probability matrix: alength= 4 w= ", ncol(pwm), 
           " nsites= 1000 E= 0"),
    ""
  )
  
  # Add PWM matrix
  for (i in 1:ncol(pwm)) {
    row <- paste(sprintf("%.6f", pwm[, i]), collapse = " ")
    lines <- c(lines, row)
  }
  
  writeLines(lines, filename)
}

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input FASTA file with sequences", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = "results",
                help = "Output directory [default: %default]", metavar = "character"),
    make_option(c("-t", "--test-chromosomes"), type = "character", default = "chr21,chr22",
                help = "Test chromosomes (comma-separated) [default: %default]", metavar = "character"),
    make_option(c("-c", "--config"), type = "character", default = NULL,
                help = "Configuration file", metavar = "character"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Chromosome Split Validation for PWM Testing Pipeline")
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input file is required.", call. = FALSE)
  }
  
  if (!file.exists(opt$input)) {
    stop("Input file does not exist: ", opt$input, call. = FALSE)
  }
  
  # Run validation
  tryCatch({
    result <- chromosome_split_validation(
      input_file = opt$input,
      output_dir = opt$output,
      test_chromosomes = opt$test_chromosomes,
      config_file = opt$config,
      verbose = opt$verbose
    )
    
    cat("\n=== Chromosome Split Validation Results ===\n")
    cat("Status:", result$validation_status, "\n")
    cat("Train/Test IC:", result$train_ic, "/", result$test_ic, "\n")
    cat("IC Correlation:", result$ic_correlation, "\n")
    cat("PWM Quality:", result$train_quality, "/", result$test_quality, "\n")
    
    if (result$validation_status == "PASS") {
      cat("\n✓ Validation PASSED - PWM shows consistent quality across chromosomes\n")
      quit(status = 0)
    } else {
      cat("\n✗ Validation FAILED - Review detailed results for issues\n")
      quit(status = 1)
    }
    
  }, error = function(e) {
    cat("Error in chromosome split validation:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
