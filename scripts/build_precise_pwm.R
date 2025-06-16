#!/usr/bin/env Rscript

# Build Precise CTCF PWM Script
# Builds high-precision Position Weight Matrices for CTCF binding prediction
# Author: CTCF Predictor Pipeline
# Usage: Rscript build_precise_pwm.R --input <sequences_file> --output <output_dir> [options]

suppressMessages({
  library(optparse)
  library(Biostrings)
  library(seqLogo)
  library(ggplot2)
  library(dplyr)
  library(jsonlite)
  library(parallel)
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input sequences file (FASTA format)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="precise_pwm_output",
              help="Output directory [default=%default]", metavar="character"),
  make_option(c("-w", "--width"), type="numeric", default=19,
              help="PWM width in base pairs [default=%default]", metavar="numeric"),
  make_option(c("--min_ic"), type="numeric", default=1.5,
              help="Minimum information content per position [default=%default]", metavar="numeric"),
  make_option(c("--pseudocount"), type="numeric", default=0.25,
              help="Pseudocount for matrix calculation [default=%default]", metavar="numeric"),
  make_option(c("--background"), type="character", default="uniform",
              help="Background model: uniform, genome, or file path [default=%default]", metavar="character"),
  make_option(c("--method"), type="character", default="enhanced",
              help="PWM building method: basic, enhanced, or iterative [default=%default]", metavar="character"),
  make_option(c("--iterations"), type="numeric", default=10,
              help="Number of iterations for iterative method [default=%default]", metavar="numeric"),
  make_option(c("--filter_quality"), action="store_true", default=FALSE,
              help="Apply sequence quality filtering"),
  make_option(c("--remove_redundant"), action="store_true", default=FALSE,
              help="Remove redundant sequences"),
  make_option(c("--similarity_threshold"), type="numeric", default=0.9,
              help="Similarity threshold for redundancy removal [default=%default]", metavar="numeric"),
  make_option(c("--cores"), type="numeric", default=1,
              help="Number of cores for parallel processing [default=%default]", metavar="numeric"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input sequences file must be specified.", call.=FALSE)
}

if (!file.exists(opt$input)) {
  stop("Input file does not exist: ", opt$input, call.=FALSE)
}

# Create output directory
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

if (opt$verbose) {
  cat("Starting precise PWM construction...\n")
  cat("Input file:", opt$input, "\n")
  cat("Output directory:", opt$output, "\n")
  cat("Method:", opt$method, "\n")
}

# Function to load and validate sequences
load_sequences <- function(input_file) {
  if (opt$verbose) cat("Loading sequences from:", input_file, "\n")
  
  # Load sequences
  sequences <- readDNAStringSet(input_file)
  
  if (length(sequences) == 0) {
    stop("No sequences found in input file.")
  }
  
  # Basic validation
  seq_lengths <- width(sequences)
  if (any(seq_lengths < opt$width)) {
    warning("Some sequences are shorter than specified PWM width. These will be filtered out.")
    sequences <- sequences[seq_lengths >= opt$width]
  }
  
  if (opt$verbose) cat("Loaded", length(sequences), "sequences\n")
  
  return(sequences)
}

# Function to filter sequences by quality
filter_sequences_by_quality <- function(sequences) {
  if (opt$verbose) cat("Filtering sequences by quality...\n")
  
  # Calculate N content
  n_content <- sapply(sequences, function(x) {
    sum(letterFrequency(x, "N")) / width(x)
  })
  
  # Filter sequences with too many Ns
  high_quality <- n_content < 0.1
  filtered_sequences <- sequences[high_quality]
  
  if (opt$verbose) cat("Filtered to", length(filtered_sequences), "high-quality sequences\n")
  
  return(filtered_sequences)
}

# Function to remove redundant sequences
remove_redundant_sequences <- function(sequences) {
  if (opt$verbose) cat("Removing redundant sequences...\n")
  
  # Calculate pairwise similarities (simplified approach)
  seq_chars <- as.character(sequences)
  unique_sequences <- unique(seq_chars)
  
  # More sophisticated similarity filtering
  if (length(unique_sequences) < length(sequences)) {
    sequences <- DNAStringSet(unique_sequences)
  }
  
  # Additional similarity-based filtering could be implemented here
  if (opt$verbose) cat("Retained", length(sequences), "non-redundant sequences\n")
  
  return(sequences)
}

# Function to get background frequencies
get_background_frequencies <- function(background_spec, sequences = NULL) {
  if (opt$verbose) cat("Setting up background model:", background_spec, "\n")
  
  if (background_spec == "uniform") {
    bg_freq <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
  } else if (background_spec == "genome") {
    # Use typical mammalian genome composition
    bg_freq <- c(A = 0.295, C = 0.205, G = 0.205, T = 0.295)
  } else if (file.exists(background_spec)) {
    # Load from file
    bg_data <- read.table(background_spec, header = TRUE, stringsAsFactors = FALSE)
    bg_freq <- setNames(bg_data$frequency, bg_data$base)
  } else {
    # Calculate from input sequences
    if (!is.null(sequences)) {
      total_freq <- letterFrequency(sequences, c("A", "C", "G", "T"), as.prob = TRUE)
      bg_freq <- colMeans(total_freq)
    } else {
      bg_freq <- c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
    }
  }
  
  return(bg_freq)
}

# Function to align sequences to motif center
align_sequences <- function(sequences, method = "center") {
  if (opt$verbose) cat("Aligning sequences...\n")
  
  # Simple center alignment - more sophisticated methods could be implemented
  seq_lengths <- width(sequences)
  target_length <- opt$width
  
  aligned_sequences <- DNAStringSet()
  
  for (i in seq_along(sequences)) {
    seq <- sequences[[i]]
    seq_len <- width(seq)
    
    if (seq_len == target_length) {
      aligned_sequences <- c(aligned_sequences, seq)
    } else if (seq_len > target_length) {
      # Extract center region
      start_pos <- ceiling((seq_len - target_length) / 2)
      end_pos <- start_pos + target_length - 1
      aligned_seq <- subseq(seq, start_pos, end_pos)
      aligned_sequences <- c(aligned_sequences, aligned_seq)
    }
  }
  
  if (opt$verbose) cat("Aligned", length(aligned_sequences), "sequences to width", target_length, "\n")
  
  return(aligned_sequences)
}

# Function to build basic PWM
build_basic_pwm <- function(sequences, background_freq, pseudocount = 0.25) {
  if (opt$verbose) cat("Building basic PWM...\n")
  
  # Convert to matrix
  seq_matrix <- consensusMatrix(sequences, as.prob = FALSE)[1:4, ]
  
  # Add pseudocounts
  pwm_counts <- seq_matrix + pseudocount
  
  # Convert to frequencies
  pwm_freq <- t(t(pwm_counts) / colSums(pwm_counts))
  
  # Convert to log-odds
  pwm_logodds <- log2(pwm_freq / background_freq)
  
  # Calculate information content
  ic <- sapply(1:ncol(pwm_freq), function(i) {
    p <- pwm_freq[, i]
    sum(p * log2(p / background_freq))
  })
  
  return(list(
    frequencies = pwm_freq,
    log_odds = pwm_logodds,
    information_content = ic,
    consensus = consensusString(DNAStringSet(apply(pwm_freq, 2, function(x) names(x)[which.max(x)])))
  ))
}

# Function to build enhanced PWM with quality filtering
build_enhanced_pwm <- function(sequences, background_freq, pseudocount = 0.25) {
  if (opt$verbose) cat("Building enhanced PWM with quality filtering...\n")
  
  # Initial PWM
  initial_pwm <- build_basic_pwm(sequences, background_freq, pseudocount)
  
  # Score sequences against initial PWM
  scores <- sapply(sequences, function(seq) {
    seq_chars <- strsplit(as.character(seq), "")[[1]]
    sum(sapply(1:length(seq_chars), function(i) {
      base <- seq_chars[i]
      if (base %in% rownames(initial_pwm$log_odds) && i <= ncol(initial_pwm$log_odds)) {
        initial_pwm$log_odds[base, i]
      } else {
        0
      }
    }))
  })
  
  # Filter sequences by score (keep top quantile)
  score_threshold <- quantile(scores, 0.25)
  high_scoring_sequences <- sequences[scores >= score_threshold]
  
  if (opt$verbose) cat("Using", length(high_scoring_sequences), "high-scoring sequences for final PWM\n")
  
  # Build final PWM
  final_pwm <- build_basic_pwm(high_scoring_sequences, background_freq, pseudocount)
  
  return(final_pwm)
}

# Function to build iterative PWM
build_iterative_pwm <- function(sequences, background_freq, pseudocount = 0.25, max_iterations = 10) {
  if (opt$verbose) cat("Building iterative PWM...\n")
  
  current_sequences <- sequences
  
  for (iteration in 1:max_iterations) {
    if (opt$verbose) cat("Iteration", iteration, "of", max_iterations, "\n")
    
    # Build PWM with current sequences
    current_pwm <- build_basic_pwm(current_sequences, background_freq, pseudocount)
    
    # Score all original sequences
    scores <- sapply(sequences, function(seq) {
      seq_chars <- strsplit(as.character(seq), "")[[1]]
      sum(sapply(1:length(seq_chars), function(i) {
        base <- seq_chars[i]
        if (base %in% rownames(current_pwm$log_odds) && i <= ncol(current_pwm$log_odds)) {
          current_pwm$log_odds[base, i]
        } else {
          0
        }
      }))
    })
    
    # Select top-scoring sequences for next iteration
    score_threshold <- quantile(scores, 1 - (0.8 - iteration * 0.05))  # Gradually more stringent
    selected_sequences <- sequences[scores >= score_threshold]
    
    # Check for convergence
    if (length(selected_sequences) == length(current_sequences)) {
      if (opt$verbose) cat("Converged at iteration", iteration, "\n")
      break
    }
    
    current_sequences <- selected_sequences
  }
  
  if (opt$verbose) cat("Final PWM built with", length(current_sequences), "sequences\n")
  
  return(current_pwm)
}

# Function to evaluate PWM quality
evaluate_pwm_quality <- function(pwm, sequences) {
  if (opt$verbose) cat("Evaluating PWM quality...\n")
  
  metrics <- list(
    width = ncol(pwm$frequencies),
    total_information_content = sum(pwm$information_content),
    mean_information_content = mean(pwm$information_content),
    max_information_content = max(pwm$information_content),
    positions_above_threshold = sum(pwm$information_content >= opt$min_ic),
    consensus_sequence = pwm$consensus
  )
  
  # Calculate sequence scores
  scores <- sapply(sequences, function(seq) {
    seq_chars <- strsplit(as.character(seq), "")[[1]]
    sum(sapply(1:min(length(seq_chars), ncol(pwm$log_odds)), function(i) {
      base <- seq_chars[i]
      if (base %in% rownames(pwm$log_odds)) {
        pwm$log_odds[base, i]
      } else {
        0
      }
    }))
  })
  
  metrics$mean_sequence_score <- mean(scores)
  metrics$score_std <- sd(scores)
  metrics$score_range <- range(scores)
  
  return(metrics)
}

# Function to create PWM visualizations
create_pwm_plots <- function(pwm, metrics, output_dir) {
  if (opt$verbose) cat("Creating PWM visualizations...\n")
  
  # Sequence logo
  if (requireNamespace("seqLogo", quietly = TRUE)) {
    png(file.path(output_dir, "sequence_logo.png"), width = 800, height = 400)
    seqLogo(pwm$frequencies)
    dev.off()
  }
  
  # Information content plot
  ic_data <- data.frame(
    position = 1:length(pwm$information_content),
    information_content = pwm$information_content
  )
  
  p1 <- ggplot(ic_data, aes(x = position, y = information_content)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = opt$min_ic, color = "red", linetype = "dashed") +
    labs(title = "Information Content per Position",
         x = "Position", y = "Information Content (bits)") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "information_content.png"), p1, width = 10, height = 6)
  
  # PWM heatmap
  pwm_df <- as.data.frame(pwm$log_odds)
  pwm_df$base <- rownames(pwm_df)
  pwm_long <- tidyr::gather(pwm_df, key = "position", value = "log_odds", -base)
  pwm_long$position <- as.numeric(gsub("V", "", pwm_long$position))
  
  p2 <- ggplot(pwm_long, aes(x = position, y = base, fill = log_odds)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Log Odds") +
    labs(title = "PWM Log-Odds Heatmap",
         x = "Position", y = "Base") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "pwm_heatmap.png"), p2, width = 10, height = 4)
  
  if (opt$verbose) cat("Plots saved to:", output_dir, "\n")
}

# Function to export PWM in various formats
export_pwm <- function(pwm, metrics, output_dir) {
  if (opt$verbose) cat("Exporting PWM in multiple formats...\n")
  
  # Save frequency matrix
  write.table(pwm$frequencies, file.path(output_dir, "pwm_frequencies.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Save log-odds matrix
  write.table(pwm$log_odds, file.path(output_dir, "pwm_log_odds.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Save information content
  write.table(data.frame(position = 1:length(pwm$information_content),
                        information_content = pwm$information_content),
              file.path(output_dir, "information_content.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save metrics as JSON
  writeLines(toJSON(metrics, pretty = TRUE),
             file.path(output_dir, "pwm_metrics.json"))
  
  # Export in MEME format
  meme_content <- paste0(
    "MEME version 4\n\n",
    "ALPHABET= ACGT\n\n",
    "strands: + -\n\n",
    "Background letter frequencies\n",
    "A 0.25 C 0.25 G 0.25 T 0.25\n\n",
    "MOTIF CTCF_precise\n",
    "letter-probability matrix: alength= 4 w= ", ncol(pwm$frequencies), " nsites= 1000 E= 0\n"
  )
  
  for (i in 1:ncol(pwm$frequencies)) {
    meme_content <- paste0(meme_content, paste(pwm$frequencies[, i], collapse = " "), "\n")
  }
  
  writeLines(meme_content, file.path(output_dir, "ctcf_precise.meme"))
  
  if (opt$verbose) cat("PWM exported in multiple formats\n")
}

# Main execution
tryCatch({
  # Load sequences
  sequences <- load_sequences(opt$input)
  
  # Apply quality filtering if requested
  if (opt$filter_quality) {
    sequences <- filter_sequences_by_quality(sequences)
  }
  
  # Remove redundant sequences if requested
  if (opt$remove_redundant) {
    sequences <- remove_redundant_sequences(sequences)
  }
  
  # Align sequences
  aligned_sequences <- align_sequences(sequences)
  
  if (length(aligned_sequences) == 0) {
    stop("No sequences remaining after processing.")
  }
  
  # Get background frequencies
  background_freq <- get_background_frequencies(opt$background, aligned_sequences)
  
  # Build PWM using specified method
  if (opt$method == "basic") {
    pwm <- build_basic_pwm(aligned_sequences, background_freq, opt$pseudocount)
  } else if (opt$method == "enhanced") {
    pwm <- build_enhanced_pwm(aligned_sequences, background_freq, opt$pseudocount)
  } else if (opt$method == "iterative") {
    pwm <- build_iterative_pwm(aligned_sequences, background_freq, opt$pseudocount, opt$iterations)
  } else {
    stop("Unknown PWM building method: ", opt$method)
  }
  
  # Evaluate PWM quality
  metrics <- evaluate_pwm_quality(pwm, aligned_sequences)
  
  # Create visualizations
  create_pwm_plots(pwm, metrics, opt$output)
  
  # Export PWM
  export_pwm(pwm, metrics, opt$output)
  
  # Print summary
  cat("\n=== Precise PWM Construction Summary ===\n")
  cat("Method:", opt$method, "\n")
  cat("Input sequences:", length(sequences), "\n")
  cat("Final sequences used:", length(aligned_sequences), "\n")
  cat("PWM width:", metrics$width, "\n")
  cat("Total information content:", round(metrics$total_information_content, 3), "bits\n")
  cat("Mean information content:", round(metrics$mean_information_content, 3), "bits\n")
  cat("Positions above IC threshold:", metrics$positions_above_threshold, "\n")
  cat("Consensus sequence:", metrics$consensus_sequence, "\n")
  cat("Mean sequence score:", round(metrics$mean_sequence_score, 3), "\n")
  cat("\nResults saved to:", opt$output, "\n")
  
}, error = function(e) {
  cat("Error in precise PWM construction:", e$message, "\n")
  quit(status = 1)
})

if (opt$verbose) cat("Precise PWM construction completed successfully!\n")
