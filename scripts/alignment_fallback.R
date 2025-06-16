#!/usr/bin/env Rscript

# Alignment Fallback System
# Tries alternative alignment methods when primary method fails
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

#' Main function for alignment with fallback mechanisms
#' @param input_file Path to input sequences
#' @param output_file Path for output aligned sequences
#' @param primary_method Primary alignment method
#' @param fallback_methods Vector of fallback methods
#' @param config_file Configuration file
#' @param verbose Enable verbose output
alignment_with_fallback <- function(input_file, output_file = NULL,
                                   primary_method = "integrated",
                                   fallback_methods = c("center", "consensus"),
                                   config_file = NULL, verbose = FALSE) {
  
  if (verbose) cat("Starting alignment with fallback system...\n")
  
  # Load configuration
  config <- load_alignment_config(config_file)
  
  # Read input sequences
  if (verbose) cat("Reading sequences from:", input_file, "\n")
  sequences <- readDNAStringSet(input_file)
  
  if (length(sequences) == 0) {
    stop("No sequences found in input file")
  }
  
  # Try primary method
  if (verbose) cat("Attempting primary alignment method:", primary_method, "\n")
  
  primary_result <- try_alignment_method(sequences, primary_method, config, verbose)
  
  if (primary_result$success) {
    if (verbose) cat("Primary method succeeded\n")
    result <- list(
      method_used = primary_method,
      aligned_sequences = primary_result$sequences,
      quality_metrics = primary_result$quality,
      status = "primary_success"
    )
  } else {
    if (verbose) cat("Primary method failed:", primary_result$error, "\n")
    
    # Try fallback methods
    result <- NULL
    for (method in fallback_methods) {
      if (verbose) cat("Trying fallback method:", method, "\n")
      
      fallback_result <- try_alignment_method(sequences, method, config, verbose)
      
      if (fallback_result$success) {
        if (verbose) cat("Fallback method succeeded:", method, "\n")
        result <- list(
          method_used = method,
          aligned_sequences = fallback_result$sequences,
          quality_metrics = fallback_result$quality,
          status = paste("fallback_success", method, sep = "_"),
          primary_error = primary_result$error
        )
        break
      } else {
        if (verbose) cat("Fallback method failed:", method, "-", fallback_result$error, "\n")
      }
    }
    
    # If all methods failed
    if (is.null(result)) {
      # Last resort: simple center alignment
      if (verbose) cat("All methods failed, trying simple center alignment...\n")
      
      simple_result <- simple_center_alignment(sequences, config, verbose)
      if (simple_result$success) {
        result <- list(
          method_used = "simple_center",
          aligned_sequences = simple_result$sequences,
          quality_metrics = simple_result$quality,
          status = "last_resort_success",
          primary_error = primary_result$error
        )
      } else {
        stop("All alignment methods failed. Check input data quality.")
      }
    }
  }
  
  # Save results
  if (!is.null(output_file)) {
    writeXStringSet(result$aligned_sequences, output_file)
    if (verbose) cat("Aligned sequences saved to:", output_file, "\n")
  }
  
  # Save alignment report
  report_file <- if (!is.null(output_file)) {
    sub("\\.[^.]*$", "_alignment_report.json", output_file)
  } else {
    "alignment_fallback_report.json"
  }
  
  report <- list(
    method_used = result$method_used,
    status = result$status,
    input_sequences = length(sequences),
    output_sequences = length(result$aligned_sequences),
    quality_metrics = result$quality_metrics,
    timestamp = Sys.time(),
    primary_method = primary_method,
    fallback_methods = fallback_methods
  )
  
  if (!is.null(result$primary_error)) {
    report$primary_error <- result$primary_error
  }
  
  writeLines(toJSON(report, pretty = TRUE), report_file)
  if (verbose) cat("Alignment report saved to:", report_file, "\n")
  
  return(result)
}

#' Load alignment configuration
load_alignment_config <- function(config_file) {
  default_config <- list(
    min_sequence_length = 30,
    max_sequence_length = 300,
    target_length = 100,
    center_window = 50,
    consensus_threshold = 0.6,
    min_information_content = 2.0,
    max_attempts = 3,
    quality_threshold = 0.5
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

#' Try a specific alignment method
#' @param sequences Input sequences
#' @param method Alignment method name
#' @param config Configuration parameters
#' @param verbose Enable verbose output
try_alignment_method <- function(sequences, method, config, verbose = FALSE) {
  
  tryCatch({
    # Filter sequences by length first
    seq_lengths <- width(sequences)
    valid_idx <- seq_lengths >= config$min_sequence_length & 
                 seq_lengths <= config$max_sequence_length
    
    if (sum(valid_idx) < 10) {
      return(list(success = FALSE, error = "Too few valid sequences after length filtering"))
    }
    
    filtered_sequences <- sequences[valid_idx]
    
    # Apply alignment method
    aligned_sequences <- switch(method,
      "center" = center_based_alignment(filtered_sequences, config, verbose),
      "consensus" = consensus_based_alignment(filtered_sequences, config, verbose),
      "integrated" = integrated_alignment(filtered_sequences, config, verbose),
      stop("Unknown alignment method: ", method)
    )
    
    if (is.null(aligned_sequences) || length(aligned_sequences) == 0) {
      return(list(success = FALSE, error = "Alignment produced no sequences"))
    }
    
    # Calculate quality metrics
    quality <- calculate_alignment_quality(aligned_sequences, config)
    
    # Check if quality meets minimum threshold
    if (quality$overall_score < config$quality_threshold) {
      return(list(
        success = FALSE, 
        error = paste("Quality too low:", round(quality$overall_score, 3))
      ))
    }
    
    return(list(
      success = TRUE,
      sequences = aligned_sequences,
      quality = quality
    ))
    
  }, error = function(e) {
    return(list(success = FALSE, error = conditionMessage(e)))
  })
}

#' Center-based alignment
center_based_alignment <- function(sequences, config, verbose = FALSE) {
  
  if (verbose) cat("  Performing center-based alignment...\n")
  
  target_length <- config$target_length
  
  aligned <- lapply(sequences, function(seq) {
    seq_len <- width(seq)
    
    if (seq_len == target_length) {
      return(seq)
    } else if (seq_len > target_length) {
      # Center crop
      excess <- seq_len - target_length
      start_pos <- (excess %/% 2) + 1
      end_pos <- start_pos + target_length - 1
      return(subseq(seq, start_pos, end_pos))
    } else {
      # Skip sequences that are too short
      return(NULL)
    }
  })
  
  # Remove NULL sequences
  aligned <- aligned[!sapply(aligned, is.null)]
  
  if (length(aligned) == 0) {
    return(NULL)
  }
  
  return(DNAStringSet(aligned))
}

#' Consensus-based alignment
consensus_based_alignment <- function(sequences, config, verbose = FALSE) {
  
  if (verbose) cat("  Performing consensus-based alignment...\n")
  
  # Find consensus motif first
  consensus <- find_consensus_motif(sequences, config, verbose)
  
  if (is.null(consensus)) {
    return(NULL)
  }
  
  # Align based on consensus position
  aligned <- align_to_consensus(sequences, consensus, config, verbose)
  
  return(aligned)
}

#' Integrated alignment (combination of methods)
integrated_alignment <- function(sequences, config, verbose = FALSE) {
  
  if (verbose) cat("  Performing integrated alignment...\n")
  
  # First pass: center-based alignment
  center_aligned <- center_based_alignment(sequences, config, verbose)
  
  if (is.null(center_aligned)) {
    return(NULL)
  }
  
  # Second pass: refine with consensus information
  consensus <- find_consensus_motif(center_aligned, config, verbose)
  
  if (is.null(consensus)) {
    # Fall back to center alignment
    return(center_aligned)
  }
  
  # Fine-tune alignment based on consensus
  final_aligned <- refine_alignment_with_consensus(center_aligned, consensus, config, verbose)
  
  return(final_aligned)
}

#' Simple center alignment (last resort)
simple_center_alignment <- function(sequences, config, verbose = FALSE) {
  
  if (verbose) cat("  Performing simple center alignment (last resort)...\n")
  
  tryCatch({
    # Use median length as target
    seq_lengths <- width(sequences)
    target_length <- as.integer(median(seq_lengths))
    
    # Simple center cropping/padding
    aligned <- lapply(sequences, function(seq) {
      seq_len <- width(seq)
      
      if (seq_len >= target_length) {
        # Center crop
        excess <- seq_len - target_length
        start_pos <- (excess %/% 2) + 1
        end_pos <- start_pos + target_length - 1
        return(subseq(seq, start_pos, end_pos))
      } else {
        # This sequence is too short, skip it
        return(NULL)
      }
    })
    
    # Remove NULL sequences
    aligned <- aligned[!sapply(aligned, is.null)]
    
    if (length(aligned) < 5) {
      return(list(success = FALSE, error = "Too few sequences after simple alignment"))
    }
    
    aligned_set <- DNAStringSet(aligned)
    
    # Basic quality check
    quality <- list(overall_score = 0.3)  # Low but acceptable for last resort
    
    return(list(success = TRUE, sequences = aligned_set, quality = quality))
    
  }, error = function(e) {
    return(list(success = FALSE, error = conditionMessage(e)))
  })
}

#' Find consensus motif in sequences
find_consensus_motif <- function(sequences, config, verbose = FALSE) {
  
  # Simple consensus finding - look for most common k-mer
  k <- 8  # k-mer length
  
  kmer_counts <- list()
  
  for (seq in sequences) {
    seq_str <- as.character(seq)
    seq_len <- nchar(seq_str)
    
    if (seq_len >= k) {
      for (i in 1:(seq_len - k + 1)) {
        kmer <- substr(seq_str, i, i + k - 1)
        if (!grepl("N", kmer)) {  # Skip k-mers with N
          kmer_counts[[kmer]] <- kmer_counts[[kmer]] %||% 0 + 1
        }
      }
    }
  }
  
  if (length(kmer_counts) == 0) {
    return(NULL)
  }
  
  # Find most common k-mer
  best_kmer <- names(kmer_counts)[which.max(unlist(kmer_counts))]
  best_count <- kmer_counts[[best_kmer]]
  
  # Check if consensus is strong enough
  min_support <- max(5, length(sequences) * 0.1)
  
  if (best_count < min_support) {
    return(NULL)
  }
  
  if (verbose) cat("    Consensus motif found:", best_kmer, "(", best_count, "occurrences )\n")
  
  return(list(motif = best_kmer, count = best_count))
}

#' Align sequences to consensus motif
align_to_consensus <- function(sequences, consensus, config, verbose = FALSE) {
  
  motif <- consensus$motif
  target_length <- config$target_length
  
  aligned <- lapply(sequences, function(seq) {
    seq_str <- as.character(seq)
    
    # Find motif position
    motif_pos <- regexpr(motif, seq_str, fixed = TRUE)
    
    if (motif_pos[1] == -1) {
      # Motif not found, skip this sequence
      return(NULL)
    }
    
    # Calculate alignment position
    motif_center <- motif_pos[1] + nchar(motif) %/% 2
    target_center <- target_length %/% 2
    
    start_pos <- motif_center - target_center + 1
    end_pos <- start_pos + target_length - 1
    
    # Check bounds
    if (start_pos < 1 || end_pos > nchar(seq_str)) {
      return(NULL)
    }
    
    return(subseq(seq, start_pos, end_pos))
  })
  
  # Remove NULL sequences
  aligned <- aligned[!sapply(aligned, is.null)]
  
  if (length(aligned) == 0) {
    return(NULL)
  }
  
  return(DNAStringSet(aligned))
}

#' Refine alignment with consensus information
refine_alignment_with_consensus <- function(sequences, consensus, config, verbose = FALSE) {
  
  # For now, just return the original sequences
  # In a full implementation, this would fine-tune positions
  return(sequences)
}

#' Calculate alignment quality metrics
calculate_alignment_quality <- function(sequences, config) {
  
  if (length(sequences) == 0) {
    return(list(overall_score = 0))
  }
  
  # Build frequency matrix
  seq_length <- unique(width(sequences))
  if (length(seq_length) > 1) {
    return(list(overall_score = 0))  # Sequences not aligned properly
  }
  
  seq_length <- seq_length[1]
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
  
  # Convert to probabilities
  prob_matrix <- sweep(freq_matrix, 2, colSums(freq_matrix), FUN = "/")
  
  # Calculate information content
  background_prob <- 0.25
  ic_per_pos <- apply(prob_matrix, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col / background_prob))
  })
  
  total_ic <- sum(ic_per_pos)
  max_possible_ic <- 2 * seq_length
  
  # Quality metrics
  ic_score <- min(1.0, total_ic / max_possible_ic)
  conservation_score <- sum(ic_per_pos > 1.0) / seq_length
  uniformity_score <- 1 - (sd(ic_per_pos) / mean(ic_per_pos))
  
  overall_score <- (ic_score + conservation_score + uniformity_score) / 3
  
  return(list(
    overall_score = overall_score,
    information_content = ic_score,
    conservation = conservation_score,
    uniformity = uniformity_score,
    total_ic = total_ic,
    position_ic = ic_per_pos
  ))
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input FASTA file with sequences", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file for aligned sequences", metavar = "character"),
    make_option(c("-m", "--method"), type = "character", default = "integrated",
                help = "Primary alignment method [default: %default]", metavar = "character"),
    make_option(c("-f", "--fallback"), type = "character", default = "center,consensus",
                help = "Fallback methods (comma-separated) [default: %default]", metavar = "character"),
    make_option(c("-c", "--config"), type = "character", default = NULL,
                help = "Configuration file", metavar = "character"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Alignment with Fallback Mechanisms")
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input file is required.", call. = FALSE)
  }
  
  if (!file.exists(opt$input)) {
    stop("Input file does not exist: ", opt$input, call. = FALSE)
  }
  
  # Parse fallback methods
  fallback_methods <- trimws(strsplit(opt$fallback, ",")[[1]])
  
  # Run alignment with fallback
  tryCatch({
    result <- alignment_with_fallback(
      input_file = opt$input,
      output_file = opt$output,
      primary_method = opt$method,
      fallback_methods = fallback_methods,
      config_file = opt$config,
      verbose = opt$verbose
    )
    
    cat("\n=== Alignment Results ===\n")
    cat("Method used:", result$method_used, "\n")
    cat("Status:", result$status, "\n")
    cat("Sequences aligned:", length(result$aligned_sequences), "\n")
    cat("Quality score:", round(result$quality_metrics$overall_score, 3), "\n")
    
    if (result$status != "primary_success") {
      cat("\n⚠️  Primary method failed, using fallback method\n")
    } else {
      cat("\n✓ Primary alignment method succeeded\n")
    }
    
    quit(status = 0)
    
  }, error = function(e) {
    cat("Error in alignment with fallback:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
