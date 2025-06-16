#!/usr/bin/env Rscript

# Handle Missing Data with Fallback Mechanisms
# Provides fallback options when input data is missing or corrupted
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

#' Main function for handling missing data
#' @param input_file Path to primary input file
#' @param fallback_files Vector of fallback file paths
#' @param output_file Output file path
#' @param generate_example Generate example data if all files are missing
#' @param verbose Enable verbose output
handle_missing_data <- function(input_file, fallback_files = NULL, 
                               output_file = NULL, generate_example = TRUE,
                               verbose = FALSE) {
  
  if (verbose) cat("Checking data availability...\n")
  
  # Check primary input file
  if (!is.null(input_file) && file.exists(input_file)) {
    if (verbose) cat("Primary input file found:", input_file, "\n")
    result <- validate_and_process_file(input_file, verbose)
    if (result$valid) {
      if (!is.null(output_file)) {
        file.copy(input_file, output_file, overwrite = TRUE)
        if (verbose) cat("Copied primary file to output:", output_file, "\n")
      }
      return(list(status = "primary", file = input_file, sequences = result$sequences))
    } else {
      if (verbose) cat("Primary file is invalid:", result$error, "\n")
    }
  } else {
    if (verbose) cat("Primary input file not found or not specified\n")
  }
  
  # Try fallback files
  if (!is.null(fallback_files)) {
    for (fallback in fallback_files) {
      if (verbose) cat("Trying fallback file:", fallback, "\n")
      
      if (file.exists(fallback)) {
        result <- validate_and_process_file(fallback, verbose)
        if (result$valid) {
          if (verbose) cat("Valid fallback file found:", fallback, "\n")
          if (!is.null(output_file)) {
            file.copy(fallback, output_file, overwrite = TRUE)
            if (verbose) cat("Copied fallback file to output:", output_file, "\n")
          }
          return(list(status = "fallback", file = fallback, sequences = result$sequences))
        } else {
          if (verbose) cat("Fallback file is invalid:", result$error, "\n")
        }
      } else {
        if (verbose) cat("Fallback file not found:", fallback, "\n")
      }
    }
  }
  
  # Generate example data as last resort
  if (generate_example) {
    if (verbose) cat("Generating example data...\n")
    
    example_file <- if (!is.null(output_file)) output_file else "example_sequences.fasta"
    sequences <- generate_example_sequences(verbose)
    
    writeXStringSet(sequences, example_file)
    if (verbose) cat("Generated example data:", example_file, "\n")
    
    return(list(status = "generated", file = example_file, sequences = sequences))
  }
  
  # All options exhausted
  stop("No valid data sources available and example generation disabled")
}

#' Validate and process a sequence file
#' @param file_path Path to sequence file
#' @param verbose Enable verbose output
validate_and_process_file <- function(file_path, verbose = FALSE) {
  tryCatch({
    # Check file size
    file_info <- file.info(file_path)
    if (is.na(file_info$size) || file_info$size == 0) {
      return(list(valid = FALSE, error = "File is empty"))
    }
    
    if (file_info$size > 1e9) {  # 1GB limit
      return(list(valid = FALSE, error = "File too large (>1GB)"))
    }
    
    # Try to read as FASTA
    sequences <- readDNAStringSet(file_path)
    
    if (length(sequences) == 0) {
      return(list(valid = FALSE, error = "No sequences found"))
    }
    
    # Basic validation
    if (length(sequences) < 10) {
      return(list(valid = FALSE, error = "Too few sequences (<10)"))
    }
    
    # Check sequence lengths
    seq_lengths <- width(sequences)
    if (all(seq_lengths < 20)) {
      return(list(valid = FALSE, error = "All sequences too short (<20bp)"))
    }
    
    if (verbose) {
      cat("  Sequences found:", length(sequences), "\n")
      cat("  Length range:", min(seq_lengths), "-", max(seq_lengths), "bp\n")
    }
    
    return(list(valid = TRUE, sequences = sequences))
    
  }, error = function(e) {
    return(list(valid = FALSE, error = paste("Read error:", conditionMessage(e))))
  })
}

#' Generate example CTCF-like sequences for testing
#' @param verbose Enable verbose output
generate_example_sequences <- function(verbose = FALSE) {
  
  if (verbose) cat("Generating CTCF-like example sequences...\n")
  
  # CTCF consensus sequence (from JASPAR)
  ctcf_core <- "CCGCGNGGNGGCAG"
  
  # Generate variations of the core motif
  n_sequences <- 1000
  sequence_length <- 100
  sequences <- character(n_sequences)
  
  for (i in 1:n_sequences) {
    # Create flanking sequences
    flank_length <- (sequence_length - nchar(ctcf_core)) %/% 2
    left_flank <- paste(sample(c("A", "C", "G", "T"), flank_length, replace = TRUE), collapse = "")
    right_flank <- paste(sample(c("A", "C", "G", "T"), flank_length, replace = TRUE), collapse = "")
    
    # Create motif with some variation
    motif <- ctcf_core
    
    # Replace N with random nucleotides
    motif <- gsub("N", sample(c("A", "C", "G", "T"), 1), motif)
    
    # Add some mutations (10% chance per position)
    motif_chars <- strsplit(motif, "")[[1]]
    for (j in seq_along(motif_chars)) {
      if (runif(1) < 0.1) {  # 10% mutation rate
        motif_chars[j] <- sample(c("A", "C", "G", "T"), 1)
      }
    }
    motif <- paste(motif_chars, collapse = "")
    
    # Combine flanks and motif
    full_sequence <- paste0(left_flank, motif, right_flank)
    
    # Ensure exact length
    if (nchar(full_sequence) > sequence_length) {
      full_sequence <- substr(full_sequence, 1, sequence_length)
    } else if (nchar(full_sequence) < sequence_length) {
      padding <- sequence_length - nchar(full_sequence)
      full_sequence <- paste0(full_sequence, paste(rep("A", padding), collapse = ""))
    }
    
    sequences[i] <- full_sequence
  }
  
  # Create DNAStringSet with proper names
  names(sequences) <- paste0("example_peak_", 1:n_sequences, "_chr", 
                            sample(1:22, n_sequences, replace = TRUE),
                            ":", sample(1000000:50000000, n_sequences, replace = TRUE),
                            "-", sample(1000000:50000000, n_sequences, replace = TRUE))
  
  dna_sequences <- DNAStringSet(sequences)
  
  if (verbose) {
    cat("Generated", length(dna_sequences), "example sequences\n")
    cat("Sequence length:", unique(width(dna_sequences)), "bp\n")
  }
  
  return(dna_sequences)
}

#' Check data availability and return recommendations
#' @param data_dir Data directory to check
#' @param required_files List of required file patterns
#' @param verbose Enable verbose output
check_data_availability <- function(data_dir = "data", 
                                   required_files = c("*.fasta", "*.bed", "*.peaks"),
                                   verbose = FALSE) {
  
  if (verbose) cat("Checking data availability in:", data_dir, "\n")
  
  if (!dir.exists(data_dir)) {
    if (verbose) cat("Data directory does not exist\n")
    return(list(
      status = "missing_directory",
      recommendations = "Create data directory and add input files"
    ))
  }
  
  # Check for required file types
  available_files <- list.files(data_dir, recursive = TRUE)
  
  if (length(available_files) == 0) {
    return(list(
      status = "empty_directory",
      recommendations = "Add input files to data directory"
    ))
  }
  
  # Check for specific file types
  file_checks <- list()
  for (pattern in required_files) {
    matching_files <- list.files(data_dir, pattern = glob2rx(pattern), 
                                recursive = TRUE, full.names = TRUE)
    file_checks[[pattern]] <- list(
      pattern = pattern,
      found = length(matching_files) > 0,
      files = matching_files,
      count = length(matching_files)
    )
  }
  
  if (verbose) {
    for (check in file_checks) {
      cat("Pattern", check$pattern, ":", check$count, "files found\n")
    }
  }
  
  # Generate recommendations
  missing_patterns <- sapply(file_checks, function(x) !x$found)
  
  if (all(missing_patterns)) {
    status <- "no_required_files"
    recommendations <- paste("Add files matching patterns:", 
                           paste(required_files, collapse = ", "))
  } else if (any(missing_patterns)) {
    status <- "partial_files"
    missing <- names(file_checks)[missing_patterns]
    recommendations <- paste("Consider adding files matching:", 
                           paste(missing, collapse = ", "))
  } else {
    status <- "complete"
    recommendations <- "All required file types found"
  }
  
  return(list(
    status = status,
    file_checks = file_checks,
    recommendations = recommendations
  ))
}

#' Create data directory structure with example files
#' @param base_dir Base directory for data structure
#' @param verbose Enable verbose output
create_data_structure <- function(base_dir = "data", verbose = FALSE) {
  
  if (verbose) cat("Creating data directory structure...\n")
  
  # Create directories
  dirs <- c(
    base_dir,
    file.path(base_dir, "raw"),
    file.path(base_dir, "processed"),
    file.path(base_dir, "reference_genome"),
    file.path(base_dir, "examples")
  )
  
  for (dir in dirs) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (verbose) cat("Created directory:", dir, "\n")
  }
  
  # Create example files
  example_sequences <- generate_example_sequences(verbose)
  example_file <- file.path(base_dir, "examples", "example_sequences.fasta")
  writeXStringSet(example_sequences, example_file)
  
  # Create README
  readme_content <- c(
    "# Data Directory Structure",
    "",
    "## Directories:",
    "- `raw/`: Original input files (FASTA, BED, etc.)",
    "- `processed/`: Intermediate processing results",
    "- `reference_genome/`: Reference genome files",
    "- `examples/`: Example data for testing",
    "",
    "## File Formats:",
    "- FASTA files: Sequence data",
    "- BED files: Genomic intervals",
    "- JSON files: Configuration and metadata",
    "",
    paste("Created:", Sys.time())
  )
  
  readme_file <- file.path(base_dir, "README.md")
  writeLines(readme_content, readme_file)
  
  if (verbose) {
    cat("Data structure created successfully\n")
    cat("Example file:", example_file, "\n")
    cat("README:", readme_file, "\n")
  }
  
  return(list(
    base_dir = base_dir,
    example_file = example_file,
    directories = dirs
  ))
}

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Primary input file", metavar = "character"),
    make_option(c("-f", "--fallback"), type = "character", default = NULL,
                help = "Comma-separated fallback files", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file path", metavar = "character"),
    make_option(c("-g", "--generate"), action = "store_true", default = TRUE,
                help = "Generate example data if needed [default: %default]"),
    make_option(c("-c", "--check"), type = "character", default = NULL,
                help = "Check data availability in directory", metavar = "character"),
    make_option(c("-s", "--setup"), type = "character", default = NULL,
                help = "Setup data directory structure", metavar = "character"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Handle Missing Data with Fallback Mechanisms")
  opt <- parse_args(opt_parser)
  
  # Handle different modes
  if (!is.null(opt$setup)) {
    # Setup mode
    result <- create_data_structure(opt$setup, opt$verbose)
    cat("Data structure created at:", result$base_dir, "\n")
    cat("Example file available:", result$example_file, "\n")
    quit(status = 0)
  }
  
  if (!is.null(opt$check)) {
    # Check mode
    result <- check_data_availability(opt$check, verbose = opt$verbose)
    cat("Data availability status:", result$status, "\n")
    cat("Recommendations:", result$recommendations, "\n")
    quit(status = 0)
  }
  
  # Handle missing data mode
  fallback_files <- NULL
  if (!is.null(opt$fallback)) {
    fallback_files <- trimws(strsplit(opt$fallback, ",")[[1]])
  }
  
  tryCatch({
    result <- handle_missing_data(
      input_file = opt$input,
      fallback_files = fallback_files,
      output_file = opt$output,
      generate_example = opt$generate,
      verbose = opt$verbose
    )
    
    cat("\n=== Data Handling Results ===\n")
    cat("Status:", result$status, "\n")
    cat("File used:", result$file, "\n")
    cat("Sequences:", length(result$sequences), "\n")
    
    if (result$status == "generated") {
      cat("\n⚠️  Using generated example data - replace with real data for production\n")
    } else if (result$status == "fallback") {
      cat("\n⚠️  Using fallback data - check primary data source\n")
    } else {
      cat("\n✓ Using primary data source\n")
    }
    
    quit(status = 0)
    
  }, error = function(e) {
    cat("Error handling missing data:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
