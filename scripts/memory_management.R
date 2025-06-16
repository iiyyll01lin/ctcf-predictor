#!/usr/bin/env Rscript

# Memory Management System
# Handles memory issues with reduced batch sizes and optimized processing
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

#' Main function for memory-aware processing
#' @param input_file Path to input sequences
#' @param output_file Path for output
#' @param operation Operation to perform ('alignment', 'pwm_build', 'validation')
#' @param memory_limit Memory limit in GB
#' @param batch_size Initial batch size
#' @param config_file Configuration file
#' @param verbose Enable verbose output
memory_aware_processing <- function(input_file, output_file = NULL,
                                   operation = "alignment", 
                                   memory_limit = 8, batch_size = 1000,
                                   config_file = NULL, verbose = FALSE) {
  
  if (verbose) cat("Starting memory-aware processing...\n")
  
  # Load configuration
  config <- load_memory_config(config_file, memory_limit)
  
  # Monitor initial memory usage
  initial_memory <- get_memory_usage()
  if (verbose) cat("Initial memory usage:", round(initial_memory$used_gb, 2), "GB\n")
  
  # Check available memory
  available_memory <- initial_memory$available_gb
  if (available_memory < 1) {
    warning("Less than 1GB memory available. Processing may be slow.")
  }
  
  # Read input file size and estimate memory requirements
  file_size <- file.size(input_file) / (1024^3)  # GB
  estimated_memory <- estimate_memory_usage(file_size, operation)
  
  if (verbose) {
    cat("Input file size:", round(file_size, 2), "GB\n")
    cat("Estimated memory needed:", round(estimated_memory, 2), "GB\n")
    cat("Available memory:", round(available_memory, 2), "GB\n")
  }
  
  # Determine processing strategy
  if (estimated_memory <= available_memory * 0.8) {
    # Process all at once
    if (verbose) cat("Processing entire dataset in memory\n")
    result <- process_full_dataset(input_file, output_file, operation, config, verbose)
  } else {
    # Use batch processing
    if (verbose) cat("Using batch processing to manage memory\n")
    optimal_batch_size <- calculate_optimal_batch_size(file_size, available_memory, 
                                                      batch_size, operation)
    result <- process_in_batches(input_file, output_file, operation, 
                               optimal_batch_size, config, verbose)
  }
  
  # Final memory report
  final_memory <- get_memory_usage()
  if (verbose) {
    cat("Final memory usage:", round(final_memory$used_gb, 2), "GB\n")
    cat("Peak memory increase:", round(final_memory$used_gb - initial_memory$used_gb, 2), "GB\n")
  }
  
  return(result)
}

#' Load memory management configuration
load_memory_config <- function(config_file, memory_limit) {
  default_config <- list(
    memory_limit_gb = memory_limit,
    safety_margin = 0.2,  # Reserve 20% of memory
    batch_size_min = 100,
    batch_size_max = 10000,
    temp_dir = tempdir(),
    enable_gc = TRUE,
    gc_frequency = 100,
    disk_cache = TRUE
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

#' Get current memory usage
get_memory_usage <- function() {
  # Try to get system memory info
  if (.Platform$OS.type == "windows") {
    # Windows memory info
    tryCatch({
      mem_info <- system("wmic OS get TotalVisibleMemorySize,FreePhysicalMemory /value", 
                        intern = TRUE)
      total_line <- grep("TotalVisibleMemorySize", mem_info, value = TRUE)
      free_line <- grep("FreePhysicalMemory", mem_info, value = TRUE)
      
      total_kb <- as.numeric(gsub("TotalVisibleMemorySize=", "", total_line))
      free_kb <- as.numeric(gsub("FreePhysicalMemory=", "", free_line))
      
      total_gb <- total_kb / (1024^2)
      free_gb <- free_kb / (1024^2)
      used_gb <- total_gb - free_gb
      
      return(list(
        total_gb = total_gb,
        used_gb = used_gb,
        free_gb = free_gb,
        available_gb = free_gb
      ))
    }, error = function(e) {
      # Fallback to R memory info
      mem_used <- as.numeric(object.size(ls(envir = .GlobalEnv))) / (1024^3)
      return(list(
        total_gb = 8,
        used_gb = mem_used,
        free_gb = 8 - mem_used,
        available_gb = 8 - mem_used
      ))
    })
  } else {
    # Unix-like systems
    tryCatch({
      mem_info <- system("free -b", intern = TRUE)
      mem_line <- mem_info[2]  # Second line has memory info
      mem_values <- as.numeric(strsplit(mem_line, "\\s+")[[1]][-1])
      
      total_gb <- mem_values[1] / (1024^3)
      used_gb <- mem_values[2] / (1024^3)
      free_gb <- mem_values[3] / (1024^3)
      available_gb <- mem_values[6] / (1024^3)  # Available column
      
      return(list(
        total_gb = total_gb,
        used_gb = used_gb,
        free_gb = free_gb,
        available_gb = available_gb
      ))
    }, error = function(e) {
      # Fallback
      mem_used <- sum(sapply(ls(envir = .GlobalEnv), function(x) {
        object.size(get(x, envir = .GlobalEnv))
      })) / (1024^3)
      return(list(
        total_gb = 8,
        used_gb = mem_used,
        free_gb = 8 - mem_used,
        available_gb = 8 - mem_used
      ))
    })
  }
}

#' Estimate memory usage for operation
estimate_memory_usage <- function(file_size_gb, operation) {
  # Memory multipliers based on operation type
  multipliers <- list(
    "alignment" = 3,      # Sequences + aligned + intermediate
    "pwm_build" = 2,      # Sequences + matrices
    "validation" = 4,     # Multiple PWMs + cross-validation
    "default" = 3
  )
  
  multiplier <- multipliers[[operation]] %||% multipliers[["default"]]
  return(file_size_gb * multiplier)
}

#' Calculate optimal batch size
calculate_optimal_batch_size <- function(file_size_gb, available_memory_gb, 
                                       initial_batch_size, operation) {
  
  # Estimate sequences per GB (rough approximation)
  sequences_per_gb <- 50000  # Assuming ~200bp sequences
  total_sequences <- file_size_gb * sequences_per_gb
  
  # Memory per sequence during processing
  memory_per_seq <- estimate_memory_usage(1/sequences_per_gb, operation)
  
  # Maximum sequences that fit in available memory
  max_sequences <- floor((available_memory_gb * 0.8) / memory_per_seq)
  
  # Optimal batch size
  optimal_size <- min(max_sequences, initial_batch_size)
  optimal_size <- max(100, optimal_size)  # Minimum batch size
  
  return(as.integer(optimal_size))
}

#' Process full dataset in memory
process_full_dataset <- function(input_file, output_file, operation, config, verbose) {
  
  if (verbose) cat("Loading full dataset...\n")
  
  # Read all sequences
  sequences <- readDNAStringSet(input_file)
  
  if (verbose) cat("Loaded", length(sequences), "sequences\n")
  
  # Monitor memory
  if (config$enable_gc) {
    gc()
  }
  
  # Process based on operation type
  result <- switch(operation,
    "alignment" = perform_alignment(sequences, config, verbose),
    "pwm_build" = perform_pwm_build(sequences, config, verbose),
    "validation" = perform_validation(sequences, config, verbose),
    stop("Unknown operation: ", operation)
  )
  
  # Save results
  if (!is.null(output_file) && !is.null(result$output)) {
    save_results(result$output, output_file, operation, verbose)
  }
  
  return(result)
}

#' Process dataset in batches
process_in_batches <- function(input_file, output_file, operation, batch_size, config, verbose) {
  
  if (verbose) cat("Processing in batches of size:", batch_size, "\n")
  
  # Count total sequences first
  temp_seq <- readDNAStringSet(input_file)
  total_sequences <- length(temp_seq)
  rm(temp_seq)
  gc()
  
  if (verbose) cat("Total sequences to process:", total_sequences, "\n")
  
  n_batches <- ceiling(total_sequences / batch_size)
  if (verbose) cat("Number of batches:", n_batches, "\n")
  
  # Process each batch
  batch_results <- list()
  temp_files <- character()
  
  for (batch_num in 1:n_batches) {
    if (verbose) cat("Processing batch", batch_num, "of", n_batches, "\n")
    
    # Read batch
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, total_sequences)
    
    batch_sequences <- readDNAStringSet(input_file)[start_idx:end_idx]
    
    if (verbose) cat("  Batch size:", length(batch_sequences), "\n")
    
    # Process batch
    batch_result <- switch(operation,
      "alignment" = perform_alignment(batch_sequences, config, verbose),
      "pwm_build" = perform_pwm_build(batch_sequences, config, verbose),
      "validation" = perform_validation(batch_sequences, config, verbose),
      stop("Unknown operation: ", operation)
    )
    
    # Save batch result to temporary file
    temp_file <- file.path(config$temp_dir, paste0("batch_", batch_num, ".rds"))
    saveRDS(batch_result, temp_file)
    temp_files <- c(temp_files, temp_file)
    
    # Clean up batch from memory
    rm(batch_sequences, batch_result)
    
    # Force garbage collection
    if (config$enable_gc && (batch_num %% config$gc_frequency == 0)) {
      gc()
      if (verbose) {
        mem_usage <- get_memory_usage()
        cat("  Memory usage after batch", batch_num, ":", round(mem_usage$used_gb, 2), "GB\n")
      }
    }
  }
  
  # Combine batch results
  if (verbose) cat("Combining batch results...\n")
  combined_result <- combine_batch_results(temp_files, operation, config, verbose)
  
  # Clean up temporary files
  file.remove(temp_files)
  
  # Save final results
  if (!is.null(output_file) && !is.null(combined_result$output)) {
    save_results(combined_result$output, output_file, operation, verbose)
  }
  
  return(combined_result)
}

#' Perform alignment operation
perform_alignment <- function(sequences, config, verbose) {
  
  if (verbose) cat("  Performing alignment...\n")
  
  # Simple center-based alignment for memory efficiency
  target_length <- 100
  
  aligned <- lapply(sequences, function(seq) {
    seq_len <- width(seq)
    if (seq_len >= target_length) {
      excess <- seq_len - target_length
      start_pos <- (excess %/% 2) + 1
      end_pos <- start_pos + target_length - 1
      return(subseq(seq, start_pos, end_pos))
    } else {
      return(NULL)
    }
  })
  
  # Remove NULL sequences
  aligned <- aligned[!sapply(aligned, is.null)]
  aligned_set <- DNAStringSet(aligned)
  
  return(list(
    output = aligned_set,
    type = "alignment",
    n_sequences = length(aligned_set)
  ))
}

#' Perform PWM building operation
perform_pwm_build <- function(sequences, config, verbose) {
  
  if (verbose) cat("  Building PWM...\n")
  
  # Ensure sequences are same length
  seq_lengths <- width(sequences)
  common_length <- as.integer(median(seq_lengths))
  
  # Filter and align sequences
  aligned_sequences <- sequences[seq_lengths == common_length]
  
  if (length(aligned_sequences) == 0) {
    return(list(output = NULL, type = "pwm", error = "No sequences of common length"))
  }
  
  # Build frequency matrix
  bases <- c("A", "C", "G", "T")
  freq_matrix <- matrix(0, nrow = 4, ncol = common_length, 
                       dimnames = list(bases, paste0("pos", 1:common_length)))
  
  # Count nucleotides
  for (pos in 1:common_length) {
    pos_chars <- as.character(subseq(aligned_sequences, pos, pos))
    for (base in bases) {
      freq_matrix[base, pos] <- sum(pos_chars == base)
    }
  }
  
  # Convert to PWM
  pwm <- sweep(freq_matrix, 2, colSums(freq_matrix), FUN = "/")
  
  return(list(
    output = pwm,
    type = "pwm",
    n_sequences = length(aligned_sequences)
  ))
}

#' Perform validation operation
perform_validation <- function(sequences, config, verbose) {
  
  if (verbose) cat("  Performing validation...\n")
  
  # Simple validation: build PWM and calculate quality
  pwm_result <- perform_pwm_build(sequences, config, verbose)
  
  if (is.null(pwm_result$output)) {
    return(list(output = NULL, type = "validation", error = "PWM build failed"))
  }
  
  pwm <- pwm_result$output
  
  # Calculate information content
  background_prob <- 0.25
  ic_per_pos <- apply(pwm, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col / background_prob))
  })
  
  total_ic <- sum(ic_per_pos)
  
  validation_result <- list(
    pwm = pwm,
    information_content = total_ic,
    position_ic = ic_per_pos,
    n_sequences = pwm_result$n_sequences
  )
  
  return(list(
    output = validation_result,
    type = "validation",
    n_sequences = pwm_result$n_sequences
  ))
}

#' Combine results from batch processing
combine_batch_results <- function(temp_files, operation, config, verbose) {
  
  if (length(temp_files) == 1) {
    # Single batch, just return it
    return(readRDS(temp_files[1]))
  }
  
  # Combine multiple batches
  if (operation == "alignment") {
    return(combine_alignment_batches(temp_files, verbose))
  } else if (operation == "pwm_build") {
    return(combine_pwm_batches(temp_files, verbose))
  } else if (operation == "validation") {
    return(combine_validation_batches(temp_files, verbose))
  }
}

#' Combine alignment batch results
combine_alignment_batches <- function(temp_files, verbose) {
  
  all_sequences <- list()
  total_sequences <- 0
  
  for (temp_file in temp_files) {
    batch_result <- readRDS(temp_file)
    if (!is.null(batch_result$output)) {
      all_sequences <- c(all_sequences, as.list(batch_result$output))
      total_sequences <- total_sequences + batch_result$n_sequences
    }
  }
  
  combined_sequences <- DNAStringSet(all_sequences)
  
  return(list(
    output = combined_sequences,
    type = "alignment",
    n_sequences = total_sequences
  ))
}

#' Combine PWM batch results
combine_pwm_batches <- function(temp_files, verbose) {
  
  # For PWM, we need to merge frequency matrices
  freq_matrices <- list()
  sequence_counts <- integer()
  
  for (temp_file in temp_files) {
    batch_result <- readRDS(temp_file)
    if (!is.null(batch_result$output)) {
      # Convert PWM back to frequency matrix (approximate)
      pwm <- batch_result$output
      n_seqs <- batch_result$n_sequences
      
      freq_matrix <- pwm * n_seqs  # Approximate frequency counts
      freq_matrices[[length(freq_matrices) + 1]] <- freq_matrix
      sequence_counts <- c(sequence_counts, n_seqs)
    }
  }
  
  if (length(freq_matrices) == 0) {
    return(list(output = NULL, type = "pwm", error = "No valid PWM batches"))
  }
  
  # Sum frequency matrices
  combined_freq <- Reduce("+", freq_matrices)
  
  # Convert back to PWM
  combined_pwm <- sweep(combined_freq, 2, colSums(combined_freq), FUN = "/")
  
  return(list(
    output = combined_pwm,
    type = "pwm",
    n_sequences = sum(sequence_counts)
  ))
}

#' Combine validation batch results
combine_validation_batches <- function(temp_files, verbose) {
  
  # For validation, take average of information content scores
  ic_scores <- numeric()
  sequence_counts <- integer()
  
  for (temp_file in temp_files) {
    batch_result <- readRDS(temp_file)
    if (!is.null(batch_result$output)) {
      ic_scores <- c(ic_scores, batch_result$output$information_content)
      sequence_counts <- c(sequence_counts, batch_result$n_sequences)
    }
  }
  
  if (length(ic_scores) == 0) {
    return(list(output = NULL, type = "validation", error = "No valid validation batches"))
  }
  
  # Weighted average by sequence count
  weighted_ic <- sum(ic_scores * sequence_counts) / sum(sequence_counts)
  
  validation_result <- list(
    information_content = weighted_ic,
    batch_count = length(temp_files),
    total_n_sequences = sum(sequence_counts)
  )
  
  return(list(
    output = validation_result,
    type = "validation",
    n_sequences = sum(sequence_counts)
  ))
}

#' Save results to output file
save_results <- function(result, output_file, operation, verbose) {
  
  if (operation == "alignment" && is(result, "DNAStringSet")) {
    writeXStringSet(result, output_file)
    if (verbose) cat("Aligned sequences saved to:", output_file, "\n")
  } else if (operation == "pwm_build" && is.matrix(result)) {
    saveRDS(result, output_file)
    if (verbose) cat("PWM saved to:", output_file, "\n")
  } else if (operation == "validation") {
    writeLines(toJSON(result, pretty = TRUE), output_file)
    if (verbose) cat("Validation results saved to:", output_file, "\n")
  } else {
    saveRDS(result, output_file)
    if (verbose) cat("Results saved to:", output_file, "\n")
  }
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file", metavar = "character"),
    make_option(c("-p", "--operation"), type = "character", default = "alignment",
                help = "Operation to perform [default: %default]", metavar = "character"),
    make_option(c("-m", "--memory"), type = "numeric", default = 8,
                help = "Memory limit in GB [default: %default]", metavar = "numeric"),
    make_option(c("-b", "--batch-size"), type = "integer", default = 1000,
                help = "Initial batch size [default: %default]", metavar = "integer"),
    make_option(c("-c", "--config"), type = "character", default = NULL,
                help = "Configuration file", metavar = "character"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Memory-Aware Processing System")
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input file is required.", call. = FALSE)
  }
  
  if (!file.exists(opt$input)) {
    stop("Input file does not exist: ", opt$input, call. = FALSE)
  }
  
  # Run memory-aware processing
  tryCatch({
    result <- memory_aware_processing(
      input_file = opt$input,
      output_file = opt$output,
      operation = opt$operation,
      memory_limit = opt$memory,
      batch_size = opt$`batch-size`,
      config_file = opt$config,
      verbose = opt$verbose
    )
    
    cat("\n=== Memory-Aware Processing Results ===\n")
    cat("Operation:", opt$operation, "\n")
    cat("Sequences processed:", result$n_sequences, "\n")
    cat("Result type:", result$type, "\n")
    
    if (!is.null(result$error)) {
      cat("Error:", result$error, "\n")
      quit(status = 1)
    } else {
      cat("\n✓ Processing completed successfully\n")
      quit(status = 0)
    }
    
  }, error = function(e) {
    cat("Error in memory-aware processing:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
