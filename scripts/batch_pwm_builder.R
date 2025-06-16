#!/usr/bin/env Rscript

# Batch PWM Builder
# Processes multiple datasets and builds PWMs in batch
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(parallel)
  library(jsonlite)
})

#' Batch PWM Builder - process multiple datasets
#' @param input_pattern File pattern for input sequences (e.g., "data/*.fasta")
#' @param output_dir Output directory for results
#' @param batch_config Configuration for batch processing
#' @param n_cores Number of cores for parallel processing
#' @param verbose Enable verbose output
batch_pwm_builder <- function(input_pattern = "data/*.fasta",
                             output_dir = "results/batch_pwm",
                             batch_config = NULL,
                             n_cores = 1,
                             verbose = FALSE) {
  
  if (verbose) cat("Batch PWM Builder starting...\n")
  
  start_time <- Sys.time()
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Find input files
  if (verbose) cat("Finding input files...\n")
  input_files <- find_input_files(input_pattern, verbose)
  
  if (length(input_files) == 0) {
    stop("No input files found matching pattern: ", input_pattern)
  }
  
  # Load batch configuration
  if (is.null(batch_config)) {
    batch_config <- get_default_batch_config()
  }
  
  # Process files
  if (verbose) cat("Processing", length(input_files), "files...\n")
  
  if (n_cores > 1 && length(input_files) > 1) {
    # Parallel processing
    if (verbose) cat("Using parallel processing with", n_cores, "cores\n")
    batch_results <- process_files_parallel(input_files, output_dir, batch_config, n_cores, verbose)
  } else {
    # Sequential processing
    batch_results <- process_files_sequential(input_files, output_dir, batch_config, verbose)
  }
  
  # Analyze batch results
  if (verbose) cat("Analyzing batch results...\n")
  batch_analysis <- analyze_batch_results(batch_results, verbose)
  
  # Generate batch report
  batch_report <- generate_batch_report(input_files, batch_results, batch_analysis, batch_config)
  
  # Save batch report
  report_file <- file.path(output_dir, "batch_report.json")
  writeLines(toJSON(batch_report, pretty = TRUE), report_file)
  
  # Save batch summary
  summary_file <- file.path(output_dir, "batch_summary.txt")
  save_batch_summary(batch_report, summary_file, verbose)
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Display summary
  display_batch_summary(batch_report, total_time, verbose)
  
  if (verbose) cat("Batch processing completed in", round(total_time, 2), "seconds\n")
  
  return(batch_report)
}

#' Find input files based on pattern
find_input_files <- function(input_pattern, verbose) {
  
  # Handle different pattern types
  if (grepl("\\*", input_pattern)) {
    # Glob pattern
    input_files <- Sys.glob(input_pattern)
  } else if (file.exists(input_pattern)) {
    # Single file
    input_files <- input_pattern
  } else if (dir.exists(input_pattern)) {
    # Directory - find FASTA files
    input_files <- list.files(input_pattern, pattern = "\\.(fasta|fa|fas)$", 
                             full.names = TRUE, recursive = FALSE)
  } else {
    input_files <- character()
  }
  
  # Filter for existing files
  input_files <- input_files[file.exists(input_files)]
  
  if (verbose) cat("  Found", length(input_files), "input files\n")
  
  return(input_files)
}

#' Get default batch configuration
get_default_batch_config <- function() {
  
  config <- list(
    # PWM building parameters
    pseudocount = 0.1,
    max_sequence_length = 200,
    min_sequence_length = 10,
    
    # Quality filtering
    max_n_content = 0.25,
    min_gc_content = 0.2,
    max_gc_content = 0.8,
    
    # Subsampling
    max_sequences_per_file = 10000,
    subsample_method = "random",
    
    # Output options
    save_individual_results = TRUE,
    save_combined_matrix = TRUE,
    generate_logos = FALSE,
    
    # Processing options
    continue_on_error = TRUE,
    cleanup_intermediate = FALSE
  )
  
  return(config)
}

#' Process files in parallel
process_files_parallel <- function(input_files, output_dir, batch_config, n_cores, verbose) {
  
  # Set up cluster
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  
  # Load libraries on workers
  clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(Biostrings))
  })
  
  # Export necessary objects
  clusterExport(cl, c("batch_config", "output_dir", "verbose"), envir = environment())
  
  # Process files in parallel
  batch_results <- parLapply(cl, input_files, function(input_file) {
    process_single_file(input_file, output_dir, batch_config, verbose)
  })
  
  # Add file names
  names(batch_results) <- basename(input_files)
  
  return(batch_results)
}

#' Process files sequentially
process_files_sequential <- function(input_files, output_dir, batch_config, verbose) {
  
  batch_results <- list()
  
  for (i in seq_along(input_files)) {
    input_file <- input_files[i]
    file_name <- basename(input_file)
    
    if (verbose) cat("  Processing file", i, "of", length(input_files), ":", file_name, "\n")
    
    # Process single file
    result <- process_single_file(input_file, output_dir, batch_config, verbose)
    batch_results[[file_name]] <- result
  }
  
  return(batch_results)
}

#' Process a single file
process_single_file <- function(input_file, output_dir, batch_config, verbose) {
  
  file_name <- basename(input_file)
  base_name <- gsub("\\.[^.]*$", "", file_name)
  
  start_time <- Sys.time()
  
  tryCatch({
    # Load sequences
    sequences <- readDNAStringSet(input_file)
    
    if (length(sequences) == 0) {
      return(list(
        file = file_name,
        status = "error",
        error = "No sequences found",
        processing_time = 0
      ))
    }
    
    # Apply filters and subsampling
    processed_sequences <- preprocess_sequences_batch(sequences, batch_config, verbose)
    
    if (length(processed_sequences) == 0) {
      return(list(
        file = file_name,
        status = "error", 
        error = "No sequences remaining after filtering",
        processing_time = 0
      ))
    }
    
    # Build PWM
    pwm_result <- build_pwm_batch(processed_sequences, batch_config, verbose)
    
    # Save individual results if requested
    if (batch_config$save_individual_results) {
      output_file <- file.path(output_dir, paste0(base_name, "_pwm.rds"))
      saveRDS(pwm_result, output_file)
    }
    
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Return result summary
    return(list(
      file = file_name,
      base_name = base_name,
      status = "success",
      n_input_sequences = length(sequences),
      n_processed_sequences = length(processed_sequences),
      total_information_content = pwm_result$total_info,
      quality_grade = assess_quality_grade(pwm_result$total_info),
      consensus = pwm_result$consensus,
      processing_time = processing_time,
      pwm_result = pwm_result,
      error = NULL
    ))
    
  }, error = function(e) {
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    if (verbose) cat("    Error processing", file_name, ":", conditionMessage(e), "\n")
    
    if (!batch_config$continue_on_error) {
      stop("Processing stopped due to error in ", file_name, ": ", conditionMessage(e))
    }
    
    return(list(
      file = file_name,
      status = "error",
      error = conditionMessage(e),
      processing_time = processing_time
    ))
  })
}

#' Preprocess sequences for batch processing
preprocess_sequences_batch <- function(sequences, batch_config, verbose) {
  
  # Length filtering
  seq_lengths <- width(sequences)
  length_filter <- (seq_lengths >= batch_config$min_sequence_length) & 
                   (seq_lengths <= batch_config$max_sequence_length)
  sequences <- sequences[length_filter]
  
  if (length(sequences) == 0) {
    return(sequences)
  }
  
  # Quality filtering
  n_content <- letterFrequency(sequences, "N") / width(sequences)
  gc_content <- letterFrequency(sequences, "GC") / width(sequences)
  
  quality_filter <- (n_content <= batch_config$max_n_content) &
                    (gc_content >= batch_config$min_gc_content) &
                    (gc_content <= batch_config$max_gc_content)
  
  sequences <- sequences[quality_filter]
  
  if (length(sequences) == 0) {
    return(sequences)
  }
  
  # Subsampling if needed
  if (length(sequences) > batch_config$max_sequences_per_file) {
    if (batch_config$subsample_method == "random") {
      sample_indices <- sample(length(sequences), batch_config$max_sequences_per_file)
      sequences <- sequences[sample_indices]
    } else {
      # Take first N sequences
      sequences <- sequences[1:batch_config$max_sequences_per_file]
    }
  }
  
  # Length standardization - use most common length
  seq_lengths <- width(sequences)
  common_length <- as.numeric(names(sort(table(seq_lengths), decreasing = TRUE))[1])
  sequences <- sequences[seq_lengths == common_length]
  
  return(sequences)
}

#' Build PWM for batch processing
build_pwm_batch <- function(sequences, batch_config, verbose) {
  
  if (length(sequences) == 0) {
    stop("No sequences provided for PWM building")
  }
  
  seq_length <- unique(width(sequences))[1]
  n_sequences <- length(sequences)
  
  # Build frequency matrix
  bases <- c("A", "C", "G", "T")
  freq_matrix <- matrix(0, nrow = 4, ncol = seq_length, 
                       dimnames = list(bases, paste0("pos", 1:seq_length)))
  
  # Count nucleotides
  for (pos in 1:seq_length) {
    pos_chars <- as.character(subseq(sequences, pos, pos))
    for (base in bases) {
      freq_matrix[base, pos] <- sum(pos_chars == base)
    }
  }
  
  # Add pseudocounts
  freq_matrix <- freq_matrix + batch_config$pseudocount
  
  # Convert to probabilities
  prob_matrix <- sweep(freq_matrix, 2, colSums(freq_matrix), FUN = "/")
  
  # Calculate information content
  background_prob <- 0.25
  ic_per_pos <- apply(prob_matrix, 2, function(col) {
    col[col == 0] <- 1e-10
    2 + sum(col * log2(col / background_prob))
  })
  
  total_ic <- sum(ic_per_pos)
  
  # Generate consensus sequence
  consensus <- paste(bases[apply(prob_matrix, 2, which.max)], collapse = "")
  
  # Return PWM result
  pwm_result <- list(
    pwm = prob_matrix,
    freq_matrix = freq_matrix,
    info_content = ic_per_pos,
    total_info = total_ic,
    consensus = consensus,
    n_sequences = n_sequences,
    method = "batch_pwm_builder"
  )
  
  return(pwm_result)
}

#' Assess quality grade
assess_quality_grade <- function(total_ic) {
  if (total_ic >= 20) {
    return("Excellent")
  } else if (total_ic >= 15) {
    return("Very Good")
  } else if (total_ic >= 10) {
    return("Good")
  } else if (total_ic >= 5) {
    return("Fair")
  } else {
    return("Poor")
  }
}

#' Analyze batch results
analyze_batch_results <- function(batch_results, verbose) {
  
  # Separate successful and failed results
  successful <- batch_results[sapply(batch_results, function(x) x$status == "success")]
  failed <- batch_results[sapply(batch_results, function(x) x$status == "error")]
  
  if (length(successful) == 0) {
    return(list(
      n_successful = 0,
      n_failed = length(failed),
      success_rate = 0,
      summary_stats = NULL,
      quality_distribution = NULL,
      failed_files = sapply(failed, function(x) x$file)
    ))
  }
  
  # Extract metrics from successful results
  ic_values <- sapply(successful, function(x) x$total_information_content)
  quality_grades <- sapply(successful, function(x) x$quality_grade)
  processing_times <- sapply(successful, function(x) x$processing_time)
  n_sequences <- sapply(successful, function(x) x$n_processed_sequences)
  
  # Summary statistics
  summary_stats <- list(
    mean_ic = mean(ic_values),
    median_ic = median(ic_values),
    sd_ic = sd(ic_values),
    min_ic = min(ic_values),
    max_ic = max(ic_values),
    mean_processing_time = mean(processing_times),
    median_processing_time = median(processing_times),
    total_sequences_processed = sum(n_sequences),
    mean_sequences_per_file = mean(n_sequences)
  )
  
  # Quality distribution
  quality_distribution <- table(quality_grades)
  
  # Best performing file
  best_idx <- which.max(ic_values)
  best_file <- successful[[best_idx]]
  
  analysis <- list(
    n_successful = length(successful),
    n_failed = length(failed),
    success_rate = length(successful) / length(batch_results),
    summary_stats = summary_stats,
    quality_distribution = quality_distribution,
    best_file = best_file,
    failed_files = sapply(failed, function(x) x$file)
  )
  
  return(analysis)
}

#' Generate batch report
generate_batch_report <- function(input_files, batch_results, batch_analysis, batch_config) {
  
  report <- list(
    metadata = list(
      timestamp = Sys.time(),
      n_input_files = length(input_files),
      batch_config = batch_config
    ),
    processing_summary = list(
      n_successful = batch_analysis$n_successful,
      n_failed = batch_analysis$n_failed,
      success_rate = batch_analysis$success_rate
    ),
    results = batch_results,
    analysis = batch_analysis
  )
  
  return(report)
}

#' Save batch summary
save_batch_summary <- function(batch_report, summary_file, verbose) {
  
  summary_lines <- c(
    "=== Batch PWM Builder Summary ===",
    paste("Timestamp:", batch_report$metadata$timestamp),
    paste("Input files:", batch_report$metadata$n_input_files),
    "",
    "Processing Results:",
    paste("  Successful:", batch_report$processing_summary$n_successful),
    paste("  Failed:", batch_report$processing_summary$n_failed), 
    paste("  Success rate:", round(batch_report$processing_summary$success_rate * 100, 1), "%"),
    ""
  )
  
  if (batch_report$processing_summary$n_successful > 0) {
    stats <- batch_report$analysis$summary_stats
    
    summary_lines <- c(summary_lines,
      "Information Content Statistics:",
      paste("  Mean:", round(stats$mean_ic, 3), "bits"),
      paste("  Median:", round(stats$median_ic, 3), "bits"),
      paste("  Range:", round(stats$min_ic, 3), "-", round(stats$max_ic, 3), "bits"),
      paste("  SD:", round(stats$sd_ic, 3), "bits"),
      "",
      "Processing Time:",
      paste("  Mean:", round(stats$mean_processing_time, 2), "seconds"),
      paste("  Median:", round(stats$median_processing_time, 2), "seconds"),
      "",
      "Quality Grade Distribution:"
    )
    
    quality_dist <- batch_report$analysis$quality_distribution
    for (grade in names(quality_dist)) {
      summary_lines <- c(summary_lines,
        paste("  ", grade, ":", quality_dist[grade]))
    }
    
    # Best file
    best_file <- batch_report$analysis$best_file
    summary_lines <- c(summary_lines,
      "",
      "Best Performing File:",
      paste("  File:", best_file$file),
      paste("  Information Content:", round(best_file$total_information_content, 3), "bits"),
      paste("  Quality Grade:", best_file$quality_grade)
    )
  }
  
  if (batch_report$processing_summary$n_failed > 0) {
    failed_files <- batch_report$analysis$failed_files
    summary_lines <- c(summary_lines,
      "",
      "Failed Files:",
      paste("  ", failed_files)
    )
  }
  
  writeLines(summary_lines, summary_file)
  
  if (verbose) cat("  Batch summary saved to:", summary_file, "\n")
}

#' Display batch summary
display_batch_summary <- function(batch_report, total_time, verbose) {
  
  cat("\n=== Batch PWM Builder Summary ===\n")
  
  cat("Input files:", batch_report$metadata$n_input_files, "\n")
  cat("Successful:", batch_report$processing_summary$n_successful, "\n")
  cat("Failed:", batch_report$processing_summary$n_failed, "\n")
  cat("Success rate:", round(batch_report$processing_summary$success_rate * 100, 1), "%\n")
  cat("Total processing time:", round(total_time, 2), "seconds\n")
  
  if (batch_report$processing_summary$n_successful > 0) {
    stats <- batch_report$analysis$summary_stats
    
    cat("\nInformation Content Statistics:\n")
    cat("  Mean:", round(stats$mean_ic, 3), "bits\n")
    cat("  Range:", round(stats$min_ic, 3), "-", round(stats$max_ic, 3), "bits\n")
    
    cat("\nQuality Grade Distribution:\n")
    quality_dist <- batch_report$analysis$quality_distribution
    for (grade in names(quality_dist)) {
      cat("  ", grade, ":", quality_dist[grade], "\n")
    }
    
    # Processing efficiency
    sequences_per_second <- stats$total_sequences_processed / total_time
    cat("\nProcessing Efficiency:", round(sequences_per_second, 1), "sequences/second\n")
    
    # Best file
    best_file <- batch_report$analysis$best_file
    cat("\nBest Result:", best_file$file, "(", round(best_file$total_information_content, 3), "bits )\n")
  }
  
  if (batch_report$processing_summary$n_failed > 0) {
    cat("\nFailed files:", paste(batch_report$analysis$failed_files, collapse = ", "), "\n")
  }
}

# Command-line interface
if (!interactive()) {
  # Simple argument parsing
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript batch_pwm_builder.R <input_pattern> [output_dir] [n_cores]\n")
    cat("Examples:\n")
    cat("  Rscript batch_pwm_builder.R 'data/*.fasta'\n")
    cat("  Rscript batch_pwm_builder.R 'data/' results/batch 4\n")
    quit(status = 1)
  }
  
  input_pattern <- args[1]
  output_dir <- if (length(args) >= 2) args[2] else "results/batch_pwm"
  n_cores <- if (length(args) >= 3) as.numeric(args[3]) else 1
  
  # Run batch processing
  tryCatch({
    cat("Starting batch PWM building...\n")
    cat("Input pattern:", input_pattern, "\n")
    cat("Output directory:", output_dir, "\n")
    cat("Cores:", n_cores, "\n\n")
    
    batch_report <- batch_pwm_builder(
      input_pattern = input_pattern,
      output_dir = output_dir,
      batch_config = NULL,
      n_cores = n_cores,
      verbose = TRUE
    )
    
    cat("Batch PWM building completed successfully.\n")
    
    # Exit with appropriate status
    if (batch_report$processing_summary$success_rate >= 0.8) {
      quit(status = 0)
    } else {
      quit(status = 1) 
    }
    
  }, error = function(e) {
    cat("Error in batch PWM building:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
