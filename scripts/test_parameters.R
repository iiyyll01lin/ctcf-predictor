#!/usr/bin/env Rscript

# Test Parameters
# Tests different parameter combinations for PWM building
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
  library(parallel)
})

#' Main function for parameter testing
#' @param sequences_file Path to input sequences file
#' @param output_dir Output directory for results
#' @param test_configs List of parameter configurations to test
#' @param n_cores Number of cores for parallel processing
#' @param verbose Enable verbose output
test_parameters <- function(sequences_file = "data/aligned_sequences.fasta",
                           output_dir = "results/parameter_testing",
                           test_configs = NULL,
                           n_cores = 1,
                           verbose = FALSE) {
  
  if (verbose) cat("Starting parameter testing...\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load sequences
  if (verbose) cat("Loading sequences...\n")
  sequences <- load_sequences(sequences_file, verbose)
  
  # Generate parameter configurations
  if (is.null(test_configs)) {
    if (verbose) cat("Generating parameter configurations...\n")
    test_configs <- generate_parameter_configs(verbose)
  }
  
  # Test each configuration
  if (verbose) cat("Testing", length(test_configs), "parameter configurations...\n")
  
  if (n_cores > 1 && length(test_configs) > 1) {
    # Parallel testing
    if (verbose) cat("Using parallel processing with", n_cores, "cores\n")
    
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl))
    
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(Biostrings)
      })
    })
    
    clusterExport(cl, c("sequences", "output_dir", "verbose"), envir = environment())
    
    test_results <- parLapply(cl, test_configs, function(config) {
      test_single_configuration(sequences, config, output_dir, verbose)
    })
    
  } else {
    # Sequential testing
    test_results <- lapply(test_configs, function(config) {
      test_single_configuration(sequences, config, output_dir, verbose)
    })
  }
  
  # Combine and analyze results
  if (verbose) cat("Analyzing parameter test results...\n")
  analysis <- analyze_parameter_results(test_results, test_configs, verbose)
  
  # Save comprehensive results
  comprehensive_results <- list(
    sequences_file = sequences_file,
    n_sequences = length(sequences),
    test_configs = test_configs,
    test_results = test_results,
    analysis = analysis,
    timestamp = Sys.time()
  )
  
  # Save results
  results_file <- file.path(output_dir, "parameter_testing_results.rds")
  saveRDS(comprehensive_results, results_file)
  
  # Save analysis summary
  analysis_file <- file.path(output_dir, "parameter_analysis.json")
  writeLines(toJSON(analysis, pretty = TRUE), analysis_file)
  
  # Display summary
  display_parameter_summary(analysis, verbose)
  
  if (verbose) cat("Parameter testing completed. Results saved to:", output_dir, "\n")
  
  return(comprehensive_results)
}

#' Load sequences from file
load_sequences <- function(sequences_file, verbose) {
  
  if (!file.exists(sequences_file)) {
    stop("Sequences file not found: ", sequences_file)
  }
  
  tryCatch({
    sequences <- readDNAStringSet(sequences_file)
    
    if (length(sequences) == 0) {
      stop("No sequences found in file")
    }
    
    if (verbose) cat("  Loaded", length(sequences), "sequences\n")
    
    return(sequences)
    
  }, error = function(e) {
    stop("Error loading sequences: ", conditionMessage(e))
  })
}

#' Generate parameter configurations to test
generate_parameter_configs <- function(verbose) {
  
  # Define parameter ranges
  subset_sizes <- c(500, 1000, 2000, 5000)
  pseudocounts <- c(0.01, 0.1, 0.5, 1.0)
  alignment_methods <- c("none", "consensus", "muscle")
  quality_filters <- c("none", "low", "medium", "high")
  
  # Generate all combinations
  configs <- expand.grid(
    subset_size = subset_sizes,
    pseudocount = pseudocounts,
    alignment_method = alignment_methods,
    quality_filter = quality_filters,
    stringsAsFactors = FALSE
  )
  
  # Convert to list of configurations
  config_list <- list()
  
  for (i in 1:nrow(configs)) {
    config <- list(
      id = i,
      name = paste0("config_", i),
      subset_size = configs$subset_size[i],
      pseudocount = configs$pseudocount[i],
      alignment_method = configs$alignment_method[i],
      quality_filter = configs$quality_filter[i]
    )
    
    config_list[[i]] <- config
  }
  
  if (verbose) cat("  Generated", length(config_list), "parameter configurations\n")
  
  return(config_list)
}

#' Test a single parameter configuration
test_single_configuration <- function(sequences, config, output_dir, verbose) {
  
  config_id <- config$id
  
  if (verbose) cat("  Testing configuration", config_id, ":", config$name, "\n")
  
  start_time <- Sys.time()
  
  tryCatch({
    # Apply configuration
    processed_sequences <- apply_configuration(sequences, config, verbose)
    
    # Build PWM
    pwm_result <- build_pwm_with_config(processed_sequences, config, verbose)
    
    # Calculate metrics
    metrics <- calculate_config_metrics(pwm_result, config, verbose)
    
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Compile results
    result <- list(
      config = config,
      pwm_result = pwm_result,
      metrics = metrics,
      processing_time = processing_time,
      status = "success",
      error = NULL
    )
    
    # Save individual result
    result_file <- file.path(output_dir, paste0("config_", config_id, "_result.rds"))
    saveRDS(result, result_file)
    
    if (verbose) cat("    Configuration", config_id, "completed in", 
                    round(processing_time, 2), "seconds\n")
    
    return(result)
    
  }, error = function(e) {
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    if (verbose) cat("    Configuration", config_id, "failed:", conditionMessage(e), "\n")
    
    # Return error result
    return(list(
      config = config,
      pwm_result = NULL,
      metrics = NULL,
      processing_time = processing_time,
      status = "error",
      error = conditionMessage(e)
    ))
  })
}

#' Apply configuration to sequences
apply_configuration <- function(sequences, config, verbose) {
  
  # Apply subset size
  if (config$subset_size < length(sequences)) {
    # Random sampling
    sample_indices <- sample(length(sequences), config$subset_size)
    sequences <- sequences[sample_indices]
  }
  
  # Apply quality filter
  if (config$quality_filter != "none") {
    sequences <- apply_quality_filter(sequences, config$quality_filter, verbose)
  }
  
  # Apply alignment method
  if (config$alignment_method != "none") {
    sequences <- apply_alignment_method(sequences, config$alignment_method, verbose)
  }
  
  return(sequences)
}

#' Apply quality filter
apply_quality_filter <- function(sequences, filter_level, verbose) {
  
  if (filter_level == "none") {
    return(sequences)
  }
  
  # Calculate sequence quality metrics
  sequence_lengths <- width(sequences)
  n_content <- letterFrequency(sequences, "N") / sequence_lengths
  gc_content <- letterFrequency(sequences, "GC") / sequence_lengths
  
  # Define filter thresholds
  thresholds <- list(
    low = list(max_n = 0.5, min_gc = 0.1, max_gc = 0.9),
    medium = list(max_n = 0.25, min_gc = 0.2, max_gc = 0.8),
    high = list(max_n = 0.1, min_gc = 0.3, max_gc = 0.7)
  )
  
  thresh <- thresholds[[filter_level]]
  
  # Apply filters
  keep <- (n_content <= thresh$max_n) & 
          (gc_content >= thresh$min_gc) & 
          (gc_content <= thresh$max_gc)
  
  filtered_sequences <- sequences[keep]
  
  if (verbose) cat("      Quality filter:", length(filtered_sequences), "of", 
                  length(sequences), "sequences retained\n")
  
  return(filtered_sequences)
}

#' Apply alignment method
apply_alignment_method <- function(sequences, method, verbose) {
  
  if (method == "none") {
    return(sequences)
  }
  
  if (method == "consensus") {
    # Simple consensus alignment
    aligned_sequences <- consensus_align_sequences(sequences, verbose)
  } else if (method == "muscle") {
    # MUSCLE alignment (simplified version)
    aligned_sequences <- muscle_align_sequences(sequences, verbose)
  } else {
    warning("Unknown alignment method: ", method)
    return(sequences)
  }
  
  return(aligned_sequences)
}

#' Simple consensus alignment
consensus_align_sequences <- function(sequences, verbose) {
  
  # Find most common length
  lengths <- width(sequences)
  common_length <- as.numeric(names(sort(table(lengths), decreasing = TRUE))[1])
  
  # Filter to common length
  aligned_sequences <- sequences[lengths == common_length]
  
  if (verbose) cat("      Consensus alignment:", length(aligned_sequences), 
                  "sequences of length", common_length, "\n")
  
  return(aligned_sequences)
}

#' MUSCLE alignment (simplified)
muscle_align_sequences <- function(sequences, verbose) {
  
  # For testing purposes, use a simple alignment approach
  # In practice, this would call external MUSCLE software
  
  # Find median length and trim/pad sequences
  lengths <- width(sequences)
  target_length <- median(lengths)
  
  aligned_sequences <- DNAStringSet(lapply(sequences, function(seq) {
    current_length <- length(seq)
    if (current_length > target_length) {
      # Trim from center
      start <- (current_length - target_length) %/% 2 + 1
      return(subseq(seq, start, start + target_length - 1))
    } else if (current_length < target_length) {
      # Pad with N's (not ideal, but for testing)
      padding <- paste(rep("N", target_length - current_length), collapse = "")
      return(DNAString(paste0(as.character(seq), padding)))
    } else {
      return(seq)
    }
  }))
  
  if (verbose) cat("      MUSCLE alignment:", length(aligned_sequences), 
                  "sequences aligned to length", target_length, "\n")
  
  return(aligned_sequences)
}

#' Build PWM with configuration
build_pwm_with_config <- function(sequences, config, verbose) {
  
  if (length(sequences) == 0) {
    stop("No sequences to build PWM")
  }
  
  # Ensure all sequences have same length
  seq_lengths <- width(sequences)
  if (length(unique(seq_lengths)) > 1) {
    # Use most common length
    common_length <- as.numeric(names(sort(table(seq_lengths), decreasing = TRUE))[1])
    sequences <- sequences[seq_lengths == common_length]
  }
  
  seq_length <- unique(width(sequences))[1]
  
  # Build frequency matrix
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
  
  # Add pseudocounts
  freq_matrix <- freq_matrix + config$pseudocount
  
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
  
  # Compile PWM result
  pwm_result <- list(
    pwm = prob_matrix,
    freq_matrix = freq_matrix,
    info_content = ic_per_pos,
    total_info = total_ic,
    consensus = consensus,
    n_sequences = length(sequences),
    pseudocount = config$pseudocount,
    config = config
  )
  
  return(pwm_result)
}

#' Calculate metrics for configuration
calculate_config_metrics <- function(pwm_result, config, verbose) {
  
  if (is.null(pwm_result)) {
    return(NULL)
  }
  
  # Basic metrics
  metrics <- list(
    total_information_content = pwm_result$total_info,
    mean_ic_per_position = mean(pwm_result$info_content),
    max_ic_per_position = max(pwm_result$info_content),
    sd_ic_per_position = sd(pwm_result$info_content),
    n_high_ic_positions = sum(pwm_result$info_content > 1.5),
    n_conserved_positions = sum(pwm_result$info_content > 1.0),
    pwm_length = length(pwm_result$info_content),
    n_sequences_used = pwm_result$n_sequences,
    consensus_length = nchar(pwm_result$consensus)
  )
  
  # Quality assessment
  if (metrics$total_information_content >= 15) {
    metrics$quality_grade <- "Excellent"
  } else if (metrics$total_information_content >= 10) {
    metrics$quality_grade <- "Good"
  } else if (metrics$total_information_content >= 5) {
    metrics$quality_grade <- "Fair"
  } else {
    metrics$quality_grade <- "Poor"
  }
  
  # Configuration efficiency (IC per sequence)
  metrics$ic_per_sequence <- metrics$total_information_content / metrics$n_sequences_used
  
  # Pattern strength
  metrics$pattern_strength <- metrics$n_high_ic_positions / metrics$pwm_length
  
  return(metrics)
}

#' Analyze parameter test results
analyze_parameter_results <- function(test_results, test_configs, verbose) {
  
  # Extract successful results
  successful_results <- test_results[sapply(test_results, function(x) x$status == "success")]
  failed_results <- test_results[sapply(test_results, function(x) x$status == "error")]
  
  if (verbose) {
    cat("  Successful configurations:", length(successful_results), "\n")
    cat("  Failed configurations:", length(failed_results), "\n")
  }
  
  if (length(successful_results) == 0) {
    return(list(
      n_successful = 0,
      n_failed = length(failed_results),
      best_config = NULL,
      parameter_effects = NULL,
      summary_stats = NULL
    ))
  }
  
  # Create results data frame
  results_df <- data.frame()
  
  for (result in successful_results) {
    config <- result$config
    metrics <- result$metrics
    
    row_data <- data.frame(
      config_id = config$id,
      subset_size = config$subset_size,
      pseudocount = config$pseudocount,
      alignment_method = config$alignment_method,
      quality_filter = config$quality_filter,
      total_ic = metrics$total_information_content,
      mean_ic = metrics$mean_ic_per_position,
      max_ic = metrics$max_ic_per_position,
      n_high_ic_pos = metrics$n_high_ic_positions,
      quality_grade = metrics$quality_grade,
      ic_per_sequence = metrics$ic_per_sequence,
      pattern_strength = metrics$pattern_strength,
      processing_time = result$processing_time,
      stringsAsFactors = FALSE
    )
    
    results_df <- rbind(results_df, row_data)
  }
  
  # Find best configuration
  best_idx <- which.max(results_df$total_ic)
  best_config <- successful_results[[best_idx]]$config
  
  # Analyze parameter effects
  parameter_effects <- analyze_parameter_effects(results_df, verbose)
  
  # Summary statistics
  summary_stats <- list(
    mean_ic = mean(results_df$total_ic),
    median_ic = median(results_df$total_ic),
    sd_ic = sd(results_df$total_ic),
    min_ic = min(results_df$total_ic),
    max_ic = max(results_df$total_ic),
    mean_processing_time = mean(results_df$processing_time),
    median_processing_time = median(results_df$processing_time)
  )
  
  # Grade distribution
  grade_distribution <- table(results_df$quality_grade)
  
  # Compile analysis
  analysis <- list(
    n_successful = length(successful_results),
    n_failed = length(failed_results),
    best_config = best_config,
    best_ic = results_df$total_ic[best_idx],
    parameter_effects = parameter_effects,
    summary_stats = summary_stats,
    grade_distribution = grade_distribution,
    results_df = results_df,
    failed_configs = sapply(failed_results, function(x) x$config$id)
  )
  
  return(analysis)
}

#' Analyze effects of individual parameters
analyze_parameter_effects <- function(results_df, verbose) {
  
  effects <- list()
  
  # Subset size effect
  if (length(unique(results_df$subset_size)) > 1) {
    size_effect <- aggregate(results_df$total_ic, by = list(subset_size = results_df$subset_size), 
                            FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
    effects$subset_size <- size_effect
  }
  
  # Pseudocount effect
  if (length(unique(results_df$pseudocount)) > 1) {
    pc_effect <- aggregate(results_df$total_ic, by = list(pseudocount = results_df$pseudocount), 
                          FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
    effects$pseudocount <- pc_effect
  }
  
  # Alignment method effect
  if (length(unique(results_df$alignment_method)) > 1) {
    align_effect <- aggregate(results_df$total_ic, by = list(alignment_method = results_df$alignment_method), 
                             FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
    effects$alignment_method <- align_effect
  }
  
  # Quality filter effect
  if (length(unique(results_df$quality_filter)) > 1) {
    quality_effect <- aggregate(results_df$total_ic, by = list(quality_filter = results_df$quality_filter), 
                               FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
    effects$quality_filter <- quality_effect
  }
  
  return(effects)
}

#' Display parameter testing summary
display_parameter_summary <- function(analysis, verbose) {
  
  cat("\n=== Parameter Testing Summary ===\n")
  
  cat("Configurations tested:", analysis$n_successful + analysis$n_failed, "\n")
  cat("Successful:", analysis$n_successful, "\n")
  cat("Failed:", analysis$n_failed, "\n")
  
  if (analysis$n_successful > 0) {
    cat("\nBest Configuration:\n")
    best_config <- analysis$best_config
    cat("  ID:", best_config$id, "\n")
    cat("  Subset size:", best_config$subset_size, "\n")
    cat("  Pseudocount:", best_config$pseudocount, "\n")
    cat("  Alignment method:", best_config$alignment_method, "\n")
    cat("  Quality filter:", best_config$quality_filter, "\n")
    cat("  Information content:", round(analysis$best_ic, 3), "bits\n")
    
    # Summary statistics
    cat("\nInformation Content Statistics:\n")
    stats <- analysis$summary_stats
    cat("  Mean:", round(stats$mean_ic, 3), "bits\n")
    cat("  Median:", round(stats$median_ic, 3), "bits\n")
    cat("  Range:", round(stats$min_ic, 3), "-", round(stats$max_ic, 3), "bits\n")
    cat("  SD:", round(stats$sd_ic, 3), "bits\n")
    
    # Quality grade distribution
    cat("\nQuality Grade Distribution:\n")
    for (grade in names(analysis$grade_distribution)) {
      count <- analysis$grade_distribution[grade]
      cat("  ", grade, ":", count, "\n")
    }
    
    # Processing time
    cat("\nProcessing Time:\n")
    cat("  Mean:", round(stats$mean_processing_time, 2), "seconds\n")
    cat("  Median:", round(stats$median_processing_time, 2), "seconds\n")
    
    # Parameter effects
    if (!is.null(analysis$parameter_effects$subset_size)) {
      cat("\nSubset Size Effect (mean IC):\n")
      size_effects <- analysis$parameter_effects$subset_size
      for (i in 1:nrow(size_effects)) {
        size <- size_effects$subset_size[i]
        mean_ic <- size_effects$x[i, "mean"]
        cat("  ", size, ":", round(mean_ic, 3), "bits\n")
      }
    }
    
    if (!is.null(analysis$parameter_effects$pseudocount)) {
      cat("\nPseudocount Effect (mean IC):\n")
      pc_effects <- analysis$parameter_effects$pseudocount
      for (i in 1:nrow(pc_effects)) {
        pc <- pc_effects$pseudocount[i]
        mean_ic <- pc_effects$x[i, "mean"]
        cat("  ", pc, ":", round(mean_ic, 3), "bits\n")
      }
    }
  }
  
  if (analysis$n_failed > 0) {
    cat("\nFailed Configuration IDs:", paste(analysis$failed_configs, collapse = ", "), "\n")
  }
}

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-s", "--sequences"), type = "character", default = "data/aligned_sequences.fasta",
                help = "Input sequences file [default: %default]", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = "results/parameter_testing",
                help = "Output directory [default: %default]", metavar = "character"),
    make_option(c("-c", "--cores"), type = "integer", default = 1,
                help = "Number of cores for parallel processing [default: %default]", metavar = "integer"),
    make_option(c("--subset-sizes"), type = "character", default = "500,1000,2000,5000",
                help = "Subset sizes to test (comma-separated)", metavar = "character"),
    make_option(c("--pseudocounts"), type = "character", default = "0.01,0.1,0.5,1.0",
                help = "Pseudocounts to test (comma-separated)", metavar = "character"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Test Different Parameter Combinations for PWM Building")
  opt <- parse_args(opt_parser)
  
  if (!file.exists(opt$sequences)) {
    stop("Sequences file does not exist: ", opt$sequences, call. = FALSE)
  }
  
  # Parse parameter lists
  subset_sizes <- as.numeric(strsplit(opt$`subset-sizes`, ",")[[1]])
  pseudocounts <- as.numeric(strsplit(opt$pseudocounts, ",")[[1]])
  
  # Generate custom configurations if provided
  custom_configs <- NULL
  if (!is.null(subset_sizes) && !is.null(pseudocounts)) {
    configs_grid <- expand.grid(
      subset_size = subset_sizes,
      pseudocount = pseudocounts,
      alignment_method = c("none", "consensus"),
      quality_filter = c("none", "medium"),
      stringsAsFactors = FALSE
    )
    
    custom_configs <- list()
    for (i in 1:nrow(configs_grid)) {
      custom_configs[[i]] <- list(
        id = i,
        name = paste0("custom_config_", i),
        subset_size = configs_grid$subset_size[i],
        pseudocount = configs_grid$pseudocount[i],
        alignment_method = configs_grid$alignment_method[i],
        quality_filter = configs_grid$quality_filter[i]
      )
    }
  }
  
  # Run parameter testing
  tryCatch({
    results <- test_parameters(
      sequences_file = opt$sequences,
      output_dir = opt$output,
      test_configs = custom_configs,
      n_cores = opt$cores,
      verbose = opt$verbose
    )
    
    cat("Parameter testing completed successfully.\n")
    
  }, error = function(e) {
    cat("Error in parameter testing:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
