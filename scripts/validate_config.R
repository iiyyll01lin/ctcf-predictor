#!/usr/bin/env Rscript
# Configuration Validation Script
# Validates CTCF pipeline configuration files

# Load required libraries
suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(optparse)
})

#' Parse command line arguments
parse_arguments <- function() {
  option_list <- list(
    make_option(c("-c", "--config"), type = "character", default = "config/pipeline.yml",
                help = "Configuration file to validate [default: %default]"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output validation report file")
  )
  
  parser <- OptionParser(
    usage = "%prog [options]",
    description = "Validate CTCF pipeline configuration files",
    option_list = option_list
  )
  
  return(parse_args(parser))
}

#' Load configuration file
load_config <- function(config_file) {
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }
  
  tryCatch({
    config <- yaml.load_file(config_file)
    return(config)
  }, error = function(e) {
    stop("Error loading configuration file: ", e$message)
  })
}

#' Validate configuration structure
validate_configuration <- function(config, config_file) {
  validation_results <- list(
    file = config_file,
    errors = c(),
    warnings = c(),
    passed = TRUE
  )
  
  # Required top-level sections
  required_sections <- c("data", "alignment", "pwm", "statistics", "performance", "output", "logging")
  
  for (section in required_sections) {
    if (!section %in% names(config)) {
      validation_results$errors <- c(validation_results$errors, 
                                   paste("Missing required section:", section))
      validation_results$passed <- FALSE
    }
  }
  
  # Validate data section
  if ("data" %in% names(config)) {
    validation_results <- validate_data_section(config$data, validation_results)
  }
  
  # Validate alignment section
  if ("alignment" %in% names(config)) {
    validation_results <- validate_alignment_section(config$alignment, validation_results)
  }
  
  # Validate PWM section
  if ("pwm" %in% names(config)) {
    validation_results <- validate_pwm_section(config$pwm, validation_results)
  }
  
  # Validate statistics section
  if ("statistics" %in% names(config)) {
    validation_results <- validate_statistics_section(config$statistics, validation_results)
  }
  
  # Validate performance section
  if ("performance" %in% names(config)) {
    validation_results <- validate_performance_section(config$performance, validation_results)
  }
  
  # Validate output section
  if ("output" %in% names(config)) {
    validation_results <- validate_output_section(config$output, validation_results)
  }
  
  return(validation_results)
}

#' Validate data section
validate_data_section <- function(data_config, validation_results) {
  
  # Check input parameters
  if ("input" %in% names(data_config)) {
    input_config <- data_config$input
    
    # Validate format
    if ("format" %in% names(input_config)) {
      valid_formats <- c("fasta", "bed", "custom")
      if (!input_config$format %in% valid_formats) {
        validation_results$errors <- c(validation_results$errors,
                                     paste("Invalid input format:", input_config$format,
                                           "Valid formats:", paste(valid_formats, collapse = ", ")))
        validation_results$passed <- FALSE
      }
    }
    
    # Validate sequence limits
    if ("max_sequences" %in% names(input_config)) {
      if (!is.numeric(input_config$max_sequences) || input_config$max_sequences <= 0) {
        validation_results$errors <- c(validation_results$errors,
                                     "max_sequences must be a positive number")
        validation_results$passed <- FALSE
      }
    }
  }
  
  # Check preprocessing parameters
  if ("preprocessing" %in% names(data_config)) {
    preprocessing_config <- data_config$preprocessing
    
    # Validate thresholds
    threshold_params <- c("quality_threshold", "complexity_threshold")
    for (param in threshold_params) {
      if (param %in% names(preprocessing_config)) {
        value <- preprocessing_config[[param]]
        if (!is.numeric(value) || value < 0 || value > 1) {
          validation_results$errors <- c(validation_results$errors,
                                       paste(param, "must be between 0 and 1"))
          validation_results$passed <- FALSE
        }
      }
    }
    
    # Validate length parameters
    length_params <- c("min_length", "max_length")
    for (param in length_params) {
      if (param %in% names(preprocessing_config)) {
        value <- preprocessing_config[[param]]
        if (!is.numeric(value) || value <= 0) {
          validation_results$errors <- c(validation_results$errors,
                                       paste(param, "must be a positive number"))
          validation_results$passed <- FALSE
        }
      }
    }
    
    # Check min_length < max_length
    if ("min_length" %in% names(preprocessing_config) && 
        "max_length" %in% names(preprocessing_config)) {
      if (preprocessing_config$min_length >= preprocessing_config$max_length) {
        validation_results$errors <- c(validation_results$errors,
                                     "min_length must be less than max_length")
        validation_results$passed <- FALSE
      }
    }
  }
  
  return(validation_results)
}

#' Validate alignment section
validate_alignment_section <- function(alignment_config, validation_results) {
  
  # Validate alignment method
  if ("method" %in% names(alignment_config)) {
    valid_methods <- c("center", "consensus", "integrated", "auto")
    if (is.character(alignment_config$method)) {
      if (!alignment_config$method %in% valid_methods) {
        validation_results$errors <- c(validation_results$errors,
                                     paste("Invalid alignment method:", alignment_config$method,
                                           "Valid methods:", paste(valid_methods, collapse = ", ")))
        validation_results$passed <- FALSE
      }
    } else if (is.list(alignment_config$method)) {
      # Adaptive method configuration
      for (rule in alignment_config$method) {
        if (!"condition" %in% names(rule) || !"value" %in% names(rule)) {
          validation_results$errors <- c(validation_results$errors,
                                       "Adaptive method rules must have 'condition' and 'value' fields")
          validation_results$passed <- FALSE
        }
      }
    }
  }
  
  # Validate numeric parameters
  numeric_params <- list(
    "center_window" = c(1, 1000),
    "consensus_threshold" = c(0, 1),
    "max_iterations" = c(1, 100),
    "convergence_threshold" = c(0, 1)
  )
  
  for (param in names(numeric_params)) {
    if (param %in% names(alignment_config)) {
      value <- alignment_config[[param]]
      if (param != "max_iterations" || value != "auto") {  # Allow "auto" for max_iterations
        range <- numeric_params[[param]]
        if (!is.numeric(value) || value < range[1] || value > range[2]) {
          validation_results$errors <- c(validation_results$errors,
                                       paste(param, "must be between", range[1], "and", range[2]))
          validation_results$passed <- FALSE
        }
      }
    }
  }
  
  return(validation_results)
}

#' Validate PWM section
validate_pwm_section <- function(pwm_config, validation_results) {
  
  # Validate pseudocount
  if ("pseudocount" %in% names(pwm_config)) {
    if (!is.numeric(pwm_config$pseudocount) || pwm_config$pseudocount < 0) {
      validation_results$errors <- c(validation_results$errors,
                                   "pseudocount must be a non-negative number")
      validation_results$passed <- FALSE
    }
  }
  
  # Validate background frequencies
  if ("background_freq" %in% names(pwm_config)) {
    bg_freq <- pwm_config$background_freq
    if (is.list(bg_freq)) {
      required_bases <- c("A", "C", "G", "T")
      for (base in required_bases) {
        if (!base %in% names(bg_freq)) {
          validation_results$errors <- c(validation_results$errors,
                                       paste("Missing background frequency for base:", base))
          validation_results$passed <- FALSE
        } else {
          freq <- bg_freq[[base]]
          if (!is.numeric(freq) || freq < 0 || freq > 1) {
            validation_results$errors <- c(validation_results$errors,
                                         paste("Background frequency for", base, "must be between 0 and 1"))
            validation_results$passed <- FALSE
          }
        }
      }
      
      # Check if frequencies sum to 1
      if (all(required_bases %in% names(bg_freq))) {
        total_freq <- sum(unlist(bg_freq[required_bases]))
        if (abs(total_freq - 1.0) > 0.01) {
          validation_results$warnings <- c(validation_results$warnings,
                                         paste("Background frequencies sum to", round(total_freq, 3), 
                                               "instead of 1.0"))
        }
      }
    }
  }
  
  # Validate information content parameters
  if ("information_content" %in% names(pwm_config)) {
    ic_config <- pwm_config$information_content
    
    if ("min_total_ic" %in% names(ic_config)) {
      if (!is.numeric(ic_config$min_total_ic) || ic_config$min_total_ic < 0) {
        validation_results$errors <- c(validation_results$errors,
                                     "min_total_ic must be a non-negative number")
        validation_results$passed <- FALSE
      }
    }
    
    if ("min_conserved_positions" %in% names(ic_config)) {
      if (!is.numeric(ic_config$min_conserved_positions) || ic_config$min_conserved_positions < 0) {
        validation_results$errors <- c(validation_results$errors,
                                     "min_conserved_positions must be a non-negative integer")
        validation_results$passed <- FALSE
      }
    }
  }
  
  return(validation_results)
}

#' Validate statistics section
validate_statistics_section <- function(stats_config, validation_results) {
  
  # Validate significance level
  if ("significance_level" %in% names(stats_config)) {
    if (!is.numeric(stats_config$significance_level) || 
        stats_config$significance_level <= 0 || stats_config$significance_level >= 1) {
      validation_results$errors <- c(validation_results$errors,
                                   "significance_level must be between 0 and 1")
      validation_results$passed <- FALSE
    }
  }
  
  # Validate correction method
  if ("multiple_testing_correction" %in% names(stats_config)) {
    valid_corrections <- c("bonferroni", "fdr", "holm", "none")
    if (!stats_config$multiple_testing_correction %in% valid_corrections) {
      validation_results$errors <- c(validation_results$errors,
                                   paste("Invalid multiple testing correction:",
                                         stats_config$multiple_testing_correction,
                                         "Valid methods:", paste(valid_corrections, collapse = ", ")))
      validation_results$passed <- FALSE
    }
  }
  
  # Validate iteration counts
  iteration_params <- c("bootstrap_iterations", "null_model_iterations")
  for (param in iteration_params) {
    if (param %in% names(stats_config)) {
      value <- stats_config[[param]]
      if (!is.numeric(value) || value < 1) {
        validation_results$errors <- c(validation_results$errors,
                                     paste(param, "must be a positive integer"))
        validation_results$passed <- FALSE
      }
    }
  }
  
  return(validation_results)
}

#' Validate performance section
validate_performance_section <- function(perf_config, validation_results) {
  
  # Validate threads
  if ("threads" %in% names(perf_config)) {
    threads <- perf_config$threads
    if (is.numeric(threads)) {
      if (threads < 1) {
        validation_results$errors <- c(validation_results$errors,
                                     "threads must be a positive integer")
        validation_results$passed <- FALSE
      }
    } else if (is.list(threads)) {
      # Adaptive threads configuration
      for (rule in threads) {
        if (!"condition" %in% names(rule) || !"value" %in% names(rule)) {
          validation_results$errors <- c(validation_results$errors,
                                       "Adaptive thread rules must have 'condition' and 'value' fields")
          validation_results$passed <- FALSE
        }
      }
    } else if (threads != "auto") {
      validation_results$errors <- c(validation_results$errors,
                                   "threads must be a positive integer, 'auto', or adaptive configuration")
      validation_results$passed <- FALSE
    }
  }
  
  # Validate memory limit
  if ("memory_limit" %in% names(perf_config)) {
    memory_limit <- perf_config$memory_limit
    if (memory_limit != "auto") {
      if (!grepl("^[0-9]+[GMK]?$", memory_limit)) {
        validation_results$errors <- c(validation_results$errors,
                                     "memory_limit must be in format like '8G', '4096M', or 'auto'")
        validation_results$passed <- FALSE
      }
    }
  }
  
  # Validate batch size
  if ("batch_size" %in% names(perf_config)) {
    batch_size <- perf_config$batch_size
    if (batch_size != "auto") {
      if (!is.numeric(batch_size) || batch_size < 1) {
        validation_results$errors <- c(validation_results$errors,
                                     "batch_size must be a positive integer or 'auto'")
        validation_results$passed <- FALSE
      }
    }
  }
  
  return(validation_results)
}

#' Validate output section
validate_output_section <- function(output_config, validation_results) {
  
  # Validate output formats
  if ("formats" %in% names(output_config)) {
    valid_formats <- c("meme", "jaspar", "transfac", "pfm", "pwm")
    formats <- output_config$formats
    
    if (is.character(formats)) {
      for (format in formats) {
        if (!format %in% valid_formats) {
          validation_results$errors <- c(validation_results$errors,
                                       paste("Invalid output format:", format,
                                             "Valid formats:", paste(valid_formats, collapse = ", ")))
          validation_results$passed <- FALSE
        }
      }
    }
  }
  
  return(validation_results)
}

#' Generate validation report
generate_report <- function(validation_results, output_file = NULL) {
  
  # Create report content
  report <- list(
    summary = list(
      file = validation_results$file,
      status = ifelse(validation_results$passed, "PASSED", "FAILED"),
      errors = length(validation_results$errors),
      warnings = length(validation_results$warnings)
    ),
    details = list(
      errors = validation_results$errors,
      warnings = validation_results$warnings
    ),
    timestamp = Sys.time()
  )
  
  # Output report
  if (!is.null(output_file)) {
    if (grepl("\\.json$", output_file)) {
      writeLines(toJSON(report, pretty = TRUE), output_file)
    } else {
      # Text format
      cat("=== Configuration Validation Report ===\n", file = output_file)
      cat("File:", validation_results$file, "\n", file = output_file, append = TRUE)
      cat("Status:", report$summary$status, "\n", file = output_file, append = TRUE)
      cat("Errors:", report$summary$errors, "\n", file = output_file, append = TRUE)
      cat("Warnings:", report$summary$warnings, "\n\n", file = output_file, append = TRUE)
      
      if (length(validation_results$errors) > 0) {
        cat("ERRORS:\n", file = output_file, append = TRUE)
        for (error in validation_results$errors) {
          cat("  -", error, "\n", file = output_file, append = TRUE)
        }
        cat("\n", file = output_file, append = TRUE)
      }
      
      if (length(validation_results$warnings) > 0) {
        cat("WARNINGS:\n", file = output_file, append = TRUE)
        for (warning in validation_results$warnings) {
          cat("  -", warning, "\n", file = output_file, append = TRUE)
        }
      }
    }
    cat("Report saved to:", output_file, "\n")
  }
  
  return(report)
}

#' Main function
main <- function() {
  args <- parse_arguments()
  
  cat("Validating configuration file:", args$config, "\n")
  
  # Load and validate configuration
  config <- load_config(args$config)
  validation_results <- validate_configuration(config, args$config)
  
  # Generate report
  report <- generate_report(validation_results, args$output)
  
  # Console output
  cat("\n=== Validation Results ===\n")
  cat("Status:", report$summary$status, "\n")
  cat("Errors:", report$summary$errors, "\n")
  cat("Warnings:", report$summary$warnings, "\n")
  
  if (args$verbose) {
    if (length(validation_results$errors) > 0) {
      cat("\nERRORS:\n")
      for (error in validation_results$errors) {
        cat("  -", error, "\n")
      }
    }
    
    if (length(validation_results$warnings) > 0) {
      cat("\nWARNINGS:\n")
      for (warning in validation_results$warnings) {
        cat("  -", warning, "\n")
      }
    }
  }
  
  # Exit with appropriate code
  if (validation_results$passed) {
    cat("\nConfiguration validation PASSED\n")
    quit(status = 0)
  } else {
    cat("\nConfiguration validation FAILED\n")
    quit(status = 1)
  }
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
