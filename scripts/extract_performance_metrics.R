#!/usr/bin/env Rscript

# Extract Performance Metrics
# Extracts and summarizes performance metrics from analysis results
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(data.table)
})

#' Main function for extracting performance metrics
#' @param results_dir Directory containing analysis results
#' @param output_file Output file for metrics summary
#' @param format Output format (csv, json, rds)
#' @param include_details Include detailed metrics
#' @param verbose Enable verbose output
extract_performance_metrics <- function(results_dir = "results", 
                                       output_file = NULL,
                                       format = "csv", 
                                       include_details = TRUE,
                                       verbose = FALSE) {
  
  if (verbose) cat("Extracting performance metrics...\n")
  
  # Initialize metrics collection
  metrics <- list(
    summary = list(),
    pwm_metrics = list(),
    method_comparison = list(),
    statistical_results = list(),
    quality_assessment = list(),
    processing_times = list()
  )
  
  # Extract PWM metrics
  if (verbose) cat("Extracting PWM metrics...\n")
  metrics$pwm_metrics <- extract_pwm_metrics(results_dir, verbose)
  
  # Extract method comparison metrics
  if (verbose) cat("Extracting method comparison metrics...\n")
  metrics$method_comparison <- extract_method_comparison_metrics(results_dir, verbose)
  
  # Extract statistical results
  if (verbose) cat("Extracting statistical results...\n")
  metrics$statistical_results <- extract_statistical_metrics(results_dir, verbose)
  
  # Extract quality assessment
  if (verbose) cat("Extracting quality assessment...\n")
  metrics$quality_assessment <- extract_quality_metrics(results_dir, verbose)
  
  # Extract processing times
  if (verbose) cat("Extracting processing times...\n")
  metrics$processing_times <- extract_processing_times(results_dir, verbose)
  
  # Generate summary
  metrics$summary <- generate_metrics_summary(metrics, verbose)
  
  # Save results
  if (!is.null(output_file)) {
    save_metrics(metrics, output_file, format, include_details, verbose)
  }
  
  # Display summary
  display_metrics_summary(metrics$summary, verbose)
  
  return(metrics)
}

#' Extract PWM-specific metrics
extract_pwm_metrics <- function(results_dir, verbose) {
  
  pwm_files <- list.files(results_dir, pattern = ".*pwm.*\\.rds$", full.names = TRUE)
  
  if (length(pwm_files) == 0) {
    if (verbose) cat("  No PWM files found\n")
    return(data.frame())
  }
  
  pwm_metrics <- data.frame()
  
  for (pwm_file in pwm_files) {
    tryCatch({
      pwm_result <- readRDS(pwm_file)
      
      # Extract method name
      method_name <- extract_method_name(pwm_file)
      
      # Extract basic metrics
      metrics_row <- data.frame(
        Method = method_name,
        File = basename(pwm_file),
        stringsAsFactors = FALSE
      )
      
      # Information content metrics
      if ("total_info" %in% names(pwm_result)) {
        metrics_row$Total_Information_Content <- pwm_result$total_info
      } else if ("info_content" %in% names(pwm_result)) {
        metrics_row$Total_Information_Content <- sum(pwm_result$info_content, na.rm = TRUE)
      } else {
        metrics_row$Total_Information_Content <- NA
      }
      
      # Position-wise metrics
      if ("info_content" %in% names(pwm_result) && is.vector(pwm_result$info_content)) {
        ic_vector <- pwm_result$info_content
        metrics_row$Max_Position_IC <- max(ic_vector, na.rm = TRUE)
        metrics_row$Mean_Position_IC <- mean(ic_vector, na.rm = TRUE)
        metrics_row$SD_Position_IC <- sd(ic_vector, na.rm = TRUE)
        metrics_row$High_IC_Positions <- sum(ic_vector > 1.5, na.rm = TRUE)
        metrics_row$Conserved_Positions <- sum(ic_vector > 1.0, na.rm = TRUE)
        metrics_row$PWM_Length <- length(ic_vector)
      }
      
      # Sequence count
      if ("n_sequences" %in% names(pwm_result)) {
        metrics_row$N_Sequences <- pwm_result$n_sequences
      }
      
      # Quality grade
      if ("quality_grade" %in% names(pwm_result)) {
        metrics_row$Quality_Grade <- pwm_result$quality_grade
      }
      
      # Consensus sequence
      if ("consensus" %in% names(pwm_result)) {
        metrics_row$Consensus_Sequence <- pwm_result$consensus
      }
      
      # Statistical significance
      if ("p_value" %in% names(pwm_result)) {
        metrics_row$P_Value <- pwm_result$p_value
        metrics_row$Significant <- pwm_result$p_value < 0.05
      }
      
      # Add to results
      pwm_metrics <- rbind(pwm_metrics, metrics_row)
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(pwm_file), ":", conditionMessage(e), "\n")
    })
  }
  
  return(pwm_metrics)
}

#' Extract method comparison metrics
extract_method_comparison_metrics <- function(results_dir, verbose) {
  
  comparison_files <- list.files(results_dir, pattern = "comparison.*\\.rds$|.*comparison.*\\.rds$", 
                                full.names = TRUE)
  
  if (length(comparison_files) == 0) {
    if (verbose) cat("  No comparison files found\n")
    return(data.frame())
  }
  
  comparison_metrics <- data.frame()
  
  for (comp_file in comparison_files) {
    tryCatch({
      comp_result <- readRDS(comp_file)
      
      if ("method_scores" %in% names(comp_result)) {
        for (method_name in names(comp_result$method_scores)) {
          method_score <- comp_result$method_scores[[method_name]]
          
          metrics_row <- data.frame(
            Method = method_name,
            Weighted_Score = method_score$weighted_score,
            IC_Component = method_score$component_scores$information_content,
            Conservation_Component = method_score$component_scores$conservation,
            Alignment_Component = method_score$component_scores$alignment_quality,
            Time_Component = method_score$component_scores$processing_time,
            Meets_Thresholds = method_score$meets_thresholds,
            Processing_Time = method_score$result$processing_time %||% NA,
            N_Sequences = method_score$result$n_sequences %||% NA,
            stringsAsFactors = FALSE
          )
          
          comparison_metrics <- rbind(comparison_metrics, metrics_row)
        }
      }
      
      # Add ranking information
      if ("ranking" %in% names(comp_result)) {
        ranking_df <- data.frame(
          Method = comp_result$ranking,
          Rank = seq_along(comp_result$ranking),
          stringsAsFactors = FALSE
        )
        
        comparison_metrics <- merge(comparison_metrics, ranking_df, by = "Method", all.x = TRUE)
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(comp_file), ":", conditionMessage(e), "\n")
    })
  }
  
  return(comparison_metrics)
}

#' Extract statistical metrics
extract_statistical_metrics <- function(results_dir, verbose) {
  
  stat_files <- list.files(results_dir, pattern = "statistical.*\\.rds$|.*stats.*\\.rds$", 
                          full.names = TRUE)
  
  if (length(stat_files) == 0) {
    if (verbose) cat("  No statistical files found\n")
    return(data.frame())
  }
  
  stat_metrics <- data.frame()
  
  for (stat_file in stat_files) {
    tryCatch({
      stat_result <- readRDS(stat_file)
      
      method_name <- extract_method_name(stat_file)
      
      metrics_row <- data.frame(
        Method = method_name,
        File = basename(stat_file),
        stringsAsFactors = FALSE
      )
      
      # P-value and significance
      if ("p_value" %in% names(stat_result)) {
        metrics_row$P_Value <- stat_result$p_value
        metrics_row$Significant <- stat_result$p_value < 0.05
        metrics_row$Log10_P_Value <- -log10(stat_result$p_value)
      }
      
      # Effect size
      if ("effect_size" %in% names(stat_result)) {
        metrics_row$Effect_Size <- stat_result$effect_size
      }
      
      # Confidence intervals
      if ("confidence_interval" %in% names(stat_result)) {
        ci <- stat_result$confidence_interval
        if (length(ci) == 2) {
          metrics_row$CI_Lower <- ci[1]
          metrics_row$CI_Upper <- ci[2]
        }
      }
      
      # Statistical test used
      if ("test_method" %in% names(stat_result)) {
        metrics_row$Test_Method <- stat_result$test_method
      }
      
      # Cross-validation results
      if ("cv_results" %in% names(stat_result)) {
        cv <- stat_result$cv_results
        if ("mean_auc" %in% names(cv)) {
          metrics_row$CV_Mean_AUC <- cv$mean_auc
          metrics_row$CV_SD_AUC <- cv$sd_auc %||% NA
        }
        if ("mean_accuracy" %in% names(cv)) {
          metrics_row$CV_Mean_Accuracy <- cv$mean_accuracy
          metrics_row$CV_SD_Accuracy <- cv$sd_accuracy %||% NA
        }
      }
      
      stat_metrics <- rbind(stat_metrics, metrics_row)
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(stat_file), ":", conditionMessage(e), "\n")
    })
  }
  
  return(stat_metrics)
}

#' Extract quality assessment metrics
extract_quality_metrics <- function(results_dir, verbose) {
  
  quality_files <- list.files(results_dir, pattern = "quality.*\\.rds$|.*quality.*\\.rds$", 
                             full.names = TRUE)
  
  if (length(quality_files) == 0) {
    if (verbose) cat("  No quality files found\n")
    return(data.frame())
  }
  
  quality_metrics <- data.frame()
  
  for (quality_file in quality_files) {
    tryCatch({
      quality_result <- readRDS(quality_file)
      
      method_name <- extract_method_name(quality_file)
      
      metrics_row <- data.frame(
        Method = method_name,
        File = basename(quality_file),
        stringsAsFactors = FALSE
      )
      
      # Assessment scores
      if ("assessment" %in% names(quality_result)) {
        assessment <- quality_result$assessment
        
        metrics_row$Information_Content_Score <- assessment$information_content %||% NA
        metrics_row$Conservation_Score <- assessment$conservation_score %||% NA
        metrics_row$Pattern_Quality <- assessment$pattern_quality %||% NA
        metrics_row$Overall_Grade <- assessment$overall_grade %||% "Unknown"
        metrics_row$Quality_Score <- assessment$quality_score %||% NA
      }
      
      # Detailed quality metrics
      if ("details" %in% names(quality_result)) {
        details <- quality_result$details
        
        if ("motif_strength" %in% names(details)) {
          metrics_row$Motif_Strength <- details$motif_strength
        }
        
        if ("positional_conservation" %in% names(details)) {
          metrics_row$Positional_Conservation <- details$positional_conservation
        }
        
        if ("background_discrimination" %in% names(details)) {
          metrics_row$Background_Discrimination <- details$background_discrimination
        }
      }
      
      quality_metrics <- rbind(quality_metrics, metrics_row)
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(quality_file), ":", conditionMessage(e), "\n")
    })
  }
  
  return(quality_metrics)
}

#' Extract processing time metrics
extract_processing_times <- function(results_dir, verbose) {
  
  # Look for log files and metadata files
  log_files <- list.files(results_dir, pattern = "\\.log$", full.names = TRUE)
  metadata_files <- list.files(results_dir, pattern = "metadata.*\\.json$", full.names = TRUE)
  
  time_metrics <- data.frame()
  
  # Extract from metadata files
  for (metadata_file in metadata_files) {
    tryCatch({
      metadata <- fromJSON(metadata_file)
      
      if ("duration_seconds" %in% names(metadata)) {
        method_name <- gsub("metadata_|_metadata", "", basename(metadata_file))
        method_name <- gsub("\\.json$", "", method_name)
        
        time_row <- data.frame(
          Method = method_name,
          Processing_Time_Seconds = metadata$duration_seconds,
          Processing_Time_Minutes = metadata$duration_seconds / 60,
          Start_Time = metadata$start_time %||% NA,
          End_Time = metadata$end_time %||% NA,
          Source = "metadata",
          stringsAsFactors = FALSE
        )
        
        time_metrics <- rbind(time_metrics, time_row)
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(metadata_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Extract from log files (simple pattern matching)
  for (log_file in log_files) {
    tryCatch({
      log_content <- readLines(log_file, warn = FALSE)
      
      # Look for timing patterns
      time_patterns <- grep("elapsed|duration|took|completed in", log_content, 
                           ignore.case = TRUE, value = TRUE)
      
      if (length(time_patterns) > 0) {
        method_name <- gsub("\\.log$", "", basename(log_file))
        
        # Extract numeric time values (simple heuristic)
        time_values <- regmatches(time_patterns, 
                                gregexpr("\\d+\\.?\\d*\\s*(second|minute|hour)", time_patterns))
        
        if (length(time_values) > 0 && length(time_values[[1]]) > 0) {
          time_str <- time_values[[1]][1]
          
          time_row <- data.frame(
            Method = method_name,
            Time_String = time_str,
            Source = "log",
            stringsAsFactors = FALSE
          )
          
          time_metrics <- rbind(time_metrics, time_row)
        }
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(log_file), ":", conditionMessage(e), "\n")
    })
  }
  
  return(time_metrics)
}

#' Generate metrics summary
generate_metrics_summary <- function(metrics, verbose) {
  
  summary <- list(
    timestamp = Sys.time(),
    total_methods = 0,
    best_performing_method = "Unknown",
    metrics_overview = list()
  )
  
  # PWM metrics summary
  if (nrow(metrics$pwm_metrics) > 0) {
    pwm_data <- metrics$pwm_metrics
    
    summary$total_methods <- nrow(pwm_data)
    
    # Best method by information content
    if ("Total_Information_Content" %in% names(pwm_data)) {
      best_ic_idx <- which.max(pwm_data$Total_Information_Content)
      if (length(best_ic_idx) > 0) {
        summary$best_performing_method <- pwm_data$Method[best_ic_idx]
        summary$best_ic_value <- pwm_data$Total_Information_Content[best_ic_idx]
      }
    }
    
    # IC statistics
    if ("Total_Information_Content" %in% names(pwm_data)) {
      summary$metrics_overview$information_content <- list(
        mean = mean(pwm_data$Total_Information_Content, na.rm = TRUE),
        median = median(pwm_data$Total_Information_Content, na.rm = TRUE),
        sd = sd(pwm_data$Total_Information_Content, na.rm = TRUE),
        min = min(pwm_data$Total_Information_Content, na.rm = TRUE),
        max = max(pwm_data$Total_Information_Content, na.rm = TRUE),
        n_methods = sum(!is.na(pwm_data$Total_Information_Content))
      )
    }
    
    # Quality grade distribution
    if ("Quality_Grade" %in% names(pwm_data)) {
      grade_table <- table(pwm_data$Quality_Grade)
      summary$metrics_overview$quality_grades <- as.list(grade_table)
    }
    
    # Significance statistics
    if ("Significant" %in% names(pwm_data)) {
      n_significant <- sum(pwm_data$Significant, na.rm = TRUE)
      summary$metrics_overview$significance <- list(
        n_significant = n_significant,
        n_total = sum(!is.na(pwm_data$Significant)),
        proportion_significant = n_significant / sum(!is.na(pwm_data$Significant))
      )
    }
  }
  
  # Method comparison summary
  if (nrow(metrics$method_comparison) > 0) {
    comp_data <- metrics$method_comparison
    
    # Best method by weighted score
    if ("Weighted_Score" %in% names(comp_data)) {
      best_score_idx <- which.max(comp_data$Weighted_Score)
      if (length(best_score_idx) > 0) {
        summary$best_method_by_score <- comp_data$Method[best_score_idx]
        summary$best_weighted_score <- comp_data$Weighted_Score[best_score_idx]
      }
    }
    
    # Methods meeting thresholds
    if ("Meets_Thresholds" %in% names(comp_data)) {
      n_meets_thresholds <- sum(comp_data$Meets_Thresholds, na.rm = TRUE)
      summary$methods_meeting_thresholds <- n_meets_thresholds
    }
  }
  
  # Processing time summary
  if (nrow(metrics$processing_times) > 0) {
    time_data <- metrics$processing_times
    
    if ("Processing_Time_Seconds" %in% names(time_data)) {
      summary$metrics_overview$processing_times <- list(
        mean_seconds = mean(time_data$Processing_Time_Seconds, na.rm = TRUE),
        median_seconds = median(time_data$Processing_Time_Seconds, na.rm = TRUE),
        min_seconds = min(time_data$Processing_Time_Seconds, na.rm = TRUE),
        max_seconds = max(time_data$Processing_Time_Seconds, na.rm = TRUE),
        total_methods_timed = sum(!is.na(time_data$Processing_Time_Seconds))
      )
    }
  }
  
  return(summary)
}

#' Save metrics to file
save_metrics <- function(metrics, output_file, format, include_details, verbose) {
  
  if (format == "csv") {
    # Save each metrics table as separate CSV
    base_name <- gsub("\\.[^.]*$", "", output_file)
    
    if (nrow(metrics$pwm_metrics) > 0) {
      write.csv(metrics$pwm_metrics, paste0(base_name, "_pwm_metrics.csv"), row.names = FALSE)
    }
    
    if (nrow(metrics$method_comparison) > 0) {
      write.csv(metrics$method_comparison, paste0(base_name, "_method_comparison.csv"), row.names = FALSE)
    }
    
    if (nrow(metrics$statistical_results) > 0) {
      write.csv(metrics$statistical_results, paste0(base_name, "_statistical_results.csv"), row.names = FALSE)
    }
    
    if (nrow(metrics$quality_assessment) > 0) {
      write.csv(metrics$quality_assessment, paste0(base_name, "_quality_assessment.csv"), row.names = FALSE)
    }
    
    if (nrow(metrics$processing_times) > 0) {
      write.csv(metrics$processing_times, paste0(base_name, "_processing_times.csv"), row.names = FALSE)
    }
    
    # Save summary as JSON
    writeLines(toJSON(metrics$summary, pretty = TRUE), paste0(base_name, "_summary.json"))
    
  } else if (format == "json") {
    writeLines(toJSON(metrics, pretty = TRUE), output_file)
    
  } else if (format == "rds") {
    saveRDS(metrics, output_file)
    
  } else {
    stop("Unsupported format: ", format)
  }
  
  if (verbose) cat("Metrics saved to:", output_file, "\n")
}

#' Display metrics summary
display_metrics_summary <- function(summary, verbose) {
  
  cat("\n=== Performance Metrics Summary ===\n")
  cat("Timestamp:", as.character(summary$timestamp), "\n")
  cat("Total Methods Analyzed:", summary$total_methods, "\n")
  
  if (!is.null(summary$best_performing_method)) {
    cat("Best Performing Method:", summary$best_performing_method)
    if (!is.null(summary$best_ic_value)) {
      cat(" (IC:", round(summary$best_ic_value, 3), "bits)")
    }
    cat("\n")
  }
  
  # Information content overview
  if (!is.null(summary$metrics_overview$information_content)) {
    ic_stats <- summary$metrics_overview$information_content
    cat("\nInformation Content Statistics:\n")
    cat("  Mean:", round(ic_stats$mean, 3), "bits\n")
    cat("  Median:", round(ic_stats$median, 3), "bits\n")
    cat("  Range:", round(ic_stats$min, 3), "-", round(ic_stats$max, 3), "bits\n")
    cat("  Methods with IC data:", ic_stats$n_methods, "\n")
  }
  
  # Quality grades
  if (!is.null(summary$metrics_overview$quality_grades)) {
    cat("\nQuality Grade Distribution:\n")
    for (grade in names(summary$metrics_overview$quality_grades)) {
      cat("  ", grade, ":", summary$metrics_overview$quality_grades[[grade]], "\n")
    }
  }
  
  # Statistical significance
  if (!is.null(summary$metrics_overview$significance)) {
    sig_stats <- summary$metrics_overview$significance
    cat("\nStatistical Significance:\n")
    cat("  Significant methods:", sig_stats$n_significant, "of", sig_stats$n_total, "\n")
    cat("  Proportion significant:", round(sig_stats$proportion_significant * 100, 1), "%\n")
  }
  
  # Processing times
  if (!is.null(summary$metrics_overview$processing_times)) {
    time_stats <- summary$metrics_overview$processing_times
    cat("\nProcessing Time Statistics:\n")
    cat("  Mean:", round(time_stats$mean_seconds / 60, 2), "minutes\n")
    cat("  Median:", round(time_stats$median_seconds / 60, 2), "minutes\n")
    cat("  Range:", round(time_stats$min_seconds / 60, 2), "-", 
        round(time_stats$max_seconds / 60, 2), "minutes\n")
  }
  
  # Method comparison
  if (!is.null(summary$best_method_by_score)) {
    cat("\nBest Method by Weighted Score:", summary$best_method_by_score)
    if (!is.null(summary$best_weighted_score)) {
      cat(" (Score:", round(summary$best_weighted_score, 3), ")")
    }
    cat("\n")
  }
  
  if (!is.null(summary$methods_meeting_thresholds)) {
    cat("Methods Meeting Quality Thresholds:", summary$methods_meeting_thresholds, "\n")
  }
}

#' Extract method name from filename
extract_method_name <- function(filename) {
  base_name <- basename(filename)
  method_name <- gsub("\\.[^.]*$", "", base_name)  # Remove extension
  method_name <- gsub("_pwm|pwm_|_results|results_", "", method_name)  # Remove common suffixes
  return(method_name)
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-r", "--results"), type = "character", default = "results",
                help = "Results directory [default: %default]", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file for metrics", metavar = "character"),
    make_option(c("-f", "--format"), type = "character", default = "csv",
                help = "Output format (csv, json, rds) [default: %default]", metavar = "character"),
    make_option(c("-d", "--details"), action = "store_true", default = FALSE,
                help = "Include detailed metrics"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Extract Performance Metrics from Analysis Results")
  opt <- parse_args(opt_parser)
  
  if (!dir.exists(opt$results)) {
    stop("Results directory does not exist: ", opt$results, call. = FALSE)
  }
  
  # Generate default output filename if not provided
  if (is.null(opt$output)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    opt$output <- file.path(opt$results, paste0("performance_metrics_", timestamp))
  }
  
  # Extract performance metrics
  tryCatch({
    metrics <- extract_performance_metrics(
      results_dir = opt$results,
      output_file = opt$output,
      format = opt$format,
      include_details = opt$details,
      verbose = opt$verbose
    )
    
    cat("Performance metrics extracted successfully.\n")
    
  }, error = function(e) {
    cat("Error extracting performance metrics:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
