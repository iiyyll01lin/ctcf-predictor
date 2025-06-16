#!/usr/bin/env Rscript

# Generate Supplementary Data
# Creates supplementary data files for research publications
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
  library(jsonlite)
  library(Biostrings)
})

#' Main function for generating supplementary data
#' @param results_dir Directory containing analysis results
#' @param output_dir Output directory for supplementary files
#' @param include_raw Include raw data files
#' @param format Output format (xlsx, csv, json)
#' @param verbose Enable verbose output
generate_supplementary_data <- function(results_dir = "results", 
                                       output_dir = "results/supplementary",
                                       include_raw = FALSE,
                                       format = "xlsx",
                                       verbose = FALSE) {
  
  if (verbose) cat("Generating supplementary data...\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Collect all supplementary data
  supplementary_data <- list()
  
  # Table S1: PWM Information Content Data
  if (verbose) cat("Generating Table S1: PWM Information Content Data...\n")
  supplementary_data$table_s1 <- generate_pwm_ic_table(results_dir, verbose)
  
  # Table S2: Method Comparison Results
  if (verbose) cat("Generating Table S2: Method Comparison Results...\n")
  supplementary_data$table_s2 <- generate_method_comparison_table(results_dir, verbose)
  
  # Table S3: Statistical Validation Results
  if (verbose) cat("Generating Table S3: Statistical Validation Results...\n")
  supplementary_data$table_s3 <- generate_statistical_table(results_dir, verbose)
  
  # Table S4: Quality Assessment Scores
  if (verbose) cat("Generating Table S4: Quality Assessment Scores...\n")
  supplementary_data$table_s4 <- generate_quality_table(results_dir, verbose)
  
  # Table S5: Sequence Statistics
  if (verbose) cat("Generating Table S5: Sequence Statistics...\n")
  supplementary_data$table_s5 <- generate_sequence_stats_table(results_dir, verbose)
  
  # Supplementary Files
  if (verbose) cat("Generating supplementary files...\n")
  supplementary_files <- generate_supplementary_files(results_dir, output_dir, include_raw, verbose)
  
  # Save supplementary data
  save_supplementary_data(supplementary_data, supplementary_files, output_dir, format, verbose)
  
  # Generate README
  generate_supplementary_readme(supplementary_data, supplementary_files, output_dir, verbose)
  
  if (verbose) cat("Supplementary data generated successfully in:", output_dir, "\n")
  
  return(list(
    tables = supplementary_data,
    files = supplementary_files
  ))
}

#' Generate PWM Information Content Table (Table S1)
generate_pwm_ic_table <- function(results_dir, verbose) {
  
  pwm_files <- list.files(results_dir, pattern = ".*pwm.*\\.rds$", full.names = TRUE)
  
  if (length(pwm_files) == 0) {
    if (verbose) cat("  No PWM files found\n")
    return(data.frame())
  }
  
  ic_table <- data.frame()
  
  for (pwm_file in pwm_files) {
    tryCatch({
      pwm_result <- readRDS(pwm_file)
      
      # Extract method name
      method_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", pwm_file)
      method_name <- gsub("_pwm|pwm_", "", method_name)
      
      # Basic information
      row_data <- data.frame(
        Method = method_name,
        File = basename(pwm_file),
        stringsAsFactors = FALSE
      )
      
      # Total information content
      if ("total_info" %in% names(pwm_result)) {
        row_data$Total_IC_bits <- round(pwm_result$total_info, 4)
      } else if ("info_content" %in% names(pwm_result)) {
        row_data$Total_IC_bits <- round(sum(pwm_result$info_content, na.rm = TRUE), 4)
      }
      
      # Position-wise information content
      if ("info_content" %in% names(pwm_result) && is.vector(pwm_result$info_content)) {
        ic_vector <- pwm_result$info_content
        
        row_data$PWM_Length <- length(ic_vector)
        row_data$Max_Position_IC <- round(max(ic_vector, na.rm = TRUE), 4)
        row_data$Mean_Position_IC <- round(mean(ic_vector, na.rm = TRUE), 4)
        row_data$Median_Position_IC <- round(median(ic_vector, na.rm = TRUE), 4)
        row_data$SD_Position_IC <- round(sd(ic_vector, na.rm = TRUE), 4)
        
        # Count high IC positions
        row_data$Positions_IC_gt_1 <- sum(ic_vector > 1.0, na.rm = TRUE)
        row_data$Positions_IC_gt_1_5 <- sum(ic_vector > 1.5, na.rm = TRUE)
        row_data$Positions_IC_gt_2 <- sum(ic_vector > 2.0, na.rm = TRUE)
        
        # Position-wise IC values (for detailed table)
        for (i in seq_along(ic_vector)) {
          col_name <- paste0("Position_", i, "_IC")
          row_data[[col_name]] <- round(ic_vector[i], 4)
        }
      }
      
      # Sequence information
      if ("n_sequences" %in% names(pwm_result)) {
        row_data$N_Sequences <- pwm_result$n_sequences
      }
      
      # Consensus sequence
      if ("consensus" %in% names(pwm_result)) {
        row_data$Consensus_Sequence <- pwm_result$consensus
      }
      
      # Quality grade
      if ("quality_grade" %in% names(pwm_result)) {
        row_data$Quality_Grade <- pwm_result$quality_grade
      }
      
      ic_table <- rbind(ic_table, row_data)
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(pwm_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Sort by total IC descending
  if ("Total_IC_bits" %in% names(ic_table)) {
    ic_table <- ic_table[order(-ic_table$Total_IC_bits), ]
    rownames(ic_table) <- NULL
  }
  
  return(ic_table)
}

#' Generate Method Comparison Table (Table S2)
generate_method_comparison_table <- function(results_dir, verbose) {
  
  comparison_files <- list.files(results_dir, pattern = "comparison.*\\.rds$|.*comparison.*\\.rds$", 
                                full.names = TRUE)
  
  if (length(comparison_files) == 0) {
    if (verbose) cat("  No comparison files found\n")
    return(data.frame())
  }
  
  comparison_table <- data.frame()
  
  for (comp_file in comparison_files) {
    tryCatch({
      comp_result <- readRDS(comp_file)
      
      if ("method_scores" %in% names(comp_result)) {
        for (method_name in names(comp_result$method_scores)) {
          method_score <- comp_result$method_scores[[method_name]]
          
          # Basic comparison data
          row_data <- data.frame(
            Method = method_name,
            Weighted_Score = round(method_score$weighted_score, 4),
            Meets_Quality_Thresholds = method_score$meets_thresholds,
            stringsAsFactors = FALSE
          )
          
          # Component scores
          if ("component_scores" %in% names(method_score)) {
            components <- method_score$component_scores
            row_data$IC_Component_Score <- round(components$information_content, 4)
            row_data$Conservation_Component_Score <- round(components$conservation, 4)
            row_data$Alignment_Component_Score <- round(components$alignment_quality, 4)
            row_data$Time_Component_Score <- round(components$processing_time, 4)
          }
          
          # Method result details
          if ("result" %in% names(method_score)) {
            result <- method_score$result
            
            row_data$N_Sequences <- result$n_sequences %||% NA
            row_data$Processing_Time_seconds <- result$processing_time %||% NA
            row_data$Processing_Time_minutes <- ifelse(!is.na(result$processing_time %||% NA), 
                                                      round((result$processing_time %||% NA) / 60, 2), NA)
            
            # Quality metrics from result
            if ("quality" %in% names(result)) {
              quality <- result$quality
              row_data$Total_Information_Content <- round(quality$information_content, 4)
              row_data$Conservation_Ratio <- round(quality$conservation_ratio, 4)
              row_data$Alignment_Quality <- round(quality$alignment_quality, 4)
              row_data$Pattern_Score <- round(quality$pattern_score, 4)
            }
          }
          
          comparison_table <- rbind(comparison_table, row_data)
        }
      }
      
      # Add ranking information
      if ("ranking" %in% names(comp_result)) {
        ranking_df <- data.frame(
          Method = comp_result$ranking,
          Overall_Rank = seq_along(comp_result$ranking),
          stringsAsFactors = FALSE
        )
        
        comparison_table <- merge(comparison_table, ranking_df, by = "Method", all.x = TRUE)
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(comp_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Sort by weighted score descending
  if ("Weighted_Score" %in% names(comparison_table)) {
    comparison_table <- comparison_table[order(-comparison_table$Weighted_Score), ]
    rownames(comparison_table) <- NULL
  }
  
  return(comparison_table)
}

#' Generate Statistical Validation Table (Table S3)
generate_statistical_table <- function(results_dir, verbose) {
  
  stat_files <- list.files(results_dir, pattern = "statistical.*\\.rds$|.*stats.*\\.rds$", 
                          full.names = TRUE)
  
  if (length(stat_files) == 0) {
    if (verbose) cat("  No statistical files found\n")
    return(data.frame())
  }
  
  stat_table <- data.frame()
  
  for (stat_file in stat_files) {
    tryCatch({
      stat_result <- readRDS(stat_file)
      
      method_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", stat_file)
      method_name <- gsub("statistical_|_stats|_statistical", "", method_name)
      
      row_data <- data.frame(
        Method = method_name,
        File = basename(stat_file),
        stringsAsFactors = FALSE
      )
      
      # Statistical test results
      if ("p_value" %in% names(stat_result)) {
        row_data$P_Value <- format(stat_result$p_value, scientific = TRUE)
        row_data$P_Value_Numeric <- stat_result$p_value
        row_data$Significant_p_lt_0_05 <- stat_result$p_value < 0.05
        row_data$Significant_p_lt_0_01 <- stat_result$p_value < 0.01
        row_data$Significant_p_lt_0_001 <- stat_result$p_value < 0.001
        row_data$Neg_Log10_P_Value <- round(-log10(stat_result$p_value), 4)
      }
      
      # Test method and parameters
      if ("test_method" %in% names(stat_result)) {
        row_data$Statistical_Test <- stat_result$test_method
      }
      
      if ("null_model" %in% names(stat_result)) {
        row_data$Null_Model <- stat_result$null_model
      }
      
      # Effect size
      if ("effect_size" %in% names(stat_result)) {
        row_data$Effect_Size <- round(stat_result$effect_size, 4)
      }
      
      # Confidence intervals
      if ("confidence_interval" %in% names(stat_result)) {
        ci <- stat_result$confidence_interval
        if (length(ci) == 2) {
          row_data$CI_Lower_95 <- round(ci[1], 4)
          row_data$CI_Upper_95 <- round(ci[2], 4)
          row_data$CI_Width <- round(ci[2] - ci[1], 4)
        }
      }
      
      # Cross-validation results
      if ("cv_results" %in% names(stat_result)) {
        cv <- stat_result$cv_results
        if ("mean_auc" %in% names(cv)) {
          row_data$CV_Mean_AUC <- round(cv$mean_auc, 4)
          row_data$CV_SD_AUC <- round(cv$sd_auc %||% NA, 4)
        }
        if ("mean_accuracy" %in% names(cv)) {
          row_data$CV_Mean_Accuracy <- round(cv$mean_accuracy, 4)
          row_data$CV_SD_Accuracy <- round(cv$sd_accuracy %||% NA, 4)
        }
        if ("n_folds" %in% names(cv)) {
          row_data$CV_N_Folds <- cv$n_folds
        }
      }
      
      # Performance metrics
      if ("performance" %in% names(stat_result)) {
        perf <- stat_result$performance
        row_data$Sensitivity <- round(perf$sensitivity %||% NA, 4)
        row_data$Specificity <- round(perf$specificity %||% NA, 4)
        row_data$Precision <- round(perf$precision %||% NA, 4)
        row_data$F1_Score <- round(perf$f1_score %||% NA, 4)
      }
      
      stat_table <- rbind(stat_table, row_data)
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(stat_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Sort by significance (most significant first)
  if ("P_Value_Numeric" %in% names(stat_table)) {
    stat_table <- stat_table[order(stat_table$P_Value_Numeric), ]
    rownames(stat_table) <- NULL
  }
  
  return(stat_table)
}

#' Generate Quality Assessment Table (Table S4)
generate_quality_table <- function(results_dir, verbose) {
  
  quality_files <- list.files(results_dir, pattern = "quality.*\\.rds$|.*quality.*\\.rds$", 
                             full.names = TRUE)
  
  if (length(quality_files) == 0) {
    if (verbose) cat("  No quality files found\n")
    return(data.frame())
  }
  
  quality_table <- data.frame()
  
  for (quality_file in quality_files) {
    tryCatch({
      quality_result <- readRDS(quality_file)
      
      method_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", quality_file)
      method_name <- gsub("quality_|_quality", "", method_name)
      
      row_data <- data.frame(
        Method = method_name,
        File = basename(quality_file),
        stringsAsFactors = FALSE
      )
      
      # Assessment scores
      if ("assessment" %in% names(quality_result)) {
        assessment <- quality_result$assessment
        
        row_data$Information_Content_Score <- round(assessment$information_content %||% NA, 4)
        row_data$Conservation_Score <- round(assessment$conservation_score %||% NA, 4)
        row_data$Pattern_Quality_Score <- round(assessment$pattern_quality %||% NA, 4)
        row_data$Overall_Quality_Score <- round(assessment$quality_score %||% NA, 4)
        row_data$Quality_Grade <- assessment$overall_grade %||% "Unknown"
      }
      
      # Detailed quality metrics
      if ("details" %in% names(quality_result)) {
        details <- quality_result$details
        
        row_data$Motif_Strength <- round(details$motif_strength %||% NA, 4)
        row_data$Positional_Conservation <- round(details$positional_conservation %||% NA, 4)
        row_data$Background_Discrimination <- round(details$background_discrimination %||% NA, 4)
        row_data$Sequence_Complexity <- round(details$sequence_complexity %||% NA, 4)
        row_data$Alignment_Consistency <- round(details$alignment_consistency %||% NA, 4)
      }
      
      # Quality flags and warnings
      if ("flags" %in% names(quality_result)) {
        flags <- quality_result$flags
        row_data$Low_IC_Warning <- "low_information_content" %in% flags
        row_data$Poor_Conservation_Warning <- "poor_conservation" %in% flags
        row_data$Alignment_Issues_Warning <- "alignment_issues" %in% flags
        row_data$N_Quality_Warnings <- length(flags)
      }
      
      # Recommendations
      if ("recommendations" %in% names(quality_result)) {
        recommendations <- quality_result$recommendations
        if (length(recommendations) > 0) {
          row_data$Recommendations <- paste(recommendations, collapse = "; ")
        }
      }
      
      quality_table <- rbind(quality_table, row_data)
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(quality_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Sort by overall quality score descending
  if ("Overall_Quality_Score" %in% names(quality_table)) {
    quality_table <- quality_table[order(-quality_table$Overall_Quality_Score, na.last = TRUE), ]
    rownames(quality_table) <- NULL
  }
  
  return(quality_table)
}

#' Generate Sequence Statistics Table (Table S5)
generate_sequence_stats_table <- function(results_dir, verbose) {
  
  # Look for sequence analysis files
  seq_files <- list.files(results_dir, pattern = "sequence.*\\.rds$|.*sequence.*\\.rds$", 
                         full.names = TRUE)
  
  # Also look for FASTA files to analyze
  fasta_files <- list.files(c(results_dir, "data"), pattern = "\\.fasta$|\\.fa$", 
                           full.names = TRUE, recursive = TRUE)
  
  seq_table <- data.frame()
  
  # Process sequence analysis files
  for (seq_file in seq_files) {
    tryCatch({
      seq_result <- readRDS(seq_file)
      
      dataset_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", seq_file)
      dataset_name <- gsub("sequence_|_sequence", "", dataset_name)
      
      if ("stats" %in% names(seq_result)) {
        stats <- seq_result$stats
        
        row_data <- data.frame(
          Dataset = dataset_name,
          Source = "analysis_file",
          N_Sequences = stats$n_sequences %||% NA,
          Mean_Length = round(stats$mean_length %||% NA, 2),
          Median_Length = stats$median_length %||% NA,
          Min_Length = stats$min_length %||% NA,
          Max_Length = stats$max_length %||% NA,
          SD_Length = round(stats$sd_length %||% NA, 2),
          Mean_GC_Content = round(stats$mean_gc_content %||% NA, 4),
          Mean_N_Content = round(stats$mean_n_content %||% NA, 4),
          stringsAsFactors = FALSE
        )
        
        seq_table <- rbind(seq_table, row_data)
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error processing", basename(seq_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Analyze FASTA files directly
  for (fasta_file in fasta_files) {
    tryCatch({
      sequences <- readDNAStringSet(fasta_file)
      
      if (length(sequences) > 0) {
        dataset_name <- gsub(".*?([^/\\\\]+)\\.[^.]*$", "\\1", fasta_file)
        
        # Calculate statistics
        lengths <- width(sequences)
        gc_content <- letterFrequency(sequences, "GC") / width(sequences)
        n_content <- letterFrequency(sequences, "N") / width(sequences)
        
        row_data <- data.frame(
          Dataset = dataset_name,
          Source = "fasta_file",
          File_Path = fasta_file,
          N_Sequences = length(sequences),
          Mean_Length = round(mean(lengths), 2),
          Median_Length = median(lengths),
          Min_Length = min(lengths),
          Max_Length = max(lengths),
          SD_Length = round(sd(lengths), 2),
          Mean_GC_Content = round(mean(gc_content), 4),
          SD_GC_Content = round(sd(gc_content), 4),
          Mean_N_Content = round(mean(n_content), 4),
          Total_Nucleotides = sum(lengths),
          stringsAsFactors = FALSE
        )
        
        seq_table <- rbind(seq_table, row_data)
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error analyzing", basename(fasta_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Sort by number of sequences descending
  if ("N_Sequences" %in% names(seq_table)) {
    seq_table <- seq_table[order(-seq_table$N_Sequences), ]
    rownames(seq_table) <- NULL
  }
  
  return(seq_table)
}

#' Generate supplementary files
generate_supplementary_files <- function(results_dir, output_dir, include_raw, verbose) {
  
  supplementary_files <- list()
  
  # File S1: Best PWM matrices in JASPAR format
  if (verbose) cat("  Creating File S1: PWM matrices...\n")
  supplementary_files$file_s1 <- create_pwm_matrices_file(results_dir, output_dir, verbose)
  
  # File S2: Sequence alignments
  if (verbose) cat("  Creating File S2: Sequence alignments...\n")
  supplementary_files$file_s2 <- create_alignments_file(results_dir, output_dir, verbose)
  
  # File S3: Method configuration files
  if (verbose) cat("  Creating File S3: Method configurations...\n")
  supplementary_files$file_s3 <- create_configurations_file(results_dir, output_dir, verbose)
  
  # File S4: Detailed statistical results
  if (verbose) cat("  Creating File S4: Statistical details...\n")
  supplementary_files$file_s4 <- create_statistical_details_file(results_dir, output_dir, verbose)
  
  # Raw data files (optional)
  if (include_raw) {
    if (verbose) cat("  Including raw data files...\n")
    supplementary_files$raw_files <- copy_raw_files(results_dir, output_dir, verbose)
  }
  
  return(supplementary_files)
}

#' Create PWM matrices file (File S1)
create_pwm_matrices_file <- function(results_dir, output_dir, verbose) {
  
  pwm_files <- list.files(results_dir, pattern = ".*pwm.*\\.rds$", full.names = TRUE)
  
  if (length(pwm_files) == 0) {
    return(NULL)
  }
  
  output_file <- file.path(output_dir, "File_S1_PWM_Matrices.txt")
  
  pwm_content <- c(
    "# Supplementary File S1: Position Weight Matrices",
    "# CTCF PWM Testing Pipeline Results",
    paste("#", Sys.time()),
    "",
    "# Format: JASPAR-style PWM matrices",
    "# Each PWM is preceded by metadata lines starting with #",
    ""
  )
  
  for (pwm_file in pwm_files) {
    tryCatch({
      pwm_result <- readRDS(pwm_file)
      
      if ("pwm" %in% names(pwm_result)) {
        method_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", pwm_file)
        method_name <- gsub("_pwm|pwm_", "", method_name)
        
        # Add metadata
        pwm_content <- c(pwm_content,
                        paste("# Method:", method_name),
                        paste("# File:", basename(pwm_file)))
        
        if ("total_info" %in% names(pwm_result)) {
          pwm_content <- c(pwm_content,
                          paste("# Total Information Content:", round(pwm_result$total_info, 4), "bits"))
        }
        
        if ("n_sequences" %in% names(pwm_result)) {
          pwm_content <- c(pwm_content,
                          paste("# Number of Sequences:", pwm_result$n_sequences))
        }
        
        if ("consensus" %in% names(pwm_result)) {
          pwm_content <- c(pwm_content,
                          paste("# Consensus Sequence:", pwm_result$consensus))
        }
        
        # Add PWM matrix
        pwm_matrix <- pwm_result$pwm
        if (is.matrix(pwm_matrix)) {
          pwm_content <- c(pwm_content, paste(">", method_name))
          
          # Ensure proper row names
          if (nrow(pwm_matrix) == 4) {
            rownames(pwm_matrix) <- c("A", "C", "G", "T")
          }
          
          for (base in rownames(pwm_matrix)) {
            row_values <- paste(round(pwm_matrix[base, ], 6), collapse = " ")
            pwm_content <- c(pwm_content, paste(base, "[", row_values, "]"))
          }
        }
        
        pwm_content <- c(pwm_content, "")
      }
      
    }, error = function(e) {
      if (verbose) cat("    Error processing", basename(pwm_file), ":", conditionMessage(e), "\n")
    })
  }
  
  writeLines(pwm_content, output_file)
  
  if (verbose) cat("    PWM matrices saved to:", basename(output_file), "\n")
  
  return(output_file)
}

#' Create alignments file (File S2)
create_alignments_file <- function(results_dir, output_dir, verbose) {
  
  # Look for alignment files
  alignment_files <- list.files(results_dir, pattern = "aligned.*\\.fasta$|.*aligned.*\\.fasta$", 
                               full.names = TRUE, recursive = TRUE)
  
  if (length(alignment_files) == 0) {
    if (verbose) cat("    No alignment files found\n")
    return(NULL)
  }
  
  output_file <- file.path(output_dir, "File_S2_Sequence_Alignments.fasta")
  
  # Combine all alignments
  all_alignments <- c(
    paste("# Supplementary File S2: Sequence Alignments"),
    paste("# CTCF PWM Testing Pipeline Results"),
    paste("#", Sys.time()),
    ""
  )
  
  for (alignment_file in alignment_files) {
    tryCatch({
      method_name <- gsub(".*?([^/\\\\]+)\\.fasta$", "\\1", alignment_file)
      
      # Add method header
      all_alignments <- c(all_alignments,
                         paste("# Method:", method_name),
                         paste("# Source:", basename(alignment_file)),
                         "")
      
      # Read and add sequences
      alignment_content <- readLines(alignment_file)
      all_alignments <- c(all_alignments, alignment_content, "")
      
    }, error = function(e) {
      if (verbose) cat("    Error processing", basename(alignment_file), ":", conditionMessage(e), "\n")
    })
  }
  
  writeLines(all_alignments, output_file)
  
  if (verbose) cat("    Alignments saved to:", basename(output_file), "\n")
  
  return(output_file)
}

#' Create configurations file (File S3)
create_configurations_file <- function(results_dir, output_dir, verbose) {
  
  # Look for configuration files
  config_files <- list.files(c(results_dir, "config", "scripts"), 
                            pattern = "\\.json$|\\.yml$|\\.yaml$", 
                            full.names = TRUE, recursive = TRUE)
  
  output_file <- file.path(output_dir, "File_S3_Method_Configurations.json")
  
  all_configs <- list(
    metadata = list(
      title = "Method Configurations",
      description = "Configuration files used in CTCF PWM Testing Pipeline",
      timestamp = Sys.time()
    ),
    configurations = list()
  )
  
  for (config_file in config_files) {
    tryCatch({
      config_name <- gsub(".*?([^/\\\\]+)\\.[^.]*$", "\\1", config_file)
      
      if (grepl("\\.json$", config_file)) {
        config_content <- fromJSON(config_file)
      } else if (grepl("\\.yml$|\\.yaml$", config_file)) {
        config_content <- yaml::read_yaml(config_file)
      } else {
        next
      }
      
      all_configs$configurations[[config_name]] <- list(
        file = basename(config_file),
        path = config_file,
        content = config_content
      )
      
    }, error = function(e) {
      if (verbose) cat("    Error processing", basename(config_file), ":", conditionMessage(e), "\n")
    })
  }
  
  writeLines(toJSON(all_configs, pretty = TRUE), output_file)
  
  if (verbose) cat("    Configurations saved to:", basename(output_file), "\n")
  
  return(output_file)
}

#' Create statistical details file (File S4)
create_statistical_details_file <- function(results_dir, output_dir, verbose) {
  
  stat_files <- list.files(results_dir, pattern = "statistical.*\\.rds$|.*stats.*\\.rds$", 
                          full.names = TRUE)
  
  if (length(stat_files) == 0) {
    if (verbose) cat("    No statistical files found\n")
    return(NULL)
  }
  
  output_file <- file.path(output_dir, "File_S4_Statistical_Details.json")
  
  statistical_details <- list(
    metadata = list(
      title = "Detailed Statistical Results",
      description = "Complete statistical analysis results from CTCF PWM Testing Pipeline",
      timestamp = Sys.time()
    ),
    results = list()
  )
  
  for (stat_file in stat_files) {
    tryCatch({
      stat_result <- readRDS(stat_file)
      method_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", stat_file)
      method_name <- gsub("statistical_|_stats|_statistical", "", method_name)
      
      statistical_details$results[[method_name]] <- list(
        file = basename(stat_file),
        results = stat_result
      )
      
    }, error = function(e) {
      if (verbose) cat("    Error processing", basename(stat_file), ":", conditionMessage(e), "\n")
    })
  }
  
  writeLines(toJSON(statistical_details, pretty = TRUE), output_file)
  
  if (verbose) cat("    Statistical details saved to:", basename(output_file), "\n")
  
  return(output_file)
}

#' Copy raw files (optional)
copy_raw_files <- function(results_dir, output_dir, verbose) {
  
  raw_dir <- file.path(output_dir, "raw_data")
  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, recursive = TRUE)
  }
  
  # Copy important result files
  important_files <- list.files(results_dir, pattern = "\\.rds$|\\.json$|\\.log$", 
                               full.names = TRUE)
  
  copied_files <- character()
  
  for (file in important_files) {
    tryCatch({
      dest_file <- file.path(raw_dir, basename(file))
      file.copy(file, dest_file, overwrite = TRUE)
      copied_files <- c(copied_files, dest_file)
      
    }, error = function(e) {
      if (verbose) cat("    Error copying", basename(file), ":", conditionMessage(e), "\n")
    })
  }
  
  if (verbose) cat("    Copied", length(copied_files), "raw files\n")
  
  return(copied_files)
}

#' Save supplementary data
save_supplementary_data <- function(supplementary_data, supplementary_files, output_dir, format, verbose) {
  
  if (format == "xlsx") {
    # Create Excel workbook with multiple sheets
    wb <- createWorkbook()
    
    for (table_name in names(supplementary_data)) {
      if (nrow(supplementary_data[[table_name]]) > 0) {
        table_title <- toupper(gsub("_", " ", table_name))
        addWorksheet(wb, table_title)
        writeData(wb, table_title, supplementary_data[[table_name]])
      }
    }
    
    excel_file <- file.path(output_dir, "Supplementary_Tables.xlsx")
    saveWorkbook(wb, excel_file, overwrite = TRUE)
    
    if (verbose) cat("  Excel file saved to:", basename(excel_file), "\n")
    
  } else if (format == "csv") {
    # Save each table as separate CSV
    for (table_name in names(supplementary_data)) {
      if (nrow(supplementary_data[[table_name]]) > 0) {
        csv_file <- file.path(output_dir, paste0(table_name, ".csv"))
        write.csv(supplementary_data[[table_name]], csv_file, row.names = FALSE)
        
        if (verbose) cat("  CSV saved:", basename(csv_file), "\n")
      }
    }
    
  } else if (format == "json") {
    # Save as JSON
    json_file <- file.path(output_dir, "supplementary_data.json")
    writeLines(toJSON(supplementary_data, pretty = TRUE), json_file)
    
    if (verbose) cat("  JSON file saved to:", basename(json_file), "\n")
  }
}

#' Generate supplementary README
generate_supplementary_readme <- function(supplementary_data, supplementary_files, output_dir, verbose) {
  
  readme_content <- c(
    "# Supplementary Data - CTCF PWM Testing Pipeline",
    "",
    paste("Generated:", Sys.time()),
    "",
    "## Contents",
    "",
    "### Supplementary Tables",
    ""
  )
  
  # Document tables
  for (table_name in names(supplementary_data)) {
    if (nrow(supplementary_data[[table_name]]) > 0) {
      table_title <- toupper(gsub("_", " ", table_name))
      n_rows <- nrow(supplementary_data[[table_name]])
      n_cols <- ncol(supplementary_data[[table_name]])
      
      readme_content <- c(readme_content,
                         paste("**", table_title, "**"),
                         paste("- Rows:", n_rows),
                         paste("- Columns:", n_cols),
                         "")
    }
  }
  
  # Document files
  readme_content <- c(readme_content,
                     "### Supplementary Files",
                     "")
  
  for (file_name in names(supplementary_files)) {
    if (!is.null(supplementary_files[[file_name]])) {
      if (is.character(supplementary_files[[file_name]])) {
        if (length(supplementary_files[[file_name]]) == 1) {
          readme_content <- c(readme_content,
                             paste("**", toupper(gsub("_", " ", file_name)), "**"),
                             paste("- File:", basename(supplementary_files[[file_name]])),
                             "")
        } else {
          readme_content <- c(readme_content,
                             paste("**", toupper(gsub("_", " ", file_name)), "**"),
                             paste("- Files:", length(supplementary_files[[file_name]])),
                             "")
        }
      }
    }
  }
  
  # Add usage instructions
  readme_content <- c(readme_content,
                     "## Usage Instructions",
                     "",
                     "1. **Excel Files**: Open with Microsoft Excel, LibreOffice Calc, or similar",
                     "2. **CSV Files**: Open with any spreadsheet software or text editor",
                     "3. **JSON Files**: Open with text editor or JSON viewer",
                     "4. **FASTA Files**: Open with sequence analysis software",
                     "5. **PWM Files**: Use with motif analysis tools (JASPAR format)",
                     "",
                     "## File Descriptions",
                     "",
                     "- **Table S1**: Position Weight Matrix information content data",
                     "- **Table S2**: Method comparison and performance results", 
                     "- **Table S3**: Statistical validation and significance testing",
                     "- **Table S4**: Quality assessment scores and metrics",
                     "- **Table S5**: Input sequence dataset statistics",
                     "- **File S1**: PWM matrices in JASPAR format",
                     "- **File S2**: Sequence alignments in FASTA format",
                     "- **File S3**: Method configuration files",
                     "- **File S4**: Detailed statistical analysis results",
                     "",
                     "## Contact",
                     "",
                     "For questions about this supplementary data, please refer to the main manuscript",
                     "or contact the corresponding author.",
                     ""
  )
  
  readme_file <- file.path(output_dir, "README.md")
  writeLines(readme_content, readme_file)
  
  if (verbose) cat("  README generated:", basename(readme_file), "\n")
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-r", "--results"), type = "character", default = "results",
                help = "Results directory [default: %default]", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = "results/supplementary",
                help = "Output directory [default: %default]", metavar = "character"),
    make_option(c("--include-raw"), action = "store_true", default = FALSE,
                help = "Include raw data files"),
    make_option(c("-f", "--format"), type = "character", default = "xlsx",
                help = "Output format (xlsx, csv, json) [default: %default]", metavar = "character"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Generate Supplementary Data for Publications")
  opt <- parse_args(opt_parser)
  
  if (!dir.exists(opt$results)) {
    stop("Results directory does not exist: ", opt$results, call. = FALSE)
  }
  
  # Generate supplementary data
  tryCatch({
    supplementary <- generate_supplementary_data(
      results_dir = opt$results,
      output_dir = opt$output,
      include_raw = opt$`include-raw`,
      format = opt$format,
      verbose = opt$verbose
    )
    
    cat("Supplementary data generated successfully.\n")
    
  }, error = function(e) {
    cat("Error generating supplementary data:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
