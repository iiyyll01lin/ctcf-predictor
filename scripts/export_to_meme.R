#!/usr/bin/env Rscript

# Export to MEME
# Exports PWM results to MEME Suite format
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
})

#' Main function to export PWM to MEME format
#' @param pwm_file Path to PWM results file
#' @param output_file Output file for MEME format
#' @param motif_name Name for the motif
#' @param background_freqs Background frequencies (A, C, G, T)
#' @param url URL for motif information
#' @param verbose Enable verbose output
export_to_meme <- function(pwm_file, output_file = NULL, motif_name = NULL,
                          background_freqs = c(0.25, 0.25, 0.25, 0.25),
                          url = NULL, verbose = FALSE) {
  
  if (verbose) cat("Exporting PWM to MEME format...\n")
  
  # Load PWM results
  if (verbose) cat("Loading PWM results...\n")
  pwm_result <- readRDS(pwm_file)
  
  # Extract PWM matrix
  if ("pwm" %in% names(pwm_result)) {
    pwm_matrix <- pwm_result$pwm
  } else if ("prob_matrix" %in% names(pwm_result)) {
    pwm_matrix <- pwm_result$prob_matrix
  } else {
    stop("No PWM matrix found in results file")
  }
  
  # Validate PWM matrix
  pwm_matrix <- validate_pwm_matrix(pwm_matrix, verbose)
  
  # Generate output filename if not provided
  if (is.null(output_file)) {
    base_name <- gsub("\\.[^.]*$", "", basename(pwm_file))
    output_file <- paste0(base_name, ".meme")
  }
  
  # Generate motif name if not provided
  if (is.null(motif_name)) {
    motif_name <- gsub("\\.[^.]*$", "", basename(pwm_file))
    motif_name <- gsub("_pwm|pwm_", "", motif_name)
    motif_name <- paste0("CTCF_", motif_name)
  }
  
  # Create MEME format content
  if (verbose) cat("Creating MEME format content...\n")
  meme_content <- create_meme_content(pwm_matrix, pwm_result, motif_name, 
                                     background_freqs, url, verbose)
  
  # Write MEME file
  if (verbose) cat("Writing MEME file...\n")
  writeLines(meme_content, output_file)
  
  # Validate output file
  if (file.exists(output_file)) {
    file_size <- file.info(output_file)$size
    if (verbose) cat("MEME file created:", output_file, "(", file_size, "bytes)\n")
  } else {
    stop("Failed to create MEME file")
  }
  
  # Create summary
  summary_info <- list(
    pwm_file = pwm_file,
    output_file = output_file,
    motif_name = motif_name,
    motif_width = ncol(pwm_matrix),
    background_freqs = background_freqs,
    information_content = pwm_result$total_info %||% sum(pwm_result$info_content %||% 0),
    timestamp = Sys.time()
  )
  
  # Display summary
  display_meme_export_summary(summary_info, verbose)
  
  return(summary_info)
}

#' Validate PWM matrix format
validate_pwm_matrix <- function(pwm_matrix, verbose) {
  
  if (!is.matrix(pwm_matrix)) {
    stop("PWM must be a matrix")
  }
  
  # Check dimensions
  if (nrow(pwm_matrix) == 4 && ncol(pwm_matrix) > 4) {
    # Correct orientation (4 bases x positions)
    if (is.null(rownames(pwm_matrix))) {
      rownames(pwm_matrix) <- c("A", "C", "G", "T")
    }
  } else if (ncol(pwm_matrix) == 4 && nrow(pwm_matrix) > 4) {
    # Transpose needed (positions x 4 bases)
    pwm_matrix <- t(pwm_matrix)
    if (is.null(rownames(pwm_matrix))) {
      rownames(pwm_matrix) <- c("A", "C", "G", "T")
    }
  } else {
    stop("PWM matrix must have 4 rows (bases) and positions as columns")
  }
  
  # Ensure proper row names
  expected_bases <- c("A", "C", "G", "T")
  if (!all(rownames(pwm_matrix) %in% expected_bases)) {
    if (nrow(pwm_matrix) == 4) {
      rownames(pwm_matrix) <- expected_bases
      if (verbose) cat("  Applied standard base names (A, C, G, T)\n")
    } else {
      stop("PWM matrix must have 4 rows for bases A, C, G, T")
    }
  }
  
  # Reorder to standard order
  pwm_matrix <- pwm_matrix[expected_bases, , drop = FALSE]
  
  # Check if probabilities sum to 1 (approximately)
  col_sums <- colSums(pwm_matrix)
  if (any(abs(col_sums - 1.0) > 0.01)) {
    warning("PWM columns do not sum to 1.0. Renormalizing...")
    pwm_matrix <- sweep(pwm_matrix, 2, col_sums, FUN = "/")
  }
  
  # Ensure positive values
  if (any(pwm_matrix < 0)) {
    warning("PWM contains negative values. Setting to small positive value...")
    pwm_matrix[pwm_matrix < 0] <- 1e-10
  }
  
  # Ensure no zero values (causes issues with log calculations)
  pwm_matrix[pwm_matrix == 0] <- 1e-10
  
  if (verbose) cat("  PWM matrix validated:", nrow(pwm_matrix), "bases x", ncol(pwm_matrix), "positions\n")
  
  return(pwm_matrix)
}

#' Create MEME format content
create_meme_content <- function(pwm_matrix, pwm_result, motif_name, background_freqs, url, verbose) {
  
  # MEME header
  meme_lines <- c(
    "MEME version 4",
    "",
    paste("ALPHABET= ACGT"),
    "",
    "strands: + -",
    "",
    paste("Background letter frequencies (from uniform background):"),
    paste("A", sprintf("%.6f", background_freqs[1]), 
          "C", sprintf("%.6f", background_freqs[2]),
          "G", sprintf("%.6f", background_freqs[3]), 
          "T", sprintf("%.6f", background_freqs[4])),
    ""
  )
  
  # Motif information
  motif_width <- ncol(pwm_matrix)
  n_sites <- pwm_result$n_sequences %||% 1000  # Default if not available
  
  # Calculate E-value (simplified estimation)
  total_ic <- pwm_result$total_info %||% sum(pwm_result$info_content %||% rep(1, motif_width))
  e_value <- 10^(-total_ic/2)  # Rough approximation
  
  meme_lines <- c(meme_lines,
    paste("MOTIF", motif_name),
    ""
  )
  
  # Add URL if provided
  if (!is.null(url)) {
    meme_lines <- c(meme_lines, paste("URL", url), "")
  }
  
  # Letter probability matrix header
  meme_lines <- c(meme_lines,
    paste("letter-probability matrix: alength= 4 w=", motif_width, 
          "nsites=", n_sites, "E=", sprintf("%.2e", e_value))
  )
  
  # Add PWM matrix data
  for (pos in 1:motif_width) {
    pos_probs <- pwm_matrix[, pos]
    prob_line <- paste(sprintf("%.6f", pos_probs), collapse = " ")
    meme_lines <- c(meme_lines, prob_line)
  }
  
  # Add blank line at end
  meme_lines <- c(meme_lines, "")
  
  return(meme_lines)
}

#' Export multiple PWMs to MEME format
export_multiple_to_meme <- function(pwm_files, output_file, motif_prefix = "CTCF",
                                   background_freqs = c(0.25, 0.25, 0.25, 0.25),
                                   verbose = FALSE) {
  
  if (verbose) cat("Exporting", length(pwm_files), "PWMs to MEME format...\n")
  
  # MEME header
  meme_lines <- c(
    "MEME version 4",
    "",
    paste("ALPHABET= ACGT"),
    "",
    "strands: + -",
    "",
    paste("Background letter frequencies (from uniform background):"),
    paste("A", sprintf("%.6f", background_freqs[1]), 
          "C", sprintf("%.6f", background_freqs[2]),
          "G", sprintf("%.6f", background_freqs[3]), 
          "T", sprintf("%.6f", background_freqs[4])),
    ""
  )
  
  # Process each PWM file
  for (i in seq_along(pwm_files)) {
    pwm_file <- pwm_files[i]
    
    if (verbose) cat("  Processing PWM", i, ":", basename(pwm_file), "\n")
    
    tryCatch({
      # Load PWM
      pwm_result <- readRDS(pwm_file)
      
      # Extract and validate matrix
      if ("pwm" %in% names(pwm_result)) {
        pwm_matrix <- pwm_result$pwm
      } else if ("prob_matrix" %in% names(pwm_result)) {
        pwm_matrix <- pwm_result$prob_matrix
      } else {
        if (verbose) cat("    Skipping - no PWM matrix found\n")
        next
      }
      
      pwm_matrix <- validate_pwm_matrix(pwm_matrix, verbose = FALSE)
      
      # Generate motif name
      base_name <- gsub("\\.[^.]*$", "", basename(pwm_file))
      base_name <- gsub("_pwm|pwm_", "", base_name)
      motif_name <- paste0(motif_prefix, "_", base_name)
      
      # Add motif to MEME content
      motif_width <- ncol(pwm_matrix)
      n_sites <- pwm_result$n_sequences %||% 1000
      total_ic <- pwm_result$total_info %||% sum(pwm_result$info_content %||% rep(1, motif_width))
      e_value <- 10^(-total_ic/2)
      
      meme_lines <- c(meme_lines,
        paste("MOTIF", motif_name),
        "",
        paste("letter-probability matrix: alength= 4 w=", motif_width, 
              "nsites=", n_sites, "E=", sprintf("%.2e", e_value))
      )
      
      # Add PWM matrix
      for (pos in 1:motif_width) {
        pos_probs <- pwm_matrix[, pos]
        prob_line <- paste(sprintf("%.6f", pos_probs), collapse = " ")
        meme_lines <- c(meme_lines, prob_line)
      }
      
      meme_lines <- c(meme_lines, "")
      
    }, error = function(e) {
      if (verbose) cat("    Error processing", basename(pwm_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Write combined MEME file
  writeLines(meme_lines, output_file)
  
  if (verbose) cat("Combined MEME file created:", output_file, "\n")
  
  return(output_file)
}

#' Create MEME database format
create_meme_database <- function(pwm_files, output_file, database_name = "CTCF_PWMs",
                                version = "1.0", description = NULL, verbose = FALSE) {
  
  if (verbose) cat("Creating MEME database format...\n")
  
  # Database header
  header_lines <- c(
    "MEME version 4",
    "",
    paste("ALPHABET= ACGT"),
    "",
    "strands: + -",
    "",
    paste("Background letter frequencies (from uniform background):"),
    "A 0.250000 C 0.250000 G 0.250000 T 0.250000",
    ""
  )
  
  if (!is.null(description)) {
    header_lines <- c(header_lines, paste("//", description), "")
  }
  
  # Process PWMs and create database
  database_lines <- header_lines
  
  for (i in seq_along(pwm_files)) {
    pwm_file <- pwm_files[i]
    
    tryCatch({
      pwm_result <- readRDS(pwm_file)
      
      if ("pwm" %in% names(pwm_result)) {
        pwm_matrix <- pwm_result$pwm
      } else {
        next
      }
      
      pwm_matrix <- validate_pwm_matrix(pwm_matrix, verbose = FALSE)
      
      # Database entry
      base_name <- gsub("\\.[^.]*$", "", basename(pwm_file))
      motif_id <- paste0(database_name, "_", sprintf("%03d", i))
      motif_name <- paste0("CTCF_", base_name)
      
      motif_width <- ncol(pwm_matrix)
      n_sites <- pwm_result$n_sequences %||% 1000
      
      database_lines <- c(database_lines,
        paste("MOTIF", motif_id, motif_name),
        "",
        paste("letter-probability matrix: alength= 4 w=", motif_width, "nsites=", n_sites)
      )
      
      # Add matrix
      for (pos in 1:motif_width) {
        pos_probs <- pwm_matrix[, pos]
        prob_line <- paste(sprintf("%.6f", pos_probs), collapse = " ")
        database_lines <- c(database_lines, prob_line)
      }
      
      database_lines <- c(database_lines, "")
      
    }, error = function(e) {
      if (verbose) cat("    Error processing", basename(pwm_file), ":", conditionMessage(e), "\n")
    })
  }
  
  # Write database file
  writeLines(database_lines, output_file)
  
  if (verbose) cat("MEME database created:", output_file, "\n")
  
  return(output_file)
}

#' Display export summary
display_meme_export_summary <- function(summary_info, verbose) {
  
  cat("\n=== MEME Export Summary ===\n")
  
  cat("Input PWM file:", basename(summary_info$pwm_file), "\n")
  cat("Output MEME file:", basename(summary_info$output_file), "\n")
  cat("Motif name:", summary_info$motif_name, "\n")
  cat("Motif width:", summary_info$motif_width, "positions\n")
  
  if (!is.null(summary_info$information_content)) {
    cat("Information content:", round(summary_info$information_content, 3), "bits\n")
  }
  
  cat("Background frequencies: A=", summary_info$background_freqs[1],
      " C=", summary_info$background_freqs[2],
      " G=", summary_info$background_freqs[3], 
      " T=", summary_info$background_freqs[4], "\n")
  
  cat("Export completed:", summary_info$timestamp, "\n")
  
  cat("\n✅ MEME file ready for use with MEME Suite tools\n")
  cat("   Compatible with: FIMO, MAST, CENTRIMO, SPAMO, etc.\n")
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-p", "--pwm"), type = "character", default = NULL,
                help = "PWM results file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output MEME file", metavar = "character"),
    make_option(c("-n", "--name"), type = "character", default = NULL,
                help = "Motif name", metavar = "character"),
    make_option(c("-b", "--background"), type = "character", default = "0.25,0.25,0.25,0.25",
                help = "Background frequencies (A,C,G,T) [default: %default]", metavar = "character"),
    make_option(c("-u", "--url"), type = "character", default = NULL,
                help = "URL for motif information", metavar = "character"),
    make_option(c("-m", "--multiple"), type = "character", default = NULL,
                help = "Directory with multiple PWM files", metavar = "character"),
    make_option(c("-d", "--database"), action = "store_true", default = FALSE,
                help = "Create MEME database format"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Export PWM Results to MEME Suite Format")
  opt <- parse_args(opt_parser)
  
  # Parse background frequencies
  background_freqs <- as.numeric(strsplit(opt$background, ",")[[1]])
  if (length(background_freqs) != 4 || abs(sum(background_freqs) - 1.0) > 0.01) {
    stop("Background frequencies must be 4 values that sum to 1.0", call. = FALSE)
  }
  
  # Multiple PWM files mode
  if (!is.null(opt$multiple)) {
    if (!dir.exists(opt$multiple)) {
      stop("Multiple PWM directory does not exist: ", opt$multiple, call. = FALSE)
    }
    
    pwm_files <- list.files(opt$multiple, pattern = ".*pwm.*\\.rds$", full.names = TRUE)
    
    if (length(pwm_files) == 0) {
      stop("No PWM files found in directory: ", opt$multiple, call. = FALSE)
    }
    
    if (is.null(opt$output)) {
      opt$output <- file.path(opt$multiple, "combined_pwms.meme")
    }
    
    # Export multiple PWMs
    if (opt$database) {
      output_file <- create_meme_database(pwm_files, opt$output, verbose = opt$verbose)
    } else {
      output_file <- export_multiple_to_meme(pwm_files, opt$output, verbose = opt$verbose)
    }
    
    cat("Multiple PWM export completed:", output_file, "\n")
    
  } else {
    # Single PWM file mode
    if (is.null(opt$pwm)) {
      print_help(opt_parser)
      stop("PWM file is required for single file mode.", call. = FALSE)
    }
    
    if (!file.exists(opt$pwm)) {
      stop("PWM file does not exist: ", opt$pwm, call. = FALSE)
    }
    
    # Generate default output filename if not provided
    if (is.null(opt$output)) {
      base_name <- gsub("\\.[^.]*$", "", basename(opt$pwm))
      opt$output <- paste0(base_name, ".meme")
    }
    
    # Export single PWM
    tryCatch({
      summary_info <- export_to_meme(
        pwm_file = opt$pwm,
        output_file = opt$output,
        motif_name = opt$name,
        background_freqs = background_freqs,
        url = opt$url,
        verbose = opt$verbose
      )
      
      cat("MEME export completed successfully.\n")
      
    }, error = function(e) {
      cat("Error exporting to MEME format:", conditionMessage(e), "\n")
      quit(status = 1)
    })
  }
}
