#!/usr/bin/env Rscript

# Validate Known CTCF Sites
# Validates PWM models against known CTCF binding sites
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

#' Main function for validating known CTCF sites
#' @param pwm_file Path to PWM results file
#' @param known_sites_file Path to known CTCF sites (BED format)
#' @param reference_genome Path to reference genome FASTA
#' @param output_file Output file for validation results
#' @param score_threshold Score threshold for binding site prediction
#' @param extend_region Extend regions around known sites (bp)
#' @param verbose Enable verbose output
validate_known_sites <- function(pwm_file, known_sites_file = NULL,
                                reference_genome = "data/reference_genome/hg38.fa",
                                output_file = NULL, score_threshold = 0.8,
                                extend_region = 50, verbose = FALSE) {
  
  if (verbose) cat("Validating PWM against known CTCF sites...\n")
  
  # Load PWM model
  if (verbose) cat("Loading PWM model...\n")
  pwm_result <- readRDS(pwm_file)
  
  # Extract PWM matrix
  if ("pwm" %in% names(pwm_result)) {
    pwm_matrix <- pwm_result$pwm
  } else if ("prob_matrix" %in% names(pwm_result)) {
    pwm_matrix <- pwm_result$prob_matrix
  } else {
    stop("No PWM matrix found in results file")
  }
  
  # Load known CTCF sites
  if (verbose) cat("Loading known CTCF sites...\n")
  known_sites <- load_known_sites(known_sites_file, verbose)
  
  # Load reference genome
  if (verbose) cat("Loading reference genome...\n")
  ref_genome <- load_reference_genome(reference_genome, verbose)
  
  # Extract sequences around known sites
  if (verbose) cat("Extracting sequences around known sites...\n")
  site_sequences <- extract_site_sequences(known_sites, ref_genome, extend_region, verbose)
  
  # Score sequences with PWM
  if (verbose) cat("Scoring sequences with PWM...\n")
  validation_results <- score_sequences_with_pwm(site_sequences, pwm_matrix, 
                                                score_threshold, verbose)
  
  # Analyze validation results
  if (verbose) cat("Analyzing validation results...\n")
  analysis <- analyze_validation_results(validation_results, verbose)
  
  # Combine results
  final_results <- list(
    pwm_file = pwm_file,
    known_sites_file = known_sites_file,
    reference_genome = reference_genome,
    parameters = list(
      score_threshold = score_threshold,
      extend_region = extend_region
    ),
    validation_results = validation_results,
    analysis = analysis,
    timestamp = Sys.time()
  )
  
  # Save results
  if (!is.null(output_file)) {
    save_validation_results(final_results, output_file, verbose)
  }
  
  # Display summary
  display_validation_summary(analysis, verbose)
  
  return(final_results)
}

#' Load known CTCF sites
load_known_sites <- function(known_sites_file, verbose) {
  
  if (is.null(known_sites_file)) {
    # Use default known sites
    known_sites_file <- "data/K562_CTCF_peaks.bed"
  }
  
  if (!file.exists(known_sites_file)) {
    stop("Known sites file not found: ", known_sites_file)
  }
  
  # Read BED format file
  bed_data <- read.table(known_sites_file, header = FALSE, sep = "\t", 
                        stringsAsFactors = FALSE, comment.char = "#")
  
  # Ensure minimum required columns
  if (ncol(bed_data) < 3) {
    stop("BED file must have at least 3 columns (chr, start, end)")
  }
  
  # Standardize column names
  colnames(bed_data)[1:3] <- c("chr", "start", "end")
  
  # Add additional columns if available
  if (ncol(bed_data) >= 4) colnames(bed_data)[4] <- "name"
  if (ncol(bed_data) >= 5) colnames(bed_data)[5] <- "score"
  if (ncol(bed_data) >= 6) colnames(bed_data)[6] <- "strand"
  
  # Filter for standard chromosomes
  standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
  bed_data <- bed_data[bed_data$chr %in% standard_chrs, ]
  
  if (verbose) cat("  Loaded", nrow(bed_data), "known CTCF sites\n")
  
  return(bed_data)
}

#' Load reference genome
load_reference_genome <- function(reference_genome, verbose) {
  
  if (!file.exists(reference_genome)) {
    stop("Reference genome file not found: ", reference_genome)
  }
  
  # Read FASTA file
  tryCatch({
    ref_genome <- readDNAStringSet(reference_genome)
    
    if (verbose) cat("  Loaded reference genome with", length(ref_genome), "sequences\n")
    
    return(ref_genome)
    
  }, error = function(e) {
    stop("Error loading reference genome: ", conditionMessage(e))
  })
}

#' Extract sequences around known sites
extract_site_sequences <- function(known_sites, ref_genome, extend_region, verbose) {
  
  site_sequences <- list()
  failed_extractions <- 0
  
  for (i in 1:nrow(known_sites)) {
    site <- known_sites[i, ]
    
    tryCatch({
      # Find chromosome in reference genome
      chr_name <- site$chr
      
      # Try different chromosome naming conventions
      chr_candidates <- c(chr_name, 
                         gsub("chr", "", chr_name), 
                         paste0("chr", gsub("chr", "", chr_name)))
      
      chr_seq <- NULL
      for (candidate in chr_candidates) {
        if (candidate %in% names(ref_genome)) {
          chr_seq <- ref_genome[[candidate]]
          break
        }
      }
      
      if (is.null(chr_seq)) {
        failed_extractions <- failed_extractions + 1
        next
      }
      
      # Calculate extended coordinates
      site_start <- max(1, site$start - extend_region)
      site_end <- min(length(chr_seq), site$end + extend_region)
      
      # Extract sequence
      site_sequence <- subseq(chr_seq, site_start, site_end)
      
      # Store sequence with metadata
      site_info <- list(
        sequence = site_sequence,
        chr = site$chr,
        original_start = site$start,
        original_end = site$end,
        extended_start = site_start,
        extended_end = site_end,
        name = site$name %||% paste0("site_", i),
        score = site$score %||% NA
      )
      
      site_sequences[[i]] <- site_info
      
    }, error = function(e) {
      failed_extractions <- failed_extractions + 1
      if (verbose) cat("  Failed to extract sequence for site", i, ":", conditionMessage(e), "\n")
    })
  }
  
  # Remove NULL entries
  site_sequences <- site_sequences[!sapply(site_sequences, is.null)]
  
  if (verbose) {
    cat("  Extracted", length(site_sequences), "sequences\n")
    if (failed_extractions > 0) {
      cat("  Failed extractions:", failed_extractions, "\n")
    }
  }
  
  return(site_sequences)
}

#' Score sequences with PWM
score_sequences_with_pwm <- function(site_sequences, pwm_matrix, score_threshold, verbose) {
  
  if (length(site_sequences) == 0) {
    return(data.frame())
  }
  
  # Ensure PWM matrix is properly formatted
  if (is.matrix(pwm_matrix)) {
    if (nrow(pwm_matrix) != 4) {
      pwm_matrix <- t(pwm_matrix)
    }
    
    if (nrow(pwm_matrix) == 4 && is.null(rownames(pwm_matrix))) {
      rownames(pwm_matrix) <- c("A", "C", "G", "T")
    }
  } else {
    stop("PWM matrix must be a matrix")
  }
  
  pwm_length <- ncol(pwm_matrix)
  
  # Score each sequence
  scoring_results <- data.frame()
  
  for (i in seq_along(site_sequences)) {
    site_info <- site_sequences[[i]]
    sequence <- site_info$sequence
    
    tryCatch({
      # Score sequence in both orientations
      scores_forward <- score_sequence_pwm(sequence, pwm_matrix)
      scores_reverse <- score_sequence_pwm(reverseComplement(sequence), pwm_matrix)
      
      # Find best score and position
      best_forward <- if (length(scores_forward) > 0) max(scores_forward, na.rm = TRUE) else -Inf
      best_reverse <- if (length(scores_reverse) > 0) max(scores_reverse, na.rm = TRUE) else -Inf
      
      if (best_forward >= best_reverse) {
        best_score <- best_forward
        best_position <- which.max(scores_forward)
        best_strand <- "+"
      } else {
        best_score <- best_reverse
        best_position <- which.max(scores_reverse)
        best_strand <- "-"
      }
      
      # Determine if site passes threshold
      passes_threshold <- best_score >= score_threshold
      
      # Create result row
      result_row <- data.frame(
        Site_Index = i,
        Site_Name = site_info$name,
        Chr = site_info$chr,
        Original_Start = site_info$original_start,
        Original_End = site_info$original_end,
        Extended_Start = site_info$extended_start,
        Extended_End = site_info$extended_end,
        Sequence_Length = length(sequence),
        Best_Score = best_score,
        Best_Position = best_position,
        Best_Strand = best_strand,
        Passes_Threshold = passes_threshold,
        Known_Score = site_info$score,
        stringsAsFactors = FALSE
      )
      
      scoring_results <- rbind(scoring_results, result_row)
      
    }, error = function(e) {
      if (verbose) cat("  Error scoring sequence", i, ":", conditionMessage(e), "\n")
    })
  }
  
  if (verbose) cat("  Scored", nrow(scoring_results), "sequences\n")
  
  return(scoring_results)
}

#' Score a single sequence with PWM
score_sequence_pwm <- function(sequence, pwm_matrix) {
  
  seq_length <- length(sequence)
  pwm_length <- ncol(pwm_matrix)
  
  if (seq_length < pwm_length) {
    return(numeric(0))
  }
  
  # Convert sequence to character
  seq_chars <- as.character(sequence)
  
  scores <- numeric(seq_length - pwm_length + 1)
  
  for (pos in 1:(seq_length - pwm_length + 1)) {
    
    # Extract subsequence
    subseq <- substr(seq_chars, pos, pos + pwm_length - 1)
    subseq_chars <- strsplit(subseq, "")[[1]]
    
    # Calculate score
    score <- 0
    valid_position <- TRUE
    
    for (i in 1:pwm_length) {
      base <- subseq_chars[i]
      if (base %in% c("A", "C", "G", "T")) {
        score <- score + log2(pwm_matrix[base, i] + 1e-10)
      } else {
        # Handle ambiguous bases
        valid_position <- FALSE
        break
      }
    }
    
    scores[pos] <- if (valid_position) score else -Inf
  }
  
  return(scores)
}

#' Analyze validation results
analyze_validation_results <- function(validation_results, verbose) {
  
  if (nrow(validation_results) == 0) {
    return(list(
      n_sites = 0,
      n_validated = 0,
      validation_rate = 0,
      mean_score = NA,
      score_distribution = list()
    ))
  }
  
  # Basic statistics
  n_sites <- nrow(validation_results)
  n_validated <- sum(validation_results$Passes_Threshold, na.rm = TRUE)
  validation_rate <- n_validated / n_sites
  
  # Score statistics
  scores <- validation_results$Best_Score
  mean_score <- mean(scores, na.rm = TRUE)
  median_score <- median(scores, na.rm = TRUE)
  sd_score <- sd(scores, na.rm = TRUE)
  
  # Score distribution
  score_quantiles <- quantile(scores, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  
  # Strand analysis
  strand_table <- table(validation_results$Best_Strand)
  
  # Performance by chromosome
  chr_performance <- aggregate(validation_results$Passes_Threshold, 
                              by = list(Chr = validation_results$Chr), 
                              FUN = function(x) c(n = length(x), validated = sum(x, na.rm = TRUE)))
  
  # Correlation with known scores if available
  score_correlation <- NA
  if ("Known_Score" %in% names(validation_results) && 
      any(!is.na(validation_results$Known_Score))) {
    known_scores <- validation_results$Known_Score
    pwm_scores <- validation_results$Best_Score
    
    valid_pairs <- !is.na(known_scores) & !is.na(pwm_scores)
    if (sum(valid_pairs) > 5) {
      score_correlation <- cor(known_scores[valid_pairs], pwm_scores[valid_pairs])
    }
  }
  
  # Compile analysis results
  analysis <- list(
    n_sites = n_sites,
    n_validated = n_validated,
    validation_rate = validation_rate,
    mean_score = mean_score,
    median_score = median_score,
    sd_score = sd_score,
    score_quantiles = score_quantiles,
    strand_distribution = strand_table,
    chr_performance = chr_performance,
    score_correlation = score_correlation,
    high_confidence_sites = sum(scores > mean_score + sd_score, na.rm = TRUE),
    low_confidence_sites = sum(scores < mean_score - sd_score, na.rm = TRUE)
  )
  
  return(analysis)
}

#' Save validation results
save_validation_results <- function(results, output_file, verbose) {
  
  if (grepl("\\.rds$", output_file)) {
    saveRDS(results, output_file)
  } else if (grepl("\\.json$", output_file)) {
    writeLines(toJSON(results, pretty = TRUE), output_file)
  } else {
    # Default to RDS
    saveRDS(results, paste0(output_file, ".rds"))
  }
  
  # Also save validation table as CSV
  if (nrow(results$validation_results) > 0) {
    csv_file <- gsub("\\.[^.]*$", "_validation_table.csv", output_file)
    write.csv(results$validation_results, csv_file, row.names = FALSE)
    
    if (verbose) cat("Validation table saved to:", basename(csv_file), "\n")
  }
  
  if (verbose) cat("Validation results saved to:", basename(output_file), "\n")
}

#' Display validation summary
display_validation_summary <- function(analysis, verbose) {
  
  cat("\n=== CTCF Site Validation Summary ===\n")
  
  cat("Total sites analyzed:", analysis$n_sites, "\n")
  cat("Sites validated:", analysis$n_validated, "\n")
  cat("Validation rate:", round(analysis$validation_rate * 100, 1), "%\n")
  
  if (!is.na(analysis$mean_score)) {
    cat("Mean PWM score:", round(analysis$mean_score, 3), "\n")
    cat("Median PWM score:", round(analysis$median_score, 3), "\n")
    cat("Score SD:", round(analysis$sd_score, 3), "\n")
  }
  
  # Score distribution
  if (length(analysis$score_quantiles) > 0) {
    cat("\nScore distribution:\n")
    for (i in seq_along(analysis$score_quantiles)) {
      quantile_name <- names(analysis$score_quantiles)[i]
      quantile_value <- analysis$score_quantiles[i]
      cat("  ", quantile_name, ": ", round(quantile_value, 3), "\n")
    }
  }
  
  # Strand preference
  if (length(analysis$strand_distribution) > 0) {
    cat("\nStrand preference:\n")
    for (strand in names(analysis$strand_distribution)) {
      count <- analysis$strand_distribution[strand]
      cat("  ", strand, ": ", count, "\n")
    }
  }
  
  # High/low confidence sites
  cat("\nConfidence analysis:\n")
  cat("  High confidence sites:", analysis$high_confidence_sites, "\n")
  cat("  Low confidence sites:", analysis$low_confidence_sites, "\n")
  
  # Correlation with known scores
  if (!is.na(analysis$score_correlation)) {
    cat("Correlation with known scores:", round(analysis$score_correlation, 3), "\n")
  }
  
  # Validation quality assessment
  if (analysis$validation_rate >= 0.8) {
    cat("\nValidation Quality: EXCELLENT (≥80% validation rate)\n")
  } else if (analysis$validation_rate >= 0.6) {
    cat("\nValidation Quality: GOOD (≥60% validation rate)\n")
  } else if (analysis$validation_rate >= 0.4) {
    cat("\nValidation Quality: MODERATE (≥40% validation rate)\n")
  } else {
    cat("\nValidation Quality: POOR (<40% validation rate)\n")
  }
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-p", "--pwm"), type = "character", default = NULL,
                help = "PWM results file", metavar = "character"),
    make_option(c("-s", "--sites"), type = "character", default = NULL,
                help = "Known CTCF sites file (BED format)", metavar = "character"),
    make_option(c("-g", "--genome"), type = "character", default = "data/reference_genome/hg38.fa",
                help = "Reference genome FASTA [default: %default]", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file for validation results", metavar = "character"),
    make_option(c("-t", "--threshold"), type = "double", default = 0.8,
                help = "Score threshold [default: %default]", metavar = "double"),
    make_option(c("-e", "--extend"), type = "integer", default = 50,
                help = "Extend regions around sites (bp) [default: %default]", metavar = "integer"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Validate PWM Against Known CTCF Sites")
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$pwm)) {
    print_help(opt_parser)
    stop("PWM file is required.", call. = FALSE)
  }
  
  if (!file.exists(opt$pwm)) {
    stop("PWM file does not exist: ", opt$pwm, call. = FALSE)
  }
  
  # Generate default output filename if not provided
  if (is.null(opt$output)) {
    base_name <- gsub("\\.[^.]*$", "", basename(opt$pwm))
    opt$output <- file.path("results", paste0(base_name, "_validation_results.rds"))
  }
  
  # Validate known sites
  tryCatch({
    results <- validate_known_sites(
      pwm_file = opt$pwm,
      known_sites_file = opt$sites,
      reference_genome = opt$genome,
      output_file = opt$output,
      score_threshold = opt$threshold,
      extend_region = opt$extend,
      verbose = opt$verbose
    )
    
    cat("Validation completed successfully.\n")
    
    # Return appropriate exit code based on validation rate
    if (results$analysis$validation_rate >= 0.6) {
      quit(status = 0)
    } else {
      cat("Warning: Low validation rate (", round(results$analysis$validation_rate * 100, 1), "%)\n")
      quit(status = 1)
    }
    
  }, error = function(e) {
    cat("Error in site validation:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
