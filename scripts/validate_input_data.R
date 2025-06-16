#!/usr/bin/env Rscript

# Input Data Validation Script
# Comprehensive validation of input sequence data for CTCF PWM pipeline
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
  library(jsonlite)
})

# Command line options
option_list <- list(
  make_option(c("--input", "-i"), type="character", default=NULL,
              help="Input FASTA file path [required]", metavar="character"),
  make_option(c("--output", "-o"), type="character", default="results/data_validation.html",
              help="Output validation report file [default: %default]", metavar="character"),
  make_option(c("--format"), type="character", default="html",
              help="Output format (html, txt, json) [default: %default]", metavar="character"),
  make_option(c("--min-length"), type="integer", default=50,
              help="Minimum sequence length [default: %default]", metavar="integer"),
  make_option(c("--max-length"), type="integer", default=500,
              help="Maximum sequence length [default: %default]", metavar="integer"),
  make_option(c("--max-n-ratio"), type="double", default=0.01,
              help="Maximum N ratio allowed [default: %default]", metavar="double"),
  make_option(c("--min-complexity"), type="double", default=1.0,
              help="Minimum sequence complexity [default: %default]", metavar="double"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Validate input sequence data for CTCF PWM pipeline")
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  cat("Error: Input file is required\n")
  print_help(opt_parser)
  quit(status=1)
}

if (!file.exists(opt$input)) {
  cat("Error: Input file does not exist:", opt$input, "\n")
  quit(status=1)
}

# Create output directory
output_dir <- dirname(opt$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (opt$verbose) {
  cat("Starting input data validation...\n")
  cat("Input file:", opt$input, "\n")
  cat("Output file:", opt$output, "\n")
}

#' Assess sequence data quality
#' @param sequences_file Path to FASTA file
#' @param min_length Minimum acceptable sequence length
#' @param max_length Maximum acceptable sequence length
#' @param max_n_ratio Maximum allowed N ratio
#' @param min_complexity Minimum sequence complexity
#' @param verbose Enable verbose output
assess_data_quality <- function(sequences_file, min_length = 50, max_length = 500,
                               max_n_ratio = 0.01, min_complexity = 1.0, verbose = FALSE) {
  
  if (verbose) cat("Loading sequences from:", sequences_file, "\n")
  
  # Load sequences
  tryCatch({
    sequences <- readDNAStringSet(sequences_file)
  }, error = function(e) {
    stop("Failed to read FASTA file: ", e$message)
  })
  
  if (verbose) cat("Loaded", length(sequences), "sequences\n")
  
  # Initialize results
  results <- list(
    total_sequences = length(sequences),
    validation_passes = list(),
    validation_failures = list(),
    quality_metrics = list(),
    recommendations = list()
  )
  
  # 1. Sequence format validation
  if (verbose) cat("Validating sequence format...\n")
  format_check <- validate_sequence_format(sequences)
  results$validation_passes$format_valid <- format_check$valid
  results$validation_failures$format_issues <- format_check$issues
  
  # 2. Length distribution analysis
  if (verbose) cat("Analyzing length distribution...\n")
  length_analysis <- analyze_length_distribution(sequences, min_length, max_length)
  results$quality_metrics$length_distribution <- length_analysis
  results$validation_passes$length_acceptable <- length_analysis$pass_rate > 0.8
  
  # 3. GC content analysis
  if (verbose) cat("Analyzing GC content...\n")
  gc_analysis <- analyze_gc_content(sequences)
  results$quality_metrics$gc_content <- gc_analysis
  results$validation_passes$gc_content_normal <- gc_analysis$is_normal
  
  # 4. Sequence complexity analysis
  if (verbose) cat("Analyzing sequence complexity...\n")
  complexity_analysis <- analyze_complexity(sequences, min_complexity)
  results$quality_metrics$complexity <- complexity_analysis
  results$validation_passes$complexity_adequate <- complexity_analysis$pass_rate > 0.9
  
  # 5. N-base content analysis
  if (verbose) cat("Analyzing N-base content...\n")
  n_content_analysis <- analyze_n_content(sequences, max_n_ratio)
  results$quality_metrics$n_content <- n_content_analysis
  results$validation_passes$n_content_acceptable <- n_content_analysis$pass_rate > 0.95
  
  # 6. Duplicate detection
  if (verbose) cat("Detecting duplicates...\n")
  duplicate_analysis <- detect_duplicates(sequences)
  results$quality_metrics$duplicates <- duplicate_analysis
  results$validation_passes$duplicate_rate_low <- duplicate_analysis$duplicate_rate < 0.1
  
  # 7. Overall quality grade
  results$quality_grade <- calculate_overall_quality_grade(results)
  
  # 8. Generate recommendations
  results$recommendations <- generate_quality_recommendations(results)
  
  return(results)
}

#' Validate sequence format
validate_sequence_format <- function(sequences) {
  valid_chars <- c("A", "T", "G", "C", "N")
  issues <- list()
  
  # Check for valid DNA characters
  seq_chars <- unique(unlist(strsplit(as.character(sequences), "")))
  invalid_chars <- setdiff(seq_chars, valid_chars)
  
  if (length(invalid_chars) > 0) {
    issues$invalid_characters <- invalid_chars
  }
  
  # Check for empty sequences
  empty_seqs <- which(width(sequences) == 0)
  if (length(empty_seqs) > 0) {
    issues$empty_sequences <- length(empty_seqs)
  }
  
  return(list(
    valid = length(issues) == 0,
    issues = issues
  ))
}

#' Analyze length distribution
analyze_length_distribution <- function(sequences, min_length, max_length) {
  lengths <- width(sequences)
  
  # Calculate statistics
  length_stats <- list(
    mean = mean(lengths),
    median = median(lengths),
    min = min(lengths),
    max = max(lengths),
    sd = sd(lengths),
    range = range(lengths)
  )
  
  # Check length criteria
  valid_lengths <- lengths >= min_length & lengths <= max_length
  pass_rate <- sum(valid_lengths) / length(lengths)
  
  # Length distribution
  length_bins <- cut(lengths, breaks = 10)
  length_distribution <- table(length_bins)
  
  return(list(
    statistics = length_stats,
    pass_rate = pass_rate,
    valid_count = sum(valid_lengths),
    invalid_count = sum(!valid_lengths),
    distribution = length_distribution
  ))
}

#' Analyze GC content
analyze_gc_content <- function(sequences) {
  # Calculate GC content for each sequence
  gc_contents <- letterFrequency(sequences, letters = "GC", as.prob = TRUE)[,1]
  
  # Statistics
  gc_stats <- list(
    mean = mean(gc_contents),
    median = median(gc_contents),
    sd = sd(gc_contents),
    range = range(gc_contents)
  )
  
  # Check if GC content is in normal range (30-70%)
  normal_gc <- gc_contents >= 0.3 & gc_contents <= 0.7
  is_normal <- sum(normal_gc) / length(gc_contents) > 0.8
  
  return(list(
    statistics = gc_stats,
    distribution = gc_contents,
    is_normal = is_normal,
    normal_rate = sum(normal_gc) / length(gc_contents)
  ))
}

#' Analyze sequence complexity
analyze_complexity <- function(sequences, min_complexity) {
  # Calculate Shannon entropy for each sequence
  complexities <- sapply(sequences, function(seq) {
    seq_str <- as.character(seq)
    chars <- table(strsplit(seq_str, "")[[1]])
    probs <- chars / sum(chars)
    -sum(probs * log2(probs))
  })
  
  # Check complexity threshold
  adequate_complexity <- complexities >= min_complexity
  pass_rate <- sum(adequate_complexity) / length(complexities)
  
  complexity_stats <- list(
    mean = mean(complexities),
    median = median(complexities),
    sd = sd(complexities),
    range = range(complexities)
  )
  
  return(list(
    statistics = complexity_stats,
    values = complexities,
    pass_rate = pass_rate,
    adequate_count = sum(adequate_complexity)
  ))
}

#' Analyze N content
analyze_n_content <- function(sequences, max_n_ratio) {
  # Calculate N ratio for each sequence
  n_ratios <- letterFrequency(sequences, letters = "N", as.prob = TRUE)[,1]
  
  # Check N ratio threshold
  acceptable_n <- n_ratios <= max_n_ratio
  pass_rate <- sum(acceptable_n) / length(n_ratios)
  
  n_stats <- list(
    mean = mean(n_ratios),
    median = median(n_ratios),
    max = max(n_ratios),
    sd = sd(n_ratios)
  )
  
  return(list(
    statistics = n_stats,
    values = n_ratios,
    pass_rate = pass_rate,
    acceptable_count = sum(acceptable_n)
  ))
}

#' Detect duplicate sequences
detect_duplicates <- function(sequences) {
  seq_strings <- as.character(sequences)
  unique_seqs <- unique(seq_strings)
  
  duplicate_rate <- 1 - (length(unique_seqs) / length(seq_strings))
  
  return(list(
    total_sequences = length(seq_strings),
    unique_sequences = length(unique_seqs),
    duplicate_count = length(seq_strings) - length(unique_seqs),
    duplicate_rate = duplicate_rate
  ))
}

#' Calculate overall quality grade
calculate_overall_quality_grade <- function(results) {
  # Count passed validations
  passes <- unlist(results$validation_passes)
  pass_count <- sum(passes, na.rm = TRUE)
  total_tests <- length(passes)
  
  pass_rate <- pass_count / total_tests
  
  if (pass_rate >= 0.9) {
    return("EXCELLENT")
  } else if (pass_rate >= 0.8) {
    return("GOOD")
  } else if (pass_rate >= 0.7) {
    return("ACCEPTABLE")
  } else {
    return("POOR")
  }
}

#' Generate quality recommendations
generate_quality_recommendations <- function(results) {
  recommendations <- list()
  
  # Length recommendations
  if (!results$validation_passes$length_acceptable) {
    recommendations$length <- "Consider filtering sequences by length to improve data quality"
  }
  
  # GC content recommendations
  if (!results$validation_passes$gc_content_normal) {
    recommendations$gc_content <- "GC content distribution is abnormal. Check for bias in sequence selection"
  }
  
  # Complexity recommendations
  if (!results$validation_passes$complexity_adequate) {
    recommendations$complexity <- "Many sequences have low complexity. Consider filtering low-complexity regions"
  }
  
  # N content recommendations
  if (!results$validation_passes$n_content_acceptable) {
    recommendations$n_content <- "High N content detected. Consider removing sequences with excessive Ns"
  }
  
  # Duplicate recommendations
  if (!results$validation_passes$duplicate_rate_low) {
    recommendations$duplicates <- "High duplicate rate detected. Consider removing redundant sequences"
  }
  
  return(recommendations)
}

#' Generate validation report
generate_validation_report <- function(validation_results, output_file, format = "html") {
  if (format == "html") {
    generate_html_report(validation_results, output_file)
  } else if (format == "json") {
    writeLines(toJSON(validation_results, pretty = TRUE), output_file)
  } else {
    generate_text_report(validation_results, output_file)
  }
}

#' Generate HTML validation report
generate_html_report <- function(results, output_file) {
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>CTCF Input Data Validation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background-color: #f0f8ff; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }
        .pass { color: green; font-weight: bold; }
        .fail { color: red; font-weight: bold; }
        .metric { margin: 10px 0; }
        .grade-excellent { color: green; font-size: 24px; font-weight: bold; }
        .grade-good { color: blue; font-size: 24px; font-weight: bold; }
        .grade-acceptable { color: orange; font-size: 24px; font-weight: bold; }
        .grade-poor { color: red; font-size: 24px; font-weight: bold; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <div class='header'>
        <h1>CTCF Input Data Validation Report</h1>
        <p>Generated: ", Sys.time(), "</p>
        <p>Total Sequences: ", results$total_sequences, "</p>
        <p>Overall Quality: <span class='grade-", tolower(results$quality_grade), "'>", results$quality_grade, "</span></p>
    </div>
    
    <div class='section'>
        <h2>Validation Summary</h2>
        <table>
            <tr><th>Test</th><th>Result</th><th>Details</th></tr>
            <tr><td>Format Validation</td><td class='", if(results$validation_passes$format_valid) "pass'>PASS" else "fail'>FAIL", "</td><td>", 
                if(length(results$validation_failures$format_issues) > 0) paste("Issues:", names(results$validation_failures$format_issues)) else "No issues", "</td></tr>
            <tr><td>Length Distribution</td><td class='", if(results$validation_passes$length_acceptable) "pass'>PASS" else "fail'>FAIL", "</td><td>Pass rate: ", 
                round(results$quality_metrics$length_distribution$pass_rate * 100, 1), "%</td></tr>
            <tr><td>GC Content</td><td class='", if(results$validation_passes$gc_content_normal) "pass'>PASS" else "fail'>FAIL", "</td><td>Normal rate: ", 
                round(results$quality_metrics$gc_content$normal_rate * 100, 1), "%</td></tr>
            <tr><td>Sequence Complexity</td><td class='", if(results$validation_passes$complexity_adequate) "pass'>PASS" else "fail'>FAIL", "</td><td>Pass rate: ", 
                round(results$quality_metrics$complexity$pass_rate * 100, 1), "%</td></tr>
            <tr><td>N Content</td><td class='", if(results$validation_passes$n_content_acceptable) "pass'>PASS" else "fail'>FAIL", "</td><td>Pass rate: ", 
                round(results$quality_metrics$n_content$pass_rate * 100, 1), "%</td></tr>
            <tr><td>Duplicate Rate</td><td class='", if(results$validation_passes$duplicate_rate_low) "pass'>PASS" else "fail'>FAIL", "</td><td>Duplicate rate: ", 
                round(results$quality_metrics$duplicates$duplicate_rate * 100, 1), "%</td></tr>
        </table>
    </div>")
    
  if (length(results$recommendations) > 0) {
    html_content <- paste0(html_content, "
    <div class='section'>
        <h2>Recommendations</h2>
        <ul>")
    
    for (rec_name in names(results$recommendations)) {
      html_content <- paste0(html_content, "<li><strong>", rec_name, ":</strong> ", results$recommendations[[rec_name]], "</li>")
    }
    
    html_content <- paste0(html_content, "
        </ul>
    </div>")
  }
  
  html_content <- paste0(html_content, "
</body>
</html>")
  
  writeLines(html_content, output_file)
}

#' Generate text validation report
generate_text_report <- function(results, output_file) {
  report_lines <- c(
    "=== CTCF Input Data Validation Report ===",
    paste("Generated:", Sys.time()),
    paste("Total Sequences:", results$total_sequences),
    paste("Overall Quality:", results$quality_grade),
    "",
    "=== Validation Results ===",
    paste("Format Validation:", if(results$validation_passes$format_valid) "PASS" else "FAIL"),
    paste("Length Distribution:", if(results$validation_passes$length_acceptable) "PASS" else "FAIL"),
    paste("GC Content:", if(results$validation_passes$gc_content_normal) "PASS" else "FAIL"),
    paste("Sequence Complexity:", if(results$validation_passes$complexity_adequate) "PASS" else "FAIL"),
    paste("N Content:", if(results$validation_passes$n_content_acceptable) "PASS" else "FAIL"),
    paste("Duplicate Rate:", if(results$validation_passes$duplicate_rate_low) "PASS" else "FAIL"),
    ""
  )
  
  if (length(results$recommendations) > 0) {
    report_lines <- c(report_lines, "=== Recommendations ===")
    for (rec_name in names(results$recommendations)) {
      report_lines <- c(report_lines, paste(rec_name, ":", results$recommendations[[rec_name]]))
    }
  }
  
  writeLines(report_lines, output_file)
}

# Main execution
if (opt$verbose) cat("Assessing data quality...\n")

validation_results <- assess_data_quality(
  sequences_file = opt$input,
  min_length = opt$`min-length`,
  max_length = opt$`max-length`,
  max_n_ratio = opt$`max-n-ratio`,
  min_complexity = opt$`min-complexity`,
  verbose = opt$verbose
)

if (opt$verbose) cat("Generating validation report...\n")

generate_validation_report(validation_results, opt$output, opt$format)

if (opt$verbose) {
  cat("Validation complete!\n")
  cat("Overall quality grade:", validation_results$quality_grade, "\n")
  cat("Report saved to:", opt$output, "\n")
}

# Exit with appropriate code
if (validation_results$quality_grade %in% c("EXCELLENT", "GOOD")) {
  quit(status = 0)
} else {
  quit(status = 1)
}
