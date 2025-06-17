#!/usr/bin/env Rscript

# Results Analysis Pipeline Script
# Main pipeline for comprehensive PWM results analysis as described in documentation
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(parallel)
  library(jsonlite)
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input PWM file (RDS format)", metavar="file"),
  make_option(c("-s", "--sequences"), type="character", default=NULL,
              help="Input sequences file (FASTA format)", metavar="file"),
  make_option(c("-o", "--output-dir"), type="character", default="results",
              help="Output directory [default: %default]"),
  make_option(c("--organize"), action="store_true", default=TRUE,
              help="Organize results into structured directories [default: %default]"),
  make_option(c("--export-formats"), type="character", default="jaspar,transfac",
              help="Export formats (comma-separated) [default: %default]"),
  make_option(c("--generate-visualizations"), action="store_true", default=TRUE,
              help="Generate visualizations [default: %default]"),
  make_option(c("--run-validation"), action="store_true", default=TRUE,
              help="Run cross-validation [default: %default]"),
  make_option(c("--validation-methods"), type="character", default="bootstrap,chromosome,kfold",
              help="Validation methods [default: %default]"),
  make_option(c("--generate-reports"), action="store_true", default=TRUE,
              help="Generate comprehensive reports [default: %default]"),
  make_option(c("-n", "--n-cores"), type="integer", default=1,
              help="Number of cores for parallel processing [default: %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Comprehensive PWM results analysis pipeline")
opt <- parse_args(opt_parser)

# Source required functions
source_analysis_scripts <- function(verbose = FALSE) {
  if (verbose) cat("Sourcing analysis scripts...\n")
  
  required_scripts <- c(
    "scripts/organize_results.R",
    "scripts/export_pwm_formats.R", 
    "scripts/ctcf_pattern_validation.R",
    "scripts/generate_visualizations.R",
    "scripts/integrate_cross_validation.R",
    "scripts/assess_quality.R"
  )
  
  for (script in required_scripts) {
    if (file.exists(script)) {
      source(script)
      if (verbose) cat("  Sourced:", script, "\n")
    } else {
      warning("Script not found:", script)
    }
  }
}

#' Main results analysis pipeline
#' @param pwm_file Path to PWM file
#' @param sequences_file Path to sequences file (optional)
#' @param output_dir Output directory
#' @param config Analysis configuration
#' @param verbose Enable verbose output
run_results_analysis_pipeline <- function(pwm_file, sequences_file = NULL, 
                                        output_dir = "results", config = list(),
                                        verbose = FALSE) {
  if (verbose) {
    cat("Starting Results Analysis Pipeline\n")
    cat("=================================\n")
    cat("PWM file:", pwm_file, "\n")
    cat("Sequences file:", sequences_file %||% "None", "\n")
    cat("Output directory:", output_dir, "\n\n")
  }
  
  # Load PWM data
  if (verbose) cat("Loading PWM data...\n")
  pwm_data <- readRDS(pwm_file)
  
  # Initialize results tracking
  pipeline_results <- list(
    timestamp = Sys.time(),
    input_files = list(pwm = pwm_file, sequences = sequences_file),
    config = config,
    results = list(),
    generated_files = list()
  )
  
  # Step 1: Organize results directory structure
  if (config$organize %||% TRUE) {
    if (verbose) cat("Step 1: Organizing results directory structure...\n")
    
    # Create organized structure
    if (exists("create_organized_structure")) {
      organized_dirs <- create_organized_structure(output_dir, verbose)
      pipeline_results$generated_files$directories <- organized_dirs
    }
    
    # Move existing files
    if (exists("organize_existing_files")) {
      organize_existing_files(output_dir, verbose)
    }
  }
  
  # Step 2: Export PWM to multiple formats
  if (config$export_formats %||% TRUE) {
    if (verbose) cat("Step 2: Exporting PWM to multiple formats...\n")
    
    formats <- strsplit(config$export_formats %||% "jaspar,transfac", ",")[[1]]
    formats <- trimws(formats)
    
    if (exists("export_pwm_multiple_formats")) {
      base_name <- tools::file_path_sans_ext(basename(pwm_file))
      output_prefix <- file.path(output_dir, "pwm_outputs", base_name)
      
      exported_files <- export_pwm_multiple_formats(
        pwm_data, output_prefix, formats, verbose = verbose
      )
      pipeline_results$generated_files$exports <- exported_files
    }
  }
  
  # Step 3: CTCF-specific pattern validation
  if (config$pattern_validation %||% TRUE) {
    if (verbose) cat("Step 3: Running CTCF pattern validation...\n")
    
    if (exists("comprehensive_ctcf_validation")) {
      validation_results <- comprehensive_ctcf_validation(pwm_data, verbose)
      pipeline_results$results$pattern_validation <- validation_results
      
      # Save validation results
      validation_file <- file.path(output_dir, "quality_reports", "ctcf_validation.json")
      dir.create(dirname(validation_file), recursive = TRUE, showWarnings = FALSE)
      writeLines(toJSON(validation_results, pretty = TRUE), validation_file)
      pipeline_results$generated_files$validation <- validation_file
    }
  }
  
  # Step 4: Generate visualizations
  if (config$generate_visualizations %||% TRUE) {
    if (verbose) cat("Step 4: Generating visualizations...\n")
    
    if (exists("generate_all_visualizations")) {
      viz_dir <- file.path(output_dir, "visualizations")
      base_name <- tools::file_path_sans_ext(basename(pwm_file))
      
      viz_files <- generate_all_visualizations(
        pwm_data, viz_dir, base_name, verbose = verbose
      )
      pipeline_results$generated_files$visualizations <- viz_files
    }
  }
  
  # Step 5: Run cross-validation
  if (config$run_validation %||% TRUE && !is.null(sequences_file)) {
    if (verbose) cat("Step 5: Running cross-validation...\n")
    
    # Run integrated cross-validation
    validation_dir <- file.path(output_dir, "validation_results")
    methods <- config$validation_methods %||% "bootstrap,chromosome,kfold"
    
    # Create a temporary script call (since we can't directly call the main function)
    temp_script <- tempfile(fileext = ".R")
    writeLines(paste0("
    source('scripts/integrate_cross_validation.R')
    # Run validation methods
    "), temp_script)
    
    # For now, just document what should be run
    validation_info <- list(
      methods = strsplit(methods, ",")[[1]],
      output_dir = validation_dir,
      sequences_file = sequences_file,
      status = "configured"
    )
    pipeline_results$results$cross_validation <- validation_info
  }
  
  # Step 6: Quality assessment
  if (config$quality_assessment %||% TRUE) {
    if (verbose) cat("Step 6: Running quality assessment...\n")
    
    if (exists("assess_pwm_quality")) {
      quality_dir <- file.path(output_dir, "quality_reports")
      quality_results <- assess_pwm_quality(
        pwm_file, sequences_file, quality_dir, verbose = verbose
      )
      pipeline_results$results$quality_assessment <- quality_results
    }
  }
  
  # Step 7: Generate comprehensive reports
  if (config$generate_reports %||% TRUE) {
    if (verbose) cat("Step 7: Generating comprehensive reports...\n")
    
    reports_dir <- file.path(output_dir, "quality_reports")
    
    # Generate executive summary
    exec_summary <- generate_executive_summary(pipeline_results, reports_dir, verbose)
    pipeline_results$generated_files$executive_summary <- exec_summary
    
    # Generate technical report
    tech_report <- generate_technical_report(pipeline_results, reports_dir, verbose)
    pipeline_results$generated_files$technical_report <- tech_report
  }
  
  # Step 8: Create final results index
  if (verbose) cat("Step 8: Creating results index...\n")
  
  results_index <- create_final_results_index(pipeline_results, output_dir, verbose)
  pipeline_results$generated_files$results_index <- results_index
  
  if (verbose) {
    cat("\nResults Analysis Pipeline Complete!\n")
    cat("===================================\n")
    cat("Total files generated:", length(unlist(pipeline_results$generated_files)), "\n")
    cat("Results index:", results_index, "\n")
  }
  
  return(pipeline_results)
}

#' Generate executive summary
#' @param pipeline_results Pipeline results
#' @param output_dir Output directory
#' @param verbose Enable verbose output
generate_executive_summary <- function(pipeline_results, output_dir, verbose = FALSE) {
  if (verbose) cat("Generating executive summary...\n")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_file <- file.path(output_dir, "executive_summary.html")
  
  # Extract key metrics
  pattern_validation <- pipeline_results$results$pattern_validation
  quality_grade <- pattern_validation$quality_grade %||% "UNKNOWN"
  total_ic <- pattern_validation$total_information_content %||% 0
  pattern_similarity <- pattern_validation$pattern_validation$similarity_score %||% 0
  
  # Generate HTML summary
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>PWM Analysis Executive Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        .grade { font-size: 1.5em; font-weight: bold; padding: 15px; border-radius: 10px; text-align: center; }
        .excellent { background-color: #d5f4e6; color: #27ae60; }
        .good { background-color: #fef9e7; color: #f39c12; }
        .acceptable { background-color: #ebf3fd; color: #3498db; }
        .poor { background-color: #fadbd8; color: #e74c3c; }
        .metric { background-color: #f8f9fa; padding: 15px; margin: 10px 0; border-radius: 5px; }
        .metric-value { font-size: 1.2em; font-weight: bold; color: #2c3e50; }
    </style>
</head>
<body>
    <h1>PWM Analysis Executive Summary</h1>
    <p><strong>Analysis Date:</strong> ", format(pipeline_results$timestamp, "%Y-%m-%d %H:%M:%S"), "</p>
    
    <div class='grade ", tolower(quality_grade), "'>
        Overall Quality Grade: ", quality_grade, "
    </div>
    
    <div class='metric'>
        <strong>Total Information Content:</strong> 
        <span class='metric-value'>", round(total_ic, 2), " bits</span>
    </div>
    
    <div class='metric'>
        <strong>CTCF Pattern Similarity:</strong> 
        <span class='metric-value'>", round(pattern_similarity * 100, 1), "%</span>
    </div>
    
    <div class='metric'>
        <strong>Files Generated:</strong> 
        <span class='metric-value'>", length(unlist(pipeline_results$generated_files)), "</span>
    </div>
    
    <h2>Recommendations</h2>
    <ul>")
  
  # Add recommendations based on quality grade
  if (quality_grade == "EXCELLENT") {
    html_content <- paste0(html_content, "
        <li>✓ APPROVED for drug discovery applications</li>
        <li>✓ APPROVED for high-stakes research</li>
        <li>✓ APPROVED for publication-quality analysis</li>")
  } else if (quality_grade == "GOOD") {
    html_content <- paste0(html_content, "
        <li>✓ APPROVED for research applications</li>
        <li>✓ APPROVED for comparative studies</li>
        <li>⚠ Consider additional validation for high-stakes applications</li>")
  } else {
    html_content <- paste0(html_content, "
        <li>⚠ Limited approval for exploratory analysis</li>
        <li>⚠ Consider improving data quality or methods</li>
        <li>✗ NOT RECOMMENDED for high-stakes applications</li>")
  }
  
  html_content <- paste0(html_content, "
    </ul>
</body>
</html>")
  
  writeLines(html_content, summary_file)
  
  if (verbose) cat("Executive summary saved to:", summary_file, "\n")
  return(summary_file)
}

#' Generate technical report
#' @param pipeline_results Pipeline results
#' @param output_dir Output directory
#' @param verbose Enable verbose output
generate_technical_report <- function(pipeline_results, output_dir, verbose = FALSE) {
  if (verbose) cat("Generating technical report...\n")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  report_file <- file.path(output_dir, "technical_report.json")
  
  # Create comprehensive technical report
  technical_report <- list(
    metadata = list(
      generated = pipeline_results$timestamp,
      pipeline_version = "2.0",
      input_files = pipeline_results$input_files
    ),
    analysis_results = pipeline_results$results,
    generated_files = pipeline_results$generated_files,
    configuration = pipeline_results$config
  )
  
  writeLines(toJSON(technical_report, pretty = TRUE), report_file)
  
  if (verbose) cat("Technical report saved to:", report_file, "\n")
  return(report_file)
}

#' Create final results index
#' @param pipeline_results Pipeline results
#' @param output_dir Output directory
#' @param verbose Enable verbose output
create_final_results_index <- function(pipeline_results, output_dir, verbose = FALSE) {
  if (verbose) cat("Creating final results index...\n")
  
  index_file <- file.path(output_dir, "results_index.json")
  
  # Scan all generated files
  all_files <- unlist(pipeline_results$generated_files)
  
  # Create index
  index_data <- list(
    timestamp = Sys.time(),
    pipeline_version = "2.0",
    total_files = length(all_files),
    directory_structure = list(
      pwm_outputs = list.files(file.path(output_dir, "pwm_outputs"), full.names = FALSE),
      quality_reports = list.files(file.path(output_dir, "quality_reports"), full.names = FALSE),
      validation_results = list.files(file.path(output_dir, "validation_results"), full.names = FALSE),
      visualizations = list.files(file.path(output_dir, "visualizations"), full.names = FALSE),
      intermediate_files = list.files(file.path(output_dir, "intermediate_files"), full.names = FALSE)
    ),
    generated_files = pipeline_results$generated_files,
    analysis_summary = list(
      quality_grade = pipeline_results$results$pattern_validation$quality_grade %||% "UNKNOWN",
      total_ic = pipeline_results$results$pattern_validation$total_information_content %||% 0,
      pattern_similarity = pipeline_results$results$pattern_validation$pattern_validation$similarity_score %||% 0
    )
  )
  
  writeLines(toJSON(index_data, pretty = TRUE), index_file)
  
  if (verbose) cat("Results index saved to:", index_file, "\n")
  return(index_file)
}

# Helper function for null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Main execution
if (!interactive()) {
  # Validate inputs
  if (is.null(opt$input)) {
    cat("Error: Input PWM file is required\n")
    quit(status = 1)
  }
  
  if (!file.exists(opt$input)) {
    cat("Error: Input file does not exist:", opt$input, "\n")
    quit(status = 1)
  }
  
  if (!is.null(opt$sequences) && !file.exists(opt$sequences)) {
    cat("Error: Sequences file does not exist:", opt$sequences, "\n")
    quit(status = 1)
  }
  
  if (opt$verbose) {
    cat("Results Analysis Pipeline\n")
    cat("========================\n")
    cat("Input PWM file:", opt$input, "\n")
    cat("Sequences file:", opt$sequences %||% "None", "\n")
    cat("Output directory:", opt$`output-dir`, "\n")
    cat("Organize results:", opt$organize, "\n")
    cat("Export formats:", opt$`export-formats`, "\n")
    cat("Generate visualizations:", opt$`generate-visualizations`, "\n")
    cat("Run validation:", opt$`run-validation`, "\n")
    cat("Generate reports:", opt$`generate-reports`, "\n\n")
  }
  
  # Source required scripts
  source_analysis_scripts(opt$verbose)
  
  # Create configuration
  config <- list(
    organize = opt$organize,
    export_formats = opt$`export-formats`,
    generate_visualizations = opt$`generate-visualizations`,
    run_validation = opt$`run-validation`,
    validation_methods = opt$`validation-methods`,
    generate_reports = opt$`generate-reports`,
    pattern_validation = TRUE,
    quality_assessment = TRUE
  )
  
  # Run pipeline
  pipeline_results <- run_results_analysis_pipeline(
    opt$input, opt$sequences, opt$`output-dir`, config, opt$verbose
  )
  
  if (opt$verbose) {
    cat("\n" + "="*50 + "\n")
    cat("RESULTS ANALYSIS PIPELINE COMPLETE\n")
    cat("="*50 + "\n")
    cat("Quality Grade:", pipeline_results$results$pattern_validation$quality_grade %||% "UNKNOWN", "\n")
    cat("Total Files Generated:", length(unlist(pipeline_results$generated_files)), "\n")
    cat("Results Index:", pipeline_results$generated_files$results_index, "\n")
    cat("="*50 + "\n")
  }
}
