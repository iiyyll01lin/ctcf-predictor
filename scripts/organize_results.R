#!/usr/bin/env Rscript

# Results Organization Script
# Organizes PWM analysis results into structured directory layout
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(optparse)
})

# Command line options
option_list <- list(
  make_option(c("-r", "--results-dir"), type="character", default="results",
              help="Results directory to organize [default: %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Organize PWM analysis results into structured directories")
opt <- parse_args(opt_parser)

# Create organized directory structure
create_organized_structure <- function(base_dir, verbose = FALSE) {
  if (verbose) cat("Creating organized directory structure in:", base_dir, "\n")
  
  # Define target directory structure
  target_dirs <- c(
    "pwm_outputs",
    "quality_reports", 
    "validation_results",
    "visualizations",
    "intermediate_files"
  )
  
  # Create directories
  for (dir in target_dirs) {
    target_path <- file.path(base_dir, dir)
    if (!dir.exists(target_path)) {
      dir.create(target_path, recursive = TRUE)
      if (verbose) cat("Created directory:", target_path, "\n")
    }
  }
  
  return(file.path(base_dir, target_dirs))
}

# Move files to organized structure
organize_existing_files <- function(base_dir, verbose = FALSE) {
  if (verbose) cat("Organizing existing files...\n")
  
  # Define file patterns and their target directories
  file_mappings <- list(
    pwm_outputs = c("*.rds", "*.meme", "*.jaspar", "*.transfac", "*pwm*metadata*.json"),
    quality_reports = c("*quality*.html", "*comparison*.html", "*assessment*.html"),
    validation_results = c("*validation*.html", "*cross_validation*.html", 
                          "*bootstrap*.html", "*statistical*.html", "*null*.html"),
    visualizations = c("*.pdf", "*.png", "*.svg", "*logo*", "*plot*", "*heatmap*"),
    intermediate_files = c("aligned_*.fa", "aligned_*.fasta", "frequency_*.txt", 
                          "processing_*.txt", "*_log.txt", "*.log")
  )
  
  # Move files based on patterns
  for (target_dir in names(file_mappings)) {
    target_path <- file.path(base_dir, target_dir)
    patterns <- file_mappings[[target_dir]]
    
    for (pattern in patterns) {
      files <- list.files(base_dir, pattern = glob2rx(pattern), full.names = TRUE)
      files <- files[!file.info(files)$isdir] # Exclude directories
      
      for (file in files) {
        if (file.exists(file)) {
          dest_file <- file.path(target_path, basename(file))
          if (!file.exists(dest_file)) {
            file.copy(file, dest_file)
            file.remove(file)
            if (verbose) cat("Moved:", basename(file), "to", target_dir, "\n")
          }
        }
      }
    }
  }
}

# Create index file for organized results
create_results_index <- function(base_dir, verbose = FALSE) {
  if (verbose) cat("Creating results index...\n")
  
  index_file <- file.path(base_dir, "results_index.json")
  
  # Scan organized directories
  index_data <- list(
    timestamp = Sys.time(),
    structure = list(),
    summary = list()
  )
  
  subdirs <- c("pwm_outputs", "quality_reports", "validation_results", 
               "visualizations", "intermediate_files")
  
  for (subdir in subdirs) {
    subdir_path <- file.path(base_dir, subdir)
    if (dir.exists(subdir_path)) {
      files <- list.files(subdir_path, full.names = FALSE)
      index_data$structure[[subdir]] <- files
      index_data$summary[[paste0(subdir, "_count")]] <- length(files)
    }
  }
  
  # Write index
  writeLines(jsonlite::toJSON(index_data, pretty = TRUE), index_file)
  if (verbose) cat("Created results index:", index_file, "\n")
  
  return(index_file)
}

# Main execution
if (!interactive()) {
  # Validate inputs
  if (!dir.exists(opt$`results-dir`)) {
    cat("Error: Results directory does not exist:", opt$`results-dir`, "\n")
    quit(status = 1)
  }
  
  if (opt$verbose) {
    cat("Results Organization Script\n")
    cat("==========================\n")
    cat("Results directory:", opt$`results-dir`, "\n\n")
  }
  
  # Create organized structure
  created_dirs <- create_organized_structure(opt$`results-dir`, opt$verbose)
  
  # Move existing files
  organize_existing_files(opt$`results-dir`, opt$verbose)
  
  # Create index
  index_file <- create_results_index(opt$`results-dir`, opt$verbose)
  
  if (opt$verbose) {
    cat("\nOrganization complete!\n")
    cat("Created directories:", paste(basename(created_dirs), collapse = ", "), "\n")
    cat("Results index:", index_file, "\n")
  }
}
