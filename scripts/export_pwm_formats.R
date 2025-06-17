#!/usr/bin/env Rscript

# PWM Format Export Functions
# Exports PWM data to various standard formats (JASPAR, TRANSFAC)
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input PWM file (RDS format)", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file path", metavar="file"),
  make_option(c("-f", "--format"), type="character", default="jaspar",
              help="Output format: jaspar, transfac [default: %default]"),
  make_option(c("-n", "--name"), type="character", default="CTCF",
              help="Motif name [default: %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Export PWM to standard formats")
opt <- parse_args(opt_parser)

#' Export PWM to JASPAR format
#' @param pwm_data PWM matrix or PWM result object
#' @param output_file Output file path
#' @param motif_name Name for the motif
#' @param motif_id Motif ID (optional)
#' @param verbose Enable verbose output
export_to_jaspar <- function(pwm_data, output_file, motif_name = "CTCF", 
                            motif_id = NULL, verbose = FALSE) {
  if (verbose) cat("Exporting to JASPAR format:", output_file, "\n")
  
  # Extract PWM matrix
  if (is.list(pwm_data) && "pwm" %in% names(pwm_data)) {
    pwm_matrix <- pwm_data$pwm
  } else if (is.matrix(pwm_data)) {
    pwm_matrix <- pwm_data
  } else {
    stop("Invalid PWM data format")
  }
  
  # Validate matrix
  if (!is.matrix(pwm_matrix) || nrow(pwm_matrix) != 4) {
    stop("PWM must be a 4xN matrix (A, C, G, T)")
  }
  
  # Ensure row names are nucleotides
  if (is.null(rownames(pwm_matrix))) {
    rownames(pwm_matrix) <- c("A", "C", "G", "T")
  }
  
  # Generate motif ID if not provided
  if (is.null(motif_id)) {
    motif_id <- paste0(gsub("[^A-Za-z0-9]", "_", motif_name), "_", 
                      format(Sys.Date(), "%Y%m%d"))
  }
  
  # Create JASPAR format content
  jaspar_content <- c(
    paste0(">", motif_id, "\t", motif_name),
    paste("A", paste(sprintf("%.3f", pwm_matrix["A",]), collapse="\t")),
    paste("C", paste(sprintf("%.3f", pwm_matrix["C",]), collapse="\t")),
    paste("G", paste(sprintf("%.3f", pwm_matrix["G",]), collapse="\t")),
    paste("T", paste(sprintf("%.3f", pwm_matrix["T",]), collapse="\t"))
  )
  
  # Write to file
  writeLines(jaspar_content, output_file)
  
  if (verbose) {
    cat("JASPAR export complete:\n")
    cat("  Motif ID:", motif_id, "\n")
    cat("  Motif Name:", motif_name, "\n")
    cat("  Matrix size:", nrow(pwm_matrix), "x", ncol(pwm_matrix), "\n")
  }
  
  return(output_file)
}

#' Export PWM to TRANSFAC format
#' @param pwm_data PWM matrix or PWM result object
#' @param output_file Output file path
#' @param motif_name Name for the motif
#' @param motif_id Motif ID (optional)
#' @param verbose Enable verbose output
export_to_transfac <- function(pwm_data, output_file, motif_name = "CTCF", 
                              motif_id = NULL, verbose = FALSE) {
  if (verbose) cat("Exporting to TRANSFAC format:", output_file, "\n")
  
  # Extract PWM matrix
  if (is.list(pwm_data) && "pwm" %in% names(pwm_data)) {
    pwm_matrix <- pwm_data$pwm
  } else if (is.matrix(pwm_data)) {
    pwm_matrix <- pwm_data
  } else {
    stop("Invalid PWM data format")
  }
  
  # Validate matrix
  if (!is.matrix(pwm_matrix) || nrow(pwm_matrix) != 4) {
    stop("PWM must be a 4xN matrix (A, C, G, T)")
  }
  
  # Ensure row names are nucleotides
  if (is.null(rownames(pwm_matrix))) {
    rownames(pwm_matrix) <- c("A", "C", "G", "T")
  }
  
  # Generate motif ID if not provided
  if (is.null(motif_id)) {
    motif_id <- paste0("M", format(Sys.Date(), "%Y%m%d"), "_001")
  }
  
  # Create TRANSFAC format content
  transfac_content <- c(
    "VV  TRANSFAC PWM Export",
    "XX",
    paste("ID", motif_id),
    paste("NA", motif_name),
    paste("DE", paste("CTCF Position Weight Matrix -", motif_name)),
    "BF  Species: Homo sapiens",
    paste("P0", paste(c("A", "C", "G", "T"), collapse="\t")),
    "",
    # Matrix positions
    sapply(1:ncol(pwm_matrix), function(i) {
      values <- sprintf("%8.3f", pwm_matrix[,i])
      paste(sprintf("%02d", i), paste(values, collapse="\t"))
    }),
    "",
    "XX",
    "//"
  )
  
  # Write to file
  writeLines(transfac_content, output_file)
  
  if (verbose) {
    cat("TRANSFAC export complete:\n")
    cat("  Motif ID:", motif_id, "\n")
    cat("  Motif Name:", motif_name, "\n")
    cat("  Matrix size:", nrow(pwm_matrix), "x", ncol(pwm_matrix), "\n")
  }
  
  return(output_file)
}

#' Export PWM to multiple formats
#' @param pwm_data PWM matrix or PWM result object
#' @param output_prefix Output file prefix
#' @param formats Vector of formats to export
#' @param motif_name Name for the motif
#' @param verbose Enable verbose output
export_pwm_multiple_formats <- function(pwm_data, output_prefix, 
                                       formats = c("jaspar", "transfac"),
                                       motif_name = "CTCF", verbose = FALSE) {
  if (verbose) cat("Exporting PWM to multiple formats...\n")
  
  exported_files <- list()
  
  for (format in formats) {
    output_file <- paste0(output_prefix, ".", format)
    
    if (format == "jaspar") {
      exported_files$jaspar <- export_to_jaspar(pwm_data, output_file, 
                                               motif_name, verbose = verbose)
    } else if (format == "transfac") {
      exported_files$transfac <- export_to_transfac(pwm_data, output_file, 
                                                   motif_name, verbose = verbose)
    } else {
      warning("Unknown format:", format)
    }
  }
  
  return(exported_files)
}

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
  
  if (is.null(opt$output)) {
    # Generate output filename based on input and format
    base_name <- tools::file_path_sans_ext(basename(opt$input))
    opt$output <- paste0(base_name, ".", opt$format)
  }
  
  if (opt$verbose) {
    cat("PWM Format Export\n")
    cat("=================\n")
    cat("Input file:", opt$input, "\n")
    cat("Output file:", opt$output, "\n")
    cat("Format:", opt$format, "\n")
    cat("Motif name:", opt$name, "\n\n")
  }
  
  # Load PWM data
  if (opt$verbose) cat("Loading PWM data...\n")
  pwm_data <- readRDS(opt$input)
  
  # Export based on format
  if (opt$format == "jaspar") {
    export_to_jaspar(pwm_data, opt$output, opt$name, verbose = opt$verbose)
  } else if (opt$format == "transfac") {
    export_to_transfac(pwm_data, opt$output, opt$name, verbose = opt$verbose)
  } else {
    cat("Error: Unknown format:", opt$format, "\n")
    quit(status = 1)
  }
  
  if (opt$verbose) {
    cat("\nExport complete!\n")
    cat("Output file:", opt$output, "\n")
  }
}
