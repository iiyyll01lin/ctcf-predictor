# Sequence Extraction Script - Peak-to-sequence conversion
# This script extracts sequences from genomic coordinates using bedtools

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Please install Biostrings: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
bed_file <- if (length(args) >= 1) args[1] else "data/K562_CTCF_peaks.bed"
genome_file <- if (length(args) >= 2) args[2] else "data/reference_genome/hg38.fa"
output_file <- if (length(args) >= 3) args[3] else "data/extracted_sequences.fasta"
extend_bp <- if (length(args) >= 4) as.numeric(args[4]) else 0

cat("Sequence Extraction Pipeline\n")
cat("============================\n")
cat("BED file:", bed_file, "\n")
cat("Genome file:", genome_file, "\n")
cat("Output file:", output_file, "\n")
cat("Extension:", extend_bp, "bp\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Functions ---

# Check if bedtools is available
check_bedtools <- function() {
  cat("Checking bedtools availability...\n")
  
  # Try to run bedtools version
  result <- tryCatch({
    system("bedtools --version", intern = TRUE, ignore.stderr = TRUE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(result)) {
    cat("Warning: bedtools not found in PATH\n")
    cat("Please install bedtools: https://bedtools.readthedocs.io/\n")
    return(FALSE)
  } else {
    cat("bedtools found:", result[1], "\n")
    return(TRUE)
  }
}

# Validate BED file format
validate_bed_file <- function(bed_file) {
  cat("Validating BED file format...\n")
  
  if (!file.exists(bed_file)) {
    stop("BED file not found: ", bed_file)
  }
  
  # Read first few lines to check format
  lines <- readLines(bed_file, n = 10)
  
  # Check if it looks like a BED file
  valid_lines <- 0
  for (line in lines) {
    if (startsWith(line, "#") || line == "") next  # Skip comments and empty lines
    
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) >= 3) {
      # Check if positions are numeric
      if (!is.na(as.numeric(fields[2])) && !is.na(as.numeric(fields[3]))) {
        valid_lines <- valid_lines + 1
      }
    }
  }
  
  if (valid_lines == 0) {
    stop("Invalid BED file format in: ", bed_file)
  }
  
  cat("BED file validation passed (", valid_lines, "valid lines checked)\n")
  return(TRUE)
}

# Count peaks in BED file
count_bed_peaks <- function(bed_file) {
  lines <- readLines(bed_file)
  # Count non-comment, non-empty lines
  peaks <- sum(!startsWith(lines, "#") & lines != "")
  return(peaks)
}

# Extract sequences using bedtools
extract_with_bedtools <- function(bed_file, genome_file, output_file, extend = 0) {
  cat("Extracting sequences with bedtools...\n")
  
  # Prepare bedtools command
  if (extend > 0) {
    # Create extended BED file first
    temp_bed <- tempfile(fileext = ".bed")
    extend_cmd <- paste("bedtools slop -i", bed_file, "-g", 
                       paste0(dirname(genome_file), "/genome.sizes"), 
                       "-b", extend, ">", temp_bed)
    
    cat("Extending regions by", extend, "bp\n")
    system(extend_cmd)
    bed_input <- temp_bed
  } else {
    bed_input <- bed_file
  }
  
  # Extract sequences
  getfasta_cmd <- paste("bedtools getfasta -fi", genome_file, 
                       "-bed", bed_input, "-fo", output_file)
  
  cat("Running command:", getfasta_cmd, "\n")
  result <- system(getfasta_cmd)
  
  # Clean up temporary file
  if (extend > 0 && file.exists(temp_bed)) {
    unlink(temp_bed)
  }
  
  if (result != 0) {
    stop("bedtools getfasta failed with exit code: ", result)
  }
  
  return(TRUE)
}

# Alternative R-based extraction (if bedtools not available)
extract_with_r <- function(bed_file, genome_file, output_file) {
  cat("Using R-based sequence extraction (slower)...\n")
  
  # This is a simplified implementation
  # For production use, consider using rtracklayer or GenomicRanges
  cat("R-based extraction not fully implemented.\n")
  cat("Please install bedtools for optimal performance.\n")
  return(FALSE)
}

# Validate extracted sequences
validate_extracted_sequences <- function(output_file) {
  cat("Validating extracted sequences...\n")
  
  if (!file.exists(output_file)) {
    stop("Output file not created: ", output_file)
  }
  
  # Try to read sequences
  tryCatch({
    sequences <- readDNAStringSet(output_file)
    cat("Successfully extracted", length(sequences), "sequences\n")
    
    if (length(sequences) > 0) {
      lengths <- width(sequences)
      cat("Sequence lengths: min =", min(lengths), "bp, max =", max(lengths), 
          "bp, mean =", round(mean(lengths), 1), "bp\n")
    }
    
    return(list(
      count = length(sequences),
      lengths = width(sequences)
    ))
  }, error = function(e) {
    stop("Failed to read extracted sequences: ", e$message)
  })
}

# --- Main Processing ---

# Validate inputs
if (!file.exists(bed_file)) {
  stop("BED file not found: ", bed_file)
}

if (!file.exists(genome_file)) {
  stop("Genome file not found: ", genome_file)
}

# Validate BED file
validate_bed_file(bed_file)

# Count input peaks
peak_count <- count_bed_peaks(bed_file)
cat("Input peaks:", peak_count, "\n\n")

# Create output directory
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Check for bedtools
has_bedtools <- check_bedtools()

if (has_bedtools) {
  # Use bedtools for extraction
  success <- extract_with_bedtools(bed_file, genome_file, output_file, extend_bp)
} else {
  # Fallback to R-based extraction
  success <- extract_with_r(bed_file, genome_file, output_file)
}

if (success) {
  # Validate extracted sequences
  extraction_stats <- validate_extracted_sequences(output_file)
  
  # Generate extraction report
  report_file <- gsub("\\.fasta$", "_extraction_report.txt", output_file)
  
  sink(report_file)
  cat("Sequence Extraction Report\n")
  cat("==========================\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("BED file:", bed_file, "\n")
  cat("Genome file:", genome_file, "\n")
  cat("Output file:", output_file, "\n")
  cat("Extension:", extend_bp, "bp\n\n")
  
  cat("Extraction Results:\n")
  cat("  Input peaks:", peak_count, "\n")
  cat("  Extracted sequences:", extraction_stats$count, "\n")
  cat("  Success rate:", round(extraction_stats$count / peak_count * 100, 1), "%\n\n")
  
  if (length(extraction_stats$lengths) > 0) {
    cat("Sequence Statistics:\n")
    cat("  Min length:", min(extraction_stats$lengths), "bp\n")
    cat("  Max length:", max(extraction_stats$lengths), "bp\n")
    cat("  Mean length:", round(mean(extraction_stats$lengths), 1), "bp\n")
    cat("  Median length:", median(extraction_stats$lengths), "bp\n")
  }
  sink()
  
  cat("\nSequence extraction completed successfully!\n")
  cat("Extracted", extraction_stats$count, "sequences saved to", output_file, "\n")
  cat("Report saved to", report_file, "\n")
} else {
  stop("Sequence extraction failed")
}
