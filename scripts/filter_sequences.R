# Sequence Quality Filtering Script - Quality control and filtering
# This script implements comprehensive sequence filtering based on quality metrics

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Please install Biostrings: if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/extracted_sequences.fasta"
output_file <- if (length(args) >= 2) args[2] else "data/filtered_sequences.fasta"
max_n_ratio <- if (length(args) >= 3) as.numeric(args[3]) else 0.01
min_complexity <- if (length(args) >= 4) as.numeric(args[4]) else 1.0
remove_repeats <- if (length(args) >= 5) as.logical(args[5]) else TRUE

cat("Sequence Quality Filtering Pipeline\n")
cat("===================================\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Max N ratio:", max_n_ratio * 100, "%\n")
cat("Min complexity:", min_complexity, "\n")
cat("Remove repeats:", remove_repeats, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Functions ---

# Calculate sequence complexity (Shannon entropy)
calculate_complexity <- function(sequence) {
  # Convert to character and calculate base frequencies
  seq_char <- strsplit(as.character(sequence), "")[[1]]
  base_table <- table(seq_char)
  base_freqs <- base_table / sum(base_table)
  
  # Calculate Shannon entropy
  entropy <- -sum(base_freqs * log2(base_freqs + 1e-10))  # Add small value to avoid log(0)
  
  return(entropy)
}

# Detect simple repeats
detect_repeats <- function(sequence, min_repeat_length = 3, max_repeat_unit = 10) {
  seq_char <- as.character(sequence)
  seq_length <- nchar(seq_char)
  
  # Check for different repeat unit sizes
  for (unit_size in 1:min(max_repeat_unit, seq_length %/% 2)) {
    for (start_pos in 1:(seq_length - unit_size * min_repeat_length + 1)) {
      unit <- substr(seq_char, start_pos, start_pos + unit_size - 1)
      
      # Count consecutive repeats
      repeat_count <- 1
      pos <- start_pos + unit_size
      
      while (pos + unit_size - 1 <= seq_length) {
        if (substr(seq_char, pos, pos + unit_size - 1) == unit) {
          repeat_count <- repeat_count + 1
          pos <- pos + unit_size
        } else {
          break
        }
      }
      
      if (repeat_count >= min_repeat_length) {
        repeat_length <- repeat_count * unit_size
        repeat_ratio <- repeat_length / seq_length
        return(list(
          has_repeat = TRUE,
          unit = unit,
          count = repeat_count,
          ratio = repeat_ratio,
          start = start_pos,
          length = repeat_length
        ))
      }
    }
  }
  
  return(list(has_repeat = FALSE, ratio = 0))
}

# Calculate N-base ratio
calculate_n_ratio <- function(sequence) {
  n_count <- letterFrequency(sequence, "N", as.prob = FALSE)
  total_length <- width(sequence)
  return(as.numeric(n_count / total_length))
}

# Filter sequences by quality criteria
filter_by_quality <- function(sequences, max_n_ratio = 0.01, min_complexity = 1.0, 
                              remove_repeats = TRUE, max_repeat_ratio = 0.5) {
  
  cat("Applying quality filters...\n")
  
  n_sequences <- length(sequences)
  passed_filters <- rep(TRUE, n_sequences)
  filter_stats <- list(
    n_filter = 0,
    complexity_filter = 0,
    repeat_filter = 0
  )
  
  for (i in 1:n_sequences) {
    if (i %% 1000 == 0) {
      cat("Processing sequence", i, "/", n_sequences, "\n")
    }
    
    seq <- sequences[i]
    
    # Filter 1: N-base ratio
    n_ratio <- calculate_n_ratio(seq)
    if (n_ratio > max_n_ratio) {
      passed_filters[i] <- FALSE
      filter_stats$n_filter <- filter_stats$n_filter + 1
      next
    }
    
    # Filter 2: Sequence complexity
    complexity <- calculate_complexity(seq)
    if (complexity < min_complexity) {
      passed_filters[i] <- FALSE
      filter_stats$complexity_filter <- filter_stats$complexity_filter + 1
      next
    }
    
    # Filter 3: Simple repeats
    if (remove_repeats) {
      repeat_info <- detect_repeats(seq)
      if (repeat_info$has_repeat && repeat_info$ratio > max_repeat_ratio) {
        passed_filters[i] <- FALSE
        filter_stats$repeat_filter <- filter_stats$repeat_filter + 1
        next
      }
    }
  }
  
  # Report filtering results
  cat("\nFiltering results:\n")
  cat("  Sequences failed N-base filter:", filter_stats$n_filter, "\n")
  cat("  Sequences failed complexity filter:", filter_stats$complexity_filter, "\n")
  cat("  Sequences failed repeat filter:", filter_stats$repeat_filter, "\n")
  cat("  Total sequences passing filters:", sum(passed_filters), "/", n_sequences, "\n")
  
  return(list(
    sequences = sequences[passed_filters],
    stats = filter_stats,
    passed_count = sum(passed_filters)
  ))
}

# Calculate quality metrics for reporting
calculate_quality_metrics <- function(sequences) {
  n_sequences <- length(sequences)
  
  if (n_sequences == 0) {
    return(list(
      count = 0,
      length_stats = NULL,
      gc_content = NULL,
      n_content = NULL,
      complexity_stats = NULL
    ))
  }
  
  # Length statistics
  lengths <- width(sequences)
  
  # GC content
  gc_content <- letterFrequency(sequences, "GC", as.prob = TRUE)
  
  # N content
  n_content <- letterFrequency(sequences, "N", as.prob = TRUE)
  
  # Complexity (sample first 1000 sequences for speed)
  sample_size <- min(1000, n_sequences)
  sample_idx <- sample(n_sequences, sample_size)
  complexity_values <- sapply(sample_idx, function(i) calculate_complexity(sequences[i]))
  
  return(list(
    count = n_sequences,
    length_stats = list(
      mean = mean(lengths),
      median = median(lengths),
      min = min(lengths),
      max = max(lengths),
      sd = sd(lengths)
    ),
    gc_content = list(
      mean = mean(gc_content),
      median = median(gc_content),
      min = min(gc_content),
      max = max(gc_content)
    ),
    n_content = list(
      mean = mean(n_content),
      max = max(n_content),
      sequences_with_n = sum(n_content > 0)
    ),
    complexity_stats = list(
      mean = mean(complexity_values),
      median = median(complexity_values),
      min = min(complexity_values),
      max = max(complexity_values),
      sample_size = sample_size
    )
  ))
}

# --- Main Processing ---

# Check input file
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# Load sequences
cat("Loading sequences from:", input_file, "\n")
raw_sequences <- readDNAStringSet(input_file)
cat("Raw sequences loaded:", length(raw_sequences), "\n\n")

if (length(raw_sequences) == 0) {
  stop("No sequences found in input file")
}

# Calculate initial quality metrics
cat("Calculating initial quality metrics...\n")
initial_metrics <- calculate_quality_metrics(raw_sequences)

cat("Initial dataset statistics:\n")
cat("  Count:", initial_metrics$count, "\n")
if (!is.null(initial_metrics$length_stats)) {
  cat("  Length: mean =", round(initial_metrics$length_stats$mean, 1), 
      "bp, range =", initial_metrics$length_stats$min, "-", initial_metrics$length_stats$max, "bp\n")
  cat("  GC content: mean =", round(initial_metrics$gc_content$mean * 100, 1), "%\n")
  cat("  N content: mean =", round(initial_metrics$n_content$mean * 100, 2), 
      "%, max =", round(initial_metrics$n_content$max * 100, 2), "%\n")
  cat("  Complexity: mean =", round(initial_metrics$complexity_stats$mean, 2), 
      ", range =", round(initial_metrics$complexity_stats$min, 2), "-", 
      round(initial_metrics$complexity_stats$max, 2), "\n\n")
}

# Apply quality filters
filtering_result <- filter_by_quality(raw_sequences, max_n_ratio, min_complexity, remove_repeats)
filtered_sequences <- filtering_result$sequences

# Check if any sequences passed filters
if (length(filtered_sequences) == 0) {
  stop("No sequences passed quality filters. Consider relaxing filter criteria.")
}

# Calculate final quality metrics
cat("\nCalculating final quality metrics...\n")
final_metrics <- calculate_quality_metrics(filtered_sequences)

cat("Final dataset statistics:\n")
cat("  Count:", final_metrics$count, "\n")
cat("  Length: mean =", round(final_metrics$length_stats$mean, 1), 
    "bp, range =", final_metrics$length_stats$min, "-", final_metrics$length_stats$max, "bp\n")
cat("  GC content: mean =", round(final_metrics$gc_content$mean * 100, 1), "%\n")
cat("  N content: mean =", round(final_metrics$n_content$mean * 100, 2), 
    "%, max =", round(final_metrics$n_content$max * 100, 2), "%\n")
cat("  Complexity: mean =", round(final_metrics$complexity_stats$mean, 2), 
    ", range =", round(final_metrics$complexity_stats$min, 2), "-", 
    round(final_metrics$complexity_stats$max, 2), "\n")

# Create output directory if needed
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save filtered sequences
cat("\nSaving filtered sequences to:", output_file, "\n")
writeXStringSet(filtered_sequences, output_file)

# Generate filtering report
report_file <- gsub("\\.fasta$", "_filtering_report.txt", output_file)
cat("Generating filtering report:", report_file, "\n")

sink(report_file)
cat("Sequence Quality Filtering Report\n")
cat("=================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n\n")

cat("Filtering Parameters:\n")
cat("  Max N ratio:", max_n_ratio * 100, "%\n")
cat("  Min complexity:", min_complexity, "\n")
cat("  Remove repeats:", remove_repeats, "\n\n")

cat("Filtering Results:\n")
cat("  Initial sequences:", initial_metrics$count, "\n")
cat("  Final sequences:", final_metrics$count, "\n")
cat("  Retention rate:", round(final_metrics$count / initial_metrics$count * 100, 1), "%\n\n")

cat("Filter Statistics:\n")
cat("  N-base filter failures:", filtering_result$stats$n_filter, "\n")
cat("  Complexity filter failures:", filtering_result$stats$complexity_filter, "\n")
cat("  Repeat filter failures:", filtering_result$stats$repeat_filter, "\n")

if (!is.null(final_metrics$length_stats)) {
  cat("\nFinal Dataset Quality:\n")
  cat("  Mean length:", round(final_metrics$length_stats$mean, 1), "bp\n")
  cat("  Length range:", final_metrics$length_stats$min, "-", final_metrics$length_stats$max, "bp\n")
  cat("  Mean GC content:", round(final_metrics$gc_content$mean * 100, 1), "%\n")
  cat("  Mean complexity:", round(final_metrics$complexity_stats$mean, 2), "\n")
}

sink()

cat("\nSequence filtering completed successfully!\n")
cat("Filtered", final_metrics$count, "sequences saved to", output_file, "\n")
