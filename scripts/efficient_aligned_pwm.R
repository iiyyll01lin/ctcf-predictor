# Efficient Aligned PWM Builder with Optimization
# This script builds optimized PWMs from aligned sequences with memory and speed optimization

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
aligned_file <- if (length(args) >= 1) args[1] else "data/aligned_sequences.fasta"
output_prefix <- if (length(args) >= 2) args[2] else "results/efficient_aligned_pwm"
batch_size <- if (length(args) >= 3) as.numeric(args[3]) else 10000
optimize_pseudocount <- if (length(args) >= 4) as.logical(args[4]) else TRUE

cat("Efficient Aligned PWM Builder\n")
cat("=============================\n")
cat("Input file:", aligned_file, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("Batch size:", batch_size, "\n")
cat("Optimize pseudocount:", optimize_pseudocount, "\n\n")

# --- Functions ---

# Calculate information content efficiently
calculate_info_content_fast <- function(freq_matrix) {
  # Vectorized calculation
  freq_matrix[freq_matrix == 0] <- 1e-10
  info_content <- 2 + colSums(freq_matrix * log2(freq_matrix))
  return(info_content)
}

# Build consensus matrix in batches for memory efficiency
build_consensus_batched <- function(sequences, batch_size = 10000) {
  total_seqs <- length(sequences)
  cat("Processing", total_seqs, "sequences in batches of", batch_size, "...\n")
  
  # Get the maximum sequence length for matrix initialization
  max_length <- max(width(sequences))
  
  # Initialize count matrix
  nucleotides <- c("A", "C", "G", "T")
  total_counts <- matrix(0, nrow = 4, ncol = max_length)
  rownames(total_counts) <- nucleotides
  
  # Process in batches
  num_batches <- ceiling(total_seqs / batch_size)
  
  for (i in 1:num_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, total_seqs)
    
    cat("Processing batch", i, "of", num_batches, 
        "(sequences", start_idx, "to", end_idx, ")...\n")
    
    batch_sequences <- sequences[start_idx:end_idx]
    batch_consensus <- consensusMatrix(batch_sequences, as.prob = FALSE)
    
    # Add batch counts to total
    for (nuc in nucleotides) {
      if (nuc %in% rownames(batch_consensus)) {
        batch_length <- ncol(batch_consensus)
        total_counts[nuc, 1:batch_length] <- total_counts[nuc, 1:batch_length] + 
                                             batch_consensus[nuc, ]
      }
    }
    
    # Clean up batch data
    rm(batch_sequences, batch_consensus)
    gc()  # Force garbage collection
  }
  
  # Trim matrix to actual length
  actual_length <- max(which(colSums(total_counts) > 0))
  total_counts <- total_counts[, 1:actual_length, drop = FALSE]
  
  return(total_counts)
}

# Optimize pseudocount using cross-validation
optimize_pseudocount_cv <- function(sequences, pseudocount_values = c(0.01, 0.05, 0.1, 0.5, 1.0), 
                                   cv_folds = 5, sample_size = 5000) {
  cat("Optimizing pseudocount using cross-validation...\n")
  
  # Sample sequences for faster optimization if dataset is large
  if (length(sequences) > sample_size) {
    sample_indices <- sample(length(sequences), sample_size)
    sequences <- sequences[sample_indices]
    cat("Using sample of", length(sequences), "sequences for optimization\n")
  }
  
  results <- data.frame(
    pseudocount = pseudocount_values,
    mean_info = numeric(length(pseudocount_values)),
    sd_info = numeric(length(pseudocount_values))
  )
  
  # Create CV folds
  fold_size <- ceiling(length(sequences) / cv_folds)
  fold_indices <- split(sample(length(sequences)), 
                       rep(1:cv_folds, each = fold_size, length.out = length(sequences)))
  
  for (p_idx in seq_along(pseudocount_values)) {
    pseudocount <- pseudocount_values[p_idx]
    fold_info_scores <- numeric(cv_folds)
    
    cat("Testing pseudocount", pseudocount, "...\n")
    
    for (fold in 1:cv_folds) {
      # Training sequences (all except current fold)
      train_indices <- unlist(fold_indices[-fold])
      train_sequences <- sequences[train_indices]
      
      # Build PWM
      consensus_matrix <- consensusMatrix(train_sequences, as.prob = FALSE)
      nucleotides <- c("A", "C", "G", "T")
      count_matrix <- consensus_matrix[nucleotides, , drop = FALSE]
      
      # Add pseudocount and normalize
      count_matrix <- count_matrix + pseudocount
      pwm <- apply(count_matrix, 2, function(x) x / sum(x))
      
      # Calculate information content
      info_content <- calculate_info_content_fast(pwm)
      fold_info_scores[fold] <- sum(info_content)
    }
    
    results$mean_info[p_idx] <- mean(fold_info_scores)
    results$sd_info[p_idx] <- sd(fold_info_scores)
  }
  
  # Find optimal pseudocount
  best_idx <- which.max(results$mean_info)
  optimal_pseudocount <- results$pseudocount[best_idx]
  
  cat("Pseudocount optimization results:\n")
  print(results)
  cat("Optimal pseudocount:", optimal_pseudocount, "\n")
  
  return(list(
    optimal_pseudocount = optimal_pseudocount,
    results = results
  ))
}

# Build efficient aligned PWM
build_efficient_aligned_pwm <- function(sequences, pseudocount = 0.1, 
                                       batch_size = 10000, optimize_pseudocount = TRUE) {
  
  if (optimize_pseudocount && length(sequences) > 1000) {
    cat("Optimizing pseudocount...\n")
    opt_result <- optimize_pseudocount_cv(sequences)
    pseudocount <- opt_result$optimal_pseudocount
    optimization_results <- opt_result$results
  } else {
    optimization_results <- NULL
  }
  
  cat("Building PWM with pseudocount:", pseudocount, "\n")
  
  # Build consensus matrix efficiently
  if (length(sequences) > batch_size) {
    count_matrix <- build_consensus_batched(sequences, batch_size)
  } else {
    consensus_matrix <- consensusMatrix(sequences, as.prob = FALSE)
    nucleotides <- c("A", "C", "G", "T")
    count_matrix <- consensus_matrix[nucleotides, , drop = FALSE]
  }
  
  cat("Count matrix dimensions:", nrow(count_matrix), "x", ncol(count_matrix), "\n")
  
  # Add pseudocounts and normalize
  count_matrix <- count_matrix + pseudocount
  pwm <- apply(count_matrix, 2, function(x) x / sum(x))
  
  # Calculate information content efficiently
  info_content <- calculate_info_content_fast(pwm)
  total_info <- sum(info_content)
  
  cat("Total information content:", round(total_info, 3), "bits\n")
  cat("Average per position:", round(total_info / ncol(pwm), 3), "bits\n")
  
  # Find conserved positions and core motif
  conserved_positions <- which(info_content > 1.0)
  highly_conserved <- which(info_content > 1.5)
  
  cat("Conserved positions (>1 bit):", length(conserved_positions), "\n")
  cat("Highly conserved positions (>1.5 bit):", length(highly_conserved), "\n")
  
  # Find core motif region
  if (length(conserved_positions) > 0) {
    core_start <- min(conserved_positions)
    core_end <- max(conserved_positions)
    core_length <- core_end - core_start + 1
    cat("Core motif region:", core_start, "-", core_end, "(length:", core_length, ")\n")
  } else {
    core_start <- core_end <- core_length <- NA
  }
  
  # Calculate quality metrics
  max_info_per_pos <- max(info_content)
  median_info_per_pos <- median(info_content)
  
  result <- list(
    pwm = pwm,
    count_matrix = count_matrix,
    info_content = info_content,
    total_info = total_info,
    conserved_positions = conserved_positions,
    highly_conserved_positions = highly_conserved,
    core_motif = list(start = core_start, end = core_end, length = core_length),
    num_sequences = length(sequences),
    pseudocount = pseudocount,
    optimization_results = optimization_results,
    quality_metrics = list(
      max_info = max_info_per_pos,
      median_info = median_info_per_pos,
      avg_info = total_info / ncol(pwm)
    ),
    method = "efficient_aligned",
    creation_time = Sys.time()
  )
  
  return(result)
}

# --- Main Execution ---

# Check input file
if (!file.exists(aligned_file)) {
  stop("Input file not found: ", aligned_file)
}

# Read aligned sequences
cat("Reading aligned sequences from:", aligned_file, "\n")
sequences <- readDNAStringSet(aligned_file)

if (length(sequences) == 0) {
  stop("No sequences found in input file")
}

cat("Loaded", length(sequences), "aligned sequences\n")

# Get memory usage before processing
initial_memory <- gc()
cat("Initial memory usage:", round(sum(initial_memory[,2]), 1), "MB\n")

# Build efficient PWM
start_time <- Sys.time()
pwm_result <- build_efficient_aligned_pwm(sequences, 
                                         batch_size = batch_size,
                                         optimize_pseudocount = optimize_pseudocount)
end_time <- Sys.time()
processing_time <- abs(as.numeric(difftime(end_time, start_time, units = "secs")))

cat("Processing time:", round(processing_time, 2), "seconds\n")

# Get final memory usage
final_memory <- gc()
cat("Final memory usage:", round(sum(final_memory[,2]), 1), "MB\n")

# Save results
output_file <- paste0(output_prefix, ".rds")
saveRDS(pwm_result, output_file)
cat("PWM saved to:", output_file, "\n")

# Save detailed report
report_file <- paste0(output_prefix, "_report.txt")
sink(report_file)
cat("=== Efficient Aligned PWM Report ===\n")
cat("Created:", as.character(Sys.time()), "\n")
cat("Input file:", aligned_file, "\n")
cat("Number of sequences:", pwm_result$num_sequences, "\n")
cat("PWM dimensions:", nrow(pwm_result$pwm), "x", ncol(pwm_result$pwm), "\n")
cat("Pseudocount:", pwm_result$pseudocount, "\n")
cat("Processing time:", round(processing_time, 2), "seconds\n")
cat("\n=== Quality Metrics ===\n")
cat("Total information content:", round(pwm_result$total_info, 3), "bits\n")
cat("Average per position:", round(pwm_result$quality_metrics$avg_info, 3), "bits\n")
cat("Maximum per position:", round(pwm_result$quality_metrics$max_info, 3), "bits\n")
cat("Median per position:", round(pwm_result$quality_metrics$median_info, 3), "bits\n")
cat("Conserved positions (>1 bit):", length(pwm_result$conserved_positions), "\n")
cat("Highly conserved positions (>1.5 bit):", length(pwm_result$highly_conserved_positions), "\n")

if (!is.na(pwm_result$core_motif$start)) {
  cat("Core motif region:", pwm_result$core_motif$start, "-", pwm_result$core_motif$end, 
      "(length:", pwm_result$core_motif$length, ")\n")
}

if (!is.null(pwm_result$optimization_results)) {
  cat("\n=== Pseudocount Optimization ===\n")
  print(pwm_result$optimization_results)
}
sink()

cat("Detailed report saved to:", report_file, "\n")

# Print final summary
cat("\n=== Final PWM Summary ===\n")
cat("Method: Efficient aligned PWM\n")
cat("Sequences processed:", pwm_result$num_sequences, "\n")
cat("PWM dimensions:", nrow(pwm_result$pwm), "x", ncol(pwm_result$pwm), "\n")
cat("Total information content:", round(pwm_result$total_info, 3), "bits\n")
cat("Average per position:", round(pwm_result$quality_metrics$avg_info, 3), "bits\n")
cat("Conserved positions:", length(pwm_result$conserved_positions), "\n")
cat("Optimal pseudocount:", pwm_result$pseudocount, "\n")
cat("Processing time:", round(processing_time, 2), "seconds\n")

cat("\nScript completed successfully!\n")