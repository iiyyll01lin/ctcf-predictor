# null model 生成腳本
# 為PWM評估系統生成各種null model基準
# Author: PWM Improvement Team
# Version: 1.0

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "data/training_sequences.fasta"
output_dir <- if (length(args) >= 2) args[2] else "results/null_models"
n_replicates <- if (length(args) >= 3) as.numeric(args[3]) else 100

cat("Null Model Generation Tool\n")
cat("==========================\n")
cat("Input sequences:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Replicates:", n_replicates, "\n\n")

# --- Create output directory ---
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Helper Functions ---

# Generate random sequences matching base composition
generate_random_sequences <- function(sequences, n_reps = 100) {
  cat("Generating random sequences with matched base composition...\n")
  
  # Calculate overall base composition
  all_seq <- paste(sequences, collapse = "")
  base_freq <- table(strsplit(all_seq, "")[[1]])
  base_freq <- base_freq[c("A", "C", "G", "T")]
  base_freq <- base_freq / sum(base_freq)
  
  cat("Original base frequencies:\n")
  print(round(base_freq, 3))
  
  # Get sequence lengths
  seq_lengths <- nchar(sequences)
  
  random_sequences <- list()
  for (rep in 1:n_reps) {
    rep_sequences <- character(length(sequences))
    for (i in seq_along(sequences)) {
      len <- seq_lengths[i]
      bases <- sample(c("A", "C", "G", "T"), len, replace = TRUE, prob = base_freq)
      rep_sequences[i] <- paste(bases, collapse = "")
    }
    random_sequences[[paste0("random_", rep)]] <- rep_sequences
  }
    return(random_sequences)
}

# Generate GC-content matched sequences
generate_gc_matched_sequences <- function(sequences, n_reps = 100) {
  cat("Generating GC-content matched sequences...\n")
  
  # Calculate GC content for each sequence
  gc_contents <- sapply(sequences, function(seq) {
    bases <- strsplit(seq, "")[[1]]
    gc_count <- sum(bases %in% c("G", "C"))
    return(gc_count / length(bases))
  })
  
  gc_matched_sequences <- list()
  for (rep in 1:n_reps) {
    rep_sequences <- character(length(sequences))
    for (i in seq_along(sequences)) {
      len <- nchar(sequences[i])
      gc_content <- gc_contents[i]
      at_prob <- (1 - gc_content) / 2
      gc_prob <- gc_content / 2
      
      bases <- sample(c("A", "T", "G", "C"), len, replace = TRUE,
                     prob = c(at_prob, at_prob, gc_prob, gc_prob))
      rep_sequences[i] <- paste(bases, collapse = "")
    }
    gc_matched_sequences[[paste0("gc_matched_", rep)]] <- rep_sequences
  }
  
  cat("Generated", n_reps, "GC-content matched sequence sets\n")
  return(gc_matched_sequences)
}

# Generate shuffled sequences preserving individual composition
generate_shuffled_sequences <- function(sequences, n_reps = 100) {
  cat("Generating shuffled sequences preserving individual base composition...\n")
  
  shuffled_sequences <- list()
  for (rep in 1:n_reps) {
    rep_sequences <- character(length(sequences))
    for (i in seq_along(sequences)) {
      bases <- strsplit(sequences[i], "")[[1]]
      shuffled_bases <- sample(bases)
      rep_sequences[i] <- paste(shuffled_bases, collapse = "")
    }
    shuffled_sequences[[paste0("shuffled_", rep)]] <- rep_sequences
  }
  
  return(shuffled_sequences)
}

# Generate position-shuffled sequences (shuffle columns)
generate_position_shuffled <- function(sequences, n_reps = 100) {
  cat("Generating position-shuffled sequences...\n")
  
  # Convert to matrix
  seq_matrix <- do.call(rbind, strsplit(sequences, ""))
  
  position_shuffled <- list()
  for (rep in 1:n_reps) {
    shuffled_matrix <- seq_matrix
    for (col in 1:ncol(seq_matrix)) {
      shuffled_matrix[, col] <- sample(seq_matrix[, col])
    }
    rep_sequences <- apply(shuffled_matrix, 1, paste, collapse = "")
    position_shuffled[[paste0("position_shuffled_", rep)]] <- rep_sequences
  }
    return(position_shuffled)
}

# Generate dinucleotide-preserving sequences
generate_dinucleotide_preserved <- function(sequences, n_reps = 100) {
  cat("Generating dinucleotide-preserving sequences...\n")
  
  # Calculate dinucleotide frequencies
  calculate_dinuc_freq <- function(seq) {
    bases <- strsplit(seq, "")[[1]]
    if (length(bases) < 2) return(NULL)
    
    dinucs <- paste0(bases[-length(bases)], bases[-1])
    freq_table <- table(dinucs)
    return(freq_table / sum(freq_table))
  }
  
  dinuc_preserved <- list()
  for (rep in 1:n_reps) {
    rep_sequences <- character(length(sequences))
    for (i in seq_along(sequences)) {
      seq <- sequences[i]
      len <- nchar(seq)
      
      if (len < 2) {
        rep_sequences[i] <- seq
        next
      }
      
      # Use Markov chain approach to preserve dinucleotide frequencies
      bases <- strsplit(seq, "")[[1]]
      transition_matrix <- matrix(0, nrow = 4, ncol = 4)
      rownames(transition_matrix) <- c("A", "C", "G", "T")
      colnames(transition_matrix) <- c("A", "C", "G", "T")
      
      for (j in 1:(length(bases) - 1)) {
        from_base <- bases[j]
        to_base <- bases[j + 1]
        transition_matrix[from_base, to_base] <- transition_matrix[from_base, to_base] + 1
      }
      
      # Normalize to probabilities
      for (j in 1:4) {
        if (sum(transition_matrix[j, ]) > 0) {
          transition_matrix[j, ] <- transition_matrix[j, ] / sum(transition_matrix[j, ])
        }
      }
      
      # Generate new sequence
      new_seq <- character(len)
      new_seq[1] <- sample(c("A", "C", "G", "T"), 1)
      
      for (j in 2:len) {
        prev_base <- new_seq[j - 1]
        if (sum(transition_matrix[prev_base, ]) > 0) {
          new_seq[j] <- sample(c("A", "C", "G", "T"), 1, 
                              prob = transition_matrix[prev_base, ])
        } else {
          new_seq[j] <- sample(c("A", "C", "G", "T"), 1)
        }
      }
      
      rep_sequences[i] <- paste(new_seq, collapse = "")
    }
    dinuc_preserved[[paste0("dinuc_preserved_", rep)]] <- rep_sequences
  }
  
  cat("Generated", n_reps, "dinucleotide-preserving sequence sets\n")
  return(dinuc_preserved)
}

# Build PWM from sequences
build_simple_pwm <- function(sequences, pseudocount = 0.01) {
  # Convert to matrix
  seq_matrix <- do.call(rbind, strsplit(sequences, ""))
  
  # Count bases at each position
  count_matrix <- matrix(0, nrow = 4, ncol = ncol(seq_matrix))
  rownames(count_matrix) <- c("A", "C", "G", "T")
  
  for (i in 1:ncol(seq_matrix)) {
    counts <- table(factor(seq_matrix[, i], levels = c("A", "C", "G", "T")))
    count_matrix[, i] <- as.numeric(counts)
  }
  
  # Add pseudocounts and normalize
  count_matrix <- count_matrix + pseudocount
  pwm <- apply(count_matrix, 2, function(x) x / sum(x))
  
  return(pwm)
}

# Calculate PWM metrics
calculate_pwm_metrics <- function(pwm) {
  # Information content per position
  info_content <- apply(pwm, 2, function(x) {
    x[x == 0] <- 1e-10
    2 + sum(x * log2(x))
  })
  
  return(list(
    total_info = sum(info_content),
    avg_info = mean(info_content),
    max_info = max(info_content),
    median_info = median(info_content),
    conserved_1bit = sum(info_content > 1.0),
    conserved_1_5bit = sum(info_content > 1.5),
    info_content = info_content
  ))
}

# --- Main Processing ---

# Load sequences
cat("Loading sequences from:", input_file, "\n")
if (file.exists(input_file)) {
  sequences <- readDNAStringSet(input_file)
  sequences <- as.character(sequences)
  cat("Loaded", length(sequences), "sequences\n")
  cat("Average length:", round(mean(nchar(sequences)), 1), "bp\n\n")
} else {
  stop("Input file not found: ", input_file)
}

# Filter sequences for consistency (same length, no N's)
cat("Filtering sequences for null model generation...\n")
valid_sequences <- sequences[
  nchar(sequences) == median(nchar(sequences)) &
  !grepl("N", sequences)
]
cat("Retained", length(valid_sequences), "valid sequences for null models\n\n")

if (length(valid_sequences) < 100) {
  warning("Few valid sequences available, reducing replicates to ", min(n_replicates, 20))
  n_replicates <- min(n_replicates, 20)
}

# Generate different types of null models
null_models <- list()

# 1. Random sequences (matched composition)
null_models$random <- generate_random_sequences(valid_sequences, n_replicates)

# 2. GC-content matched sequences
null_models$gc_matched <- generate_gc_matched_sequences(valid_sequences, n_replicates)

# 3. Shuffled sequences (individual composition preserved)
null_models$shuffled <- generate_shuffled_sequences(valid_sequences, n_replicates)

# 4. Position-shuffled sequences
null_models$position_shuffled <- generate_position_shuffled(valid_sequences, n_replicates)

# 5. Dinucleotide-preserving sequences
null_models$dinuc_preserved <- generate_dinucleotide_preserved(valid_sequences, n_replicates)

# Build PWMs for all null models and calculate metrics
cat("\nBuilding PWMs for null models...\n")
null_pwm_metrics <- list()

for (model_type in names(null_models)) {
  cat("Processing", model_type, "models...\n")
  type_metrics <- list()
  
  for (rep_name in names(null_models[[model_type]])) {
    pwm <- build_simple_pwm(null_models[[model_type]][[rep_name]])
    metrics <- calculate_pwm_metrics(pwm)
    type_metrics[[rep_name]] <- metrics
  }
  
  null_pwm_metrics[[model_type]] <- type_metrics
}

# Summarize null model statistics
cat("\nSummarizing null model statistics...\n")
null_summary <- list()

for (model_type in names(null_pwm_metrics)) {
  metrics_df <- do.call(rbind, lapply(null_pwm_metrics[[model_type]], function(x) {
    data.frame(
      total_info = x$total_info,
      avg_info = x$avg_info,
      max_info = x$max_info,
      median_info = x$median_info,
      conserved_1bit = x$conserved_1bit,
      conserved_1_5bit = x$conserved_1_5bit
    )
  }))
  
  summary_stats <- data.frame(
    metric = c("total_info", "avg_info", "max_info", "median_info", "conserved_1bit", "conserved_1_5bit"),
    mean = apply(metrics_df, 2, mean),
    sd = apply(metrics_df, 2, sd),
    median = apply(metrics_df, 2, median),
    q25 = apply(metrics_df, 2, quantile, 0.25),
    q75 = apply(metrics_df, 2, quantile, 0.75),
    min = apply(metrics_df, 2, min),
    max = apply(metrics_df, 2, max)
  )
  
  null_summary[[model_type]] <- list(
    raw_metrics = metrics_df,
    summary_stats = summary_stats
  )
  
  cat(model_type, "summary:\n")
  cat("  Total info - Mean:", round(summary_stats$mean[1], 3), 
      "± SD:", round(summary_stats$sd[1], 3), "\n")
  cat("  Conserved (>1bit) - Mean:", round(summary_stats$mean[5], 1), 
      "± SD:", round(summary_stats$sd[5], 1), "\n")
}

# Save results
cat("\nSaving null model results...\n")

# Save raw null model sequences
for (model_type in names(null_models)) {
  for (rep_name in names(null_models[[model_type]])) {
    output_file <- file.path(output_dir, paste0(model_type, "_", rep_name, ".fasta"))
    writeXStringSet(DNAStringSet(null_models[[model_type]][[rep_name]]), output_file)
  }
}

# Save PWM metrics
saveRDS(null_pwm_metrics, file.path(output_dir, "null_pwm_metrics.rds"))
saveRDS(null_summary, file.path(output_dir, "null_summary_statistics.rds"))

# Save summary report
report_file <- file.path(output_dir, "null_models_report.txt")
sink(report_file)
cat("=== Null Model Generation Report ===\n")
cat("Generated:", Sys.time(), "\n")
cat("Input file:", input_file, "\n")
cat("Original sequences:", length(sequences), "\n")
cat("Valid sequences used:", length(valid_sequences), "\n")
cat("Replicates per model:", n_replicates, "\n\n")

cat("=== Null Model Statistics ===\n")
for (model_type in names(null_summary)) {
  cat("\n", toupper(model_type), "NULL MODEL:\n")
  cat("Description:", switch(model_type,
    "random" = "Random sequences with matched overall base composition",
    "shuffled" = "Sequences with shuffled order, preserving individual composition",
    "position_shuffled" = "Position-wise shuffled sequences"
  ), "\n")
  
  stats <- null_summary[[model_type]]$summary_stats
  cat("Total Information Content:\n")
  cat("  Mean ± SD:", round(stats$mean[1], 3), "±", round(stats$sd[1], 3), "bits\n")
  cat("  Range:", round(stats$min[1], 3), "-", round(stats$max[1], 3), "bits\n")
  cat("  95% CI:", round(stats$mean[1] - 1.96*stats$sd[1], 3), "-", 
      round(stats$mean[1] + 1.96*stats$sd[1], 3), "bits\n")
  
  cat("Conserved Positions (>1 bit):\n")
  cat("  Mean ± SD:", round(stats$mean[5], 1), "±", round(stats$sd[5], 1), "\n")
  cat("  Range:", stats$min[5], "-", stats$max[5], "\n")
  
  cat("Average Information per Position:\n")
  cat("  Mean ± SD:", round(stats$mean[2], 4), "±", round(stats$sd[2], 4), "bits\n")
}

cat("\n=== Files Generated ===\n")
cat("Null model sequences:", length(null_models) * n_replicates, "FASTA files\n")
cat("PWM metrics:", file.path(output_dir, "null_pwm_metrics.rds"), "\n")
cat("Summary statistics:", file.path(output_dir, "null_summary_statistics.rds"), "\n")
cat("This report:", report_file, "\n")

sink()

cat("Null model generation completed successfully!\n")
cat("Report saved to:", report_file, "\n")
cat("Summary statistics saved to:", file.path(output_dir, "null_summary_statistics.rds"), "\n")
cat("\nGenerated", length(null_models) * n_replicates, "null model sequence sets\n")
