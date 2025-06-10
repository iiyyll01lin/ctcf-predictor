# PWM Comparison and Evaluation Tool
# This script compares PWMs built with different methods and provides comprehensive analysis

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results"
output_file <- if (length(args) >= 2) args[2] else "results/pwm_comparison_report.html"

cat("PWM Comparison and Evaluation Tool\n")
cat("==================================\n")
cat("Results directory:", results_dir, "\n")
cat("Output report:", output_file, "\n\n")

# --- Functions ---

# Load all PWM files
load_all_pwms <- function(results_dir) {
  cat("Loading PWM files from:", results_dir, "\n")
  
  # Find all RDS files
  rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  pwm_files <- rds_files[grepl("pwm", tolower(basename(rds_files)))]
  
  cat("Found", length(pwm_files), "PWM files\n")
  
  pwms <- list()
  
  for (file in pwm_files) {
    filename <- basename(file)
    cat("Loading:", filename, "\n")
    
    tryCatch({
      data <- readRDS(file)
      
      # Handle different data structures
      if (is.list(data) && "pwm" %in% names(data)) {
        # Single PWM result
        pwms[[filename]] <- data
      } else if (is.list(data) && all(sapply(data, function(x) is.list(x) && "pwm" %in% names(x)))) {
        # Multiple PWMs in one file
        for (i in seq_along(data)) {
          name <- paste0(filename, "_", names(data)[i])
          pwms[[name]] <- data[[i]]
        }
      } else if (is.matrix(data)) {
        # Raw PWM matrix
        pwms[[filename]] <- list(
          pwm = data,
          method = "unknown",
          creation_time = file.info(file)$mtime
        )
      }
    }, error = function(e) {
      cat("Error loading", filename, ":", e$message, "\n")
    })
  }
  
  cat("Successfully loaded", length(pwms), "PWMs\n")
  return(pwms)
}

# Calculate comprehensive PWM metrics
calculate_pwm_metrics <- function(pwm_data) {
  if (is.null(pwm_data$pwm)) {
    return(NULL)
  }
  
  pwm <- pwm_data$pwm
  
  # Basic dimensions
  n_positions <- ncol(pwm)
  
  # Information content calculation
  info_content <- apply(pwm, 2, function(x) {
    x[x == 0] <- 1e-10
    2 + sum(x * log2(x))
  })
  
  total_info <- sum(info_content)
  avg_info <- total_info / n_positions
  max_info <- max(info_content)
  min_info <- min(info_content)
  median_info <- median(info_content)
  
  # Conserved positions
  conserved_1bit <- sum(info_content > 1.0)
  conserved_1_5bit <- sum(info_content > 1.5)
  conserved_0_5bit <- sum(info_content > 0.5)
  
  # Core motif detection
  core_positions <- which(info_content > 1.0)
  if (length(core_positions) > 0) {
    core_start <- min(core_positions)
    core_end <- max(core_positions)
    core_length <- core_end - core_start + 1
    core_info <- sum(info_content[core_positions])
  } else {
    core_start <- core_end <- core_length <- core_info <- NA
  }
  
  # Positional bias analysis
  # Check if motif is centered or biased towards ends
  center_pos <- round(n_positions / 2)
  center_window <- max(1, center_pos - 5):min(n_positions, center_pos + 5)
  center_info <- sum(info_content[center_window])
  
  edge_window_left <- 1:min(10, n_positions)
  edge_window_right <- max(1, n_positions - 9):n_positions
  edge_info <- sum(info_content[edge_window_left]) + sum(info_content[edge_window_right])
  
  position_bias <- if (center_info > 0) edge_info / center_info else Inf
  
  # Sequence complexity metrics
  overall_gc <- (sum(pwm["G", ]) + sum(pwm["C", ])) / n_positions
  position_gc <- colSums(pwm[c("G", "C"), ])
  gc_variance <- var(position_gc)
  
  # Most informative nucleotide per position
  max_nucleotides <- apply(pwm, 2, function(x) names(x)[which.max(x)])
  consensus_sequence <- paste(max_nucleotides, collapse = "")
  
  return(list(
    # Basic metrics
    n_positions = n_positions,
    total_info = total_info,
    avg_info = avg_info,
    max_info = max_info,
    min_info = min_info,
    median_info = median_info,
    
    # Conservation metrics
    conserved_0_5bit = conserved_0_5bit,
    conserved_1bit = conserved_1bit,
    conserved_1_5bit = conserved_1_5bit,
    
    # Core motif
    core_start = core_start,
    core_end = core_end,
    core_length = core_length,
    core_info = core_info,
    
    # Position bias
    position_bias = position_bias,
    center_info = center_info,
    edge_info = edge_info,
    
    # Composition
    overall_gc = overall_gc,
    gc_variance = gc_variance,
    consensus_sequence = consensus_sequence,
    
    # Method info
    method = pwm_data$method %||% "unknown",
    num_sequences = pwm_data$num_sequences %||% NA,
    pseudocount = pwm_data$pseudocount %||% NA,
    creation_time = pwm_data$creation_time %||% NA
  ))
}

# Compare PWMs
compare_pwms <- function(pwms) {
  cat("Calculating metrics for", length(pwms), "PWMs...\n")
  
  metrics_list <- list()
  
  for (name in names(pwms)) {
    metrics <- calculate_pwm_metrics(pwms[[name]])
    if (!is.null(metrics)) {
      metrics_list[[name]] <- metrics
    }
  }
  
  # Convert to data frame for easier analysis
  metrics_df <- do.call(rbind, lapply(names(metrics_list), function(name) {
    m <- metrics_list[[name]]
    data.frame(
      name = name,
      method = m$method,
      n_positions = m$n_positions,
      total_info = m$total_info,
      avg_info = m$avg_info,
      max_info = m$max_info,
      median_info = m$median_info,
      conserved_1bit = m$conserved_1bit,
      conserved_1_5bit = m$conserved_1_5bit,
      core_length = m$core_length %||% 0,
      core_info = m$core_info %||% 0,
      position_bias = m$position_bias,
      overall_gc = m$overall_gc,
      num_sequences = m$num_sequences %||% NA,
      pseudocount = m$pseudocount %||% NA,
      stringsAsFactors = FALSE
    )
  }))
  
  return(list(
    metrics_df = metrics_df,
    detailed_metrics = metrics_list
  ))
}

# Generate HTML report
generate_html_report <- function(comparison_results, output_file) {
  cat("Generating HTML report...\n")
  
  metrics_df <- comparison_results$metrics_df
  detailed_metrics <- comparison_results$detailed_metrics
  
  # Create HTML content
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>PWM Comparison Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .metric-high { background-color: #d4edda; }
        .metric-medium { background-color: #fff3cd; }
        .metric-low { background-color: #f8d7da; }
        .summary { background-color: #e9ecef; padding: 15px; margin: 20px 0; border-radius: 5px; }
        .consensus { font-family: monospace; font-size: 14px; background-color: #f8f9fa; padding: 5px; }
    </style>
</head>
<body>
    <h1>PWM Comparison Report</h1>
    <p>Generated on: ", Sys.time(), "</p>
    
    <div class='summary'>
        <h2>Summary</h2>
        <p><strong>Total PWMs analyzed:</strong> ", nrow(metrics_df), "</p>
        <p><strong>Best overall information content:</strong> ", round(max(metrics_df$total_info, na.rm = TRUE), 3), " bits</p>
        <p><strong>Best average per position:</strong> ", round(max(metrics_df$avg_info, na.rm = TRUE), 3), " bits</p>
        <p><strong>Most conserved positions:</strong> ", max(metrics_df$conserved_1bit, na.rm = TRUE), "</p>
    </div>
    
    <h2>Overall Comparison</h2>
    <table>
        <tr>
            <th>PWM Name</th>
            <th>Method</th>
            <th>Positions</th>
            <th>Total Info (bits)</th>
            <th>Avg Info (bits)</th>
            <th>Conserved (>1 bit)</th>
            <th>Core Length</th>
            <th>Sequences</th>
        </tr>")
  
  # Add rows for each PWM
  for (i in 1:nrow(metrics_df)) {
    row <- metrics_df[i, ]
    
    # Color coding based on information content
    total_info_class <- if (row$total_info > quantile(metrics_df$total_info, 0.75, na.rm = TRUE)) "metric-high"
                       else if (row$total_info > quantile(metrics_df$total_info, 0.25, na.rm = TRUE)) "metric-medium"
                       else "metric-low"
    
    html_content <- paste0(html_content, "
        <tr>
            <td>", row$name, "</td>
            <td>", row$method, "</td>
            <td>", row$n_positions, "</td>
            <td class='", total_info_class, "'>", round(row$total_info, 3), "</td>
            <td>", round(row$avg_info, 3), "</td>
            <td>", row$conserved_1bit, "</td>
            <td>", ifelse(is.na(row$core_length), 0, row$core_length), "</td>
            <td>", ifelse(is.na(row$num_sequences), "N/A", row$num_sequences), "</td>
        </tr>")
  }
  
  html_content <- paste0(html_content, "
    </table>
    
    <h2>Detailed Analysis</h2>")
  
  # Top performers section
  best_total <- metrics_df[which.max(metrics_df$total_info), ]
  best_avg <- metrics_df[which.max(metrics_df$avg_info), ]
  most_conserved <- metrics_df[which.max(metrics_df$conserved_1bit), ]
  
  html_content <- paste0(html_content, "
    <h3>Top Performers</h3>
    <ul>
        <li><strong>Best Total Information:</strong> ", best_total$name, " (", round(best_total$total_info, 3), " bits)</li>
        <li><strong>Best Average per Position:</strong> ", best_avg$name, " (", round(best_avg$avg_info, 3), " bits)</li>
        <li><strong>Most Conserved Positions:</strong> ", most_conserved$name, " (", most_conserved$conserved_1bit, " positions)</li>
    </ul>")
  
  # Consensus sequences
  html_content <- paste0(html_content, "
    <h3>Consensus Sequences</h3>
    <table>
        <tr><th>PWM Name</th><th>Consensus Sequence</th></tr>")
  
  for (name in names(detailed_metrics)) {
    consensus <- detailed_metrics[[name]]$consensus_sequence
    if (!is.null(consensus)) {
      html_content <- paste0(html_content, "
        <tr>
            <td>", name, "</td>
            <td class='consensus'>", consensus, "</td>
        </tr>")
    }
  }
  
  html_content <- paste0(html_content, "
    </table>
    
    <h2>Quality Assessment Guidelines</h2>
    <div class='summary'>
        <h3>Information Content Interpretation:</h3>
        <ul>
            <li><strong>Excellent:</strong> >15 bits total, >0.3 bits average per position</li>
            <li><strong>Good:</strong> 10-15 bits total, 0.2-0.3 bits average per position</li>
            <li><strong>Poor:</strong> <10 bits total, <0.2 bits average per position</li>
        </ul>
        
        <h3>Conserved Positions:</h3>
        <ul>
            <li><strong>Strong motif:</strong> >5 positions with >1 bit</li>
            <li><strong>Moderate motif:</strong> 2-5 positions with >1 bit</li>
            <li><strong>Weak motif:</strong> <2 positions with >1 bit</li>
        </ul>
    </div>
    
</body>
</html>")
  
  # Write HTML file
  writeLines(html_content, output_file)
  cat("HTML report saved to:", output_file, "\n")
}

# Generate text summary
generate_text_summary <- function(comparison_results) {
  metrics_df <- comparison_results$metrics_df
  
  cat("\n=== PWM Comparison Summary ===\n")
  cat("Total PWMs analyzed:", nrow(metrics_df), "\n")
  
  # Overall statistics
  cat("\nOverall Statistics:\n")
  cat("Total Information Content:\n")
  cat("  Mean:", round(mean(metrics_df$total_info, na.rm = TRUE), 3), "bits\n")
  cat("  Median:", round(median(metrics_df$total_info, na.rm = TRUE), 3), "bits\n")
  cat("  Range:", round(min(metrics_df$total_info, na.rm = TRUE), 3), "-", 
      round(max(metrics_df$total_info, na.rm = TRUE), 3), "bits\n")
  
  cat("Average Information per Position:\n")
  cat("  Mean:", round(mean(metrics_df$avg_info, na.rm = TRUE), 3), "bits\n")
  cat("  Median:", round(median(metrics_df$avg_info, na.rm = TRUE), 3), "bits\n")
  cat("  Range:", round(min(metrics_df$avg_info, na.rm = TRUE), 3), "-", 
      round(max(metrics_df$avg_info, na.rm = TRUE), 3), "bits\n")
  
  # Top performers
  cat("\nTop Performers:\n")
  best_total <- metrics_df[which.max(metrics_df$total_info), ]
  cat("Best Total Information:", best_total$name, "(", round(best_total$total_info, 3), "bits)\n")
  
  best_avg <- metrics_df[which.max(metrics_df$avg_info), ]
  cat("Best Average per Position:", best_avg$name, "(", round(best_avg$avg_info, 3), "bits)\n")
  
  most_conserved <- metrics_df[which.max(metrics_df$conserved_1bit), ]
  cat("Most Conserved Positions:", most_conserved$name, "(", most_conserved$conserved_1bit, "positions)\n")
  
  # Method comparison
  cat("\nMethod Comparison:\n")
  method_summary <- aggregate(cbind(total_info, avg_info, conserved_1bit) ~ method, 
                             data = metrics_df, FUN = mean, na.rm = TRUE)
  for (i in 1:nrow(method_summary)) {
    method <- method_summary[i, ]
    cat(sprintf("%-20s: %.3f total, %.3f avg, %.1f conserved\n", 
               method$method, method$total_info, method$avg_info, method$conserved_1bit))
  }
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Main Execution ---

# Check if results directory exists
if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

# Load all PWMs
pwms <- load_all_pwms(results_dir)

if (length(pwms) == 0) {
  stop("No PWM files found in results directory")
}

# Compare PWMs
comparison_results <- compare_pwms(pwms)

# Generate text summary
generate_text_summary(comparison_results)

# Generate HTML report
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

generate_html_report(comparison_results, output_file)

# Save detailed comparison data
data_file <- sub("\\.html$", "_data.rds", output_file)
saveRDS(comparison_results, data_file)
cat("Detailed comparison data saved to:", data_file, "\n")

cat("\nComparison completed successfully!\n")
cat("View the HTML report at:", output_file, "\n")
