# 統計顯著性測試腳本
# 比較真實PWM與null model PWM的性能，計算p-value和效應量
# Author: PWM Improvement Team
# Version: 1.0

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results"
null_dir <- if (length(args) >= 2) args[2] else "results/null_models"
output_file <- if (length(args) >= 3) args[3] else "results/statistical_significance_report.html"

cat("Statistical Significance Testing Tool\n")
cat("====================================\n")
cat("Results directory:", results_dir, "\n")
cat("Null models directory:", null_dir, "\n")
cat("Output report:", output_file, "\n\n")

# --- Helper Functions ---

# Load PWM metrics
load_pwm_metrics <- function(pwm_file) {
  if (!file.exists(pwm_file)) {
    return(NULL)
  }
  
  tryCatch({
    data <- readRDS(pwm_file)
    
    # Handle different data structures
    if (is.list(data) && "pwm" %in% names(data)) {
      pwm <- data$pwm
    } else if (is.matrix(data)) {
      pwm <- data
    } else {
      return(NULL)
    }
    
    # Calculate metrics
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
      n_positions = ncol(pwm)
    ))
  }, error = function(e) {
    cat("Error loading", pwm_file, ":", e$message, "\n")
    return(NULL)
  })
}

# Calculate p-value using empirical distribution
calculate_empirical_pvalue <- function(observed_value, null_distribution, alternative = "greater") {
  if (alternative == "greater") {
    p_value <- sum(null_distribution >= observed_value) / length(null_distribution)
  } else if (alternative == "less") {
    p_value <- sum(null_distribution <= observed_value) / length(null_distribution)
  } else { # two-sided
    null_mean <- mean(null_distribution)
    if (observed_value >= null_mean) {
      p_value <- 2 * sum(null_distribution >= observed_value) / length(null_distribution)
    } else {
      p_value <- 2 * sum(null_distribution <= observed_value) / length(null_distribution)
    }
  }
  
  # Ensure p-value is not exactly 0 (for reporting)
  if (p_value == 0) {
    p_value <- 1 / length(null_distribution)
  }
  
  return(p_value)
}

# Calculate effect size (Cohen's d)
calculate_cohens_d <- function(observed_value, null_distribution) {
  null_mean <- mean(null_distribution)
  null_sd <- sd(null_distribution)
  
  if (null_sd == 0) {
    return(Inf)
  }
  
  cohens_d <- (observed_value - null_mean) / null_sd
  return(cohens_d)
}

# Interpret effect size
interpret_effect_size <- function(d) {
  abs_d <- abs(d)
  if (abs_d < 0.2) {
    return("negligible")
  } else if (abs_d < 0.5) {
    return("small")
  } else if (abs_d < 0.8) {
    return("medium")
  } else if (abs_d < 1.2) {
    return("large")
  } else {
    return("very large")
  }
}

# Format p-value for display
format_pvalue <- function(p) {
  if (p < 0.001) {
    return(sprintf("< 0.001"))
  } else if (p < 0.01) {
    return(sprintf("< 0.01"))
  } else {
    return(sprintf("%.3f", p))
  }
}

# --- Main Processing ---

# Load null model statistics
null_stats_file <- file.path(null_dir, "null_summary_statistics.rds")
if (!file.exists(null_stats_file)) {
  stop("Null model statistics not found. Please run generate_null_models.R first.")
}

null_summary <- readRDS(null_stats_file)
cat("Loaded null model statistics for", length(null_summary), "model types\n")

# Load real PWM files
cat("Loading real PWM files...\n")
pwm_files <- list.files(results_dir, pattern = ".*pwm.*\\.rds$", full.names = TRUE)
pwm_files <- pwm_files[!grepl("null", pwm_files)]

cat("Found", length(pwm_files), "real PWM files\n")

# Calculate metrics for real PWMs
real_pwm_metrics <- list()
for (file in pwm_files) {
  filename <- basename(file)
  metrics <- load_pwm_metrics(file)
  if (!is.null(metrics)) {
    real_pwm_metrics[[filename]] <- metrics
  }
}

cat("Successfully loaded", length(real_pwm_metrics), "real PWMs\n\n")

# Perform statistical tests
cat("Performing statistical significance tests...\n")

statistical_results <- list()

for (pwm_name in names(real_pwm_metrics)) {
  cat("Testing", pwm_name, "...\n")
  real_metrics <- real_pwm_metrics[[pwm_name]]
  
  pwm_results <- list()
  
  for (null_type in names(null_summary)) {
    null_data <- null_summary[[null_type]]$raw_metrics
    
    # Test total information content
    total_info_pvalue <- calculate_empirical_pvalue(
      real_metrics$total_info, 
      null_data$total_info, 
      "greater"
    )
    total_info_effect <- calculate_cohens_d(
      real_metrics$total_info, 
      null_data$total_info
    )
    
    # Test conserved positions
    conserved_pvalue <- calculate_empirical_pvalue(
      real_metrics$conserved_1bit, 
      null_data$conserved_1bit, 
      "greater"
    )
    conserved_effect <- calculate_cohens_d(
      real_metrics$conserved_1bit, 
      null_data$conserved_1bit
    )
    
    # Test average information
    avg_info_pvalue <- calculate_empirical_pvalue(
      real_metrics$avg_info, 
      null_data$avg_info, 
      "greater"
    )
    avg_info_effect <- calculate_cohens_d(
      real_metrics$avg_info, 
      null_data$avg_info
    )
    
    pwm_results[[null_type]] <- list(
      total_info = list(
        observed = real_metrics$total_info,
        null_mean = mean(null_data$total_info),
        null_sd = sd(null_data$total_info),
        p_value = total_info_pvalue,
        effect_size = total_info_effect,
        interpretation = interpret_effect_size(total_info_effect)
      ),
      conserved_positions = list(
        observed = real_metrics$conserved_1bit,
        null_mean = mean(null_data$conserved_1bit),
        null_sd = sd(null_data$conserved_1bit),
        p_value = conserved_pvalue,
        effect_size = conserved_effect,
        interpretation = interpret_effect_size(conserved_effect)
      ),
      avg_info = list(
        observed = real_metrics$avg_info,
        null_mean = mean(null_data$avg_info),
        null_sd = sd(null_data$avg_info),
        p_value = avg_info_pvalue,
        effect_size = avg_info_effect,
        interpretation = interpret_effect_size(avg_info_effect)
      )
    )
  }
  
  statistical_results[[pwm_name]] <- pwm_results
}

# Generate HTML report
cat("Generating HTML statistical report...\n")

html_content <- c(
  "<!DOCTYPE html>",
  "<html>",
  "<head>",
  "    <title>PWM Statistical Significance Report</title>",
  "    <style>",
  "        body { font-family: Arial, sans-serif; margin: 40px; }",
  "        h1 { color: #2c3e50; }",
  "        h2 { color: #34495e; border-bottom: 2px solid #ecf0f1; }",
  "        h3 { color: #7f8c8d; }",
  "        table { border-collapse: collapse; width: 100%; margin: 20px 0; }",
  "        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
  "        th { background-color: #f2f2f2; font-weight: bold; }",
  "        .significant { background-color: #d5f4e6; }",
  "        .not-significant { background-color: #fadbd8; }",
  "        .effect-large { color: #e74c3c; font-weight: bold; }",
  "        .effect-medium { color: #f39c12; }",
  "        .effect-small { color: #95a5a6; }",
  "        .summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin: 20px 0; }",
  "        .metric-box { display: inline-block; margin: 10px; padding: 10px; border: 1px solid #ddd; border-radius: 5px; }",
  "    </style>",
  "</head>",
  "<body>",
  "    <h1>PWM Statistical Significance Report</h1>",
  paste0("    <p>Generated on: ", Sys.time(), "</p>"),
  "",
  "    <div class='summary'>",
  "        <h2>Summary</h2>",
  paste0("        <p><strong>PWMs analyzed:</strong> ", length(statistical_results), "</p>"),
  paste0("        <p><strong>Null model types:</strong> ", length(null_summary), "</p>"),
  "        <p><strong>Metrics tested:</strong> Total information content, Conserved positions, Average information content</p>",
  "        <p><strong>Significance threshold:</strong> p < 0.05</p>",
  "    </div>"
)

# Add overview table
html_content <- c(html_content,
  "    <h2>Statistical Significance Overview</h2>",
  "    <table>",
  "        <tr>",
  "            <th>PWM Name</th>",
  "            <th>Null Model</th>",
  "            <th>Total Info</th>",
  "            <th>Conserved Pos</th>",
  "            <th>Avg Info</th>",
  "        </tr>"
)

for (pwm_name in names(statistical_results)) {
  for (null_type in names(statistical_results[[pwm_name]])) {
    results <- statistical_results[[pwm_name]][[null_type]]
    
    total_info_sig <- if (results$total_info$p_value < 0.05) "significant" else "not-significant"
    conserved_sig <- if (results$conserved_positions$p_value < 0.05) "significant" else "not-significant"
    avg_info_sig <- if (results$avg_info$p_value < 0.05) "significant" else "not-significant"
    
    html_content <- c(html_content,
      "        <tr>",
      paste0("            <td>", pwm_name, "</td>"),
      paste0("            <td>", null_type, "</td>"),
      paste0("            <td class='", total_info_sig, "'>p = ", format_pvalue(results$total_info$p_value), "</td>"),
      paste0("            <td class='", conserved_sig, "'>p = ", format_pvalue(results$conserved_positions$p_value), "</td>"),
      paste0("            <td class='", avg_info_sig, "'>p = ", format_pvalue(results$avg_info$p_value), "</td>"),
      "        </tr>"
    )
  }
}

html_content <- c(html_content, "    </table>")

# Add detailed results for each PWM
for (pwm_name in names(statistical_results)) {
  html_content <- c(html_content,
    paste0("    <h2>", pwm_name, "</h2>"),
    "    <table>",
    "        <tr>",
    "            <th>Metric</th>",
    "            <th>Null Model</th>",
    "            <th>Observed</th>",
    "            <th>Null Mean ± SD</th>",
    "            <th>p-value</th>",
    "            <th>Effect Size (d)</th>",
    "            <th>Interpretation</th>",
    "        </tr>"
  )
  
  for (null_type in names(statistical_results[[pwm_name]])) {
    results <- statistical_results[[pwm_name]][[null_type]]
    
    for (metric in c("total_info", "conserved_positions", "avg_info")) {
      metric_data <- results[[metric]]
      
      effect_class <- switch(metric_data$interpretation,
        "negligible" = "effect-small",
        "small" = "effect-small",
        "medium" = "effect-medium",
        "large" = "effect-large",
        "very large" = "effect-large",
        "effect-small"
      )
      
      sig_class <- if (metric_data$p_value < 0.05) "significant" else "not-significant"
      
      metric_label <- switch(metric,
        "total_info" = "Total Information",
        "conserved_positions" = "Conserved Positions",
        "avg_info" = "Average Information"
      )
      
      html_content <- c(html_content,
        "        <tr>",
        paste0("            <td>", metric_label, "</td>"),
        paste0("            <td>", null_type, "</td>"),
        paste0("            <td>", round(metric_data$observed, 3), "</td>"),
        paste0("            <td>", round(metric_data$null_mean, 3), " ± ", round(metric_data$null_sd, 3), "</td>"),
        paste0("            <td class='", sig_class, "'>", format_pvalue(metric_data$p_value), "</td>"),
        paste0("            <td class='", effect_class, "'>", round(metric_data$effect_size, 2), "</td>"),
        paste0("            <td>", metric_data$interpretation, "</td>"),
        "        </tr>"
      )
    }
  }
  
  html_content <- c(html_content, "    </table>")
}

# Add methodology section
html_content <- c(html_content,
  "    <h2>Methodology</h2>",
  "    <h3>Null Models</h3>",
  "    <ul>",
  "        <li><strong>Random:</strong> Random sequences with matched overall base composition</li>",
  "        <li><strong>Shuffled:</strong> Sequences with shuffled order, preserving individual composition</li>",
  "        <li><strong>Position Shuffled:</strong> Position-wise shuffled sequences</li>",
  "    </ul>",
  "    <h3>Statistical Tests</h3>",
  "    <ul>",
  "        <li><strong>p-value:</strong> Empirical p-value based on null distribution</li>",
  "        <li><strong>Effect Size:</strong> Cohen's d = (observed - null_mean) / null_sd</li>",
  "        <li><strong>Significance threshold:</strong> p < 0.05</li>",
  "    </ul>",
  "    <h3>Effect Size Interpretation</h3>",
  "    <ul>",
  "        <li><strong>Negligible:</strong> |d| < 0.2</li>",
  "        <li><strong>Small:</strong> 0.2 ≤ |d| < 0.5</li>",
  "        <li><strong>Medium:</strong> 0.5 ≤ |d| < 0.8</li>",
  "        <li><strong>Large:</strong> 0.8 ≤ |d| < 1.2</li>",
  "        <li><strong>Very Large:</strong> |d| ≥ 1.2</li>",
  "    </ul>",
  "</body>",
  "</html>"
)

# Write HTML file
writeLines(html_content, output_file)

# Save statistical results as RDS
saveRDS(statistical_results, gsub("\\.html$", "_data.rds", output_file))

cat("Statistical significance testing completed!\n")
cat("HTML report saved to:", output_file, "\n")
cat("Raw results saved to:", gsub("\\.html$", "_data.rds", output_file), "\n")

# Print summary to console
cat("\n=== SUMMARY OF RESULTS ===\n")
for (pwm_name in names(statistical_results)) {
  cat("\n", pwm_name, ":\n")
  for (null_type in names(statistical_results[[pwm_name]])) {
    results <- statistical_results[[pwm_name]][[null_type]]
    
    cat("  vs", null_type, "null model:\n")
    cat("    Total info: p =", format_pvalue(results$total_info$p_value), 
        ", d =", round(results$total_info$effect_size, 2), 
        "(", results$total_info$interpretation, ")\n")
    cat("    Conserved pos: p =", format_pvalue(results$conserved_positions$p_value),
        ", d =", round(results$conserved_positions$effect_size, 2),
        "(", results$conserved_positions$interpretation, ")\n")
  }
}
