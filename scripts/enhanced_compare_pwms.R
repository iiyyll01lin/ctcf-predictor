# 增強的PWM比較工具 - 整合null model統計分析
# 這個腳本擴展了原有的PWM比較功能，加入null model基準和統計顯著性測試
# Author: PWM Improvement Team  
# Version: 2.0 (Null Model Integration)

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results"
output_file <- if (length(args) >= 2) args[2] else "results/enhanced_pwm_comparison_report.html"
null_dir <- if (length(args) >= 3) args[3] else "results/null_models"

cat("Enhanced PWM Comparison Tool with Null Model Integration\n")
cat("========================================================\n")
cat("Results directory:", results_dir, "\n")
cat("Output report:", output_file, "\n")
cat("Null models directory:", null_dir, "\n\n")

# --- Load original functions ---
# Simple and robust approach: try common locations for compare_pwms.R
compare_pwms_paths <- c(
  "scripts/compare_pwms.R",         # From project root
  "./compare_pwms.R",               # Same directory as this script
  "../compare_pwms.R",              # Parent directory
  file.path(dirname(getwd()), "scripts", "compare_pwms.R")  # Alternative path
)

compare_pwms_file <- NULL
for (path in compare_pwms_paths) {
  if (file.exists(path)) {
    compare_pwms_file <- path
    break
  }
}

if (is.null(compare_pwms_file)) {
  stop("Could not find compare_pwms.R. Please ensure it exists in the scripts directory.")
}

source(compare_pwms_file)
cat("Successfully loaded compare_pwms.R from:", compare_pwms_file, "\n\n")

# --- Additional Helper Functions for Null Model Integration ---

# Load null model summary statistics
load_null_statistics <- function(null_dir) {
  null_stats_file <- file.path(null_dir, "null_summary_statistics.rds")
  
  if (!file.exists(null_stats_file)) {
    cat("Warning: Null model statistics not found. Generating simplified report without statistical testing.\n")
    return(NULL)
  }
  
  tryCatch({
    null_summary <- readRDS(null_stats_file)
    cat("Loaded null model statistics for", length(null_summary), "model types\n")
    return(null_summary)
  }, error = function(e) {
    cat("Error loading null model statistics:", e$message, "\n")
    return(NULL)
  })
}

# Calculate statistical significance for a single PWM
calculate_pwm_significance <- function(pwm_metrics, null_summary) {
  if (is.null(null_summary)) {
    return(NULL)
  }
  
  significance_results <- list()
  
  for (null_type in names(null_summary)) {
    null_data <- null_summary[[null_type]]$raw_metrics
    
    # Calculate p-values and effect sizes
    total_info_pvalue <- calculate_empirical_pvalue(
      pwm_metrics$total_info, null_data$total_info, "greater"
    )
    total_info_effect <- calculate_cohens_d(
      pwm_metrics$total_info, null_data$total_info
    )
    
    conserved_pvalue <- calculate_empirical_pvalue(
      pwm_metrics$conserved_1bit, null_data$conserved_1bit, "greater"
    )
    conserved_effect <- calculate_cohens_d(
      pwm_metrics$conserved_1bit, null_data$conserved_1bit
    )
    
    avg_info_pvalue <- calculate_empirical_pvalue(
      pwm_metrics$avg_info, null_data$avg_info, "greater"
    )
    avg_info_effect <- calculate_cohens_d(
      pwm_metrics$avg_info, null_data$avg_info
    )
    
    significance_results[[null_type]] <- list(
      total_info = list(
        p_value = total_info_pvalue,
        effect_size = total_info_effect,
        is_significant = total_info_pvalue < 0.05
      ),
      conserved_positions = list(
        p_value = conserved_pvalue,
        effect_size = conserved_effect,
        is_significant = conserved_pvalue < 0.05
      ),
      avg_info = list(
        p_value = avg_info_pvalue,
        effect_size = avg_info_effect,
        is_significant = avg_info_pvalue < 0.05
      )
    )
  }
  
  return(significance_results)
}

# Helper functions from statistical_significance_test.R
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
  
  if (p_value == 0) {
    p_value <- 1 / length(null_distribution)
  }
  
  return(p_value)
}

calculate_cohens_d <- function(observed_value, null_distribution) {
  null_mean <- mean(null_distribution)
  null_sd <- sd(null_distribution)
  
  if (null_sd == 0) {
    return(Inf)
  }
  
  cohens_d <- (observed_value - null_mean) / null_sd
  return(cohens_d)
}

format_pvalue <- function(p) {
  if (p < 0.001) {
    return("< 0.001")
  } else if (p < 0.01) {
    return("< 0.01")
  } else {
    return(sprintf("%.3f", p))
  }
}

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

# Generate enhanced HTML report with null model analysis
generate_enhanced_html_report <- function(comparison_results, null_summary, significance_results, output_file) {
  cat("Generating enhanced HTML report with null model analysis...\n")
  
  metrics_df <- comparison_results$metrics_df
  detailed_metrics <- comparison_results$detailed_metrics
  
  # Start HTML
  html_content <- paste0("
<!DOCTYPE html>
<html>
<head>
    <title>Enhanced PWM Comparison Report with Null Model Analysis</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; border-bottom: 2px solid #ecf0f1; }
        h3 { color: #7f8c8d; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; font-weight: bold; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .metric-high { background-color: #d5f4e6; }
        .metric-medium { background-color: #fff3cd; }
        .metric-low { background-color: #fadbd8; }
        .significant { background-color: #d5f4e6; font-weight: bold; }
        .not-significant { background-color: #fadbd8; }
        .effect-large { color: #e74c3c; font-weight: bold; }
        .effect-medium { color: #f39c12; }
        .effect-small { color: #95a5a6; }
        .summary { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin: 20px 0; }
        .null-summary { background-color: #e9ecef; padding: 15px; border-radius: 5px; margin: 20px 0; }
        .consensus { font-family: monospace; font-size: 14px; background-color: #f8f9fa; padding: 5px; }
        .warning { background-color: #fff3cd; padding: 10px; border-left: 4px solid #ffc107; margin: 10px 0; }
    </style>
</head>
<body>
    <h1>Enhanced PWM Comparison Report with Null Model Analysis</h1>
    <p>Generated on: ", Sys.time(), "</p>")
  
  # Summary section
  html_content <- paste0(html_content, "
    <div class='summary'>
        <h2>Summary</h2>
        <p><strong>PWMs analyzed:</strong> ", nrow(metrics_df), "</p>
        <p><strong>Best overall information content:</strong> ", round(max(metrics_df$total_info, na.rm = TRUE), 3), " bits</p>
        <p><strong>Best average per position:</strong> ", round(max(metrics_df$avg_info, na.rm = TRUE), 3), " bits</p>
        <p><strong>Most conserved positions:</strong> ", max(metrics_df$conserved_1bit, na.rm = TRUE), "</p>")
  
  if (!is.null(null_summary)) {
    html_content <- paste0(html_content, "
        <p><strong>Null model types:</strong> ", length(null_summary), "</p>
        <p><strong>Statistical significance threshold:</strong> p < 0.05</p>")
  } else {
    html_content <- paste0(html_content, "
        <div class='warning'>
            <strong>Note:</strong> Null model statistics not available. Statistical significance tests were not performed.
        </div>")
  }
  
  html_content <- paste0(html_content, "
    </div>")
  
  # Null model baselines (if available)
  if (!is.null(null_summary)) {
    html_content <- paste0(html_content, "
    <div class='null-summary'>
        <h2>Null Model Baselines</h2>
        <table>
            <tr>
                <th>Null Model Type</th>
                <th>Description</th>
                <th>Total Info (Mean ± SD)</th>
                <th>Conserved Pos (Mean ± SD)</th>
                <th>Avg Info (Mean ± SD)</th>
            </tr>")
    
    for (null_type in names(null_summary)) {
      stats <- null_summary[[null_type]]$summary_stats
      description <- switch(null_type,
        "random" = "Random sequences with matched base composition",
        "shuffled" = "Shuffled sequences preserving individual composition",
        "position_shuffled" = "Position-wise shuffled sequences",
        "Unknown null model type"
      )
      
      html_content <- paste0(html_content, "
            <tr>
                <td><strong>", null_type, "</strong></td>
                <td>", description, "</td>
                <td>", round(stats$mean[1], 3), " ± ", round(stats$sd[1], 3), "</td>
                <td>", round(stats$mean[5], 1), " ± ", round(stats$sd[5], 1), "</td>
                <td>", round(stats$mean[2], 4), " ± ", round(stats$sd[2], 4), "</td>
            </tr>")
    }
    
    html_content <- paste0(html_content, "
        </table>
    </div>")
  }
  
  # Main comparison table with significance
  html_content <- paste0(html_content, "
    <h2>PWM Performance vs Null Models</h2>
    <table>
        <tr>
            <th rowspan='2'>PWM Name</th>
            <th rowspan='2'>Method</th>
            <th rowspan='2'>Positions</th>
            <th rowspan='2'>Total Info (bits)</th>
            <th rowspan='2'>Conserved (>1 bit)</th>
            <th rowspan='2'>Avg Info (bits)</th>")
  
  if (!is.null(null_summary)) {
    html_content <- paste0(html_content, "
            <th colspan='", length(null_summary), "'>Statistical Significance (vs Null Models)</th>")
  }
  
  html_content <- paste0(html_content, "
        </tr>
        <tr>")
  
  if (!is.null(null_summary)) {
    for (null_type in names(null_summary)) {
      html_content <- paste0(html_content, "
            <th>vs ", null_type, "</th>")
    }
  }
  
  html_content <- paste0(html_content, "
        </tr>")
    # Add data rows
  for (i in seq_len(nrow(metrics_df))) {
    row <- metrics_df[i, ]
    pwm_name <- row$name
    
    # Color coding
    total_info_class <- if (row$total_info > quantile(metrics_df$total_info, 0.75, na.rm = TRUE)) "metric-high"
                       else if (row$total_info > quantile(metrics_df$total_info, 0.25, na.rm = TRUE)) "metric-medium"
                       else "metric-low"
    
    html_content <- paste0(html_content, "
        <tr>
            <td><strong>", pwm_name, "</strong></td>
            <td>", row$method, "</td>
            <td>", row$n_positions, "</td>
            <td class='", total_info_class, "'>", round(row$total_info, 3), "</td>
            <td>", row$conserved_1bit, "</td>
            <td>", round(row$avg_info, 3), "</td>")
    
    # Add significance columns
    if (!is.null(null_summary) && !is.null(significance_results[[pwm_name]])) {
      for (null_type in names(null_summary)) {
        sig_data <- significance_results[[pwm_name]][[null_type]]
        
        if (!is.null(sig_data)) {
          # Determine overall significance (any metric significant)
          is_significant <- any(
            sig_data$total_info$is_significant,
            sig_data$conserved_positions$is_significant,
            sig_data$avg_info$is_significant
          )
          
          sig_class <- if (is_significant) "significant" else "not-significant"
          
          # Create significance summary
          sig_summary <- paste0(
            "Total: p=", format_pvalue(sig_data$total_info$p_value), "<br>",
            "Conserved: p=", format_pvalue(sig_data$conserved_positions$p_value), "<br>",
            "Avg: p=", format_pvalue(sig_data$avg_info$p_value)
          )
          
          html_content <- paste0(html_content, "
            <td class='", sig_class, "'>", sig_summary, "</td>")
        } else {
          html_content <- paste0(html_content, "
            <td>N/A</td>")
        }
      }
    }
    
    html_content <- paste0(html_content, "
        </tr>")
  }
  
  html_content <- paste0(html_content, "
    </table>")
  
  # Top performers with significance analysis
  best_total <- metrics_df[which.max(metrics_df$total_info), ]
  best_avg <- metrics_df[which.max(metrics_df$avg_info), ]
  most_conserved <- metrics_df[which.max(metrics_df$conserved_1bit), ]
  
  html_content <- paste0(html_content, "
    <h2>Top Performers with Statistical Analysis</h2>
    <div class='summary'>
        <h3>Performance Leaders</h3>
        <ul>
            <li><strong>Best Total Information:</strong> ", best_total$name, " (", round(best_total$total_info, 3), " bits)</li>
            <li><strong>Best Average per Position:</strong> ", best_avg$name, " (", round(best_avg$avg_info, 3), " bits)</li>
            <li><strong>Most Conserved Positions:</strong> ", most_conserved$name, " (", most_conserved$conserved_1bit, " positions)</li>
        </ul>")
  
  if (!is.null(null_summary) && !is.null(significance_results)) {
    html_content <- paste0(html_content, "
        <h3>Statistical Significance Summary</h3>
        <p>PWMs showing significant improvement over null models:</p>
        <ul>")
    
    for (pwm_name in names(significance_results)) {
      pwm_sig <- significance_results[[pwm_name]]
      
      if (!is.null(pwm_sig)) {
        significant_metrics <- c()
        
        for (null_type in names(pwm_sig)) {
          sig_data <- pwm_sig[[null_type]]
          if (any(sig_data$total_info$is_significant, 
                  sig_data$conserved_positions$is_significant,
                  sig_data$avg_info$is_significant)) {
            significant_metrics <- c(significant_metrics, null_type)
          }
        }
        
        if (length(significant_metrics) > 0) {
          html_content <- paste0(html_content, "
            <li><strong>", pwm_name, ":</strong> Significant vs ", paste(significant_metrics, collapse = ", "), "</li>")
        }
      }
    }
    
    html_content <- paste0(html_content, "
        </ul>")
  }
  
  html_content <- paste0(html_content, "
    </div>")
  
  # Add consensus sequences table
  html_content <- paste0(html_content, "
    <h2>Consensus Sequences</h2>
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
    </table>")
  
  # Quality assessment guidelines
  html_content <- paste0(html_content, "
    <h2>Quality Assessment Guidelines</h2>
    <div class='summary'>
        <h3>Information Content Standards for CTCF PWMs:</h3>
        <ul>
            <li><strong>Excellent:</strong> >15 bits total, >0.05 bits average per position, 2-5 conserved positions</li>
            <li><strong>Good:</strong> 10-15 bits total, 0.03-0.05 bits average per position, 1-2 conserved positions</li>
            <li><strong>Fair:</strong> 5-10 bits total, 0.02-0.03 bits average per position, 0-1 conserved positions</li>
            <li><strong>Poor:</strong> <5 bits total, <0.02 bits average per position, 0 conserved positions</li>
        </ul>
        
        <h3>Statistical Significance Interpretation:</h3>
        <ul>
            <li><strong>p < 0.001:</strong> Highly significant improvement over null model</li>
            <li><strong>p < 0.01:</strong> Significant improvement over null model</li>
            <li><strong>p < 0.05:</strong> Marginally significant improvement over null model</li>
            <li><strong>p ≥ 0.05:</strong> No significant improvement over null model</li>
        </ul>
        
        <h3>Effect Size Interpretation (Cohen's d):</h3>
        <ul>
            <li><strong>|d| ≥ 1.2:</strong> Very large effect</li>
            <li><strong>0.8 ≤ |d| < 1.2:</strong> Large effect</li>
            <li><strong>0.5 ≤ |d| < 0.8:</strong> Medium effect</li>
            <li><strong>0.2 ≤ |d| < 0.5:</strong> Small effect</li>
            <li><strong>|d| < 0.2:</strong> Negligible effect</li>
        </ul>
    </div>")
  
  # Methodology section
  html_content <- paste0(html_content, "
    <h2>Methodology</h2>
    <div class='summary'>
        <h3>PWM Quality Metrics:</h3>
        <ul>
            <li><strong>Total Information Content:</strong> Sum of information content across all positions</li>
            <li><strong>Average Information Content:</strong> Mean information content per position</li>
            <li><strong>Conserved Positions:</strong> Number of positions with >1 bit information content</li>
        </ul>
        
        <h3>Null Models:</h3>
        <ul>
            <li><strong>Random:</strong> Random sequences matching overall base composition</li>
            <li><strong>Shuffled:</strong> Original sequences with shuffled base order</li>
            <li><strong>Position Shuffled:</strong> Sequences with position-wise shuffled bases</li>
        </ul>
        
        <h3>Statistical Testing:</h3>
        <ul>
            <li><strong>P-values:</strong> Empirical p-values based on null distribution comparisons</li>
            <li><strong>Effect sizes:</strong> Cohen's d calculated as (observed - null_mean) / null_sd</li>
            <li><strong>Significance threshold:</strong> p < 0.05</li>
        </ul>
    </div>
    
</body>
</html>")
  
  # Write HTML file
  writeLines(html_content, output_file)
  cat("Enhanced HTML report saved to:", output_file, "\n")
}

# --- Main Execution ---

# Check directories
if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

# Load null model statistics
null_summary <- load_null_statistics(null_dir)

# Load and compare PWMs (using original functions)
pwms <- load_all_pwms(results_dir)

if (length(pwms) == 0) {
  stop("No PWM files found in results directory")
}

comparison_results <- compare_pwms(pwms)

# Calculate statistical significance for each PWM
significance_results <- list()
if (!is.null(null_summary)) {
  cat("Calculating statistical significance...\n")
  
  for (pwm_name in names(comparison_results$detailed_metrics)) {
    pwm_metrics <- comparison_results$detailed_metrics[[pwm_name]]
    significance_results[[pwm_name]] <- calculate_pwm_significance(pwm_metrics, null_summary)
  }
  
  cat("Statistical significance testing completed\n")
}

# Generate enhanced reports
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Generate enhanced HTML report
generate_enhanced_html_report(comparison_results, null_summary, significance_results, output_file)

# Save enhanced comparison data
data_file <- sub("\\.html$", "_data.rds", output_file)
enhanced_results <- list(
  comparison_results = comparison_results,
  null_summary = null_summary,
  significance_results = significance_results,
  generation_time = Sys.time()
)
saveRDS(enhanced_results, data_file)

# Generate text summary with significance
cat("\n=== Enhanced PWM Comparison Summary ===\n")
generate_text_summary(comparison_results)

if (!is.null(null_summary)) {
  cat("\n=== Statistical Significance Summary ===\n")
  cat("Null model types:", length(null_summary), "\n")
  
  # Count significant PWMs
  significant_pwms <- 0
  for (pwm_name in names(significance_results)) {
    pwm_sig <- significance_results[[pwm_name]]
    if (!is.null(pwm_sig)) {
      for (null_type in names(pwm_sig)) {
        sig_data <- pwm_sig[[null_type]]
        if (any(sig_data$total_info$is_significant,
                sig_data$conserved_positions$is_significant,
                sig_data$avg_info$is_significant)) {
          significant_pwms <- significant_pwms + 1
          break
        }
      }
    }
  }
  
  cat("PWMs with significant improvement:", significant_pwms, "out of", length(significance_results), "\n")
  
  # Report top performers with significance
  for (pwm_name in names(significance_results)) {
    pwm_sig <- significance_results[[pwm_name]]
    if (!is.null(pwm_sig)) {
      cat("\n", pwm_name, ":\n")
      for (null_type in names(pwm_sig)) {
        sig_data <- pwm_sig[[null_type]]
        cat("  vs", null_type, ":\n")
        cat("    Total info: p =", format_pvalue(sig_data$total_info$p_value), 
            ifelse(sig_data$total_info$is_significant, " *", ""), "\n")
        cat("    Conserved: p =", format_pvalue(sig_data$conserved_positions$p_value),
            ifelse(sig_data$conserved_positions$is_significant, " *", ""), "\n")
      }
    }
  }
}

cat("\nEnhanced comparison completed successfully!\n")
cat("View the enhanced HTML report at:", output_file, "\n")
cat("Enhanced data saved to:", data_file, "\n")
