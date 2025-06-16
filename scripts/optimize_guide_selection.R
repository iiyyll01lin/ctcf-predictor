#!/usr/bin/env Rscript

# CRISPR Guide Selection Optimization Script
# Optimizes guide RNA selection based on multiple criteria including off-target risk
# Author: CTCF Predictor Pipeline
# Usage: Rscript optimize_guide_selection.R --guides <guides_file> --off_targets <off_targets_file> --output <output_dir> [options]

suppressMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(corrplot)
  library(pheatmap)
  library(RColorBrewer)
})

# Command line options
option_list <- list(
  make_option(c("-g", "--guides"), type="character", default=NULL,
              help="Input guide RNAs file (TSV format)", metavar="character"),
  make_option(c("-t", "--off_targets"), type="character", default=NULL,
              help="Off-target predictions file (TSV format)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="optimized_guides",
              help="Output directory [default=%default]", metavar="character"),
  make_option(c("-r", "--risk_assessment"), type="character", default=NULL,
              help="Risk assessment file (optional)", metavar="character"),
  make_option(c("--max_guides"), type="numeric", default=3,
              help="Maximum guides to select per CTCF site [default=%default]", metavar="numeric"),
  make_option(c("--min_score"), type="numeric", default=0.6,
              help="Minimum guide score threshold [default=%default]", metavar="numeric"),
  make_option(c("--max_off_targets"), type="numeric", default=5,
              help="Maximum allowed off-targets per guide [default=%default]", metavar="numeric"),
  make_option(c("--max_risk"), type="character", default="medium",
              help="Maximum allowed risk level: low, medium, high [default=%default]", metavar="character"),
  make_option(c("--optimization_method"), type="character", default="balanced",
              help="Optimization method: balanced, specificity, activity [default=%default]", metavar="character"),
  make_option(c("--weight_activity"), type="numeric", default=0.4,
              help="Weight for guide activity score [default=%default]", metavar="numeric"),
  make_option(c("--weight_specificity"), type="numeric", default=0.4,
              help="Weight for specificity score [default=%default]", metavar="numeric"),
  make_option(c("--weight_coverage"), type="numeric", default=0.2,
              help="Weight for target coverage [default=%default]", metavar="numeric"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$guides)) {
  print_help(opt_parser)
  stop("Input guides file must be specified.", call.=FALSE)
}

if (!file.exists(opt$guides)) {
  stop("Guides file does not exist: ", opt$guides, call.=FALSE)
}

# Create output directory
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

if (opt$verbose) {
  cat("Starting guide selection optimization...\n")
  cat("Guides file:", opt$guides, "\n")
  cat("Off-targets file:", opt$off_targets, "\n")
  cat("Output directory:", opt$output, "\n")
  cat("Optimization method:", opt$optimization_method, "\n")
}

# Function to load guide data
load_guide_data <- function(guides_file) {
  if (opt$verbose) cat("Loading guide data...\n")
  
  guides <- read.table(guides_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Validate required columns
  required_cols <- c("guide_id", "guide_sequence", "score")
  missing_cols <- setdiff(required_cols, colnames(guides))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in guides file: ", paste(missing_cols, collapse = ", "))
  }
  
  if (opt$verbose) cat("Loaded", nrow(guides), "guide RNAs\n")
  
  return(guides)
}

# Function to load off-target data
load_off_target_data <- function(off_targets_file) {
  if (is.null(off_targets_file) || !file.exists(off_targets_file)) {
    if (opt$verbose) cat("No off-target data provided\n")
    return(data.frame())
  }
  
  if (opt$verbose) cat("Loading off-target data...\n")
  
  off_targets <- read.table(off_targets_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  if (opt$verbose) cat("Loaded", nrow(off_targets), "off-target predictions\n")
  
  return(off_targets)
}

# Function to load risk assessment data
load_risk_data <- function(risk_file) {
  if (is.null(risk_file) || !file.exists(risk_file)) {
    if (opt$verbose) cat("No risk assessment data provided\n")
    return(data.frame())
  }
  
  if (opt$verbose) cat("Loading risk assessment data...\n")
  
  risk_data <- read.table(risk_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  if (opt$verbose) cat("Loaded risk data for", nrow(risk_data), "guides\n")
  
  return(risk_data)
}

# Function to calculate specificity scores
calculate_specificity_scores <- function(guides, off_targets) {
  if (opt$verbose) cat("Calculating specificity scores...\n")
  
  specificity_scores <- sapply(guides$guide_id, function(guide_id) {
    guide_off_targets <- off_targets[off_targets$guide_id == guide_id, ]
    
    if (nrow(guide_off_targets) == 0) {
      return(1.0)  # Perfect specificity if no off-targets
    }
    
    # Count different types of off-targets
    high_score_count <- sum(guide_off_targets$score > 0.8, na.rm = TRUE)
    medium_score_count <- sum(guide_off_targets$score > 0.6 & guide_off_targets$score <= 0.8, na.rm = TRUE)
    low_score_count <- sum(guide_off_targets$score <= 0.6, na.rm = TRUE)
    
    # Calculate specificity penalty
    penalty <- (high_score_count * 0.5) + (medium_score_count * 0.3) + (low_score_count * 0.1)
    
    # Specificity score (higher is better)
    specificity <- max(0, 1 - (penalty / 10))
    
    return(specificity)
  })
  
  return(specificity_scores)
}

# Function to calculate activity scores (enhanced guide scores)
calculate_activity_scores <- function(guides) {
  if (opt$verbose) cat("Calculating activity scores...\n")
  
  # Use existing guide scores as base activity
  activity_scores <- guides$score
  
  # Enhance with additional factors if available
  if ("gc_content" %in% colnames(guides)) {
    # Optimal GC content around 50%
    gc_bonus <- 1 - abs(guides$gc_content - 0.5) * 2
    activity_scores <- activity_scores * (0.8 + gc_bonus * 0.2)
  }
  
  # Avoid homopolymer runs (if sequence available)
  if ("guide_sequence" %in% colnames(guides)) {
    homo_penalty <- sapply(guides$guide_sequence, function(seq) {
      if (grepl("AAAA|TTTT|GGGG|CCCC", seq)) {
        return(0.9)  # 10% penalty
      } else {
        return(1.0)
      }
    })
    activity_scores <- activity_scores * homo_penalty
  }
  
  return(activity_scores)
}

# Function to calculate coverage scores
calculate_coverage_scores <- function(guides) {
  if (opt$verbose) cat("Calculating coverage scores...\n")
  
  # If guides target multiple sites, calculate coverage
  if ("ctcf_site" %in% colnames(guides)) {
    coverage_scores <- sapply(guides$guide_id, function(guide_id) {
      guide_row <- guides[guides$guide_id == guide_id, ]
      ctcf_site <- guide_row$ctcf_site
      
      # Count guides for this CTCF site
      site_guides <- sum(guides$ctcf_site == ctcf_site)
      
      # Diminishing returns for many guides per site
      coverage_score <- min(1.0, site_guides / 3)
      
      return(coverage_score)
    })
  } else {
    # All guides have equal coverage if no site information
    coverage_scores <- rep(1.0, nrow(guides))
  }
  
  return(coverage_scores)
}

# Function to apply risk filtering
apply_risk_filtering <- function(guides, risk_data, max_risk_level = "medium") {
  if (nrow(risk_data) == 0) {
    if (opt$verbose) cat("No risk filtering applied (no risk data)\n")
    return(guides)
  }
  
  if (opt$verbose) cat("Applying risk filtering...\n")
  
  # Define risk level hierarchy
  risk_levels <- c("low" = 1, "medium" = 2, "high" = 3)
  max_risk_num <- risk_levels[max_risk_level]
  
  # Filter guides by risk level
  guides_with_risk <- merge(guides, risk_data[, c("guide_id", "risk_level", "total_off_targets")], 
                           by = "guide_id", all.x = TRUE)
  
  # Assign default risk for guides without risk assessment
  guides_with_risk$risk_level[is.na(guides_with_risk$risk_level)] <- "medium"
  guides_with_risk$total_off_targets[is.na(guides_with_risk$total_off_targets)] <- 0
  
  # Filter by risk level
  acceptable_guides <- guides_with_risk[risk_levels[guides_with_risk$risk_level] <= max_risk_num, ]
  
  # Additional filtering by off-target count
  acceptable_guides <- acceptable_guides[acceptable_guides$total_off_targets <= opt$max_off_targets, ]
  
  if (opt$verbose) cat("Retained", nrow(acceptable_guides), "guides after risk filtering\n")
  
  return(acceptable_guides)
}

# Function to calculate composite optimization scores
calculate_optimization_scores <- function(guides, specificity_scores, activity_scores, coverage_scores) {
  if (opt$verbose) cat("Calculating optimization scores...\n")
  
  # Normalize weights to sum to 1
  total_weight <- opt$weight_activity + opt$weight_specificity + opt$weight_coverage
  w_activity <- opt$weight_activity / total_weight
  w_specificity <- opt$weight_specificity / total_weight
  w_coverage <- opt$weight_coverage / total_weight
  
  # Calculate composite scores based on method
  if (opt$optimization_method == "specificity") {
    # Prioritize specificity
    composite_scores <- specificity_scores * 0.7 + activity_scores * 0.2 + coverage_scores * 0.1
  } else if (opt$optimization_method == "activity") {
    # Prioritize activity
    composite_scores <- activity_scores * 0.7 + specificity_scores * 0.2 + coverage_scores * 0.1
  } else {
    # Balanced approach
    composite_scores <- (activity_scores * w_activity + 
                        specificity_scores * w_specificity + 
                        coverage_scores * w_coverage)
  }
  
  return(composite_scores)
}

# Function to select optimal guides
select_optimal_guides <- function(guides, composite_scores) {
  if (opt$verbose) cat("Selecting optimal guides...\n")
  
  # Add optimization scores to guides
  guides$specificity_score <- calculate_specificity_scores(guides, data.frame())
  guides$activity_score <- calculate_activity_scores(guides)
  guides$coverage_score <- calculate_coverage_scores(guides)
  guides$optimization_score <- composite_scores
  
  # Filter by minimum score
  high_quality_guides <- guides[guides$score >= opt$min_score, ]
  
  if (nrow(high_quality_guides) == 0) {
    warning("No guides meet minimum score threshold")
    return(data.frame())
  }
  
  # Select guides per CTCF site if site information available
  if ("ctcf_site" %in% colnames(high_quality_guides)) {
    selected_guides <- data.frame()
    
    for (site in unique(high_quality_guides$ctcf_site)) {
      site_guides <- high_quality_guides[high_quality_guides$ctcf_site == site, ]
      site_guides <- site_guides[order(site_guides$optimization_score, decreasing = TRUE), ]
      
      # Select top guides for this site
      n_select <- min(opt$max_guides, nrow(site_guides))
      selected_site_guides <- site_guides[1:n_select, ]
      
      selected_guides <- rbind(selected_guides, selected_site_guides)
    }
  } else {
    # Global selection if no site information
    high_quality_guides <- high_quality_guides[order(high_quality_guides$optimization_score, decreasing = TRUE), ]
    n_select <- min(opt$max_guides * 10, nrow(high_quality_guides))  # Assume multiple sites
    selected_guides <- high_quality_guides[1:n_select, ]
  }
  
  if (opt$verbose) cat("Selected", nrow(selected_guides), "optimal guides\n")
  
  return(selected_guides)
}

# Function to generate optimization summary
generate_optimization_summary <- function(original_guides, selected_guides, off_targets) {
  if (opt$verbose) cat("Generating optimization summary...\n")
  
  summary_stats <- list(
    input_guides = nrow(original_guides),
    selected_guides = nrow(selected_guides),
    selection_rate = nrow(selected_guides) / nrow(original_guides),
    optimization_method = opt$optimization_method,
    filters_applied = list(
      min_score = opt$min_score,
      max_risk = opt$max_risk,
      max_off_targets = opt$max_off_targets
    )
  )
  
  if (nrow(selected_guides) > 0) {
    summary_stats$selected_guide_stats <- list(
      mean_score = mean(selected_guides$score),
      mean_optimization_score = mean(selected_guides$optimization_score),
      mean_activity_score = mean(selected_guides$activity_score),
      mean_specificity_score = mean(selected_guides$specificity_score),
      score_range = range(selected_guides$score),
      optimization_score_range = range(selected_guides$optimization_score)
    )
    
    # CTCF site coverage if available
    if ("ctcf_site" %in% colnames(selected_guides)) {
      summary_stats$site_coverage <- list(
        total_sites_targeted = length(unique(selected_guides$ctcf_site)),
        guides_per_site = table(selected_guides$ctcf_site),
        mean_guides_per_site = mean(table(selected_guides$ctcf_site))
      )
    }
  }
  
  return(summary_stats)
}

# Function to create optimization visualizations
create_optimization_plots <- function(guides, selected_guides, output_dir) {
  if (opt$verbose) cat("Creating optimization visualizations...\n")
  
  if (nrow(selected_guides) == 0) {
    if (opt$verbose) cat("No selected guides to plot\n")
    return()
  }
  
  # Score comparison: all vs selected
  score_comparison <- data.frame(
    score = c(guides$score, selected_guides$score),
    type = c(rep("All Guides", nrow(guides)), rep("Selected Guides", nrow(selected_guides)))
  )
  
  p1 <- ggplot(score_comparison, aes(x = score, fill = type)) +
    geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
    labs(title = "Guide Score Distribution: All vs Selected",
         x = "Guide Score", y = "Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2", name = "Guide Set")
  
  ggsave(file.path(output_dir, "score_comparison.png"), p1, width = 10, height = 6)
  
  # Optimization score vs individual components
  if (all(c("activity_score", "specificity_score") %in% colnames(selected_guides))) {
    p2 <- ggplot(selected_guides, aes(x = activity_score, y = specificity_score, 
                                     size = optimization_score, color = score)) +
      geom_point(alpha = 0.7) +
      labs(title = "Activity vs Specificity for Selected Guides",
           x = "Activity Score", y = "Specificity Score") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red", name = "Guide Score") +
      scale_size_continuous(name = "Optimization Score")
    
    ggsave(file.path(output_dir, "activity_vs_specificity.png"), p2, width = 10, height = 8)
  }
  
  # Guides per site (if applicable)
  if ("ctcf_site" %in% colnames(selected_guides)) {
    guides_per_site <- as.data.frame(table(selected_guides$ctcf_site))
    colnames(guides_per_site) <- c("Site", "Count")
    
    p3 <- ggplot(guides_per_site, aes(x = Count)) +
      geom_histogram(bins = 10, fill = "lightblue", alpha = 0.7) +
      labs(title = "Distribution of Selected Guides per CTCF Site",
           x = "Number of Guides", y = "Number of Sites") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "guides_per_site.png"), p3, width = 8, height = 6)
  }
  
  # Optimization score ranking
  if (nrow(selected_guides) > 0) {
    selected_guides$rank <- rank(-selected_guides$optimization_score, ties.method = "min")
    
    p4 <- ggplot(selected_guides, aes(x = rank, y = optimization_score)) +
      geom_point(color = "darkgreen", alpha = 0.7) +
      geom_smooth(method = "loess", se = FALSE, color = "red") +
      labs(title = "Optimization Score by Rank",
           x = "Rank", y = "Optimization Score") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "optimization_ranking.png"), p4, width = 8, height = 6)
  }
  
  if (opt$verbose) cat("Plots saved to:", output_dir, "\n")
}

# Main execution
tryCatch({
  # Load data
  guides <- load_guide_data(opt$guides)
  off_targets <- load_off_target_data(opt$off_targets)
  risk_data <- load_risk_data(opt$risk_assessment)
  
  # Apply risk filtering
  filtered_guides <- apply_risk_filtering(guides, risk_data, opt$max_risk)
  
  if (nrow(filtered_guides) == 0) {
    stop("No guides remain after risk filtering.")
  }
  
  # Calculate scoring components
  specificity_scores <- calculate_specificity_scores(filtered_guides, off_targets)
  activity_scores <- calculate_activity_scores(filtered_guides)
  coverage_scores <- calculate_coverage_scores(filtered_guides)
  
  # Calculate optimization scores
  optimization_scores <- calculate_optimization_scores(filtered_guides, specificity_scores, 
                                                      activity_scores, coverage_scores)
  
  # Select optimal guides
  selected_guides <- select_optimal_guides(filtered_guides, optimization_scores)
  
  # Generate summary
  summary_stats <- generate_optimization_summary(guides, selected_guides, off_targets)
  
  # Save results
  if (nrow(selected_guides) > 0) {
    write.table(selected_guides, file.path(opt$output, "optimized_guides.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save top-ranked guides separately
    top_guides <- selected_guides[order(selected_guides$optimization_score, decreasing = TRUE), ]
    top_10 <- top_guides[1:min(10, nrow(top_guides)), ]
    
    write.table(top_10, file.path(opt$output, "top_10_guides.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Save summary as JSON
  writeLines(toJSON(summary_stats, pretty = TRUE),
             file.path(opt$output, "optimization_summary.json"))
  
  # Create visualizations
  create_optimization_plots(guides, selected_guides, opt$output)
  
  # Print summary
  cat("\n=== Guide Selection Optimization Summary ===\n")
  cat("Input guides:", summary_stats$input_guides, "\n")
  cat("Selected guides:", summary_stats$selected_guides, "\n")
  cat("Selection rate:", round(summary_stats$selection_rate * 100, 1), "%\n")
  cat("Optimization method:", summary_stats$optimization_method, "\n")
  
  if (nrow(selected_guides) > 0) {
    cat("Mean guide score:", round(summary_stats$selected_guide_stats$mean_score, 3), "\n")
    cat("Mean optimization score:", round(summary_stats$selected_guide_stats$mean_optimization_score, 3), "\n")
    
    if (!is.null(summary_stats$site_coverage)) {
      cat("Sites targeted:", summary_stats$site_coverage$total_sites_targeted, "\n")
      cat("Mean guides per site:", round(summary_stats$site_coverage$mean_guides_per_site, 2), "\n")
    }
  }
  
  cat("\nResults saved to:", opt$output, "\n")
  
}, error = function(e) {
  cat("Error in guide selection optimization:", e$message, "\n")
  quit(status = 1)
})

if (opt$verbose) cat("Guide selection optimization completed successfully!\n")
