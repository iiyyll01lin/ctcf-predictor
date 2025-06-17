#!/usr/bin/env Rscript

# PWM Visualization Generation Script
# Creates publication-quality visualizations for PWM analysis
# Author: CTCF PWM Testing Pipeline Team
# Date: 2025-06-16

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input PWM file (RDS format)", metavar="file"),
  make_option(c("-o", "--output-dir"), type="character", default="results/visualizations",
              help="Output directory for visualizations [default: %default]"),
  make_option(c("-f", "--format"), type="character", default="pdf",
              help="Output format: pdf, png, svg [default: %default]"),
  make_option(c("--width"), type="numeric", default=10,
              help="Figure width in inches [default: %default]"),
  make_option(c("--height"), type="numeric", default=6,
              help="Figure height in inches [default: %default]"),
  make_option(c("--dpi"), type="numeric", default=300,
              help="Resolution for raster formats [default: %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Enable verbose output")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Generate PWM visualizations")
opt <- parse_args(opt_parser)

#' Generate sequence logo visualization
#' @param pwm_matrix PWM matrix (4 x N)
#' @param output_file Output file path
#' @param title Plot title
#' @param width Figure width
#' @param height Figure height
#' @param format Output format
#' @param dpi Resolution
#' @param verbose Enable verbose output
generate_sequence_logo <- function(pwm_matrix, output_file, title = "CTCF Sequence Logo",
                                  width = 10, height = 6, format = "pdf", dpi = 300,
                                  verbose = FALSE) {
  if (verbose) cat("Generating sequence logo:", output_file, "\n")
  
  # Ensure row names are nucleotides
  if (is.null(rownames(pwm_matrix))) {
    rownames(pwm_matrix) <- c("A", "C", "G", "T")
  }
  
  # Create sequence logo data
  positions <- 1:ncol(pwm_matrix)
  logo_data <- data.frame()
  
  for (i in positions) {
    # Calculate information content for position
    col <- pwm_matrix[, i]
    col[col == 0] <- 1e-10  # Avoid log(0)
    ic <- 2 + sum(col * log2(col))
    
    # Create data for each nucleotide
    for (nuc in rownames(pwm_matrix)) {
      height <- col[nuc] * ic
      if (height > 0.01) {  # Only include significant heights
        logo_data <- rbind(logo_data, data.frame(
          position = i,
          nucleotide = nuc,
          height = height,
          frequency = col[nuc],
          ic = ic
        ))
      }
    }
  }
  
  # Create color scheme
  nuc_colors <- c("A" = "#FF6B6B", "C" = "#4ECDC4", "G" = "#45B7D1", "T" = "#96CEB4")
  
  # Create plot
  p <- ggplot(logo_data, aes(x = position, y = height, fill = nucleotide)) +
    geom_col(position = "stack", width = 0.8) +
    scale_fill_manual(values = nuc_colors) +
    scale_x_continuous(breaks = 1:max(positions), labels = 1:max(positions)) +
    labs(
      title = title,
      subtitle = paste("Total Information Content:", round(sum(unique(logo_data$ic)), 2), "bits"),
      x = "Position",
      y = "Information Content (bits)",
      fill = "Nucleotide"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # Save plot
  if (format == "pdf") {
    ggsave(output_file, p, width = width, height = height, device = "pdf")
  } else if (format == "png") {
    ggsave(output_file, p, width = width, height = height, dpi = dpi, device = "png")
  } else if (format == "svg") {
    ggsave(output_file, p, width = width, height = height, device = "svg")
  }
  
  if (verbose) cat("Sequence logo saved to:", output_file, "\n")
  return(output_file)
}

#' Generate information content profile plot
#' @param pwm_matrix PWM matrix (4 x N)
#' @param output_file Output file path
#' @param title Plot title
#' @param width Figure width
#' @param height Figure height
#' @param format Output format
#' @param dpi Resolution
#' @param verbose Enable verbose output
plot_ic_profile <- function(pwm_matrix, output_file, title = "Information Content Profile",
                           width = 10, height = 6, format = "pdf", dpi = 300,
                           verbose = FALSE) {
  if (verbose) cat("Generating IC profile plot:", output_file, "\n")
  
  # Calculate information content for each position
  ic_profile <- apply(pwm_matrix, 2, function(col) {
    col[col == 0] <- 1e-10  # Avoid log(0)
    2 + sum(col * log2(col))
  })
  
  positions <- 1:length(ic_profile)
  
  # Create data frame
  ic_data <- data.frame(
    position = positions,
    information_content = ic_profile
  )
  
  # Add expected CTCF pattern indicators
  expected_high_ic <- c(1, 2, 3, 4, 8, 9, 13, 14, 15)
  ic_data$expected_high <- ic_data$position %in% expected_high_ic
  
  # Create plot
  p <- ggplot(ic_data, aes(x = position, y = information_content)) +
    geom_col(aes(fill = expected_high), width = 0.8, alpha = 0.8) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_manual(
      values = c("FALSE" = "#BDC3C7", "TRUE" = "#3498DB"),
      labels = c("Variable", "Expected High IC"),
      name = "Position Type"
    ) +
    scale_x_continuous(breaks = positions, labels = positions) +
    labs(
      title = title,
      subtitle = paste("Total IC:", round(sum(ic_profile), 2), "bits,",
                      "Conserved positions (IC>1.0):", sum(ic_profile > 1.0)),
      x = "Position",
      y = "Information Content (bits)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = max(positions) * 0.8, y = max(ic_profile) * 0.9,
             label = "Threshold: 1.0 bits", color = "red", size = 3)
  
  # Save plot
  if (format == "pdf") {
    ggsave(output_file, p, width = width, height = height, device = "pdf")
  } else if (format == "png") {
    ggsave(output_file, p, width = width, height = height, dpi = dpi, device = "png")
  } else if (format == "svg") {
    ggsave(output_file, p, width = width, height = height, device = "svg")
  }
  
  if (verbose) cat("IC profile plot saved to:", output_file, "\n")
  return(output_file)
}

#' Generate conservation heatmap
#' @param pwm_matrix PWM matrix (4 x N)
#' @param output_file Output file path
#' @param title Plot title
#' @param width Figure width
#' @param height Figure height
#' @param format Output format
#' @param dpi Resolution
#' @param verbose Enable verbose output
generate_conservation_heatmap <- function(pwm_matrix, output_file, 
                                        title = "Conservation Pattern Heatmap",
                                        width = 10, height = 6, format = "pdf", dpi = 300,
                                        verbose = FALSE) {
  if (verbose) cat("Generating conservation heatmap:", output_file, "\n")
  
  # Ensure row names
  if (is.null(rownames(pwm_matrix))) {
    rownames(pwm_matrix) <- c("A", "C", "G", "T")
  }
  
  # Prepare data for heatmap
  heatmap_data <- data.frame()
  for (i in 1:ncol(pwm_matrix)) {
    for (j in 1:nrow(pwm_matrix)) {
      heatmap_data <- rbind(heatmap_data, data.frame(
        position = i,
        nucleotide = rownames(pwm_matrix)[j],
        frequency = pwm_matrix[j, i]
      ))
    }
  }
  
  # Reorder nucleotides for better visualization
  heatmap_data$nucleotide <- factor(heatmap_data$nucleotide, levels = c("T", "G", "C", "A"))
  
  # Create heatmap
  p <- ggplot(heatmap_data, aes(x = position, y = nucleotide, fill = frequency)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradientn(
      colors = c("#F8F9FA", "#3498DB", "#E74C3C"),
      values = c(0, 0.5, 1),
      name = "Frequency"
    ) +
    scale_x_continuous(breaks = 1:ncol(pwm_matrix), labels = 1:ncol(pwm_matrix)) +
    labs(
      title = title,
      subtitle = "Nucleotide frequencies across binding site positions",
      x = "Position",
      y = "Nucleotide"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid = element_blank()
    ) +
    coord_fixed(ratio = 1)
  
  # Add text annotations for high frequencies
  high_freq_data <- heatmap_data[heatmap_data$frequency > 0.6, ]
  if (nrow(high_freq_data) > 0) {
    p <- p + geom_text(data = high_freq_data, 
                      aes(label = sprintf("%.2f", frequency)), 
                      color = "white", size = 3, fontface = "bold")
  }
  
  # Save plot
  if (format == "pdf") {
    ggsave(output_file, p, width = width, height = height, device = "pdf")
  } else if (format == "png") {
    ggsave(output_file, p, width = width, height = height, dpi = dpi, device = "png")
  } else if (format == "svg") {
    ggsave(output_file, p, width = width, height = height, device = "svg")
  }
  
  if (verbose) cat("Conservation heatmap saved to:", output_file, "\n")
  return(output_file)
}

#' Generate alignment quality comparison plot
#' @param before_ic Information content before alignment
#' @param after_ic Information content after alignment
#' @param output_file Output file path
#' @param width Figure width
#' @param height Figure height
#' @param format Output format
#' @param dpi Resolution
#' @param verbose Enable verbose output
plot_alignment_quality <- function(before_ic, after_ic, output_file,
                                  width = 10, height = 6, format = "pdf", dpi = 300,
                                  verbose = FALSE) {
  if (verbose) cat("Generating alignment quality plot:", output_file, "\n")
  
  # Create comparison data
  positions <- 1:max(length(before_ic), length(after_ic))
  
  # Pad shorter vector with zeros
  if (length(before_ic) < length(positions)) {
    before_ic <- c(before_ic, rep(0, length(positions) - length(before_ic)))
  }
  if (length(after_ic) < length(positions)) {
    after_ic <- c(after_ic, rep(0, length(positions) - length(after_ic)))
  }
  
  comparison_data <- data.frame(
    position = rep(positions, 2),
    information_content = c(before_ic, after_ic),
    condition = rep(c("Before Alignment", "After Alignment"), each = length(positions))
  )
  
  # Create plot
  p <- ggplot(comparison_data, aes(x = position, y = information_content, fill = condition)) +
    geom_col(position = "dodge", width = 0.8, alpha = 0.8) +
    scale_fill_manual(
      values = c("Before Alignment" = "#E74C3C", "After Alignment" = "#27AE60"),
      name = "Condition"
    ) +
    scale_x_continuous(breaks = positions, labels = positions) +
    labs(
      title = "Alignment Quality Improvement",
      subtitle = paste("Total IC improvement:", 
                      round(sum(after_ic) - sum(before_ic), 2), "bits"),
      x = "Position",
      y = "Information Content (bits)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  # Save plot
  if (format == "pdf") {
    ggsave(output_file, p, width = width, height = height, device = "pdf")
  } else if (format == "png") {
    ggsave(output_file, p, width = width, height = height, dpi = dpi, device = "png")
  } else if (format == "svg") {
    ggsave(output_file, p, width = width, height = height, device = "svg")
  }
  
  if (verbose) cat("Alignment quality plot saved to:", output_file, "\n")
  return(output_file)
}

#' Generate all visualizations for a PWM
#' @param pwm_data PWM data object
#' @param output_dir Output directory
#' @param base_name Base name for output files
#' @param format Output format
#' @param width Figure width
#' @param height Figure height
#' @param dpi Resolution
#' @param verbose Enable verbose output
generate_all_visualizations <- function(pwm_data, output_dir = "results/visualizations",
                                       base_name = "ctcf_pwm", format = "pdf",
                                       width = 10, height = 6, dpi = 300,
                                       verbose = FALSE) {
  if (verbose) cat("Generating all PWM visualizations...\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created output directory:", output_dir, "\n")
  }
  
  # Extract PWM matrix
  if (is.list(pwm_data) && "pwm" %in% names(pwm_data)) {
    pwm_matrix <- pwm_data$pwm
  } else if (is.matrix(pwm_data)) {
    pwm_matrix <- pwm_data
  } else {
    stop("Invalid PWM data format")
  }
  
  generated_files <- list()
  
  # 1. Sequence logo
  logo_file <- file.path(output_dir, paste0(base_name, "_sequence_logo.", format))
  generated_files$sequence_logo <- generate_sequence_logo(
    pwm_matrix, logo_file, width = width, height = height, 
    format = format, dpi = dpi, verbose = verbose
  )
  
  # 2. IC profile
  ic_file <- file.path(output_dir, paste0(base_name, "_ic_profile.", format))
  generated_files$ic_profile <- plot_ic_profile(
    pwm_matrix, ic_file, width = width, height = height,
    format = format, dpi = dpi, verbose = verbose
  )
  
  # 3. Conservation heatmap
  heatmap_file <- file.path(output_dir, paste0(base_name, "_conservation_heatmap.", format))
  generated_files$conservation_heatmap <- generate_conservation_heatmap(
    pwm_matrix, heatmap_file, width = width, height = height,
    format = format, dpi = dpi, verbose = verbose
  )
  
  if (verbose) {
    cat("All visualizations generated:\n")
    for (name in names(generated_files)) {
      cat("  ", name, ":", generated_files[[name]], "\n")
    }
  }
  
  return(generated_files)
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
  
  if (opt$verbose) {
    cat("PWM Visualization Generator\n")
    cat("==========================\n")
    cat("Input file:", opt$input, "\n")
    cat("Output directory:", opt$`output-dir`, "\n")
    cat("Format:", opt$format, "\n")
    cat("Dimensions:", opt$width, "x", opt$height, "\n\n")
  }
  
  # Load PWM data
  if (opt$verbose) cat("Loading PWM data...\n")
  pwm_data <- readRDS(opt$input)
  
  # Generate visualizations
  base_name <- tools::file_path_sans_ext(basename(opt$input))
  generated_files <- generate_all_visualizations(
    pwm_data, opt$`output-dir`, base_name, opt$format,
    opt$width, opt$height, opt$dpi, opt$verbose
  )
  
  if (opt$verbose) {
    cat("\nVisualization generation complete!\n")
    cat("Generated", length(generated_files), "visualization files\n")
  }
}
