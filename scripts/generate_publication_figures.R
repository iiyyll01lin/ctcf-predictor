#!/usr/bin/env Rscript

# Generate Publication-Quality Figures
# Creates high-quality figures for research publications
# Author: CTCF PWM Testing Pipeline Team

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  library(viridis)
  library(ggseqlogo)
  library(optparse)
  library(Biostrings)
})

#' Main function for generating publication figures
#' @param results_dir Directory containing analysis results
#' @param output_dir Output directory for figures
#' @param format Output format (png, pdf, svg)
#' @param dpi Resolution for raster formats
#' @param width Figure width in inches
#' @param height Figure height in inches
#' @param verbose Enable verbose output
generate_publication_figures <- function(results_dir = "results", 
                                       output_dir = "results/figures",
                                       format = "png", dpi = 300,
                                       width = 8, height = 6, verbose = FALSE) {
  
  if (verbose) cat("Generating publication figures...\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set up plotting parameters
  setup_publication_theme()
  
  # Generate individual figures
  figures <- list()
  
  # Figure 1: PWM Information Content Comparison
  if (verbose) cat("Generating PWM comparison figure...\n")
  figures$pwm_comparison <- generate_pwm_comparison_figure(results_dir, verbose)
  
  # Figure 2: Sequence Logo Plots
  if (verbose) cat("Generating sequence logo plots...\n")
  figures$sequence_logos <- generate_sequence_logo_figure(results_dir, verbose)
  
  # Figure 3: Statistical Validation
  if (verbose) cat("Generating statistical validation figure...\n")
  figures$statistical_validation <- generate_statistical_figure(results_dir, verbose)
  
  # Figure 4: Method Performance Comparison
  if (verbose) cat("Generating method comparison figure...\n")
  figures$method_comparison <- generate_method_comparison_figure(results_dir, verbose)
  
  # Figure 5: Quality Assessment
  if (verbose) cat("Generating quality assessment figure...\n")
  figures$quality_assessment <- generate_quality_figure(results_dir, verbose)
  
  # Save all figures
  save_publication_figures(figures, output_dir, format, dpi, width, height, verbose)
  
  # Generate combined figure
  if (verbose) cat("Generating combined figure...\n")
  combined_figure <- create_combined_figure(figures)
  save_figure(combined_figure, file.path(output_dir, paste0("combined_analysis.", format)),
              width = width * 2, height = height * 2, dpi = dpi)
  
  if (verbose) cat("Publication figures saved to:", output_dir, "\n")
  
  return(figures)
}

#' Set up publication-quality theme
setup_publication_theme <- function() {
  theme_set(theme_minimal() +
    theme(
      text = element_text(size = 12, family = "Arial"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "lightgray", color = "black")
    ))
}

#' Generate PWM comparison figure
generate_pwm_comparison_figure <- function(results_dir, verbose) {
  
  # Load PWM results
  pwm_files <- list.files(results_dir, pattern = ".*pwm.*\\.rds$", full.names = TRUE)
  
  if (length(pwm_files) == 0) {
    if (verbose) cat("  No PWM files found for comparison\n")
    return(NULL)
  }
  
  # Extract information content data
  ic_data <- data.frame()
  
  for (pwm_file in pwm_files) {
    tryCatch({
      pwm_result <- readRDS(pwm_file)
      
      # Extract method name from filename
      method_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", pwm_file)
      method_name <- gsub("_pwm|pwm_", "", method_name)
      
      # Get information content
      if ("total_info" %in% names(pwm_result)) {
        total_ic <- pwm_result$total_info
      } else if ("info_content" %in% names(pwm_result)) {
        total_ic <- sum(pwm_result$info_content, na.rm = TRUE)
      } else {
        next
      }
      
      # Get position-wise IC if available
      if ("info_content" %in% names(pwm_result) && is.vector(pwm_result$info_content)) {
        pos_ic <- pwm_result$info_content
        for (i in seq_along(pos_ic)) {
          ic_data <- rbind(ic_data, data.frame(
            Method = method_name,
            Position = i,
            Information_Content = pos_ic[i],
            Total_IC = total_ic,
            stringsAsFactors = FALSE
          ))
        }
      } else {
        # If no position data, create single entry
        ic_data <- rbind(ic_data, data.frame(
          Method = method_name,
          Position = 1,
          Information_Content = total_ic,
          Total_IC = total_ic,
          stringsAsFactors = FALSE
        ))
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error loading", basename(pwm_file), ":", conditionMessage(e), "\n")
    })
  }
  
  if (nrow(ic_data) == 0) {
    if (verbose) cat("  No valid IC data found\n")
    return(NULL)
  }
  
  # Create comparison plots
  p1 <- ggplot(ic_data, aes(x = Position, y = Information_Content, color = Method)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    labs(title = "Position-wise Information Content",
         x = "Position", y = "Information Content (bits)") +
    scale_color_viridis_d() +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", alpha = 0.7) +
    theme(legend.position = "bottom")
  
  # Total IC comparison
  total_ic_summary <- unique(ic_data[, c("Method", "Total_IC")])
  
  p2 <- ggplot(total_ic_summary, aes(x = reorder(Method, Total_IC), y = Total_IC, fill = Method)) +
    geom_col() +
    geom_text(aes(label = round(Total_IC, 2)), hjust = -0.1) +
    coord_flip() +
    labs(title = "Total Information Content Comparison",
         x = "Method", y = "Total Information Content (bits)") +
    scale_fill_viridis_d() +
    theme(legend.position = "none")
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, ncol = 2, 
                               top = textGrob("PWM Information Content Analysis", 
                                            gp = gpar(fontsize = 16, fontface = "bold")))
  
  return(combined_plot)
}

#' Generate sequence logo figure
generate_sequence_logo_figure <- function(results_dir, verbose) {
  
  # Find best PWM file
  pwm_files <- list.files(results_dir, pattern = "subset_pwm_size.*\\.rds$|best.*pwm.*\\.rds$", 
                         full.names = TRUE)
  
  if (length(pwm_files) == 0) {
    pwm_files <- list.files(results_dir, pattern = ".*pwm.*\\.rds$", full.names = TRUE)
  }
  
  if (length(pwm_files) == 0) {
    if (verbose) cat("  No PWM files found for sequence logo\n")
    return(NULL)
  }
  
  # Use first available PWM
  pwm_file <- pwm_files[1]
  
  tryCatch({
    pwm_result <- readRDS(pwm_file)
    
    # Extract PWM matrix
    if ("pwm" %in% names(pwm_result)) {
      pwm_matrix <- pwm_result$pwm
    } else if ("prob_matrix" %in% names(pwm_result)) {
      pwm_matrix <- pwm_result$prob_matrix
    } else {
      if (verbose) cat("  No PWM matrix found in", basename(pwm_file), "\n")
      return(NULL)
    }
    
    # Ensure proper format for ggseqlogo
    if (is.matrix(pwm_matrix)) {
      if (nrow(pwm_matrix) == 4) {
        rownames(pwm_matrix) <- c("A", "C", "G", "T")
      }
      pwm_matrix <- t(pwm_matrix)  # ggseqlogo expects positions as rows
    }
    
    # Create sequence logo
    logo_plot <- ggseqlogo(pwm_matrix, method = "prob") +
      labs(title = "CTCF Binding Site Sequence Logo",
           subtitle = paste("Source:", basename(pwm_file)),
           x = "Position", y = "Probability") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold")
      )
    
    return(logo_plot)
    
  }, error = function(e) {
    if (verbose) cat("  Error creating sequence logo:", conditionMessage(e), "\n")
    return(NULL)
  })
}

#' Generate statistical validation figure
generate_statistical_figure <- function(results_dir, verbose) {
  
  # Look for statistical analysis results
  stat_files <- list.files(results_dir, pattern = "statistical.*\\.rds$|.*stats.*\\.rds$", 
                          full.names = TRUE)
  
  if (length(stat_files) == 0) {
    if (verbose) cat("  No statistical files found\n")
    return(NULL)
  }
  
  # Load statistical results
  stat_data <- data.frame()
  
  for (stat_file in stat_files) {
    tryCatch({
      stat_result <- readRDS(stat_file)
      
      # Extract p-values and other statistics
      if ("p_value" %in% names(stat_result)) {
        method_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", stat_file)
        
        stat_data <- rbind(stat_data, data.frame(
          Method = method_name,
          P_Value = stat_result$p_value,
          Significant = stat_result$p_value < 0.05,
          stringsAsFactors = FALSE
        ))
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error loading", basename(stat_file), ":", conditionMessage(e), "\n")
    })
  }
  
  if (nrow(stat_data) == 0) {
    if (verbose) cat("  No statistical data found\n")
    return(NULL)
  }
  
  # Create statistical significance plot
  stat_plot <- ggplot(stat_data, aes(x = reorder(Method, -log10(P_Value)), 
                                    y = -log10(P_Value), 
                                    fill = Significant)) +
    geom_col() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "darkblue"),
                     name = "Significant\n(p < 0.05)") +
    labs(title = "Statistical Significance of PWM Models",
         x = "Method", y = "-log10(p-value)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(stat_plot)
}

#' Generate method comparison figure
generate_method_comparison_figure <- function(results_dir, verbose) {
  
  # Look for comparison results
  comparison_files <- list.files(results_dir, pattern = "comparison.*\\.rds$|.*comparison.*\\.rds$", 
                                full.names = TRUE)
  
  if (length(comparison_files) == 0) {
    if (verbose) cat("  No comparison files found\n")
    return(NULL)
  }
  
  # Load comparison data
  comparison_result <- readRDS(comparison_files[1])
  
  if (!("method_scores" %in% names(comparison_result))) {
    if (verbose) cat("  No method scores found in comparison file\n")
    return(NULL)
  }
  
  # Extract method performance data
  method_data <- data.frame()
  
  for (method_name in names(comparison_result$method_scores)) {
    method_score <- comparison_result$method_scores[[method_name]]
    
    method_data <- rbind(method_data, data.frame(
      Method = method_name,
      Weighted_Score = method_score$weighted_score,
      IC_Score = method_score$component_scores$information_content,
      Conservation_Score = method_score$component_scores$conservation,
      Alignment_Score = method_score$component_scores$alignment_quality,
      Time_Score = method_score$component_scores$processing_time,
      Meets_Thresholds = method_score$meets_thresholds,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create radar chart data
  radar_data <- method_data[, c("Method", "IC_Score", "Conservation_Score", 
                               "Alignment_Score", "Time_Score")]
  
  # Reshape for plotting
  radar_long <- reshape2::melt(radar_data, id.vars = "Method", 
                              variable.name = "Metric", value.name = "Score")
  
  # Create comparison plot
  comparison_plot <- ggplot(radar_long, aes(x = Metric, y = Score, color = Method, group = Method)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Method Performance Comparison",
         x = "Performance Metric", y = "Normalized Score") +
    scale_color_viridis_d() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(comparison_plot)
}

#' Generate quality assessment figure
generate_quality_figure <- function(results_dir, verbose) {
  
  # Look for quality assessment files
  quality_files <- list.files(results_dir, pattern = "quality.*\\.rds$|.*quality.*\\.rds$", 
                             full.names = TRUE)
  
  if (length(quality_files) == 0) {
    if (verbose) cat("  No quality files found\n")
    return(NULL)
  }
  
  # Load quality data
  quality_data <- data.frame()
  
  for (quality_file in quality_files) {
    tryCatch({
      quality_result <- readRDS(quality_file)
      
      method_name <- gsub(".*?([^/\\\\]+)\\.rds$", "\\1", quality_file)
      
      # Extract quality metrics
      if ("assessment" %in% names(quality_result)) {
        assessment <- quality_result$assessment
        
        quality_data <- rbind(quality_data, data.frame(
          Method = method_name,
          Information_Content = assessment$information_content %||% NA,
          Conservation_Score = assessment$conservation_score %||% NA,
          Pattern_Quality = assessment$pattern_quality %||% NA,
          Overall_Grade = assessment$overall_grade %||% "Unknown",
          stringsAsFactors = FALSE
        ))
      }
      
    }, error = function(e) {
      if (verbose) cat("  Error loading", basename(quality_file), ":", conditionMessage(e), "\n")
    })
  }
  
  if (nrow(quality_data) == 0) {
    if (verbose) cat("  No quality data found\n")
    return(NULL)
  }
  
  # Create quality assessment plot
  quality_long <- reshape2::melt(quality_data[, c("Method", "Information_Content", 
                                                 "Conservation_Score", "Pattern_Quality")],
                                id.vars = "Method", variable.name = "Quality_Metric", 
                                value.name = "Score")
  
  quality_long <- na.omit(quality_long)
  
  quality_plot <- ggplot(quality_long, aes(x = Method, y = Score, fill = Quality_Metric)) +
    geom_col(position = "dodge") +
    labs(title = "PWM Quality Assessment",
         x = "Method", y = "Quality Score") +
    scale_fill_brewer(type = "qual", palette = "Set2", name = "Quality\nMetric") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(quality_plot)
}

#' Save publication figures
save_publication_figures <- function(figures, output_dir, format, dpi, width, height, verbose) {
  
  for (figure_name in names(figures)) {
    if (!is.null(figures[[figure_name]])) {
      filename <- file.path(output_dir, paste0(figure_name, ".", format))
      save_figure(figures[[figure_name]], filename, width, height, dpi)
      
      if (verbose) cat("  Saved:", basename(filename), "\n")
    }
  }
}

#' Save individual figure
save_figure <- function(plot, filename, width, height, dpi) {
  
  if (grepl("\\.png$", filename)) {
    ggsave(filename, plot, width = width, height = height, dpi = dpi, device = "png")
  } else if (grepl("\\.pdf$", filename)) {
    ggsave(filename, plot, width = width, height = height, device = "pdf")
  } else if (grepl("\\.svg$", filename)) {
    ggsave(filename, plot, width = width, height = height, device = "svg")
  } else {
    # Default to PNG
    ggsave(paste0(filename, ".png"), plot, width = width, height = height, dpi = dpi, device = "png")
  }
}

#' Create combined figure
create_combined_figure <- function(figures) {
  
  valid_figures <- figures[!sapply(figures, is.null)]
  
  if (length(valid_figures) == 0) {
    return(NULL)
  }
  
  # Arrange figures in grid
  if (length(valid_figures) >= 4) {
    combined <- grid.arrange(grobs = valid_figures[1:4], ncol = 2, nrow = 2,
                           top = textGrob("CTCF PWM Analysis - Complete Results", 
                                        gp = gpar(fontsize = 18, fontface = "bold")))
  } else {
    combined <- grid.arrange(grobs = valid_figures, ncol = 2,
                           top = textGrob("CTCF PWM Analysis - Complete Results", 
                                        gp = gpar(fontsize = 18, fontface = "bold")))
  }
  
  return(combined)
}

# Utility function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Command-line interface
if (!interactive()) {
  option_list <- list(
    make_option(c("-r", "--results"), type = "character", default = "results",
                help = "Results directory [default: %default]", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = "results/figures",
                help = "Output directory [default: %default]", metavar = "character"),
    make_option(c("-f", "--format"), type = "character", default = "png",
                help = "Output format (png, pdf, svg) [default: %default]", metavar = "character"),
    make_option(c("-d", "--dpi"), type = "integer", default = 300,
                help = "Resolution for raster formats [default: %default]", metavar = "integer"),
    make_option(c("-w", "--width"), type = "double", default = 8,
                help = "Figure width in inches [default: %default]", metavar = "double"),
    make_option(c("--height"), type = "double", default = 6,
                help = "Figure height in inches [default: %default]", metavar = "double"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Enable verbose output")
  )
  
  opt_parser <- OptionParser(option_list = option_list,
                            description = "Generate Publication-Quality Figures")
  opt <- parse_args(opt_parser)
  
  if (!dir.exists(opt$results)) {
    stop("Results directory does not exist: ", opt$results, call. = FALSE)
  }
  
  # Generate publication figures
  tryCatch({
    figures <- generate_publication_figures(
      results_dir = opt$results,
      output_dir = opt$output,
      format = opt$format,
      dpi = opt$dpi,
      width = opt$width,
      height = opt$height,
      verbose = opt$verbose
    )
    
    cat("Publication figures generated successfully.\n")
    
  }, error = function(e) {
    cat("Error generating publication figures:", conditionMessage(e), "\n")
    quit(status = 1)
  })
}
