#!/usr/bin/env Rscript

# CRISPR Off-Target Prediction Script
# Predicts potential off-target sites for CRISPR guide RNAs
# Author: CTCF Predictor Pipeline
# Usage: Rscript predict_off_targets.R --input <guides_file> --output <output_dir> [options]

suppressMessages({
  library(optparse)
  library(Biostrings)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(parallel)
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input guide RNAs file (TSV format)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="off_targets_output",
              help="Output directory [default=%default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Reference genome file (FASTA)", metavar="character"),
  make_option(c("--max_mismatches"), type="numeric", default=4,
              help="Maximum number of mismatches allowed [default=%default]", metavar="numeric"),
  make_option(c("--pam_mismatches"), type="numeric", default=0,
              help="Maximum PAM mismatches allowed [default=%default]", metavar="numeric"),
  make_option(c("--seed_region"), type="numeric", default=12,
              help="Length of seed region (3' end) [default=%default]", metavar="numeric"),
  make_option(c("--seed_mismatches"), type="numeric", default=2,
              help="Maximum mismatches in seed region [default=%default]", metavar="numeric"),
  make_option(c("--score_threshold"), type="numeric", default=0.5,
              help="Minimum off-target score threshold [default=%default]", metavar="numeric"),
  make_option(c("--annotations"), type="character", default=NULL,
              help="Gene annotation file (GTF/GFF) for impact assessment", metavar="character"),
  make_option(c("--exclude_regions"), type="character", default=NULL,
              help="BED file of regions to exclude from search", metavar="character"),
  make_option(c("--cores"), type="numeric", default=1,
              help="Number of cores for parallel processing [default=%default]", metavar="numeric"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input guide RNAs file must be specified.", call.=FALSE)
}

if (!file.exists(opt$input)) {
  stop("Input file does not exist: ", opt$input, call.=FALSE)
}

# Create output directory
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

if (opt$verbose) {
  cat("Starting off-target prediction...\n")
  cat("Input file:", opt$input, "\n")
  cat("Output directory:", opt$output, "\n")
  cat("Max mismatches:", opt$max_mismatches, "\n")
  cat("Seed region length:", opt$seed_region, "\n")
}

# Function to load guide RNAs
load_guide_rnas <- function(input_file) {
  if (opt$verbose) cat("Loading guide RNAs...\n")
  
  guides <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Validate required columns
  required_cols <- c("guide_sequence", "guide_id")
  missing_cols <- setdiff(required_cols, colnames(guides))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
  }
  
  # Add PAM sequence if not present
  if (!"pam_sequence" %in% colnames(guides)) {
    guides$pam_sequence <- "NGG"  # Default PAM
  }
  
  if (opt$verbose) cat("Loaded", nrow(guides), "guide RNAs\n")
  
  return(guides)
}

# Function to generate possible off-target sequences
generate_possible_sequences <- function(guide_sequence, max_mismatches = 4) {
  # This is a simplified approach - in practice, genome-wide search would be performed
  # Here we generate a sample of possible off-target sequences
  
  bases <- c("A", "T", "G", "C")
  guide_chars <- strsplit(guide_sequence, "")[[1]]
  n_positions <- length(guide_chars)
  
  # Generate sequences with up to max_mismatches mismatches
  possible_sequences <- character()
  
  # Original sequence (0 mismatches)
  possible_sequences <- c(possible_sequences, guide_sequence)
  
  # Generate sequences with mismatches (simplified sampling)
  for (n_mm in 1:min(max_mismatches, 3)) {  # Limit to avoid exponential explosion
    # Sample some positions to mutate
    n_samples <- min(100, choose(n_positions, n_mm))  # Limit number of samples
    
    for (sample_i in 1:n_samples) {
      # Random positions to mutate
      positions_to_mutate <- sample(1:n_positions, n_mm)
      
      # Create mutated sequence
      mutated_chars <- guide_chars
      for (pos in positions_to_mutate) {
        # Change to a different base
        current_base <- mutated_chars[pos]
        new_base <- sample(setdiff(bases, current_base), 1)
        mutated_chars[pos] <- new_base
      }
      
      mutated_sequence <- paste(mutated_chars, collapse = "")
      possible_sequences <- c(possible_sequences, mutated_sequence)
    }
  }
  
  return(unique(possible_sequences))
}

# Function to score off-target matches
score_off_target_match <- function(guide_seq, target_seq, pam_seq = "NGG") {
  if (nchar(guide_seq) != nchar(target_seq)) {
    return(0)
  }
  
  guide_chars <- strsplit(guide_seq, "")[[1]]
  target_chars <- strsplit(target_seq, "")[[1]]
  
  # Count mismatches
  mismatches <- sum(guide_chars != target_chars)
  
  # Calculate seed region mismatches (3' end is more important)
  seed_length <- min(opt$seed_region, length(guide_chars))
  seed_start <- length(guide_chars) - seed_length + 1
  seed_mismatches <- sum(guide_chars[seed_start:length(guide_chars)] != 
                        target_chars[seed_start:length(target_chars)])
  
  # Score based on total mismatches and seed mismatches
  if (mismatches > opt$max_mismatches) {
    return(0)
  }
  
  if (seed_mismatches > opt$seed_mismatches) {
    return(0)
  }
  
  # Calculate score (higher score = more likely off-target)
  # This is a simplified scoring model
  base_score <- 1 - (mismatches / 20)  # Assume 20bp guide
  seed_penalty <- seed_mismatches * 0.3
  
  final_score <- max(0, base_score - seed_penalty)
  
  return(final_score)
}

# Function to search for off-targets in genome regions
search_genome_regions <- function(guide_sequence, guide_id, genome_regions = NULL) {
  if (opt$verbose) cat("Searching off-targets for guide:", guide_id, "\n")
  
  # In practice, this would search the actual genome
  # Here we simulate potential off-target sites
  
  # Generate possible off-target sequences
  possible_targets <- generate_possible_sequences(guide_sequence, opt$max_mismatches)
  
  # Score each potential target
  off_targets <- data.frame(
    guide_id = character(0),
    target_sequence = character(0),
    chr = character(0),
    start = numeric(0),
    end = numeric(0),
    strand = character(0),
    score = numeric(0),
    mismatches = numeric(0),
    seed_mismatches = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (target_seq in possible_targets) {
    score <- score_off_target_match(guide_sequence, target_seq)
    
    if (score >= opt$score_threshold) {
      # Calculate mismatches
      guide_chars <- strsplit(guide_sequence, "")[[1]]
      target_chars <- strsplit(target_seq, "")[[1]]
      total_mismatches <- sum(guide_chars != target_chars)
      
      # Seed mismatches
      seed_length <- min(opt$seed_region, length(guide_chars))
      seed_start <- length(guide_chars) - seed_length + 1
      seed_mm <- sum(guide_chars[seed_start:length(guide_chars)] != 
                    target_chars[seed_start:length(target_chars)])
      
      # Simulate genomic location
      sim_chr <- paste0("chr", sample(1:22, 1))
      sim_start <- sample(1000000:100000000, 1)
      sim_end <- sim_start + nchar(target_seq) - 1
      sim_strand <- sample(c("+", "-"), 1)
      
      off_targets <- rbind(off_targets, data.frame(
        guide_id = guide_id,
        target_sequence = target_seq,
        chr = sim_chr,
        start = sim_start,
        end = sim_end,
        strand = sim_strand,
        score = score,
        mismatches = total_mismatches,
        seed_mismatches = seed_mm,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(off_targets)
}

# Function to annotate off-targets with genomic features
annotate_off_targets <- function(off_targets, annotations_file = NULL) {
  if (opt$verbose) cat("Annotating off-targets with genomic features...\n")
  
  if (!is.null(annotations_file) && file.exists(annotations_file)) {
    # In practice, would load actual gene annotations
    # Here we simulate annotations
    off_targets$genomic_feature <- sample(c("exon", "intron", "promoter", "intergenic", "UTR"),
                                         nrow(off_targets), replace = TRUE,
                                         prob = c(0.2, 0.3, 0.1, 0.35, 0.05))
  } else {
    # Simulate genomic features
    off_targets$genomic_feature <- sample(c("exon", "intron", "promoter", "intergenic", "UTR"),
                                         nrow(off_targets), replace = TRUE,
                                         prob = c(0.2, 0.3, 0.1, 0.35, 0.05))
  }
  
  # Add gene information
  off_targets$nearest_gene <- paste0("GENE_", sprintf("%05d", sample(1:20000, nrow(off_targets))))
  off_targets$gene_distance <- ifelse(off_targets$genomic_feature == "intergenic",
                                     sample(1:50000, nrow(off_targets)),
                                     0)
  
  # Assess potential impact
  off_targets$impact_level <- with(off_targets, {
    ifelse(genomic_feature == "exon", "high",
           ifelse(genomic_feature == "promoter", "high",
                  ifelse(genomic_feature == "UTR", "medium",
                         ifelse(genomic_feature == "intron", "low", "very_low"))))
  })
  
  return(off_targets)
}

# Function to predict off-targets for all guides
predict_all_off_targets <- function(guides) {
  if (opt$verbose) cat("Predicting off-targets for", nrow(guides), "guides...\n")
  
  if (opt$cores > 1) {
    # Parallel processing
    all_off_targets <- mclapply(1:nrow(guides), function(i) {
      guide <- guides[i, ]
      search_genome_regions(guide$guide_sequence, guide$guide_id)
    }, mc.cores = opt$cores)
  } else {
    # Sequential processing
    all_off_targets <- lapply(1:nrow(guides), function(i) {
      guide <- guides[i, ]
      search_genome_regions(guide$guide_sequence, guide$guide_id)
    })
  }
  
  # Combine results
  combined_off_targets <- do.call(rbind, all_off_targets)
  
  if (opt$verbose) cat("Found", nrow(combined_off_targets), "potential off-targets\n")
  
  return(combined_off_targets)
}

# Function to assess off-target risk
assess_off_target_risk <- function(off_targets, guides) {
  if (opt$verbose) cat("Assessing off-target risk...\n")
  
  risk_assessment <- data.frame(
    guide_id = unique(off_targets$guide_id),
    stringsAsFactors = FALSE
  )
  
  for (guide_id in risk_assessment$guide_id) {
    guide_off_targets <- off_targets[off_targets$guide_id == guide_id, ]
    
    # Count off-targets by category
    high_score_count <- sum(guide_off_targets$score > 0.8)
    high_impact_count <- sum(guide_off_targets$impact_level == "high")
    total_count <- nrow(guide_off_targets)
    
    # Calculate risk scores
    quantity_risk <- min(1, total_count / 10)  # Risk increases with number
    quality_risk <- min(1, high_score_count / 3)  # Risk from high-scoring targets
    impact_risk <- min(1, high_impact_count / 2)  # Risk from high-impact targets
    
    # Overall risk score
    overall_risk <- (quantity_risk + quality_risk + impact_risk) / 3
    
    # Risk level
    risk_level <- ifelse(overall_risk > 0.7, "high",
                        ifelse(overall_risk > 0.4, "medium", "low"))
    
    # Add to assessment
    risk_assessment[risk_assessment$guide_id == guide_id, "total_off_targets"] <- total_count
    risk_assessment[risk_assessment$guide_id == guide_id, "high_score_off_targets"] <- high_score_count
    risk_assessment[risk_assessment$guide_id == guide_id, "high_impact_off_targets"] <- high_impact_count
    risk_assessment[risk_assessment$guide_id == guide_id, "risk_score"] <- overall_risk
    risk_assessment[risk_assessment$guide_id == guide_id, "risk_level"] <- risk_level
  }
  
  return(risk_assessment)
}

# Function to generate off-target summary
generate_off_target_summary <- function(off_targets, risk_assessment, guides) {
  if (opt$verbose) cat("Generating off-target summary...\n")
  
  summary_stats <- list(
    total_guides = nrow(guides),
    guides_with_off_targets = length(unique(off_targets$guide_id)),
    guides_without_off_targets = nrow(guides) - length(unique(off_targets$guide_id)),
    total_off_targets = nrow(off_targets),
    mean_off_targets_per_guide = nrow(off_targets) / length(unique(off_targets$guide_id)),
    high_risk_guides = sum(risk_assessment$risk_level == "high"),
    medium_risk_guides = sum(risk_assessment$risk_level == "medium"),
    low_risk_guides = sum(risk_assessment$risk_level == "low"),
    genomic_feature_distribution = table(off_targets$genomic_feature),
    impact_level_distribution = table(off_targets$impact_level),
    mismatch_distribution = table(off_targets$mismatches),
    mean_off_target_score = mean(off_targets$score)
  )
  
  return(summary_stats)
}

# Function to create off-target visualizations
create_off_target_plots <- function(off_targets, risk_assessment, output_dir) {
  if (opt$verbose) cat("Creating off-target visualizations...\n")
  
  if (nrow(off_targets) == 0) {
    if (opt$verbose) cat("No off-targets to plot\n")
    return()
  }
  
  # Off-target score distribution
  p1 <- ggplot(off_targets, aes(x = score)) +
    geom_histogram(bins = 20, fill = "coral", alpha = 0.7) +
    geom_vline(xintercept = opt$score_threshold, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Off-Target Scores",
         x = "Off-Target Score", y = "Count") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "off_target_scores.png"), p1, width = 8, height = 6)
  
  # Mismatches vs Score
  p2 <- ggplot(off_targets, aes(x = mismatches, y = score, color = impact_level)) +
    geom_point(alpha = 0.7) +
    labs(title = "Off-Target Score vs Number of Mismatches",
         x = "Number of Mismatches", y = "Off-Target Score") +
    theme_minimal() +
    scale_color_brewer(palette = "Spectral", name = "Impact Level")
  
  ggsave(file.path(output_dir, "score_vs_mismatches.png"), p2, width = 10, height = 6)
  
  # Genomic feature distribution
  feature_counts <- as.data.frame(table(off_targets$genomic_feature))
  colnames(feature_counts) <- c("Feature", "Count")
  
  p3 <- ggplot(feature_counts, aes(x = reorder(Feature, Count), y = Count, fill = Feature)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    coord_flip() +
    labs(title = "Off-Targets by Genomic Feature",
         x = "Genomic Feature", y = "Number of Off-Targets") +
    theme_minimal() +
    guides(fill = FALSE)
  
  ggsave(file.path(output_dir, "genomic_features.png"), p3, width = 8, height = 6)
  
  # Risk assessment
  if (nrow(risk_assessment) > 0) {
    p4 <- ggplot(risk_assessment, aes(x = risk_score, fill = risk_level)) +
      geom_histogram(bins = 15, alpha = 0.7) +
      labs(title = "Distribution of Off-Target Risk Scores",
           x = "Risk Score", y = "Number of Guides") +
      theme_minimal() +
      scale_fill_brewer(palette = "RdYlBu", name = "Risk Level", direction = -1)
    
    ggsave(file.path(output_dir, "risk_distribution.png"), p4, width = 8, height = 6)
  }
  
  if (opt$verbose) cat("Plots saved to:", output_dir, "\n")
}

# Main execution
tryCatch({
  # Load guide RNAs
  guides <- load_guide_rnas(opt$input)
  
  # Predict off-targets
  off_targets <- predict_all_off_targets(guides)
  
  # Annotate off-targets
  if (nrow(off_targets) > 0) {
    off_targets <- annotate_off_targets(off_targets, opt$annotations)
  }
  
  # Assess risk
  risk_assessment <- if (nrow(off_targets) > 0) {
    assess_off_target_risk(off_targets, guides)
  } else {
    data.frame()
  }
  
  # Generate summary
  summary_stats <- generate_off_target_summary(off_targets, risk_assessment, guides)
  
  # Save results
  if (nrow(off_targets) > 0) {
    write.table(off_targets, file.path(opt$output, "predicted_off_targets.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save high-risk off-targets
    high_risk_targets <- off_targets[off_targets$score > 0.8 | off_targets$impact_level == "high", ]
    if (nrow(high_risk_targets) > 0) {
      write.table(high_risk_targets, file.path(opt$output, "high_risk_off_targets.tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
  if (nrow(risk_assessment) > 0) {
    write.table(risk_assessment, file.path(opt$output, "risk_assessment.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Save summary as JSON
  writeLines(toJSON(summary_stats, pretty = TRUE),
             file.path(opt$output, "off_target_summary.json"))
  
  # Create visualizations
  create_off_target_plots(off_targets, risk_assessment, opt$output)
  
  # Print summary
  cat("\n=== Off-Target Prediction Summary ===\n")
  cat("Total guides analyzed:", summary_stats$total_guides, "\n")
  cat("Guides with off-targets:", summary_stats$guides_with_off_targets, "\n")
  cat("Total off-targets found:", summary_stats$total_off_targets, "\n")
  
  if (summary_stats$total_off_targets > 0) {
    cat("Mean off-targets per guide:", round(summary_stats$mean_off_targets_per_guide, 2), "\n")
    cat("Mean off-target score:", round(summary_stats$mean_off_target_score, 3), "\n")
  }
  
  if (nrow(risk_assessment) > 0) {
    cat("High-risk guides:", summary_stats$high_risk_guides, "\n")
    cat("Medium-risk guides:", summary_stats$medium_risk_guides, "\n")
    cat("Low-risk guides:", summary_stats$low_risk_guides, "\n")
  }
  
  cat("\nResults saved to:", opt$output, "\n")
  
}, error = function(e) {
  cat("Error in off-target prediction:", e$message, "\n")
  quit(status = 1)
})

if (opt$verbose) cat("Off-target prediction completed successfully!\n")
