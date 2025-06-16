#!/usr/bin/env Rscript

# CRISPR Guide Design for CTCF Sites Script
# Designs CRISPR guide RNAs targeting CTCF binding sites
# Author: CTCF Predictor Pipeline
# Usage: Rscript design_crispr_guides.R --input <ctcf_sites> --output <output_dir> [options]

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
              help="Input CTCF binding sites file (BED or TSV)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="crispr_guides_output",
              help="Output directory [default=%default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Reference genome file (FASTA)", metavar="character"),
  make_option(c("--pam"), type="character", default="NGG",
              help="PAM sequence for guide design [default=%default]", metavar="character"),
  make_option(c("--guide_length"), type="numeric", default=20,
              help="Guide RNA length [default=%default]", metavar="numeric"),
  make_option(c("--extend_region"), type="numeric", default=200,
              help="Extend CTCF site region for guide search [default=%default]", metavar="numeric"),
  make_option(c("--max_guides"), type="numeric", default=5,
              help="Maximum guides per CTCF site [default=%default]", metavar="numeric"),
  make_option(c("--gc_min"), type="numeric", default=0.2,
              help="Minimum GC content [default=%default]", metavar="numeric"),
  make_option(c("--gc_max"), type="numeric", default=0.8,
              help="Maximum GC content [default=%default]", metavar="numeric"),
  make_option(c("--filter_repeats"), action="store_true", default=FALSE,
              help="Filter guides overlapping repetitive elements"),
  make_option(c("--score_threshold"), type="numeric", default=0.6,
              help="Minimum guide score threshold [default=%default]", metavar="numeric"),
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
  stop("Input CTCF sites file must be specified.", call.=FALSE)
}

if (!file.exists(opt$input)) {
  stop("Input file does not exist: ", opt$input, call.=FALSE)
}

# Create output directory
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

if (opt$verbose) {
  cat("Starting CRISPR guide design...\n")
  cat("Input file:", opt$input, "\n")
  cat("Output directory:", opt$output, "\n")
  cat("PAM sequence:", opt$pam, "\n")
  cat("Guide length:", opt$guide_length, "\n")
}

# Function to load CTCF binding sites
load_ctcf_sites <- function(input_file) {
  if (opt$verbose) cat("Loading CTCF binding sites...\n")
  
  # Try different file formats
  if (grepl("\\.bed$", input_file, ignore.case = TRUE)) {
    # BED format
    sites <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("chr", "start", "end", "name", "score", "strand"))
  } else {
    # Tab-delimited format
    sites <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  }
  
  # Ensure required columns exist
  if (!all(c("chr", "start", "end") %in% colnames(sites))) {
    stop("Input file must contain chr, start, and end columns.")
  }
  
  # Add site IDs if not present
  if (!"name" %in% colnames(sites)) {
    sites$name <- paste0("CTCF_site_", seq_len(nrow(sites)))
  }
  
  if (opt$verbose) cat("Loaded", nrow(sites), "CTCF binding sites\n")
  
  return(sites)
}

# Function to load reference genome (simplified - in practice would use proper genome loading)
load_reference_genome <- function(genome_file) {
  if (opt$verbose) cat("Loading reference genome...\n")
  
  if (!is.null(genome_file) && file.exists(genome_file)) {
    # In a real implementation, this would load the actual genome
    # For this example, we'll simulate genome sequences
    genome <- list()
    return(genome)
  } else {
    # Simulate genome for demonstration
    if (opt$verbose) cat("Using simulated genome sequences\n")
    return(NULL)
  }
}

# Function to generate sequence for a genomic region (simplified)
get_genomic_sequence <- function(chr, start, end, genome = NULL) {
  # In practice, this would extract actual genomic sequence
  # For demonstration, generate random sequence
  seq_length <- end - start + 1
  random_seq <- paste(sample(c("A", "T", "G", "C"), seq_length, replace = TRUE), collapse = "")
  return(DNAString(random_seq))
}

# Function to find PAM sites in sequence
find_pam_sites <- function(sequence, pam_pattern = "NGG") {
  if (opt$verbose && nchar(sequence) > 1000) cat("Searching for PAM sites in", nchar(sequence), "bp sequence...\n")
  
  # Convert PAM pattern to regex
  pam_regex <- gsub("N", "[ATGC]", pam_pattern)
  
  # Find forward PAM sites
  forward_matches <- gregexpr(pam_regex, as.character(sequence))[[1]]
  forward_positions <- if (forward_matches[1] != -1) forward_matches else integer(0)
  
  # Find reverse complement PAM sites
  rev_sequence <- as.character(reverseComplement(sequence))
  reverse_matches <- gregexpr(pam_regex, rev_sequence)[[1]]
  reverse_positions <- if (reverse_matches[1] != -1) {
    nchar(sequence) - reverse_matches - nchar(pam_pattern) + 2
  } else {
    integer(0)
  }
  
  # Combine positions with strand information
  pam_sites <- data.frame(
    position = c(forward_positions, reverse_positions),
    strand = c(rep("+", length(forward_positions)), rep("-", length(reverse_positions))),
    stringsAsFactors = FALSE
  )
  
  return(pam_sites)
}

# Function to extract guide sequences
extract_guide_sequences <- function(sequence, pam_sites, guide_length = 20) {
  guides <- data.frame(
    position = integer(0),
    strand = character(0),
    guide_sequence = character(0),
    pam_sequence = character(0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(pam_sites))) {
    pos <- pam_sites$position[i]
    strand <- pam_sites$strand[i]
    
    if (strand == "+") {
      # Forward strand: guide is upstream of PAM
      guide_start <- pos - guide_length
      guide_end <- pos - 1
      pam_start <- pos
      pam_end <- pos + 2
      
      if (guide_start > 0 && pam_end <= nchar(sequence)) {
        guide_seq <- as.character(subseq(sequence, guide_start, guide_end))
        pam_seq <- as.character(subseq(sequence, pam_start, pam_end))
        
        guides <- rbind(guides, data.frame(
          position = guide_start,
          strand = strand,
          guide_sequence = guide_seq,
          pam_sequence = pam_seq,
          stringsAsFactors = FALSE
        ))
      }
    } else {
      # Reverse strand: guide is downstream of PAM (in reverse complement)
      guide_start <- pos + 3
      guide_end <- pos + guide_length + 2
      pam_start <- pos
      pam_end <- pos + 2
      
      if (guide_end <= nchar(sequence) && pam_start > 0) {
        guide_seq <- as.character(reverseComplement(subseq(sequence, guide_start, guide_end)))
        pam_seq <- as.character(reverseComplement(subseq(sequence, pam_start, pam_end)))
        
        guides <- rbind(guides, data.frame(
          position = guide_start,
          strand = strand,
          guide_sequence = guide_seq,
          pam_sequence = pam_seq,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(guides)
}

# Function to calculate GC content
calculate_gc_content <- function(sequence) {
  gc_count <- sum(strsplit(sequence, "")[[1]] %in% c("G", "C"))
  return(gc_count / nchar(sequence))
}

# Function to score guide RNAs
score_guide_rnas <- function(guides) {
  if (opt$verbose) cat("Scoring guide RNAs...\n")
  
  scores <- sapply(guides$guide_sequence, function(guide) {
    # GC content score
    gc_content <- calculate_gc_content(guide)
    gc_score <- if (gc_content >= opt$gc_min && gc_content <= opt$gc_max) {
      1 - abs(gc_content - 0.5) * 2  # Optimal around 50%
    } else {
      0.2
    }
    
    # Avoid homopolymer runs
    homo_penalty <- 0
    for (base in c("A", "T", "G", "C")) {
      if (grepl(paste0(base, "{4,}"), guide)) {
        homo_penalty <- homo_penalty + 0.3
      }
    }
    
    # Avoid secondary structure (simplified)
    structure_penalty <- 0
    if (grepl("AAAA|TTTT", guide)) {
      structure_penalty <- structure_penalty + 0.2
    }
    
    # Calculate final score
    final_score <- pmax(0, gc_score - homo_penalty - structure_penalty)
    return(final_score)
  })
  
  guides$score <- scores
  guides$gc_content <- sapply(guides$guide_sequence, calculate_gc_content)
  
  return(guides)
}

# Function to filter and rank guides
filter_and_rank_guides <- function(guides, max_guides = 5) {
  if (opt$verbose) cat("Filtering and ranking guides...\n")
  
  # Filter by score threshold
  high_quality_guides <- guides[guides$score >= opt$score_threshold, ]
  
  # Remove guides with extreme GC content
  high_quality_guides <- high_quality_guides[
    high_quality_guides$gc_content >= opt$gc_min & 
    high_quality_guides$gc_content <= opt$gc_max, ]
  
  # Rank by score
  high_quality_guides <- high_quality_guides[order(high_quality_guides$score, decreasing = TRUE), ]
  
  # Limit number of guides
  if (nrow(high_quality_guides) > max_guides) {
    high_quality_guides <- high_quality_guides[1:max_guides, ]
  }
  
  return(high_quality_guides)
}

# Function to design guides for a single CTCF site
design_guides_for_site <- function(site, genome = NULL) {
  if (opt$verbose) cat("Designing guides for site:", site$name, "\n")
  
  # Extend region around CTCF site
  extended_start <- max(1, site$start - opt$extend_region)
  extended_end <- site$end + opt$extend_region
  
  # Get genomic sequence
  sequence <- get_genomic_sequence(site$chr, extended_start, extended_end, genome)
  
  # Find PAM sites
  pam_sites <- find_pam_sites(sequence, opt$pam)
  
  if (nrow(pam_sites) == 0) {
    if (opt$verbose) cat("No PAM sites found for site:", site$name, "\n")
    return(data.frame())
  }
  
  # Extract guide sequences
  guides <- extract_guide_sequences(sequence, pam_sites, opt$guide_length)
  
  if (nrow(guides) == 0) {
    if (opt$verbose) cat("No valid guides found for site:", site$name, "\n")
    return(data.frame())
  }
  
  # Score guides
  guides <- score_guide_rnas(guides)
  
  # Filter and rank
  final_guides <- filter_and_rank_guides(guides, opt$max_guides)
  
  # Add site information
  if (nrow(final_guides) > 0) {
    final_guides$ctcf_site <- site$name
    final_guides$ctcf_chr <- site$chr
    final_guides$ctcf_start <- site$start
    final_guides$ctcf_end <- site$end
    final_guides$genomic_position <- extended_start + final_guides$position - 1
    final_guides$guide_id <- paste0(site$name, "_guide_", seq_len(nrow(final_guides)))
  }
  
  return(final_guides)
}

# Function to design guides for all sites
design_all_guides <- function(sites, genome = NULL) {
  if (opt$verbose) cat("Designing guides for", nrow(sites), "CTCF sites...\n")
  
  if (opt$cores > 1) {
    # Parallel processing
    all_guides <- mclapply(1:nrow(sites), function(i) {
      design_guides_for_site(sites[i, ], genome)
    }, mc.cores = opt$cores)
  } else {
    # Sequential processing
    all_guides <- lapply(1:nrow(sites), function(i) {
      design_guides_for_site(sites[i, ], genome)
    })
  }
  
  # Combine results
  combined_guides <- do.call(rbind, all_guides)
  
  if (opt$verbose) cat("Designed", nrow(combined_guides), "total guides\n")
  
  return(combined_guides)
}

# Function to generate guide design summary
generate_guide_summary <- function(guides, sites) {
  if (opt$verbose) cat("Generating guide design summary...\n")
  
  summary_stats <- list(
    total_ctcf_sites = nrow(sites),
    total_guides_designed = nrow(guides),
    sites_with_guides = length(unique(guides$ctcf_site)),
    sites_without_guides = nrow(sites) - length(unique(guides$ctcf_site)),
    mean_guides_per_site = nrow(guides) / length(unique(guides$ctcf_site)),
    mean_guide_score = mean(guides$score),
    mean_gc_content = mean(guides$gc_content),
    score_distribution = summary(guides$score),
    gc_distribution = summary(guides$gc_content),
    strand_distribution = table(guides$strand)
  )
  
  return(summary_stats)
}

# Function to create guide design visualizations
create_guide_plots <- function(guides, summary_stats, output_dir) {
  if (opt$verbose) cat("Creating guide design visualizations...\n")
  
  if (nrow(guides) == 0) {
    if (opt$verbose) cat("No guides to plot\n")
    return()
  }
  
  # Score distribution
  p1 <- ggplot(guides, aes(x = score)) +
    geom_histogram(bins = 30, fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = opt$score_threshold, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Guide RNA Scores",
         x = "Guide Score", y = "Count") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "guide_scores.png"), p1, width = 8, height = 6)
  
  # GC content distribution
  p2 <- ggplot(guides, aes(x = gc_content)) +
    geom_histogram(bins = 20, fill = "lightgreen", alpha = 0.7) +
    geom_vline(xintercept = c(opt$gc_min, opt$gc_max), color = "red", linetype = "dashed") +
    labs(title = "Distribution of Guide RNA GC Content",
         x = "GC Content", y = "Count") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "gc_content.png"), p2, width = 8, height = 6)
  
  # Score vs GC content
  p3 <- ggplot(guides, aes(x = gc_content, y = score, color = strand)) +
    geom_point(alpha = 0.7) +
    labs(title = "Guide Score vs GC Content",
         x = "GC Content", y = "Guide Score") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  ggsave(file.path(output_dir, "score_vs_gc.png"), p3, width = 8, height = 6)
  
  # Guides per site
  guides_per_site <- as.data.frame(table(guides$ctcf_site))
  colnames(guides_per_site) <- c("Site", "Count")
  
  p4 <- ggplot(guides_per_site, aes(x = Count)) +
    geom_histogram(bins = 10, fill = "orange", alpha = 0.7) +
    labs(title = "Number of Guides per CTCF Site",
         x = "Number of Guides", y = "Number of Sites") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "guides_per_site.png"), p4, width = 8, height = 6)
  
  if (opt$verbose) cat("Plots saved to:", output_dir, "\n")
}

# Main execution
tryCatch({
  # Load CTCF sites
  sites <- load_ctcf_sites(opt$input)
  
  # Load reference genome
  genome <- load_reference_genome(opt$genome)
  
  # Design guides for all sites
  guides <- design_all_guides(sites, genome)
  
  if (nrow(guides) == 0) {
    warning("No guides were designed for any CTCF sites.")
  }
  
  # Generate summary
  summary_stats <- generate_guide_summary(guides, sites)
  
  # Save results
  if (nrow(guides) > 0) {
    write.table(guides, file.path(opt$output, "designed_guides.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save high-scoring guides separately
    high_scoring <- guides[guides$score >= 0.8, ]
    if (nrow(high_scoring) > 0) {
      write.table(high_scoring, file.path(opt$output, "high_scoring_guides.tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
  # Save summary as JSON
  writeLines(toJSON(summary_stats, pretty = TRUE),
             file.path(opt$output, "guide_design_summary.json"))
  
  # Create visualizations
  if (nrow(guides) > 0) {
    create_guide_plots(guides, summary_stats, opt$output)
  }
  
  # Print summary
  cat("\n=== CRISPR Guide Design Summary ===\n")
  cat("Total CTCF sites:", summary_stats$total_ctcf_sites, "\n")
  cat("Total guides designed:", summary_stats$total_guides_designed, "\n")
  cat("Sites with guides:", summary_stats$sites_with_guides, "\n")
  cat("Sites without guides:", summary_stats$sites_without_guides, "\n")
  
  if (summary_stats$total_guides_designed > 0) {
    cat("Mean guides per site:", round(summary_stats$mean_guides_per_site, 2), "\n")
    cat("Mean guide score:", round(summary_stats$mean_guide_score, 3), "\n")
    cat("Mean GC content:", round(summary_stats$mean_gc_content, 3), "\n")
  }
  
  cat("\nResults saved to:", opt$output, "\n")
  
}, error = function(e) {
  cat("Error in CRISPR guide design:", e$message, "\n")
  quit(status = 1)
})

if (opt$verbose) cat("CRISPR guide design completed successfully!\n")
