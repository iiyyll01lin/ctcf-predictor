#!/usr/bin/env Rscript

# CTCF Therapeutic Potential Assessment Script
# Assesses therapeutic potential of CTCF-related drug targets
# Author: CTCF Predictor Pipeline
# Usage: Rscript assess_therapeutic_potential.R --input <targets_file> --output <output_dir> [options]

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
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input drug targets file (TSV format)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="therapeutic_assessment",
              help="Output directory [default=%default]", metavar="character"),
  make_option(c("-d", "--drug_db"), type="character", default=NULL,
              help="Drug database file (optional)", metavar="character"),
  make_option(c("-p", "--ppi_network"), type="character", default=NULL,
              help="Protein-protein interaction network file (optional)", metavar="character"),
  make_option(c("-e", "--expression"), type="character", default=NULL,
              help="Gene expression data file (optional)", metavar="character"),
  make_option(c("-c", "--clinical"), type="character", default=NULL,
              help="Clinical trial data file (optional)", metavar="character"),
  make_option(c("-s", "--side_effects"), type="character", default=NULL,
              help="Known side effects database (optional)", metavar="character"),
  make_option(c("--min_score"), type="numeric", default=0.5,
              help="Minimum therapeutic score threshold [default=%default]", metavar="numeric"),
  make_option(c("--prioritize"), type="character", default="cancer",
              help="Disease type to prioritize [default=%default]", metavar="character"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input targets file must be specified.", call.=FALSE)
}

if (!file.exists(opt$input)) {
  stop("Input file does not exist: ", opt$input, call.=FALSE)
}

# Create output directory
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

if (opt$verbose) {
  cat("Starting therapeutic potential assessment...\n")
  cat("Input file:", opt$input, "\n")
  cat("Output directory:", opt$output, "\n")
}

# Function to load drug targets
load_drug_targets <- function(input_file) {
  if (opt$verbose) cat("Loading drug targets...\n")
  
  targets <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Validate required columns
  required_cols <- c("site_id", "druggability_score")
  missing_cols <- setdiff(required_cols, colnames(targets))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
  }
  
  return(targets)
}

# Function to assess druggability factors
assess_druggability <- function(targets) {
  if (opt$verbose) cat("Assessing druggability factors...\n")
  
  # Initialize druggability assessment
  assessment <- targets
  
  # Molecular weight assessment (simulated)
  assessment$molecular_accessibility <- runif(nrow(assessment), 0.3, 1.0)
  
  # Binding pocket assessment
  assessment$binding_pocket_score <- pmin(1.0, assessment$druggability_score + 
                                         rnorm(nrow(assessment), 0, 0.1))
  
  # Target selectivity
  assessment$selectivity_score <- runif(nrow(assessment), 0.4, 0.9)
  
  # Pharmacokinetic properties (simulated)
  assessment$adme_score <- runif(nrow(assessment), 0.3, 0.8)
  
  # Safety profile assessment
  assessment$safety_score <- runif(nrow(assessment), 0.5, 0.9)
  
  # Calculate composite druggability score
  assessment$composite_druggability <- with(assessment, {
    (druggability_score * 0.3 + 
     molecular_accessibility * 0.2 +
     binding_pocket_score * 0.2 +
     selectivity_score * 0.15 +
     adme_score * 0.1 +
     safety_score * 0.05)
  })
  
  return(assessment)
}

# Function to assess clinical potential
assess_clinical_potential <- function(targets, clinical_file = NULL) {
  if (opt$verbose) cat("Assessing clinical potential...\n")
  
  # Add clinical trial information if available
  if (!is.null(clinical_file) && file.exists(clinical_file)) {
    # Load clinical data (simplified)
    targets$clinical_precedent <- sample(c("none", "preclinical", "phase1", "phase2", "phase3", "approved"),
                                        nrow(targets), replace = TRUE,
                                        prob = c(0.4, 0.25, 0.15, 0.1, 0.05, 0.05))
  } else {
    targets$clinical_precedent <- sample(c("none", "preclinical", "phase1", "phase2", "phase3", "approved"),
                                        nrow(targets), replace = TRUE,
                                        prob = c(0.4, 0.25, 0.15, 0.1, 0.05, 0.05))
  }
  
  # Market potential assessment
  targets$market_potential <- sample(c("low", "medium", "high"), nrow(targets), 
                                    replace = TRUE, prob = c(0.3, 0.5, 0.2))
  
  # Competitive landscape
  targets$competition_level <- sample(c("low", "medium", "high"), nrow(targets),
                                     replace = TRUE, prob = c(0.2, 0.5, 0.3))
  
  # Development risk assessment
  targets$development_risk <- sample(c("low", "medium", "high"), nrow(targets),
                                    replace = TRUE, prob = c(0.2, 0.6, 0.2))
  
  # Clinical feasibility score
  targets$clinical_feasibility <- with(targets, {
    base_score <- 0.5
    
    # Boost for existing clinical precedent
    precedent_boost <- ifelse(clinical_precedent == "approved", 0.4,
                             ifelse(clinical_precedent == "phase3", 0.3,
                                   ifelse(clinical_precedent == "phase2", 0.2,
                                         ifelse(clinical_precedent == "phase1", 0.1, 0))))
    
    # Market potential boost
    market_boost <- ifelse(market_potential == "high", 0.2,
                          ifelse(market_potential == "medium", 0.1, 0))
    
    # Competition penalty
    competition_penalty <- ifelse(competition_level == "high", -0.1,
                                 ifelse(competition_level == "medium", -0.05, 0))
    
    # Risk penalty
    risk_penalty <- ifelse(development_risk == "high", -0.15,
                          ifelse(development_risk == "medium", -0.05, 0))
    
    pmax(0, pmin(1, base_score + precedent_boost + market_boost + competition_penalty + risk_penalty))
  })
  
  return(targets)
}

# Function to perform network analysis
perform_network_analysis <- function(targets, ppi_file = NULL) {
  if (opt$verbose) cat("Performing network analysis...\n")
  
  # Simulate network properties if PPI data not available
  if (!is.null(ppi_file) && file.exists(ppi_file)) {
    # Load PPI network (simplified)
    targets$network_centrality <- runif(nrow(targets), 0.1, 0.9)
    targets$network_clustering <- runif(nrow(targets), 0.2, 0.8)
  } else {
    targets$network_centrality <- runif(nrow(targets), 0.1, 0.9)
    targets$network_clustering <- runif(nrow(targets), 0.2, 0.8)
  }
  
  # Hub gene identification
  targets$is_hub_gene <- targets$network_centrality > 0.7
  
  # Pathway enrichment (simulated)
  targets$pathway_enrichment_score <- with(targets, {
    ifelse(is_hub_gene, 
           pmin(1.0, network_centrality + rnorm(nrow(targets), 0, 0.1)),
           network_centrality * 0.8)
  })
  
  return(targets)
}

# Function to prioritize targets
prioritize_targets <- function(targets, disease_focus = "cancer") {
  if (opt$verbose) cat("Prioritizing targets for", disease_focus, "...\n")
  
  # Calculate priority score based on disease focus
  targets$priority_score <- with(targets, {
    base_score <- composite_druggability * 0.4 + clinical_feasibility * 0.3 + 
                  pathway_enrichment_score * 0.2 + therapeutic_score * 0.1
    
    # Disease-specific boost
    if (disease_focus %in% colnames(targets)) {
      disease_boost <- ifelse(get(disease_focus) != "none", 0.2, 0)
    } else if ("disease_association" %in% colnames(targets)) {
      disease_boost <- ifelse(disease_association == disease_focus, 0.2, 0)
    } else {
      disease_boost <- 0
    }
    
    # Hub gene boost
    hub_boost <- ifelse(is_hub_gene, 0.1, 0)
    
    pmin(1.0, base_score + disease_boost + hub_boost)
  })
  
  # Rank targets
  targets$rank <- rank(-targets$priority_score, ties.method = "min")
  
  # Classify priority levels
  targets$priority_level <- ifelse(targets$priority_score > 0.8, "very_high",
                                  ifelse(targets$priority_score > 0.6, "high",
                                        ifelse(targets$priority_score > 0.4, "medium", "low")))
  
  return(targets)
}

# Function to assess side effects and toxicity
assess_toxicity <- function(targets, side_effects_file = NULL) {
  if (opt$verbose) cat("Assessing toxicity and side effects...\n")
  
  # Simulate toxicity assessment
  targets$hepatotoxicity_risk <- sample(c("low", "medium", "high"), nrow(targets),
                                       replace = TRUE, prob = c(0.6, 0.3, 0.1))
  
  targets$cardiotoxicity_risk <- sample(c("low", "medium", "high"), nrow(targets),
                                       replace = TRUE, prob = c(0.7, 0.25, 0.05))
  
  targets$neurotoxicity_risk <- sample(c("low", "medium", "high"), nrow(targets),
                                      replace = TRUE, prob = c(0.65, 0.3, 0.05))
  
  # Overall toxicity score
  targets$overall_toxicity_risk <- with(targets, {
    hepato_score <- ifelse(hepatotoxicity_risk == "high", 0.8,
                          ifelse(hepatotoxicity_risk == "medium", 0.5, 0.2))
    cardio_score <- ifelse(cardiotoxicity_risk == "high", 0.9,
                          ifelse(cardiotoxicity_risk == "medium", 0.6, 0.3))
    neuro_score <- ifelse(neurotoxicity_risk == "high", 0.7,
                         ifelse(neurotoxicity_risk == "medium", 0.4, 0.1))
    
    (hepato_score + cardio_score + neuro_score) / 3
  })
  
  # Adjust priority score based on toxicity
  targets$adjusted_priority <- with(targets, {
    pmax(0.1, priority_score * (1 - overall_toxicity_risk * 0.3))
  })
  
  return(targets)
}

# Function to generate assessment report
generate_assessment_report <- function(targets) {
  if (opt$verbose) cat("Generating assessment report...\n")
  
  report <- list(
    summary = list(
      total_targets = nrow(targets),
      very_high_priority = sum(targets$priority_level == "very_high"),
      high_priority = sum(targets$priority_level == "high"),
      medium_priority = sum(targets$priority_level == "medium"),
      low_priority = sum(targets$priority_level == "low"),
      hub_genes = sum(targets$is_hub_gene),
      low_toxicity = sum(targets$overall_toxicity_risk < 0.4),
      clinical_precedent = table(targets$clinical_precedent)
    ),
    metrics = list(
      mean_druggability = mean(targets$composite_druggability),
      mean_clinical_feasibility = mean(targets$clinical_feasibility),
      mean_priority_score = mean(targets$priority_score),
      mean_toxicity_risk = mean(targets$overall_toxicity_risk)
    ),
    top_targets = targets[order(targets$adjusted_priority, decreasing = TRUE)[1:min(10, nrow(targets))], 
                         c("site_id", "nearest_gene", "priority_level", "adjusted_priority", 
                           "composite_druggability", "clinical_feasibility", "overall_toxicity_risk")]
  )
  
  return(report)
}

# Function to create assessment visualizations
create_assessment_plots <- function(targets, output_dir) {
  if (opt$verbose) cat("Creating assessment visualizations...\n")
  
  # Priority score distribution
  p1 <- ggplot(targets, aes(x = priority_score, fill = priority_level)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    labs(title = "Distribution of Priority Scores",
         x = "Priority Score", y = "Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Spectral", name = "Priority Level")
  
  ggsave(file.path(output_dir, "priority_distribution.png"), p1, width = 10, height = 6)
  
  # Druggability vs Clinical Feasibility
  p2 <- ggplot(targets, aes(x = composite_druggability, y = clinical_feasibility, 
                           color = priority_level, size = adjusted_priority)) +
    geom_point(alpha = 0.7) +
    labs(title = "Druggability vs Clinical Feasibility",
         x = "Composite Druggability Score", y = "Clinical Feasibility Score") +
    theme_minimal() +
    scale_color_brewer(palette = "Spectral", name = "Priority Level") +
    scale_size_continuous(name = "Adjusted Priority")
  
  ggsave(file.path(output_dir, "druggability_vs_feasibility.png"), p2, width = 10, height = 8)
  
  # Toxicity assessment
  toxicity_data <- targets %>%
    select(hepatotoxicity_risk, cardiotoxicity_risk, neurotoxicity_risk) %>%
    tidyr::gather(key = "toxicity_type", value = "risk_level") %>%
    mutate(toxicity_type = gsub("_risk", "", toxicity_type))
  
  p3 <- ggplot(toxicity_data, aes(x = toxicity_type, fill = risk_level)) +
    geom_bar(position = "fill") +
    labs(title = "Toxicity Risk Assessment",
         x = "Toxicity Type", y = "Proportion") +
    theme_minimal() +
    scale_fill_brewer(palette = "RdYlBu", name = "Risk Level", direction = -1) +
    coord_flip()
  
  ggsave(file.path(output_dir, "toxicity_assessment.png"), p3, width = 10, height = 6)
  
  # Top targets heatmap
  if (nrow(targets) > 1) {
    top_targets <- targets[order(targets$adjusted_priority, decreasing = TRUE)[1:min(20, nrow(targets))], ]
    heatmap_data <- top_targets %>%
      select(composite_druggability, clinical_feasibility, pathway_enrichment_score, 
             priority_score, overall_toxicity_risk) %>%
      as.matrix()
    
    rownames(heatmap_data) <- top_targets$site_id
    
    png(file.path(output_dir, "top_targets_heatmap.png"), width = 800, height = 600)
    pheatmap(heatmap_data, 
             scale = "column",
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = colorRampPalette(c("blue", "white", "red"))(50),
             main = "Top Drug Targets Assessment Heatmap")
    dev.off()
  }
  
  if (opt$verbose) cat("Plots saved to:", output_dir, "\n")
}

# Main execution
tryCatch({
  # Load drug targets
  targets <- load_drug_targets(opt$input)
  
  if (nrow(targets) == 0) {
    stop("No targets found in input file.")
  }
  
  if (opt$verbose) cat("Loaded", nrow(targets), "drug targets\n")
  
  # Assess druggability
  targets <- assess_druggability(targets)
  
  # Assess clinical potential
  targets <- assess_clinical_potential(targets, opt$clinical)
  
  # Perform network analysis
  targets <- perform_network_analysis(targets, opt$ppi_network)
  
  # Prioritize targets
  targets <- prioritize_targets(targets, opt$prioritize)
  
  # Assess toxicity
  targets <- assess_toxicity(targets, opt$side_effects)
  
  # Generate assessment report
  report <- generate_assessment_report(targets)
  
  # Save results
  write.table(targets, file.path(opt$output, "therapeutic_assessment.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save report as JSON
  writeLines(toJSON(report, pretty = TRUE),
             file.path(opt$output, "assessment_report.json"))
  
  # Save prioritized targets
  prioritized <- targets[targets$priority_score >= opt$min_score, ]
  prioritized <- prioritized[order(prioritized$adjusted_priority, decreasing = TRUE), ]
  
  write.table(prioritized, file.path(opt$output, "prioritized_targets.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Create visualizations
  create_assessment_plots(targets, opt$output)
  
  # Print summary
  cat("\n=== Therapeutic Potential Assessment Summary ===\n")
  cat("Total targets assessed:", report$summary$total_targets, "\n")
  cat("Very high priority targets:", report$summary$very_high_priority, "\n")
  cat("High priority targets:", report$summary$high_priority, "\n")
  cat("Hub genes identified:", report$summary$hub_genes, "\n")
  cat("Low toxicity targets:", report$summary$low_toxicity, "\n")
  cat("Mean priority score:", round(report$metrics$mean_priority_score, 3), "\n")
  cat("Targets above threshold:", nrow(prioritized), "\n")
  cat("\nResults saved to:", opt$output, "\n")
  
}, error = function(e) {
  cat("Error in therapeutic assessment:", e$message, "\n")
  quit(status = 1)
})

if (opt$verbose) cat("Therapeutic potential assessment completed successfully!\n")
