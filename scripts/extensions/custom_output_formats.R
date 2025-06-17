# Custom Output Formats Extension
# Provides specialized output formats for different research needs

# Source the plugin manager
if (file.exists("scripts/core/plugin_manager.R")) {
  source("scripts/core/plugin_manager.R")
}

# Research-specific format export
export_research_format <- function(pwm, validation_results, config = list()) {
  cat("Exporting in research format...\n")
  
  # Parameters
  include_metadata <- if (is.null(config$include_metadata)) TRUE else config$include_metadata
  include_statistics <- if (is.null(config$include_statistics)) TRUE else config$include_statistics
  output_format <- if (is.null(config$format)) "json" else config$format
  
  # Build research export structure
  research_export <- list(
    # Core PWM data
    pwm = list(
      matrix = pwm,
      dimensions = dim(pwm),
      information_content = calculate_information_content(pwm),
      consensus_sequence = derive_consensus_sequence(pwm)
    ),
    
    # Validation results
    validation = validation_results
  )
  
  # Add metadata if requested
  if (include_metadata) {
    research_export$metadata <- list(
      creation_date = Sys.time(),
      pipeline_version = "2.0.0",
      export_format = "research",
      parameters_used = config
    )
  }
  
  # Add detailed statistics if requested
  if (include_statistics) {
    research_export$statistics <- calculate_detailed_statistics(pwm, validation_results)
  }
  
  # Format output
  if (output_format == "json") {
    # Convert to JSON-friendly format
    research_export <- make_json_compatible(research_export)
    
    # Convert to JSON string
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      output_string <- jsonlite::toJSON(research_export, auto_unbox = TRUE, pretty = TRUE)
    } else {
      # Fallback: create simple text representation
      output_string <- create_text_representation(research_export)
    }
  } else if (output_format == "xml") {
    output_string <- create_xml_representation(research_export)
  } else {
    # Default: structured text format
    output_string <- create_text_representation(research_export)
  }
  
  return(list(
    content = output_string,
    format = output_format,
    export_method = "research_format",
    export_time = Sys.time()
  ))
}

# Database export format (e.g., for JASPAR, TRANSFAC)
export_database_format <- function(pwm, metadata, config = list()) {
  cat("Exporting in database format...\n")
  
  # Parameters
  database_type <- if (is.null(config$database_type)) "jaspar" else config$database_type
  include_accession <- if (is.null(config$include_accession)) TRUE else config$include_accession
  
  if (database_type == "jaspar") {
    database_entry <- create_jaspar_format(pwm, metadata, config)
  } else if (database_type == "transfac") {
    database_entry <- create_transfac_format(pwm, metadata, config)
  } else if (database_type == "meme") {
    database_entry <- create_meme_format(pwm, metadata, config)
  } else {
    stop("Unsupported database type: ", database_type)
  }
  
  return(database_entry)
}

# Helper function: Calculate information content
calculate_information_content <- function(pwm) {
  ic_values <- numeric(ncol(pwm))
  
  for (pos in seq_len(ncol(pwm))) {
    probs <- pwm[, pos]
    probs <- probs[probs > 0]  # Remove zero probabilities
    
    # Information content formula: IC = 2 + sum(p * log2(p))
    ic <- 2 + sum(probs * log2(probs))
    ic_values[pos] <- max(0, ic)  # Ensure non-negative
  }
  
  return(list(
    per_position = ic_values,
    total = sum(ic_values),
    mean = mean(ic_values),
    max = max(ic_values),
    min = min(ic_values)
  ))
}

# Helper function: Derive consensus sequence
derive_consensus_sequence <- function(pwm) {
  consensus <- character(ncol(pwm))
  
  for (pos in seq_len(ncol(pwm))) {
    # Find nucleotide with highest probability
    max_idx <- which.max(pwm[, pos])
    consensus[pos] <- rownames(pwm)[max_idx]
  }
  
  return(paste(consensus, collapse = ""))
}

# Helper function: Calculate detailed statistics
calculate_detailed_statistics <- function(pwm, validation_results) {
  statistics <- list(
    # PWM statistics
    pwm_stats = list(
      nucleotide_frequencies = rowMeans(pwm),
      position_entropies = apply(pwm, 2, function(col) {
        probs <- col[col > 0]
        -sum(probs * log2(probs))
      }),
      gc_content = mean(pwm["G", ] + pwm["C", ]),
      at_content = mean(pwm["A", ] + pwm["T", ])
    ),
    
    # Validation statistics
    validation_stats = list(
      summary = validation_results,
      quality_scores = extract_quality_scores(validation_results),
      confidence_intervals = extract_confidence_intervals(validation_results)
    )
  )
  
  return(statistics)
}

# Helper function: Make data structure JSON-compatible
make_json_compatible <- function(data) {
  if (is.matrix(data)) {
    # Convert matrix to list of lists
    return(lapply(seq_len(nrow(data)), function(i) {
      row_data <- as.list(data[i, ])
      names(row_data) <- colnames(data)
      return(row_data)
    }))
  } else if (is.list(data)) {
    # Recursively process list elements
    return(lapply(data, make_json_compatible))
  } else {
    return(data)
  }
}

# Helper function: Create text representation
create_text_representation <- function(data) {
  output_lines <- c()
  
  output_lines <- c(output_lines, "# Research Format Export")
  output_lines <- c(output_lines, paste("# Generated:", Sys.time()))
  output_lines <- c(output_lines, "")
  
  # PWM section
  if (!is.null(data$pwm)) {
    output_lines <- c(output_lines, "## Position Weight Matrix")
    
    if (!is.null(data$pwm$matrix)) {
      pwm <- data$pwm$matrix
      output_lines <- c(output_lines, "Matrix:")
      for (i in seq_len(nrow(pwm))) {
        row_str <- paste(rownames(pwm)[i], ":", paste(round(pwm[i, ], 4), collapse = " "))
        output_lines <- c(output_lines, row_str)
      }
    }
    
    if (!is.null(data$pwm$consensus_sequence)) {
      output_lines <- c(output_lines, paste("Consensus:", data$pwm$consensus_sequence))
    }
    
    output_lines <- c(output_lines, "")
  }
  
  # Validation section
  if (!is.null(data$validation)) {
    output_lines <- c(output_lines, "## Validation Results")
    output_lines <- c(output_lines, paste("Results:", paste(names(data$validation), collapse = ", ")))
    output_lines <- c(output_lines, "")
  }
  
  return(paste(output_lines, collapse = "\n"))
}

# Helper function: Create XML representation
create_xml_representation <- function(data) {
  xml_lines <- c()
  
  xml_lines <- c(xml_lines, '<?xml version="1.0" encoding="UTF-8"?>')
  xml_lines <- c(xml_lines, '<research_export>')
  
  if (!is.null(data$pwm)) {
    xml_lines <- c(xml_lines, '  <pwm>')
    
    if (!is.null(data$pwm$matrix)) {
      pwm <- data$pwm$matrix
      xml_lines <- c(xml_lines, '    <matrix>')
      
      for (i in seq_len(nrow(pwm))) {
        row_xml <- paste0('      <row nucleotide="', rownames(pwm)[i], '">')
        for (j in seq_len(ncol(pwm))) {
          row_xml <- paste0(row_xml, '<pos', j, '>', round(pwm[i, j], 4), '</pos', j, '>')
        }
        row_xml <- paste0(row_xml, '</row>')
        xml_lines <- c(xml_lines, row_xml)
      }
      
      xml_lines <- c(xml_lines, '    </matrix>')
    }
    
    if (!is.null(data$pwm$consensus_sequence)) {
      xml_lines <- c(xml_lines, paste0('    <consensus>', data$pwm$consensus_sequence, '</consensus>'))
    }
    
    xml_lines <- c(xml_lines, '  </pwm>')
  }
  
  xml_lines <- c(xml_lines, '</research_export>')
  
  return(paste(xml_lines, collapse = "\n"))
}

# Helper function: Create JASPAR format
create_jaspar_format <- function(pwm, metadata, config) {
  jaspar_lines <- c()
  
  # Header
  if (!is.null(metadata$name)) {
    jaspar_lines <- c(jaspar_lines, paste0(">", metadata$name))
  } else {
    jaspar_lines <- c(jaspar_lines, ">CTCF_PWM")
  }
  
  # Matrix in JASPAR format
  for (nucleotide in c("A", "C", "G", "T")) {
    if (nucleotide %in% rownames(pwm)) {
      values <- round(pwm[nucleotide, ] * 1000)  # Convert to counts
      line <- paste0(nucleotide, " [", paste(values, collapse = " "), "]")
      jaspar_lines <- c(jaspar_lines, line)
    }
  }
  
  return(list(
    content = paste(jaspar_lines, collapse = "\n"),
    format = "jaspar",
    export_method = "database_format"
  ))
}

# Helper function: Create TRANSFAC format
create_transfac_format <- function(pwm, metadata, config) {
  transfac_lines <- c()
  
  # Header
  transfac_lines <- c(transfac_lines, "VV  TRANSFAC MATRIX")
  transfac_lines <- c(transfac_lines, "//")
  
  if (!is.null(metadata$accession)) {
    transfac_lines <- c(transfac_lines, paste("AC", metadata$accession))
  }
  
  if (!is.null(metadata$name)) {
    transfac_lines <- c(transfac_lines, paste("ID", metadata$name))
  }
  
  # Matrix
  transfac_lines <- c(transfac_lines, "P0      A      C      G      T")
  
  for (pos in seq_len(ncol(pwm))) {
    values <- round(pwm[, pos] * 1000)  # Convert to counts
    line <- sprintf("%02d    %6d %6d %6d %6d", pos, 
                   values["A"], values["C"], values["G"], values["T"])
    transfac_lines <- c(transfac_lines, line)
  }
  
  transfac_lines <- c(transfac_lines, "//")
  
  return(list(
    content = paste(transfac_lines, collapse = "\n"),
    format = "transfac",
    export_method = "database_format"
  ))
}

# Helper function: Create MEME format
create_meme_format <- function(pwm, metadata, config) {
  meme_lines <- c()
  
  # Header
  meme_lines <- c(meme_lines, "MEME version 4")
  meme_lines <- c(meme_lines, "")
  meme_lines <- c(meme_lines, "ALPHABET= ACGT")
  meme_lines <- c(meme_lines, "")
  meme_lines <- c(meme_lines, "strands: + -")
  meme_lines <- c(meme_lines, "")
  meme_lines <- c(meme_lines, "Background letter frequencies (from uniform background):")
  meme_lines <- c(meme_lines, "A 0.25 C 0.25 G 0.25 T 0.25")
  meme_lines <- c(meme_lines, "")
  
  # Motif
  motif_name <- if (!is.null(metadata$name)) metadata$name else "CTCF_MOTIF"
  meme_lines <- c(meme_lines, paste("MOTIF", motif_name))
  meme_lines <- c(meme_lines, "")
  
  # Matrix
  meme_lines <- c(meme_lines, paste("letter-probability matrix: alength= 4 w=", ncol(pwm), "nsites= 100 E= 1e-10"))
  
  for (pos in seq_len(ncol(pwm))) {
    values <- pwm[, pos]
    line <- sprintf("  %.6f  %.6f  %.6f  %.6f", 
                   values["A"], values["C"], values["G"], values["T"])
    meme_lines <- c(meme_lines, line)
  }
  
  meme_lines <- c(meme_lines, "")
  
  return(list(
    content = paste(meme_lines, collapse = "\n"),
    format = "meme",
    export_method = "database_format"
  ))
}

# Helper functions for validation data extraction
extract_quality_scores <- function(validation_results) {
  if (is.list(validation_results)) {
    quality_scores <- list()
    
    for (name in names(validation_results)) {
      result <- validation_results[[name]]
      if (is.list(result) && !is.null(result$score)) {
        quality_scores[[name]] <- result$score
      }
    }
    
    return(quality_scores)
  }
  
  return(list())
}

extract_confidence_intervals <- function(validation_results) {
  if (is.list(validation_results)) {
    confidence_intervals <- list()
    
    for (name in names(validation_results)) {
      result <- validation_results[[name]]
      if (is.list(result) && !is.null(result$confidence_interval)) {
        confidence_intervals[[name]] <- result$confidence_interval
      }
    }
    
    return(confidence_intervals)
  }
  
  return(list())
}

# Register output formats
register_output_format("research", export_research_format)
register_output_format("database", export_database_format)

cat("Custom output formats loaded and registered.\n")
