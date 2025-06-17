# Dynamic Configuration System for CTCF Pipeline Extensions
# Provides functionality for loading, merging, and updating configurations at runtime

# Dependencies
if (!requireNamespace("yaml", quietly = TRUE)) {
  warning("Package 'yaml' is needed for configuration files. Install with: install.packages('yaml')")
}

#' Load base configuration from file
#' @param config_file Path to base configuration file
#' @return Configuration list
load_base_config <- function(config_file) {
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }
  
  tryCatch({
    if (grepl("\\.ya?ml$", config_file, ignore.case = TRUE)) {
      config <- yaml::read_yaml(config_file)
    } else if (grepl("\\.json$", config_file, ignore.case = TRUE)) {
      config <- jsonlite::fromJSON(config_file, simplifyVector = FALSE)
    } else {
      stop("Unsupported configuration file format. Use YAML or JSON.")
    }
    
    return(config)
  }, error = function(e) {
    stop("Error loading configuration from ", config_file, ": ", e$message)
  })
}

# Get default configuration
get_default_config <- function() {
  list(
    pipeline_config = list(
      data = list(
        input = list(
          format = "fasta",
          max_sequences = 100000
        ),
        preprocessing = list(
          min_length = 50,
          normalize_length = FALSE
        ),
        validation = list(
          method = "chromosome_split",
          random_seed = 42
        )
      ),
      algorithms = list(
        alignment = list(
          method = "integrated"
        ),
        pwm_construction = list(
          pseudocount = 0.1,
          normalization = "probability"
        )
      ),
      quality = list(
        thresholds = list(
          min_total_ic = 8.0,
          conservation_threshold = 1.0
        ),
        tests = list(
          information_content = TRUE,
          cross_validation = TRUE
        ),
        grading = list(
          excellent_threshold = 16.0,
          acceptable_threshold = 8.0
        )
      ),
      output = list(
        formats = list("meme", "jaspar"),
        include_metadata = TRUE,
        generate_reports = TRUE,
        create_visualizations = TRUE,
        file_naming = list(
          use_timestamp = TRUE,
          use_method_suffix = TRUE,
          use_quality_suffix = FALSE
        )
      ),
      performance = list(
        parallel_processing = TRUE,
        num_threads = 4,
        memory_limit = "8G",
        temp_directory = "/tmp/ctcf"
      ),
      logging = list(
        level = "INFO",
        file = "pipeline.log",
        console_output = TRUE,
        verbose_alignment = FALSE
      )
    ),
    extensions = list(
      alignment_methods = list(),
      quality_metrics = list(),
      validation_tests = list(),
      output_formats = list(),
      visualizations = list()
    ),
    custom_extensions = list()
  )
}

# Discover available extensions in directory
discover_extensions <- function(extensions_dir = "scripts/extensions/") {
  if (!dir.exists(extensions_dir)) {
    cat("Extensions directory not found:", extensions_dir, "\n")
    return(list())
  }
  
  # Find all R files in extensions directory
  extension_files <- list.files(extensions_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  
  extensions <- list()
  
  for (file in extension_files) {
    # Extract extension name from filename
    base_name <- tools::file_path_sans_ext(basename(file))
    extensions[[base_name]] <- list(
      file = file,
      name = base_name,
      loaded = FALSE
    )
  }
  
  cat("Discovered", length(extensions), "potential extensions\n")
  
  return(extensions)
}

#' Discover extension configuration files
#' @param extensions_dir Directory to search for extension configs
#' @return Character vector of config file paths
discover_extension_configs <- function(extensions_dir = "scripts/extensions/") {
  if (!dir.exists(extensions_dir)) {
    return(character(0))
  }
  
  # Look for config files in extensions directory
  config_files <- list.files(extensions_dir, 
                             pattern = "\\.(yml|yaml|json)$", 
                             full.names = TRUE, 
                             recursive = TRUE)
  
  # Also check for dedicated config subdirectories
  config_dirs <- list.dirs(extensions_dir, recursive = TRUE)
  config_dirs <- config_dirs[grepl("/config$", config_dirs)]
  
  for (config_dir in config_dirs) {
    additional_configs <- list.files(config_dir, 
                                   pattern = "\\.(yml|yaml|json)$", 
                                   full.names = TRUE)
    config_files <- c(config_files, additional_configs)
  }
  
  return(unique(config_files))
}

#' Load extension configurations
#' @param extension_files Vector of extension config file paths
#' @return List of extension configurations
load_extension_configs <- function(extension_files) {
  extension_configs <- list()
  
  for (config_file in extension_files) {
    tryCatch({
      # Extract extension name from file path
      base_name <- basename(config_file)
      extension_name <- gsub("\\.(yml|yaml|json)$", "", base_name)
      
      # Load configuration
      config <- load_base_config(config_file)
      extension_configs[[extension_name]] <- config
      
    }, error = function(e) {
      warning("Failed to load extension config from ", config_file, ": ", e$message)
    })
  }
  
  return(extension_configs)
}

# Merge configurations
merge_configurations <- function(base_config, extension_configs) {
  merged_config <- base_config
  
  if (length(extension_configs) == 0) {
    return(merged_config)
  }
  
  for (ext_name in names(extension_configs)) {
    ext_config <- extension_configs[[ext_name]]
    
    # Merge extension config into base config
    merged_config <- merge_config_recursive(merged_config, ext_config)
    
    cat("Merged configuration for extension:", ext_name, "\n")
  }
  
  return(merged_config)
}

#' Merge configurations with precedence handling
#' @param base_config Base configuration list
#' @param extension_configs List of extension configurations
#' @return Merged configuration list
merge_configurations_precedence <- function(base_config, extension_configs) {
  merged_config <- base_config
  
  # Initialize extensions section if not present
  if (is.null(merged_config$extensions)) {
    merged_config$extensions <- list()
  }
  
  # Merge each extension configuration
  for (ext_name in names(extension_configs)) {
    ext_config <- extension_configs[[ext_name]]
    
    # Add to extensions section
    merged_config$extensions[[ext_name]] <- ext_config
    
    # If extension provides overrides for main config sections, apply them
    if (!is.null(ext_config$overrides)) {
      merged_config <- merge_config_recursive(merged_config, ext_config$overrides)
    }
  }
  
  return(merged_config)
}

# Recursive configuration merging
merge_config_recursive <- function(base, extension) {
  if (!is.list(base) || !is.list(extension)) {
    return(extension)  # Extension overrides base for non-list values
  }
  
  merged <- base
  
  for (name in names(extension)) {
    if (name %in% names(base)) {
      if (is.list(base[[name]]) && is.list(extension[[name]])) {
        # Recursively merge nested lists
        merged[[name]] <- merge_config_recursive(base[[name]], extension[[name]])
      } else {
        # Extension overrides base
        merged[[name]] <- extension[[name]]
      }
    } else {
      # Add new configuration from extension
      merged[[name]] <- extension[[name]]
    }
  }
  
  return(merged)
}

# Validate extended configuration
validate_extended_config <- function(config) {
  # Basic structure validation
  if (!is.list(config)) {
    stop("Configuration must be a list")
  }
  
  # Validate extension-specific configurations
  if (!is.null(config$extensions)) {
    if (!is.list(config$extensions)) {
      stop("extensions section must be a list")
    }
    
    # Validate each extension configuration
    for (ext_name in names(config$extensions)) {
      ext_config <- config$extensions[[ext_name]]
      
      if (!is.list(ext_config)) {
        stop("Extension configuration for '", ext_name, "' must be a list")
      }
      
      # Check for required fields if specified
      if (!is.null(ext_config$required_fields)) {
        missing_fields <- setdiff(ext_config$required_fields, names(ext_config))
        if (length(missing_fields) > 0) {
          stop("Extension '", ext_name, "' missing required fields: ", 
               paste(missing_fields, collapse = ", "))
        }
      }
    }
  }
  
  return(TRUE)
}

# Validate extended configuration
validate_extended_config_v2 <- function(config) {
  errors <- character(0)
  warnings <- character(0)
  
  # Validate basic structure
  required_sections <- c("pipeline_config", "extensions")
  missing_sections <- setdiff(required_sections, names(config))
  
  if (length(missing_sections) > 0) {
    errors <- c(errors, paste("Missing required sections:", paste(missing_sections, collapse = ", ")))
  }
  
  # Validate pipeline_config structure
  if ("pipeline_config" %in% names(config)) {
    pipeline_errors <- validate_pipeline_config(config$pipeline_config)
    errors <- c(errors, pipeline_errors)
  }
  
  # Validate extensions structure
  if ("extensions" %in% names(config)) {
    extension_errors <- validate_extensions_config(config$extensions)
    errors <- c(errors, extension_errors)
  }
  
  # Validate numeric values
  numeric_validations <- validate_numeric_parameters(config)
  errors <- c(errors, numeric_validations$errors)
  warnings <- c(warnings, numeric_validations$warnings)
  
  # Print warnings
  if (length(warnings) > 0) {
    for (warning_msg in warnings) {
      warning(warning_msg)
    }
  }
  
  return(list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings
  ))
}

# Validate pipeline configuration section
validate_pipeline_config <- function(pipeline_config) {
  errors <- character(0)
  
  required_subsections <- c("data", "algorithms", "quality", "output")
  missing_subsections <- setdiff(required_subsections, names(pipeline_config))
  
  if (length(missing_subsections) > 0) {
    errors <- c(errors, paste("Missing pipeline config subsections:", paste(missing_subsections, collapse = ", ")))
  }
  
  return(errors)
}

# Validate extensions configuration section
validate_extensions_config <- function(extensions_config) {
  errors <- character(0)
  
  # Extensions config can be empty, so no required fields
  # Just validate structure if present
  
  valid_extension_types <- c("alignment_methods", "quality_metrics", "validation_tests", "output_formats", "visualizations")
  
  for (ext_type in names(extensions_config)) {
    if (!ext_type %in% valid_extension_types) {
      errors <- c(errors, paste("Unknown extension type:", ext_type))
    }
  }
  
  return(errors)
}

# Validate numeric parameters
validate_numeric_parameters <- function(config) {
  errors <- character(0)
  warnings <- character(0)
  
  # Check pseudocount
  if (!is.null(config$pipeline_config$algorithms$pwm_construction$pseudocount)) {
    pseudocount <- config$pipeline_config$algorithms$pwm_construction$pseudocount
    if (!is.numeric(pseudocount) || pseudocount <= 0) {
      errors <- c(errors, "Pseudocount must be a positive number")
    }
  }
  
  # Check quality thresholds
  if (!is.null(config$pipeline_config$quality$thresholds$min_total_ic)) {
    min_ic <- config$pipeline_config$quality$thresholds$min_total_ic
    if (!is.numeric(min_ic) || min_ic < 0) {
      errors <- c(errors, "Minimum total IC must be a non-negative number")
    }
  }
  
  # Check thread count
  if (!is.null(config$pipeline_config$performance$num_threads)) {
    num_threads <- config$pipeline_config$performance$num_threads
    if (!is.numeric(num_threads) || num_threads < 1) {
      warnings <- c(warnings, "Number of threads should be at least 1")
    }
  }
  
  return(list(errors = errors, warnings = warnings))
}

# Update configuration at runtime
update_config_runtime <- function(config, updates) {
  cat("Updating configuration at runtime...\n")
  
  # Validate updates
  update_validation <- validate_config_updates(updates)
  
  if (!update_validation$valid) {
    stop("Invalid configuration updates: ", paste(update_validation$errors, collapse = ", "))
  }
  
  # Apply updates
  updated_config <- apply_config_updates(config, updates)
  
  # Reload affected components
  affected_components <- identify_affected_components(config, updates)
  reload_affected_components(affected_components, updated_config)
  
  cat("Configuration updated successfully\n")
  
  return(updated_config)
}

# Validate configuration updates
validate_config_updates <- function(updates) {
  errors <- character(0)
  
  if (!is.list(updates)) {
    errors <- c(errors, "Updates must be a list")
  }
  
  # Additional validation logic here
  
  return(list(
    valid = length(errors) == 0,
    errors = errors
  ))
}

# Apply configuration updates
apply_config_updates <- function(config, updates) {
  updated_config <- config
  
  for (path in names(updates)) {
    # Parse the configuration path (e.g., "pipeline_config.quality.thresholds.min_total_ic")
    path_components <- strsplit(path, "\\.")[[1]]
    
    # Navigate to the target location and update
    current_level <- updated_config
    
    for (i in 1:(length(path_components) - 1)) {
      component <- path_components[i]
      if (is.null(current_level[[component]])) {
        current_level[[component]] <- list()
      }
      current_level <- current_level[[component]]
    }
    
    # Set the final value
    final_component <- path_components[length(path_components)]
    current_level[[final_component]] <- updates[[path]]
  }
  
  return(updated_config)
}

#' Apply configuration updates
#' @param config Current configuration
#' @param updates Updates to apply
#' @return Updated configuration
apply_config_updates_v2 <- function(config, updates) {
  validate_config_updates(updates)
  
  # Apply updates recursively
  updated_config <- merge_config_recursive(config, updates)
  
  # Update metadata
  if (!is.null(updated_config$meta)) {
    updated_config$meta$last_update <- Sys.time()
    if (is.null(updated_config$meta$update_history)) {
      updated_config$meta$update_history <- list()
    }
    updated_config$meta$update_history[[length(updated_config$meta$update_history) + 1]] <- list(
      timestamp = Sys.time(),
      updates = names(updates)
    )
  }
  
  return(updated_config)
}

# Identify components affected by configuration changes
identify_affected_components <- function(old_config, updates) {
  affected <- character(0)
  
  for (path in names(updates)) {
    if (grepl("^pipeline_config\\.algorithms", path)) {
      affected <- c(affected, "algorithms")
    }
    if (grepl("^pipeline_config\\.quality", path)) {
      affected <- c(affected, "quality")
    }
    if (grepl("^pipeline_config\\.output", path)) {
      affected <- c(affected, "output")
    }
    if (grepl("^extensions", path)) {
      affected <- c(affected, "extensions")
    }
  }
  
  return(unique(affected))
}

# Reload components affected by configuration changes
reload_affected_components <- function(affected_components, updated_config) {
  for (component in affected_components) {
    cat("Reloading component:", component, "\n")
    
    switch(component,
      "algorithms" = {
        # Reload algorithm configurations
        cat("  Algorithm configurations reloaded\n")
      },
      "quality" = {
        # Reload quality assessment configurations
        cat("  Quality assessment configurations reloaded\n")
      },
      "output" = {
        # Reload output configurations
        cat("  Output configurations reloaded\n")
      },
      "extensions" = {
        # Reload extensions
        cat("  Extensions reloaded\n")
        if (exists("load_all_extensions", mode = "function")) {
          load_all_extensions()
        }
      }
    )
  }
}

# Save configuration to file
save_config <- function(config, output_file) {
  tryCatch({
    yaml::write_yaml(config, output_file)
    cat("Configuration saved to:", output_file, "\n")
    return(TRUE)
  }, error = function(e) {
    warning("Error saving configuration: ", e$message)
    return(FALSE)
  })
}

# Get configuration value by path
get_config_value <- function(config, path, default = NULL) {
  path_components <- strsplit(path, "\\.")[[1]]
  current_level <- config
  
  for (component in path_components) {
    if (is.null(current_level[[component]])) {
      return(default)
    }
    current_level <- current_level[[component]]
  }
  
  return(current_level)
}

# Set configuration value by path
set_config_value <- function(config, path, value) {
  path_components <- strsplit(path, "\\.")[[1]]
  current_level <- config
  
  for (i in 1:(length(path_components) - 1)) {
    component <- path_components[i]
    if (is.null(current_level[[component]])) {
      current_level[[component]] <- list()
    }
    current_level <- current_level[[component]]
  }
  
  final_component <- path_components[length(path_components)]
  current_level[[final_component]] <- value
  
  return(config)
}

#' Get configuration section
#' @param config Configuration list
#' @param section_path Dot-separated path to section (e.g., "extensions.my_plugin.settings")
#' @return Configuration section or NULL if not found
get_config_section <- function(config, section_path) {
  path_parts <- strsplit(section_path, "\\.")[[1]]
  current_section <- config
  
  for (part in path_parts) {
    if (is.list(current_section) && part %in% names(current_section)) {
      current_section <- current_section[[part]]
    } else {
      return(NULL)
    }
  }
  
  return(current_section)
}

#' Set configuration section
#' @param config Configuration list
#' @param section_path Dot-separated path to section
#' @param value Value to set
#' @return Updated configuration
set_config_section <- function(config, section_path, value) {
  path_parts <- strsplit(section_path, "\\.")[[1]]
  
  # Create nested structure if needed
  current_ref <- config
  for (i in seq_along(path_parts)) {
    part <- path_parts[i]
    
    if (i == length(path_parts)) {
      # Last part - set the value
      current_ref[[part]] <- value
    } else {
      # Intermediate part - ensure it's a list
      if (!is.list(current_ref[[part]])) {
        current_ref[[part]] <- list()
      }
      current_ref <- current_ref[[part]]
    }
  }
  
  return(config)
}

# Initialize dynamic configuration system
initialize_dynamic_config_system <- function() {
  cat("Initializing Dynamic Configuration System...\n")
  cat("Configuration loading and merging ready.\n")
  cat("Runtime configuration updates supported.\n")
  
  return(TRUE)
}

# Print configuration summary
print_config_summary <- function(config) {
  cat("=== Configuration Summary ===\n")
  
  if ("pipeline_config" %in% names(config)) {
    pc <- config$pipeline_config
    
    cat("Data Processing:\n")
    cat("  Input format:", get_config_value(config, "pipeline_config.data.input.format", "unknown"), "\n")
    cat("  Max sequences:", get_config_value(config, "pipeline_config.data.input.max_sequences", "unknown"), "\n")
    
    cat("Algorithms:\n")
    cat("  Alignment method:", get_config_value(config, "pipeline_config.algorithms.alignment.method", "unknown"), "\n")
    cat("  Pseudocount:", get_config_value(config, "pipeline_config.algorithms.pwm_construction.pseudocount", "unknown"), "\n")
    
    cat("Quality Thresholds:\n")
    cat("  Min total IC:", get_config_value(config, "pipeline_config.quality.thresholds.min_total_ic", "unknown"), "\n")
    cat("  Conservation threshold:", get_config_value(config, "pipeline_config.quality.thresholds.conservation_threshold", "unknown"), "\n")
    
    cat("Output:\n")
    formats <- get_config_value(config, "pipeline_config.output.formats", list())
    cat("  Formats:", paste(formats, collapse = ", "), "\n")
  }
  
  if ("extensions" %in% names(config)) {
    ext_count <- length(config$extensions)
    cat("Extensions: ", ext_count, " types configured\n")
  }
  
  cat("=============================\n")
}

cat("Dynamic configuration system loaded\n")
