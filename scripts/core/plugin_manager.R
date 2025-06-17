# Plugin Management System for CTCF Pipeline Extensions
# Provides core functionality for plugin registration, loading, and execution

# Global plugin registry
PLUGIN_REGISTRY <- new.env()

# Global extension method registries
ALIGNMENT_METHODS <- new.env()
QUALITY_METRICS <- new.env()
VALIDATION_TESTS <- new.env()
OUTPUT_FORMATS <- new.env()
VISUALIZATION_METHODS <- new.env()

#' Validate plugin interface
#' @param plugin_interface List containing plugin components
#' @return TRUE if valid, stops with error if invalid
validate_plugin_interface <- function(plugin_interface) {
  required_components <- c("info", "initialize", "execute", "cleanup")
  
  if (!is.list(plugin_interface)) {
    stop("Plugin interface must be a list")
  }
  
  missing_components <- setdiff(required_components, names(plugin_interface))
  if (length(missing_components) > 0) {
    stop("Invalid plugin interface. Missing required components: ", 
         paste(missing_components, collapse = ", "))
  }
  
  # Validate info structure
  if (!is.list(plugin_interface$info)) {
    stop("Plugin info must be a list")
  }
  
  required_info_fields <- c("name", "version", "description", "author")
  missing_info <- setdiff(required_info_fields, names(plugin_interface$info))
  if (length(missing_info) > 0) {
    stop("Plugin info missing required fields: ", 
         paste(missing_info, collapse = ", "))
  }
  
  return(TRUE)
}

#' Check version compatibility
#' @param plugin_info Plugin information list
#' @return TRUE if compatible, stops with error if incompatible
check_version_compatibility <- function(plugin_info) {
  if (!is.null(plugin_info$pipeline_version_min)) {
    # For now, we'll assume version compatibility
    # In a real implementation, you'd compare version numbers
    current_version <- "2.0.0"  # This should come from a config file
    min_version <- plugin_info$pipeline_version_min
    
    # Simple version check (would need proper semantic versioning in production)
    if (min_version > current_version) {
      stop("Plugin requires pipeline version ", min_version, 
           " but current version is ", current_version)
    }
  }
  
  return(TRUE)
}

#' Register plugin
#' @param plugin_name Character string with plugin name
#' @param plugin_interface List containing plugin components
#' @return TRUE if successful
register_plugin <- function(plugin_name, plugin_interface) {
  # Validate plugin interface
  validate_plugin_interface(plugin_interface)
  
  # Check version compatibility
  check_version_compatibility(plugin_interface$info)
  
  # Check for name conflicts
  if (exists(plugin_name, envir = PLUGIN_REGISTRY)) {
    warning("Plugin '", plugin_name, "' is already registered. Overwriting.")
  }
  
  # Register plugin
  assign(plugin_name, plugin_interface, envir = PLUGIN_REGISTRY)
  
  cat("Plugin registered:", plugin_name, "v", plugin_interface$info$version, "\n")
  return(TRUE)
}

#' Load plugin from file
#' @param plugin_path Path to plugin R file
#' @return Plugin name if successful
load_plugin <- function(plugin_path) {
  if (!file.exists(plugin_path)) {
    stop("Plugin file not found: ", plugin_path)
  }
  
  # Create isolated environment for plugin loading
  plugin_env <- new.env()
  
  # Source plugin file in isolated environment
  tryCatch({
    source(plugin_path, local = plugin_env)
  }, error = function(e) {
    stop("Error loading plugin from ", plugin_path, ": ", e$message)
  })
  
  # Check if export_plugin_interface function exists
  if (!exists("export_plugin_interface", envir = plugin_env)) {
    stop("Plugin file must contain export_plugin_interface() function")
  }
  
  # Get plugin interface
  export_func <- get("export_plugin_interface", envir = plugin_env)
  plugin_interface <- export_func()
  
  # Register plugin
  plugin_name <- plugin_interface$info$name
  register_plugin(plugin_name, plugin_interface)
  
  return(plugin_name)
}

#' Execute plugin
#' @param plugin_name Name of registered plugin
#' @param input_data Input data for plugin
#' @param config Configuration list
#' @return Plugin execution result
execute_plugin <- function(plugin_name, input_data, config = list()) {
  if (!exists(plugin_name, envir = PLUGIN_REGISTRY)) {
    stop("Plugin not found: ", plugin_name)
  }
  
  plugin <- get(plugin_name, envir = PLUGIN_REGISTRY)
  
  # Initialize plugin if needed
  if (!is.null(plugin$initialize)) {
    plugin$initialize(config)
  }
  
  # Execute plugin
  result <- plugin$execute(input_data, config)
  
  return(result)
}

#' List registered plugins
#' @return Character vector of plugin names
list_plugins <- function() {
  return(ls(envir = PLUGIN_REGISTRY))
}

#' Get plugin information
#' @param plugin_name Name of plugin
#' @return Plugin information list
get_plugin_info <- function(plugin_name) {
  if (!exists(plugin_name, envir = PLUGIN_REGISTRY)) {
    stop("Plugin not found: ", plugin_name)
  }
  
  plugin <- get(plugin_name, envir = PLUGIN_REGISTRY)
  return(plugin$info)
}

#' Unregister plugin
#' @param plugin_name Name of plugin to unregister
#' @return TRUE if successful
unregister_plugin <- function(plugin_name) {
  if (!exists(plugin_name, envir = PLUGIN_REGISTRY)) {
    warning("Plugin not found: ", plugin_name)
    return(FALSE)
  }
  
  plugin <- get(plugin_name, envir = PLUGIN_REGISTRY)
  
  # Run cleanup if available
  if (!is.null(plugin$cleanup)) {
    plugin$cleanup()
  }
  
  # Remove from registry
  rm(list = plugin_name, envir = PLUGIN_REGISTRY)
  
  cat("Plugin unregistered:", plugin_name, "\n")
  return(TRUE)
}

# Extension method registration functions

#' Register alignment method
#' @param method_name Name of the alignment method
#' @param method_function Function implementing the method
#' @return TRUE if successful
register_alignment_method <- function(method_name, method_function) {
  if (!is.function(method_function)) {
    stop("method_function must be a function")
  }
  
  assign(method_name, method_function, envir = ALIGNMENT_METHODS)
  cat("Alignment method registered:", method_name, "\n")
  return(TRUE)
}

#' Register quality metric
#' @param metric_name Name of the quality metric
#' @param metric_function Function implementing the metric
#' @return TRUE if successful
register_quality_metric <- function(metric_name, metric_function) {
  if (!is.function(metric_function)) {
    stop("metric_function must be a function")
  }
  
  assign(metric_name, metric_function, envir = QUALITY_METRICS)
  cat("Quality metric registered:", metric_name, "\n")
  return(TRUE)
}

#' Register validation test
#' @param test_name Name of the validation test
#' @param test_function Function implementing the test
#' @return TRUE if successful
register_validation_test <- function(test_name, test_function) {
  if (!is.function(test_function)) {
    stop("test_function must be a function")
  }
  
  assign(test_name, test_function, envir = VALIDATION_TESTS)
  cat("Validation test registered:", test_name, "\n")
  return(TRUE)
}

#' Register output format
#' @param format_name Name of the output format
#' @param format_function Function implementing the format
#' @return TRUE if successful
register_output_format <- function(format_name, format_function) {
  if (!is.function(format_function)) {
    stop("format_function must be a function")
  }
  
  assign(format_name, format_function, envir = OUTPUT_FORMATS)
  cat("Output format registered:", format_name, "\n")
  return(TRUE)
}

#' Register visualization method
#' @param viz_name Name of the visualization method
#' @param viz_function Function implementing the visualization
#' @return TRUE if successful
register_visualization_method <- function(viz_name, viz_function) {
  if (!is.function(viz_function)) {
    stop("viz_function must be a function")
  }
  
  assign(viz_name, viz_function, envir = VISUALIZATION_METHODS)
  cat("Visualization method registered:", viz_name, "\n")
  return(TRUE)
}

# Helper functions for accessing registered methods

#' Get registered alignment methods
#' @return Character vector of method names
get_alignment_methods <- function() {
  return(ls(envir = ALIGNMENT_METHODS))
}

#' Get registered quality metrics
#' @return Character vector of metric names
get_quality_metrics <- function() {
  return(ls(envir = QUALITY_METRICS))
}

#' Get registered validation tests
#' @return Character vector of test names
get_validation_tests <- function() {
  return(ls(envir = VALIDATION_TESTS))
}

#' Get registered output formats
#' @return Character vector of format names
get_output_formats <- function() {
  return(ls(envir = OUTPUT_FORMATS))
}

#' Get registered visualization methods
#' @return Character vector of visualization names
get_visualization_methods <- function() {
  return(ls(envir = VISUALIZATION_METHODS))
}

#' Execute registered alignment method
#' @param method_name Name of alignment method
#' @param sequences Input sequences
#' @param config Configuration list
#' @return Alignment result
execute_alignment_method <- function(method_name, sequences, config = list()) {
  if (!exists(method_name, envir = ALIGNMENT_METHODS)) {
    stop("Alignment method not found: ", method_name)
  }
  
  method_func <- get(method_name, envir = ALIGNMENT_METHODS)
  return(method_func(sequences, config))
}

#' Execute registered quality metric
#' @param metric_name Name of quality metric
#' @param pwm PWM matrix
#' @param sequences Input sequences
#' @param config Configuration list
#' @return Quality metric result
execute_quality_metric <- function(metric_name, pwm, sequences, config = list()) {
  if (!exists(metric_name, envir = QUALITY_METRICS)) {
    stop("Quality metric not found: ", metric_name)
  }
  
  metric_func <- get(metric_name, envir = QUALITY_METRICS)
  return(metric_func(pwm, sequences, config))
}

#' Execute registered validation test
#' @param test_name Name of validation test
#' @param pwm PWM matrix
#' @param sequences Input sequences
#' @param config Configuration list
#' @return Validation test result
execute_validation_test <- function(test_name, pwm, sequences, config = list()) {
  if (!exists(test_name, envir = VALIDATION_TESTS)) {
    stop("Validation test not found: ", test_name)
  }
  
  test_func <- get(test_name, envir = VALIDATION_TESTS)
  return(test_func(pwm, sequences, config))
}

#' Execute registered output format
#' @param format_name Name of output format
#' @param pwm PWM matrix
#' @param validation_results Validation results
#' @param config Configuration list
#' @return Output format result
execute_output_format <- function(format_name, pwm, validation_results, config = list()) {
  if (!exists(format_name, envir = OUTPUT_FORMATS)) {
    stop("Output format not found: ", format_name)
  }
  
  format_func <- get(format_name, envir = OUTPUT_FORMATS)
  return(format_func(pwm, validation_results, config))
}

#' Execute registered visualization method
#' @param viz_name Name of visualization method
#' @param pwm PWM matrix
#' @param validation_results Validation results
#' @param config Configuration list
#' @return Visualization result
execute_visualization_method <- function(viz_name, pwm, validation_results, config = list()) {
  if (!exists(viz_name, envir = VISUALIZATION_METHODS)) {
    stop("Visualization method not found: ", viz_name)
  }
  
  viz_func <- get(viz_name, envir = VISUALIZATION_METHODS)
  return(viz_func(pwm, validation_results, config))
}

# Auto-discovery and loading functions

#' Discover available extensions in directory
#' @param extensions_dir Directory to search for extensions
#' @return Character vector of extension file paths
discover_extensions <- function(extensions_dir = "scripts/extensions/") {
  if (!dir.exists(extensions_dir)) {
    return(character(0))
  }
  
  extension_files <- list.files(extensions_dir, pattern = "\\.R$", 
                               full.names = TRUE, recursive = TRUE)
  return(extension_files)
}

#' Load all extensions from directory
#' @param extensions_dir Directory containing extensions
#' @return Character vector of loaded plugin names
load_all_extensions <- function(extensions_dir = "scripts/extensions/") {
  extension_files <- discover_extensions(extensions_dir)
  loaded_plugins <- character(0)
  
  for (file in extension_files) {
    tryCatch({
      plugin_name <- load_plugin(file)
      loaded_plugins <- c(loaded_plugins, plugin_name)
    }, error = function(e) {
      warning("Failed to load extension from ", file, ": ", e$message)
    })
  }
  
  return(loaded_plugins)
}

# Initialize extension system
initialize_extension_system <- function() {
  cat("Initializing CTCF Pipeline Extension System...\n")
  cat("Plugin registry ready.\n")
  cat("Extension method registries ready.\n")
  
  # Auto-load extensions if directory exists
  if (dir.exists("scripts/extensions/")) {
    loaded_plugins <- load_all_extensions()
    if (length(loaded_plugins) > 0) {
      cat("Auto-loaded plugins:", paste(loaded_plugins, collapse = ", "), "\n")
    }
  }
  
  return(TRUE)
}
