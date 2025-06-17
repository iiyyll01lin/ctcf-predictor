# CTCF Pipeline Extension Installation Framework
# Provides functionality for installing, managing, and removing extensions

# Source required dependencies
if (file.exists("scripts/core/plugin_manager.R")) {
  source("scripts/core/plugin_manager.R")
}

if (file.exists("scripts/core/dynamic_config.R")) {
  source("scripts/core/dynamic_config.R")
}

# Extension installation functions

#' Install CTCF pipeline extension
#' @param extension_path Path to extension directory or archive
#' @param target_dir Target directory for installation
#' @param force Force installation (overwrite existing)
#' @return TRUE if successful
install_ctcf_extension <- function(extension_path, target_dir = "scripts/extensions/", force = FALSE) {
  cat("Installing CTCF pipeline extension...\n")
  
  # Validate inputs
  if (!file.exists(extension_path) && !dir.exists(extension_path)) {
    stop("Extension path not found: ", extension_path)
  }
  
  # Create target directory if it doesn't exist
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
    cat("Created target directory:", target_dir, "\n")
  }
  
  # Determine installation type
  if (dir.exists(extension_path)) {
    # Directory installation
    install_from_directory(extension_path, target_dir, force)
  } else if (grepl("\\.(zip|tar\\.gz|tgz)$", extension_path, ignore.case = TRUE)) {
    # Archive installation
    install_from_archive(extension_path, target_dir, force)
  } else {
    stop("Unsupported extension format. Use directory or archive (zip, tar.gz)")
  }
  
  cat("Extension installed successfully!\n")
  return(TRUE)
}

#' Install extension from directory
#' @param source_dir Source directory
#' @param target_dir Target directory
#' @param force Force installation
#' @return TRUE if successful
install_from_directory <- function(source_dir, target_dir, force) {
  # Get extension name from directory
  extension_name <- basename(source_dir)
  target_path <- file.path(target_dir, extension_name)
  
  # Check if extension already exists
  if (dir.exists(target_path) && !force) {
    stop("Extension already exists: ", target_path, ". Use force=TRUE to overwrite.")
  }
  
  # Validate extension structure
  validate_extension_structure(source_dir)
  
  # Copy extension files
  if (dir.exists(target_path)) {
    unlink(target_path, recursive = TRUE)
  }
  
  file.copy(source_dir, target_dir, recursive = TRUE)
  
  # Install dependencies if specified
  install_extension_dependencies(target_path)
  
  # Validate installation
  validate_extension_installation(target_path)
  
  cat("Extension installed from directory:", source_dir, "\n")
  return(TRUE)
}

#' Install extension from archive
#' @param archive_path Path to archive file
#' @param target_dir Target directory
#' @param force Force installation
#' @return TRUE if successful
install_from_archive <- function(archive_path, target_dir, force) {
  # Create temporary extraction directory
  temp_dir <- tempfile("extension_install_")
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  
  dir.create(temp_dir)
  
  # Extract archive
  if (grepl("\\.zip$", archive_path, ignore.case = TRUE)) {
    utils::unzip(archive_path, exdir = temp_dir)
  } else if (grepl("\\.(tar\\.gz|tgz)$", archive_path, ignore.case = TRUE)) {
    utils::untar(archive_path, exdir = temp_dir)
  }
  
  # Find the extension directory (should be the only top-level directory)
  extracted_contents <- list.dirs(temp_dir, full.names = TRUE, recursive = FALSE)
  
  if (length(extracted_contents) != 1) {
    stop("Archive should contain exactly one top-level directory")
  }
  
  extension_dir <- extracted_contents[1]
  
  # Install from extracted directory
  install_from_directory(extension_dir, target_dir, force)
  
  cat("Extension installed from archive:", archive_path, "\n")
  return(TRUE)
}

#' Validate extension structure
#' @param extension_dir Extension directory path
#' @return TRUE if valid, stops with error if invalid
validate_extension_structure <- function(extension_dir) {
  # Check for required files
  required_files <- c("DESCRIPTION", "R/")
  
  for (file in required_files) {
    file_path <- file.path(extension_dir, file)
    if (!file.exists(file_path) && !dir.exists(file_path)) {
      stop("Extension missing required file/directory: ", file)
    }
  }
  
  # Validate DESCRIPTION file
  desc_file <- file.path(extension_dir, "DESCRIPTION")
  if (file.exists(desc_file)) {
    validate_description_file(desc_file)
  }
  
  return(TRUE)
}

#' Validate DESCRIPTION file
#' @param desc_file Path to DESCRIPTION file
#' @return TRUE if valid, stops with error if invalid
validate_description_file <- function(desc_file) {
  # Read DESCRIPTION file
  desc_content <- readLines(desc_file)
  
  # Parse key-value pairs
  desc_data <- list()
  for (line in desc_content) {
    if (grepl(":", line)) {
      parts <- strsplit(line, ":", fixed = TRUE)[[1]]
      if (length(parts) >= 2) {
        key <- trimws(parts[1])
        value <- trimws(paste(parts[-1], collapse = ":"))
        desc_data[[key]] <- value
      }
    }
  }
  
  # Check required fields
  required_fields <- c("Package", "Version", "Description", "Author")
  missing_fields <- setdiff(required_fields, names(desc_data))
  
  if (length(missing_fields) > 0) {
    stop("DESCRIPTION file missing required fields: ", paste(missing_fields, collapse = ", "))
  }
  
  cat("Extension metadata validated:\n")
  cat("  Package:", desc_data$Package, "\n")
  cat("  Version:", desc_data$Version, "\n")
  cat("  Author:", desc_data$Author, "\n")
  
  return(TRUE)
}

#' Install extension dependencies
#' @param extension_dir Extension directory
#' @return TRUE if successful
install_extension_dependencies <- function(extension_dir) {
  desc_file <- file.path(extension_dir, "DESCRIPTION")
  
  if (!file.exists(desc_file)) {
    return(TRUE)  # No dependencies to install
  }
  
  # Read DESCRIPTION file
  desc_content <- readLines(desc_file)
  
  # Find Dependencies line
  deps_line <- grep("^Dependencies:", desc_content, value = TRUE)
  
  if (length(deps_line) > 0) {
    # Parse dependencies
    deps_str <- sub("^Dependencies:\\s*", "", deps_line[1])
    dependencies <- trimws(strsplit(deps_str, ",")[[1]])
    
    if (length(dependencies) > 0 && dependencies[1] != "") {
      cat("Installing dependencies:", paste(dependencies, collapse = ", "), "\n")
      
      for (dep in dependencies) {
        if (!requireNamespace(dep, quietly = TRUE)) {
          cat("Installing package:", dep, "\n")
          tryCatch({
            utils::install.packages(dep, repos = "https://cran.r-project.org")
          }, error = function(e) {
            warning("Failed to install dependency: ", dep, ". Error: ", e$message)
          })
        } else {
          cat("Dependency already installed:", dep, "\n")
        }
      }
    }
  }
  
  return(TRUE)
}

#' Validate extension installation
#' @param extension_dir Extension directory
#' @return TRUE if valid installation
validate_extension_installation <- function(extension_dir) {
  # Check if R files can be sourced without errors
  r_dir <- file.path(extension_dir, "R")
  
  if (dir.exists(r_dir)) {
    r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
    
    for (r_file in r_files) {
      tryCatch({
        source(r_file, local = new.env())
        cat("Validated R file:", basename(r_file), "\n")
      }, error = function(e) {
        warning("Error in R file ", basename(r_file), ": ", e$message)
      })
    }
  }
  
  return(TRUE)
}

#' Uninstall extension
#' @param extension_name Name of extension to uninstall
#' @param extensions_dir Extensions directory
#' @return TRUE if successful
uninstall_ctcf_extension <- function(extension_name, extensions_dir = "scripts/extensions/") {
  cat("Uninstalling extension:", extension_name, "\n")
  
  extension_path <- file.path(extensions_dir, extension_name)
  
  if (!dir.exists(extension_path)) {
    stop("Extension not found: ", extension_path)
  }
  
  # Unregister extension components if loaded
  if (exists("unregister_plugin")) {
    tryCatch({
      unregister_plugin(extension_name)
    }, error = function(e) {
      # Extension might not be loaded
    })
  }
  
  # Remove extension directory
  unlink(extension_path, recursive = TRUE)
  
  cat("Extension uninstalled successfully!\n")
  return(TRUE)
}

#' List installed extensions
#' @param extensions_dir Extensions directory
#' @return Character vector of extension names
list_installed_extensions <- function(extensions_dir = "scripts/extensions/") {
  if (!dir.exists(extensions_dir)) {
    return(character(0))
  }
  
  # Get all subdirectories
  extensions <- list.dirs(extensions_dir, full.names = FALSE, recursive = FALSE)
  
  # Filter out hidden directories and files
  extensions <- extensions[!grepl("^\\.", extensions)]
  
  return(extensions)
}

#' Get extension information
#' @param extension_name Extension name
#' @param extensions_dir Extensions directory
#' @return List with extension information
get_extension_info <- function(extension_name, extensions_dir = "scripts/extensions/") {
  extension_path <- file.path(extensions_dir, extension_name)
  
  if (!dir.exists(extension_path)) {
    stop("Extension not found: ", extension_path)
  }
  
  # Read DESCRIPTION file if it exists
  desc_file <- file.path(extension_path, "DESCRIPTION")
  info <- list(
    name = extension_name,
    path = extension_path,
    installed = TRUE
  )
  
  if (file.exists(desc_file)) {
    desc_content <- readLines(desc_file)
    
    # Parse description
    for (line in desc_content) {
      if (grepl(":", line)) {
        parts <- strsplit(line, ":", fixed = TRUE)[[1]]
        if (length(parts) >= 2) {
          key <- trimws(parts[1])
          value <- trimws(paste(parts[-1], collapse = ":"))
          info[[tolower(key)]] <- value
        }
      }
    }
  }
  
  # Check if extension is loaded
  if (exists("list_plugins")) {
    loaded_plugins <- list_plugins()
    info$loaded <- extension_name %in% loaded_plugins
  } else {
    info$loaded <- FALSE
  }
  
  return(info)
}

#' Update extension
#' @param extension_name Extension name
#' @param extension_path New extension path
#' @param extensions_dir Extensions directory
#' @return TRUE if successful
update_ctcf_extension <- function(extension_name, extension_path, extensions_dir = "scripts/extensions/") {
  cat("Updating extension:", extension_name, "\n")
  
  # Backup existing extension
  existing_path <- file.path(extensions_dir, extension_name)
  if (dir.exists(existing_path)) {
    backup_path <- paste0(existing_path, "_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    file.rename(existing_path, backup_path)
    cat("Backed up existing extension to:", backup_path, "\n")
  }
  
  # Install new version
  tryCatch({
    install_ctcf_extension(extension_path, extensions_dir, force = TRUE)
    
    # Remove backup if installation successful
    if (dir.exists(backup_path)) {
      unlink(backup_path, recursive = TRUE)
    }
    
    cat("Extension updated successfully!\n")
    return(TRUE)
    
  }, error = function(e) {
    # Restore backup if installation failed
    if (exists("backup_path") && dir.exists(backup_path)) {
      if (dir.exists(existing_path)) {
        unlink(existing_path, recursive = TRUE)
      }
      file.rename(backup_path, existing_path)
      cat("Installation failed. Restored backup.\n")
    }
    
    stop("Extension update failed: ", e$message)
  })
}

#' Create extension template
#' @param extension_name Name for new extension
#' @param target_dir Directory to create extension in
#' @return TRUE if successful
create_extension_template <- function(extension_name, target_dir = "scripts/extensions/") {
  extension_path <- file.path(target_dir, extension_name)
  
  if (dir.exists(extension_path)) {
    stop("Extension directory already exists: ", extension_path)
  }
  
  # Create extension directory structure
  dir.create(extension_path, recursive = TRUE)
  dir.create(file.path(extension_path, "R"))
  dir.create(file.path(extension_path, "inst", "config"), recursive = TRUE)
  dir.create(file.path(extension_path, "tests", "testthat"), recursive = TRUE)
  dir.create(file.path(extension_path, "data"))
  
  # Create DESCRIPTION file
  desc_content <- paste0(
    "Package: ", extension_name, "\n",
    "Version: 1.0.0\n",
    "Description: Custom extension for CTCF PWM Testing Pipeline\n",
    "Author: Your Name\n",
    "Maintainer: Your Name <your.email@example.com>\n",
    "Dependencies: \n",
    "License: MIT\n"
  )
  
  writeLines(desc_content, file.path(extension_path, "DESCRIPTION"))
  
  # Create basic R file
  r_content <- paste0(
    "# ", extension_name, " Extension\n",
    "# Custom functionality for CTCF PWM Testing Pipeline\n\n",
    "# Your extension code here\n"
  )
  
  writeLines(r_content, file.path(extension_path, "R", paste0(extension_name, ".R")))
  
  # Create README
  readme_content <- paste0(
    "# ", extension_name, "\n\n",
    "Custom extension for CTCF PWM Testing Pipeline.\n\n",
    "## Installation\n\n",
    "```r\n",
    "install_ctcf_extension('", extension_path, "')\n",
    "```\n\n",
    "## Usage\n\n",
    "Describe how to use your extension here.\n"
  )
  
  writeLines(readme_content, file.path(extension_path, "README.md"))
  
  cat("Extension template created:", extension_path, "\n")
  return(TRUE)
}

# Initialize extension installer
cat("CTCF Pipeline Extension Installation Framework loaded.\n")
cat("Available functions:\n")
cat("  - install_ctcf_extension()\n")
cat("  - uninstall_ctcf_extension()\n")
cat("  - list_installed_extensions()\n")
cat("  - get_extension_info()\n")
cat("  - update_ctcf_extension()\n")
cat("  - create_extension_template()\n")
