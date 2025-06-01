#!/bin/bash

# Script to install all required dependencies for the CTCF Binding Site Prediction Pipeline
# This script provides an alternative to using Docker

# Check if running as root or with sudo
if [ "$EUID" -ne 0 ]; then
  echo "Please run as root or with sudo"
  exit 1
fi

echo "Installing system dependencies..."

# Detect OS
if [ -f /etc/debian_version ]; then
  # Debian/Ubuntu
  apt-get update
  # Filter comments and empty lines from requirements.txt
  DEPS=$(grep -v "^#" requirements.txt | grep -v "^$")
  apt-get install -y $DEPS
elif [ -f /etc/redhat-release ]; then
  # RHEL/CentOS/Fedora
  yum update -y
  # Note: package names might differ between distros
  yum install -y wget curl bedtools openssl-devel libcurl-devel libxml2-devel
else
  echo "Unsupported OS. Please install dependencies manually according to requirements.txt"
  exit 1
fi

echo "Installing R dependencies..."

# Check if R is installed
if ! command -v R &> /dev/null; then
  echo "R is not installed. Please install R first."
  exit 1
fi

# Install R packages
R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
R -e "BiocManager::install('Biostrings')"
R -e "install.packages(c('pROC', 'jsonlite'))"

echo "All dependencies installed successfully!"
echo "You can now run the CTCF Binding Site Prediction Pipeline scripts directly."
