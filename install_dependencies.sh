#!/bin/bash

# Script to install all required dependencies for the CTCF Binding Site Prediction Pipeline
# This script provides an alternative to using Docker and includes proxy fallback support

# Check if running as root or with sudo
if [ "$EUID" -ne 0 ]; then
  echo "Please run as root or with sudo"
  exit 1
fi

echo "=== Installing CTCF Predictor Dependencies ==="

# Check proxy connectivity before installation
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "$SCRIPT_DIR/check-proxy.sh" ]; then
  if "$SCRIPT_DIR/check-proxy.sh"; then
    echo "Using proxy for package installation"
    PROXY_AVAILABLE=true
  else
    echo "Using direct connection for package installation"
    PROXY_AVAILABLE=false
  fi
else
  echo "Proxy detection script not found, assuming direct connection"
  PROXY_AVAILABLE=false
fi

echo "Installing system dependencies..."

# Detect OS
if [ -f /etc/debian_version ]; then
  # Debian/Ubuntu
  if [ "$PROXY_AVAILABLE" == true ]; then
    # Configure apt proxy if proxy is available
    echo 'Acquire::http::Proxy "http://10.6.254.210:3128";' > /etc/apt/apt.conf.d/proxy.conf
    echo 'Acquire::https::Proxy "http://10.6.254.210:3128";' >> /etc/apt/apt.conf.d/proxy.conf
  fi
  
  apt-get update
  # Filter comments and empty lines from requirements.txt
  DEPS=$(grep -v "^#" requirements.txt | grep -v "^$")
  apt-get install -y $DEPS
  
  # Clean up proxy configuration if it was set
  if [ "$PROXY_AVAILABLE" == true ]; then
    rm -f /etc/apt/apt.conf.d/proxy.conf
  fi
  
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

echo "✓ All dependencies installed successfully!"
echo "You can now run the CTCF Binding Site Prediction Pipeline scripts directly."
