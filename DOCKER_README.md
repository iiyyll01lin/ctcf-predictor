# Docker Environment for CTCF Binding Site Prediction Pipeline

This document explains how to use Docker to run the CTCF Binding Site Prediction Pipeline without installing dependencies directly on your system.

## Prerequisites

- [Docker](https://www.docker.com/get-started) installed on your system
- Docker Compose (usually included with Docker Desktop)

## Building the Docker Image

To build the Docker image containing all necessary dependencies:

```bash
# Navigate to the project directory
cd /path/to/ctcf-predictor

# Build the Docker image
docker build -t ctcf-predictor .
```

This will create a Docker image with all the required dependencies:
- R with BiocManager, Biostrings, pROC, and jsonlite packages
- bedtools for sequence extraction
- wget and curl for downloading data

## Proxy Configuration

If you're working in a corporate environment that requires a proxy server, the project includes built-in proxy configuration:

```bash
# Configure proxy settings for your environment
# Edit proxy-setup.sh to match your proxy server details
./proxy-setup.sh
```

The Docker environment is preconfigured to use the proxy settings in:

- `docker-compose.yml` - For Docker Compose runs
- `run-in-docker.sh` - For the helper script execution
- `set-proxy.sh` - Applied within the container during runtime

You can modify the proxy settings in these files to match your network requirements.

## Requirements Files

The project includes two requirements files:

1. `r-requirements.txt` - Lists all required R packages for the pipeline
2. `requirements.txt` - Lists system dependencies (bedtools, wget, curl)

These files are used when building the Docker image, ensuring all necessary dependencies are installed.

### Manual Installation (Without Docker)

If you prefer not to use Docker, you can manually install the dependencies using the provided script:

```bash
# Make the script executable
chmod +x install_dependencies.sh

# Run the installation script with sudo
sudo ./install_dependencies.sh
```

Alternatively, you can install dependencies manually:

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y $(cat requirements.txt | grep -v "^#")

# Install R packages
R -e "install.packages('BiocManager')"
R -e "BiocManager::install('Biostrings')"
R -e "install.packages(c('pROC', 'jsonlite'))"
```

## Running the Pipeline in Docker

There are three ways to run the pipeline inside Docker:

### 1. Using the run-in-docker.sh Helper Script (Recommended)

The simplest approach is to use the provided helper script:

```bash
# Make the script executable
chmod +x run-in-docker.sh

# Run a bash script
./run-in-docker.sh scripts/download_data.sh -d

# Run an R script
./run-in-docker.sh scripts/build_pwm.R

# Or with Rscript explicitly
./run-in-docker.sh Rscript scripts/optimize_threshold.R data/test_sequences.fasta results/generated_pwm.rds
```

The helper script includes these features:

- Automatically builds the Docker image if it doesn't exist
- Mounts your current directory to the container
- Configures proxy settings if needed
- Handles both R and Bash scripts appropriately

### 2. Using Docker Compose

This approach gives you an interactive shell within the container:

```bash
# Start the container with an interactive shell
docker-compose run --rm ctcf-predictor

# Now you're inside the container and can run commands directly
./scripts/download_data.sh -d
Rscript scripts/build_pwm.R

# Exit the container when done
exit
```

### 3. Using Docker Run Directly

You can also run commands directly with Docker:

```bash
# For bash scripts:
docker run --rm -it -v "$(pwd):/app" ctcf-predictor:latest /bin/bash -c "./scripts/download_data.sh -d"

# For R scripts:
docker run --rm -it -v "$(pwd):/app" ctcf-predictor:latest Rscript scripts/build_pwm.R
```

## Example Workflow

Here's a complete workflow example using the Docker environment:

```bash
# Build the Docker image
docker build -t ctcf-predictor .

# Make the helper script executable
chmod +x run-in-docker.sh

# Download data in demo mode (chr21 only)
./run-in-docker.sh scripts/download_data.sh -d

# Alternative: Download full genome data for real analysis
# ./run-in-docker.sh scripts/download_data.sh

# Alternative: Force re-download of data
# ./run-in-docker.sh scripts/download_data.sh -f

# Preprocess the sequences
./run-in-docker.sh Rscript scripts/preprocess_sequences.R data/extracted_sequences.fasta data/preprocessed_sequences.fasta

# Prepare datasets
./run-in-docker.sh Rscript scripts/prepare_datasets.R

# Build PWM model
./run-in-docker.sh Rscript scripts/build_pwm.R

# Evaluate models
./run-in-docker.sh Rscript scripts/evaluate_models.R

# Optimize threshold
./run-in-docker.sh Rscript scripts/optimize_threshold.R data/test_sequences.fasta results/generated_pwm.rds

# Run prediction
./run-in-docker.sh Rscript scripts/predict_ctcf.R data/extracted_sequences.fasta results/predictions.tsv results/generated_pwm.rds 5.0
```

## Data Persistence

When using Docker, all your data and results are saved to your local file system, not inside the container. The container mounts your current directory, so any files created or modified are saved to your local machine.

## Disk Space Requirements

Be aware of the following file sizes when using this pipeline:

- Full hg38 reference genome: ~3.1 GB
- Chr21 only (demo mode): ~46 MB
- ENCODE CTCF peaks file: ~2.7 MB
- Extracted sequences: ~8.8 MB

Ensure you have at least 4 GB of free space for the full analysis mode, or 100 MB for the demo mode.

## Troubleshooting

If you encounter permission issues with file access:

```bash
# Fix permissions on scripts
chmod -R +x scripts/

# If needed, fix ownership of files created by Docker
sudo chown -R $(id -u):$(id -g) data/ results/
```

If you need additional R packages:

1. Add them to the Dockerfile
2. Rebuild the image: `docker build -t ctcf-predictor .`

## Network Issues

If you're having trouble with network connectivity inside the container:

1. Check your proxy settings in the `set-proxy.sh` and `docker-compose.yml` files
2. Ensure your firewall allows Docker to access the internet
3. Try adding `--network host` to your Docker run command for troubleshooting
