# Dockerfile for CTCF Binding Site Prediction Pipeline
# This Dockerfile creates an environment with all the necessary dependencies
# for running the CTCF prediction pipeline scripts.

# Use a specific version of r-base to avoid dependency issues
FROM r-base:4.3.2

# Set working directory
WORKDIR /app

# Copy requirements files
COPY requirements.txt r-requirements.txt /app/

# Install system dependencies from requirements.txt
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      wget \
      curl \
      bedtools \
      ca-certificates \
      libxml2-dev \
      libssl-dev \
      # Use libcurl4-gnutls-dev instead of openssl to avoid dependency issues
      libcurl4-gnutls-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
# Read R packages from r-requirements.txt, filter comments and install
RUN R -e "install.packages('BiocManager')" \
    && R -e "BiocManager::install(c('Biostrings'))" \
    && R -e "install.packages(c('pROC', 'jsonlite'))"

# Comment explaining why we're not using a more elegant approach with the requirements file
# We install packages explicitly rather than parsing the requirements file
# because R doesn't have a built-in equivalent to pip install -r requirements.txt

# Copy scripts to the container (optional - can be mounted instead)
# COPY scripts/ /app/scripts/

# Set default command to bash with entrypoint to allow both bash and R scripts
ENTRYPOINT ["/bin/bash"]
CMD ["--help"]
