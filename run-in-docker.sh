#!/bin/bash

# This script helps run the CTCF Predictor pipeline scripts within the Docker container
# It automatically detects proxy availability and configures the container accordingly

if [ $# -lt 1 ]; then
  echo "Usage: $0 <script_name> [script_args...]"
  echo "Example: $0 download_data.sh -d"
  echo "Example: $0 Rscript scripts/build_pwm.R"
  exit 1
fi

SCRIPT="$1"
shift

if [[ "$SCRIPT" == *.R ]]; then
  # For R scripts, prepend Rscript if not already provided
  if [[ "$SCRIPT" != Rscript* ]]; then
    SCRIPT="Rscript $SCRIPT"
  fi
fi

# Check proxy connectivity and set environment accordingly
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if "$SCRIPT_DIR/check-proxy.sh" > /dev/null 2>&1; then
  echo "Using proxy configuration for container execution"
  PROXY_ENV="-e https_proxy=http://10.6.254.210:3128 -e http_proxy=http://10.6.254.210:3128 -e HTTP_PROXY=http://10.6.254.210:3128 -e HTTPS_PROXY=http://10.6.254.210:3128"
  PROXY_SETUP="source /app/set-proxy.sh &&"
else
  echo "Using direct connection for container execution"
  PROXY_ENV=""
  PROXY_SETUP=""
fi

# Check if the image exists first
if ! docker image inspect ctcf-predictor:latest &>/dev/null; then
  echo "Docker image 'ctcf-predictor:latest' not found."
  echo "Building the image now..."
  docker build -t ctcf-predictor .
fi

# Run the container with the current directory mounted
docker run --rm -it \
  --network host \
  $PROXY_ENV \
  -v "$(pwd):/app" \
  -w /app \
  ctcf-predictor:latest \
  -c "$PROXY_SETUP $SCRIPT $*"
