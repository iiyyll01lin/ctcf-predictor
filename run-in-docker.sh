#!/bin/bash

# This script helps run the CTCF Predictor pipeline scripts within the Docker container
# It automatically mounts the current directory to the container and runs the specified script

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

# Run the container with the current directory mounted
# Check if the image exists first
if ! docker image inspect ctcf-predictor:latest &>/dev/null; then
  echo "Docker image 'ctcf-predictor:latest' not found."
  echo "Building the image now..."
  docker build -t ctcf-predictor .
fi

# Run the container with the current directory mounted
docker run --rm -it \
  --network host \
  -e https_proxy=http://10.6.254.210:3128 \
  -e http_proxy=http://10.6.254.210:3128 \
  -e HTTP_PROXY=http://10.6.254.210:3128 \
  -e HTTPS_PROXY=http://10.6.254.210:3128 \
  -v "$(pwd):/app" \
  -w /app \
  ctcf-predictor:latest \
  -c "source /app/set-proxy.sh && $SCRIPT $*"
