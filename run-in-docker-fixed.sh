#!/bin/bash

SCRIPT="$1"
shift

# Simple Docker execution without double bash
docker run --rm -it \
  --network host \
  -v "$(pwd):/app" \
  -w /app \
  ctcf-predictor:latest \
  $SCRIPT $*
