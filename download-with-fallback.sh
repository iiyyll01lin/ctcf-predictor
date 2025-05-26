#!/bin/bash

# Enhanced download script wrapper with proxy fallback
# This script wraps the original download_data.sh with automatic proxy detection

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ORIGINAL_SCRIPT="$SCRIPT_DIR/scripts/download_data.sh"

echo "=== CTCF Data Download with Proxy Fallback ==="

# Check if original download script exists
if [ ! -f "$ORIGINAL_SCRIPT" ]; then
    echo "Error: Original download script not found at $ORIGINAL_SCRIPT"
    exit 1
fi

# Check proxy connectivity and configure environment
if "$SCRIPT_DIR/check-proxy.sh"; then
    echo "Using proxy for data download"
    # Source the proxy configuration to set environment variables
    source "$SCRIPT_DIR/check-proxy.sh" > /dev/null 2>&1
else
    echo "Using direct connection for data download"
    # Ensure no proxy variables are set
    unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY
fi

# Pass all arguments to the original download script
echo "Starting data download..."
"$ORIGINAL_SCRIPT" "$@"

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "✓ Data download completed successfully"
else
    echo "✗ Data download failed with exit code $exit_code"
fi

exit $exit_code
