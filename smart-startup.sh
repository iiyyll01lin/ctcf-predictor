#!/bin/bash

# Smart startup script with automatic proxy detection
# This script detects proxy availability and starts the application accordingly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== CTCF Predictor Smart Startup ==="

# Check proxy and configure environment
if "$SCRIPT_DIR/check-proxy.sh"; then
    echo "Using proxy configuration for Docker environment"
    COMPOSE_FILE="docker-compose.yml"
else
    echo "Using direct connection for Docker environment"
    COMPOSE_FILE="docker-compose-fallback.yml"
fi

# Start the application
echo "Starting CTCF Predictor with $COMPOSE_FILE..."

if [ -f "$COMPOSE_FILE" ]; then
    docker-compose -f "$COMPOSE_FILE" up -d
    echo "✓ Application started successfully"
else
    echo "⚠ Compose file $COMPOSE_FILE not found, using default configuration"
    docker-compose up -d
fi

echo "=== Startup Complete ==="
