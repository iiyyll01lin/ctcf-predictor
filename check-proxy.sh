#!/bin/bash

# Proxy detection and configuration script with fallback
# This script checks if the proxy server is accessible and configures the environment accordingly

PROXY_HOST="10.6.254.210"
PROXY_PORT="3128"
TIMEOUT=5

# Function to check if proxy is accessible
check_proxy() {
    echo "Checking proxy connectivity to $PROXY_HOST:$PROXY_PORT..."
    
    # Test proxy connection using timeout and /dev/tcp
    if timeout $TIMEOUT bash -c "</dev/tcp/$PROXY_HOST/$PROXY_PORT" 2>/dev/null; then
        echo "✓ Proxy is accessible"
        return 0
    else
        echo "✗ Proxy is not accessible"
        return 1
    fi
}

# Function to set proxy configuration
set_proxy() {
    echo "Setting proxy configuration..."
    export http_proxy="http://$PROXY_HOST:$PROXY_PORT"
    export https_proxy="http://$PROXY_HOST:$PROXY_PORT"
    export HTTP_PROXY="http://$PROXY_HOST:$PROXY_PORT"
    export HTTPS_PROXY="http://$PROXY_HOST:$PROXY_PORT"
    export no_proxy="localhost,127.0.0.1,::1"
    export NO_PROXY="localhost,127.0.0.1,::1"
    
    # Configure git to use proxy
    git config --global http.proxy "http://$PROXY_HOST:$PROXY_PORT" 2>/dev/null || true
    git config --global https.proxy "http://$PROXY_HOST:$PROXY_PORT" 2>/dev/null || true
    
    # Configure wget
    mkdir -p ~/.wget 2>/dev/null || true
    cat > ~/.wgetrc << EOL
http_proxy = http://$PROXY_HOST:$PROXY_PORT
https_proxy = http://$PROXY_HOST:$PROXY_PORT
use_proxy = on
EOL

    # Configure curl
    mkdir -p ~/.curl 2>/dev/null || true
    cat > ~/.curlrc << EOL
proxy = http://$PROXY_HOST:$PROXY_PORT
EOL
}

# Function to unset proxy configuration
unset_proxy() {
    echo "Removing proxy configuration..."
    unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY no_proxy NO_PROXY
    
    # Remove git proxy config
    git config --global --unset http.proxy 2>/dev/null || true
    git config --global --unset https.proxy 2>/dev/null || true
    
    # Remove wget and curl proxy config
    rm -f ~/.wgetrc ~/.curlrc 2>/dev/null || true
}

# Main logic
echo "=== Proxy Detection and Configuration ==="

if check_proxy; then
    set_proxy
    echo "✓ Environment configured with proxy"
    echo "Proxy URL: http://$PROXY_HOST:$PROXY_PORT"
    exit 0
else
    unset_proxy
    echo "✓ Environment configured for direct connection"
    echo "Using direct internet connection (no proxy)"
    exit 1
fi
