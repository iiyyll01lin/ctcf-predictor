#!/bin/bash

# Proxy configuration for corporate networks
PROXY_SERVER="10.6.254.210"
PROXY_PORT="3128"
PROXY_URL="http://$PROXY_SERVER:$PROXY_PORT"

# Set environment variables
export http_proxy=$PROXY_URL
export https_proxy=$PROXY_URL
export HTTP_PROXY=$PROXY_URL
export HTTPS_PROXY=$PROXY_URL
export no_proxy="localhost,127.0.0.1"
export NO_PROXY="localhost,127.0.0.1"

# Configure git to use proxy
git config --global http.proxy $PROXY_URL
git config --global https.proxy $PROXY_URL

# Configure apt to use proxy (if on Debian/Ubuntu)
if command -v apt-get &> /dev/null; then
    echo "Configuring apt to use proxy..."
    cat > /etc/apt/apt.conf.d/proxy.conf << EOL
Acquire::http::Proxy "$PROXY_URL";
Acquire::https::Proxy "$PROXY_URL";
EOL
fi

echo "Proxy environment has been set up with $PROXY_URL"