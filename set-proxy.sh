#!/bin/bash

# Set up proxy environment variables
export http_proxy=http://10.6.254.210:3128
export https_proxy=http://10.6.254.210:3128
export HTTP_PROXY=http://10.6.254.210:3128
export HTTPS_PROXY=http://10.6.254.210:3128

# Configure wget to use proxy
mkdir -p ~/.wget
cat > ~/.wgetrc << EOL
http_proxy = http://10.6.254.210:3128
https_proxy = http://10.6.254.210:3128
use_proxy = on
EOL

# Configure curl to use proxy
mkdir -p ~/.curl
cat > ~/.curlrc << EOL
proxy = http://10.6.254.210:3128
EOL

echo "Proxy settings configured"
