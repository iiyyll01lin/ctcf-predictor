# Proxy Fallback System for CTCF Predictor

This document describes the automatic proxy detection and fallback system implemented in the CTCF Predictor project.

## Overview

The fallback system automatically detects whether the proxy server (10.6.254.210:3128) is accessible and configures the environment accordingly. If the proxy is not available, the system seamlessly switches to direct internet connection.

## Components

### 1. Core Scripts

#### `check-proxy.sh`
Main proxy detection script that:
- Tests connectivity to the proxy server
- Configures or removes proxy settings based on availability
- Sets environment variables for various tools (git, wget, curl)
- Returns exit code 0 if proxy is available, 1 if not

#### `smart-startup.sh`
Intelligent startup script that:
- Detects proxy availability
- Chooses appropriate Docker Compose configuration
- Starts the application with the correct network settings

#### `proxy_detector.py`
Python module for proxy detection that provides:
- Proxy connectivity testing
- Environment configuration for Python applications
- Integration with requests library
- pip proxy configuration

### 2. Docker Configurations

#### `docker-compose.yml` (Modified)
Enhanced version that uses environment variables for proxy settings:
```yaml
environment:
  - HTTP_PROXY=${HTTP_PROXY:-}
  - HTTPS_PROXY=${HTTPS_PROXY:-}
```

#### `docker-compose-fallback.yml`
Fallback configuration for direct connection with empty proxy variables.

### 3. Enhanced Scripts

#### `run-in-docker.sh` (Modified)
Now includes automatic proxy detection and sets container environment accordingly.

#### `install_dependencies.sh` (Modified)
Enhanced to detect proxy availability and configure package managers appropriately.

#### `download-with-fallback.sh`
Wrapper for the data download script that handles proxy fallback automatically.

## Usage

### Basic Usage

1. **Start the application:**
   ```bash
   ./smart-startup.sh
   ```

2. **Run scripts in Docker with automatic proxy detection:**
   ```bash
   ./run-in-docker.sh scripts/download_data.sh -d
   ./run-in-docker.sh Rscript scripts/build_pwm.R
   ```

3. **Download data with fallback:**
   ```bash
   ./download-with-fallback.sh -d
   ```

### Manual Proxy Detection

1. **Check proxy status:**
   ```bash
   ./check-proxy.sh
   ```

2. **Python proxy detection:**
   ```bash
   python3 proxy_detector.py --status
   python3 proxy_detector.py --configure
   ```

### Environment Configuration

The system automatically sets or unsets these environment variables:
- `HTTP_PROXY` / `http_proxy`
- `HTTPS_PROXY` / `https_proxy`
- `NO_PROXY` / `no_proxy`

It also configures:
- Git proxy settings
- wget proxy settings (`.wgetrc`)
- curl proxy settings (`.curlrc`)

## How It Works

### Proxy Detection Logic

1. **Connectivity Test:** Uses TCP socket connection to test proxy server accessibility
2. **Timeout Handling:** 5-second timeout to avoid hanging on unavailable proxies
3. **Environment Setup:** Automatically configures all necessary proxy settings
4. **Fallback Mode:** Removes proxy configuration when proxy is unavailable

### Decision Flow

```
Start Application
       ↓
Check Proxy (10.6.254.210:3128)
       ↓
   Available?
    ↙     ↘
  Yes      No
   ↓        ↓
Set Proxy  Direct
Settings  Connection
   ↓        ↓
Use docker-compose.yml  Use docker-compose-fallback.yml
   ↓        ↓
Start Application with appropriate configuration
```

## VS Code Integration

The system is designed to work with VS Code's proxy settings. You can set:

```json
{
  "http.proxy": "",
  "http.proxySupport": "fallback"
}
```

This allows VS Code to automatically handle proxy detection.

## Troubleshooting

### Common Issues

1. **Proxy detection fails but proxy is available:**
   - Check if port 3128 is blocked by firewall
   - Verify proxy server address (10.6.254.210:3128)
   - Test manually: `telnet 10.6.254.210 3128`

2. **Docker build fails with proxy:**
   - Ensure proxy environment variables are set before building
   - Use `docker build --build-arg HTTP_PROXY=http://10.6.254.210:3128 .`

3. **Downloads fail in direct mode:**
   - Check internet connectivity
   - Verify DNS resolution
   - Test with: `curl -I https://www.google.com`

### Manual Override

To force proxy usage even if detection fails:
```bash
export HTTP_PROXY=http://10.6.254.210:3128
export HTTPS_PROXY=http://10.6.254.210:3128
```

To force direct connection even if proxy is available:
```bash
unset HTTP_PROXY HTTPS_PROXY http_proxy https_proxy
```

## Benefits

1. **Automatic Detection:** No manual configuration needed
2. **Seamless Fallback:** Works in both corporate and home environments
3. **Comprehensive Coverage:** Handles Docker, system tools, and applications
4. **Error Recovery:** Graceful handling of proxy failures
5. **Development Friendly:** Easy testing and debugging

## Testing

To test the fallback system:

1. **Test with proxy available:**
   ```bash
   ./check-proxy.sh && echo "Proxy mode works"
   ```

2. **Test with proxy blocked (simulate):**
   ```bash
   # Block proxy temporarily (requires admin)
   sudo iptables -A OUTPUT -p tcp --dport 3128 -j DROP
   ./check-proxy.sh || echo "Fallback mode works"
   sudo iptables -D OUTPUT -p tcp --dport 3128 -j DROP
   ```

3. **Test Docker configurations:**
   ```bash
   # Test proxy mode
   docker-compose -f docker-compose.yml config
   
   # Test fallback mode
   docker-compose -f docker-compose-fallback.yml config
   ```

This fallback system ensures that the CTCF Predictor works reliably in any network environment, automatically adapting to available connectivity options.
