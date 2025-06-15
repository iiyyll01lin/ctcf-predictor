# Docker Setup

> **🐳 Complete Docker environment setup for cross-platform reproducibility**

## Why Docker?

Docker provides a **consistent, isolated environment** that ensures the pipeline works identically across Linux, Windows, and macOS systems. This eliminates "works on my machine" problems and ensures reproducible results.

### Benefits of Docker Mode
- ✅ **Consistent Environment**: All dependencies included
- ✅ **Cross-platform Support**: Linux/Windows/macOS compatibility
- ✅ **Isolated Execution**: No system conflicts
- ✅ **Proxy Auto-detection**: Network transparency
- ✅ **Easy Updates**: Pull latest container versions

### Considerations
- ⚠️ **Slower Startup**: Container initialization overhead
- ⚠️ **Higher Memory Usage**: ~2GB container overhead
- ⚠️ **Docker Dependency**: Requires Docker installation

## Docker Installation

### Linux (Ubuntu/Debian)
```bash
# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker

# Install Docker Compose
sudo curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

# Verify installation
docker --version
docker-compose --version
```

### Windows 10/11
1. **Enable WSL2**:
   ```powershell
   # Run as Administrator
   wsl --install
   wsl --set-default-version 2
   ```

2. **Install Docker Desktop**:
   - Download from [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop)
   - Run installer and restart system
   - Enable WSL2 integration in Docker Desktop settings

3. **Verify Installation**:
   ```powershell
   docker --version
   docker-compose --version
   ```

### macOS
1. **Install Docker Desktop**:
   - Download from [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop)
   - Drag Docker.app to Applications folder
   - Launch Docker Desktop

2. **Verify Installation**:
   ```bash
   docker --version
   docker-compose --version
   ```

## Pipeline Docker Setup

### Automatic Setup (Recommended)
```bash
# Clone repository
git clone https://github.com/organization/ctcf-predictor.git
cd ctcf-predictor

# Smart startup with automatic proxy detection
./smart-startup.sh
```

The `smart-startup.sh` script automatically:
- Detects network proxy settings
- Builds Docker containers if needed
- Configures environment variables
- Starts the pipeline environment

### Manual Docker Build
```bash
# Build the Docker image
docker build -t ctcf-predictor:latest .

# Alternative: Build with Docker Compose
docker-compose build
```

### Verify Docker Environment
```bash
# Test Docker container
docker run --rm ctcf-predictor:latest Rscript --version

# Test R packages
docker run --rm ctcf-predictor:latest Rscript -e "library(Biostrings); cat('✅ Biostrings OK\\n')"

# Test complete environment
./run-in-docker.sh test_chromosome_split.R demo
```

## Container Architecture

### Base Image
- **Base**: `r-base:4.3.0`
- **OS**: Debian 11 (Bullseye)
- **Architecture**: Multi-platform (amd64, arm64)

### Installed Components
```dockerfile
# R packages
Biostrings 2.66.0       # DNA sequence analysis
seqinr 4.2-30          # Sequence I/O
ggplot2 3.4.2          # Visualization
pROC 1.18.0            # ROC analysis
jsonlite 1.8.4         # Configuration handling

# System tools
bedtools 2.30.0        # Genomic intervals
curl 7.68.0            # Data download
wget 1.21              # Alternative download
unzip 6.0              # Archive extraction
```

### Container Specifications
- **Size**: ~1.5GB compressed, ~4GB uncompressed
- **Memory**: 2GB base + analysis memory
- **CPU**: All available cores utilized
- **Network**: Proxy-aware configuration

## Running the Pipeline

### Interactive Mode
```bash
# Start interactive container
docker run -it --rm \
  -v "$(pwd):/app" \
  -w /app \
  ctcf-predictor:latest bash

# Run commands inside container
Rscript scripts/simple_aligned_pwm.R
```

### Single Command Mode
```bash
# Run specific script
./run-in-docker.sh scripts/simple_aligned_pwm.R

# Run with arguments
./run-in-docker.sh scripts/build_subset_pwm.R 1000 2000 5000
```

### Background Processing
```bash
# Start long-running analysis
./run-in-docker.sh test_pwm_improvements_with_null_analysis.sh &

# Monitor progress
docker logs -f $(docker ps -q --filter ancestor=ctcf-predictor:latest)
```

## Configuration and Customization

### Environment Variables
```bash
# Memory limits
export CTCF_MEMORY_LIMIT="8G"
export CTCF_SWAP_LIMIT="12G"

# Processing options  
export CTCF_THREADS="4"
export CTCF_BATCH_SIZE="1000"

# Run with custom settings
./run-in-docker.sh scripts/efficient_aligned_pwm.R
```

### Volume Mounting
```bash
# Mount custom data directory
docker run --rm -it \
  -v "/path/to/data:/app/data" \
  -v "/path/to/results:/app/results" \
  ctcf-predictor:latest bash
```

### Custom Docker Compose
```yaml
# docker-compose.custom.yml
version: '3.8'
services:
  ctcf-pipeline:
    build: .
    volumes:
      - ./data:/app/data
      - ./results:/app/results
      - ./custom_scripts:/app/custom_scripts
    environment:
      - CTCF_MEMORY_LIMIT=16G
      - CTCF_THREADS=8
    networks:
      - ctcf-network
networks:
  ctcf-network:
    driver: bridge
```

## Proxy Configuration

### Automatic Proxy Detection
```bash
# Check proxy settings
./check-proxy.sh

# Expected output:
# ✅ Proxy detected: http://proxy.company.com:8080
# ✅ Docker configured for proxy
# ✅ Container proxy settings applied
```

### Manual Proxy Setup
```bash
# Set proxy environment variables
export HTTP_PROXY="http://proxy.company.com:8080"
export HTTPS_PROXY="http://proxy.company.com:8080"
export NO_PROXY="localhost,127.0.0.1,.local"

# Configure Docker daemon
sudo mkdir -p /etc/systemd/system/docker.service.d
sudo tee /etc/systemd/system/docker.service.d/proxy.conf <<EOF
[Service]
Environment="HTTP_PROXY=$HTTP_PROXY"
Environment="HTTPS_PROXY=$HTTPS_PROXY"
Environment="NO_PROXY=$NO_PROXY"
EOF

# Restart Docker
sudo systemctl daemon-reload
sudo systemctl restart docker
```

### Container Proxy Configuration
```bash
# Build with proxy settings
docker build \
  --build-arg HTTP_PROXY=$HTTP_PROXY \
  --build-arg HTTPS_PROXY=$HTTPS_PROXY \
  -t ctcf-predictor:latest .
```

## Performance Optimization

### Memory Management
```bash
# Set container memory limits
docker run --memory="8g" --memory-swap="12g" \
  ctcf-predictor:latest \
  Rscript scripts/build_pwm_robust.R
```

### CPU Optimization
```bash
# Use all available cores
docker run --cpus="0" \
  ctcf-predictor:latest \
  Rscript scripts/efficient_aligned_pwm.R

# Limit CPU usage
docker run --cpus="4.0" \
  ctcf-predictor:latest \
  Rscript scripts/simple_aligned_pwm.R
```

### Storage Optimization
```bash
# Use tmpfs for temporary files
docker run --tmpfs /tmp:rw,size=2g \
  ctcf-predictor:latest \
  ./test_pipeline_chromosome_split.sh
```

## Troubleshooting

### Common Docker Issues

**Docker daemon not running**:
```bash
# Linux
sudo systemctl start docker

# Windows/macOS
# Start Docker Desktop application
```

**Permission denied**:
```bash
# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker

# Or run with sudo (not recommended)
sudo docker run --rm ctcf-predictor:latest Rscript --version
```

**Out of memory**:
```bash
# Increase Docker memory limits
# Docker Desktop: Settings → Resources → Memory
# Or use command line limits
docker run --memory="16g" ctcf-predictor:latest
```

**Network connectivity issues**:
```bash
# Check proxy configuration
./check-proxy.sh --verbose

# Test network connectivity
docker run --rm ctcf-predictor:latest ping -c 3 google.com

# Check DNS resolution
docker run --rm ctcf-predictor:latest nslookup github.com
```

### Container Debugging

**Access container logs**:
```bash
# View container logs
docker logs $(docker ps -q --filter ancestor=ctcf-predictor:latest)

# Follow logs in real-time
docker logs -f $(docker ps -q --filter ancestor=ctcf-predictor:latest)
```

**Interactive debugging**:
```bash
# Start container with shell
docker run -it --rm --entrypoint bash ctcf-predictor:latest

# Inspect container filesystem
docker run --rm ctcf-predictor:latest ls -la /app/scripts/
```

**Performance monitoring**:
```bash
# Monitor container resource usage
docker stats $(docker ps -q --filter ancestor=ctcf-predictor:latest)
```

## Multi-Platform Support

### ARM64 Support (Apple Silicon, ARM servers)
```bash
# Build for ARM64
docker buildx build --platform linux/arm64 -t ctcf-predictor:arm64 .

# Run on ARM64
docker run --rm ctcf-predictor:arm64 Rscript --version
```

### Multi-architecture Build
```bash
# Build for multiple platforms
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  -t ctcf-predictor:multi \
  --push .
```

## Production Deployment

### Docker Swarm
```bash
# Initialize swarm
docker swarm init

# Deploy stack
docker stack deploy -c docker-compose.yml ctcf-stack
```

### Kubernetes
```yaml
# ctcf-deployment.yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: ctcf-pipeline
spec:
  replicas: 1
  selector:
    matchLabels:
      app: ctcf-pipeline
  template:
    metadata:
      labels:
        app: ctcf-pipeline
    spec:
      containers:
      - name: ctcf-pipeline
        image: ctcf-predictor:latest
        resources:
          requests:
            memory: "8Gi"
            cpu: "4"
          limits:
            memory: "16Gi"
            cpu: "8"
```

## Next Steps

### Verify Setup
1. **Test Installation**: Run `./smart-startup.sh` successfully
2. **Quick Test**: Execute `./run-in-docker.sh test_chromosome_split.R demo`
3. **Full Pipeline**: Try `./test_pwm_improvements_with_null_analysis.sh`

### Explore Pipeline
- **[User Guide](10-user-guide.md)**: Learn to operate the pipeline
- **[System Architecture](07-system-architecture.md)**: Understand the design
- **[Configuration Guide](09-configuration.md)**: Customize parameters

---

**Success Indicator**: When you see PWM files generated in the `results/` directory with information content > 8 bits, your Docker environment is properly configured and ready for production use.
