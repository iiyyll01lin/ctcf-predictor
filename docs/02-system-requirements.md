# System Requirements

> **🖥️ Hardware, software, and environment requirements for optimal pipeline performance**

## Hardware Requirements

### Minimum Requirements (Demo Mode)
- **CPU**: 2 cores, x86_64 architecture
- **Memory**: 4GB RAM 
- **Storage**: 1GB free disk space
- **Network**: Stable internet connection (500MB download)

### Recommended Configuration (Production)
- **CPU**: 4+ cores, modern x86_64 processor
- **Memory**: 8-16GB RAM
- **Storage**: 5-10GB free disk space (SSD preferred)
- **Network**: Stable broadband (3GB initial download)

### High-Performance Configuration (Research)
- **CPU**: 8+ cores, multi-threaded processing
- **Memory**: 16-32GB RAM
- **Storage**: 10-20GB free SSD space
- **Network**: High-speed connection for large datasets

## Software Dependencies

### Docker Mode (Recommended)

#### Required
- **Docker**: Version 20.0+ 
- **Docker Compose**: Version 1.27+
- **Git**: For repository cloning

#### Platform-Specific Requirements

**Linux (Ubuntu/Debian)**:
```bash
# Docker installation
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER
```

**Windows 10/11**:
- Docker Desktop with WSL2 backend
- Windows Subsystem for Linux 2 (WSL2)
- PowerShell 5.1+ or PowerShell Core 7+

**macOS**:
- Docker Desktop for Mac
- macOS 10.15+ (Catalina or later)
- Xcode Command Line Tools

### Local Installation Mode

#### Core R Environment
- **R**: Version 4.0.0+ (4.3.0+ recommended)
- **Rscript**: Command-line R execution
- **Memory**: R configured with sufficient memory limits

#### Required R Packages
```r
# Core bioinformatics packages
install.packages("Biostrings")      # DNA sequence analysis
install.packages("seqinr")         # Sequence input/output
install.packages("GenomicRanges")  # Genomic coordinate handling

# Statistical and visualization packages
install.packages("ggplot2")        # Data visualization
install.packages("pROC")           # ROC curve analysis
install.packages("jsonlite")       # JSON configuration handling
```

#### System Tools
```bash
# Linux/macOS
sudo apt-get install curl wget unzip bedtools

# Or via conda
conda install -c bioconda bedtools
```

**Windows**:
- Install via chocolatey or download binaries
- Git Bash or PowerShell with appropriate tools

## Operating System Support

### Fully Supported
- **Ubuntu**: 18.04 LTS, 20.04 LTS, 22.04 LTS
- **Debian**: 9 (Stretch), 10 (Buster), 11 (Bullseye)
- **CentOS**: 7, 8 (Stream)
- **Red Hat Enterprise Linux**: 7, 8, 9
- **macOS**: 10.15 (Catalina), 11 (Big Sur), 12 (Monterey), 13 (Ventura)
- **Windows**: 10 (version 1903+), 11 with WSL2

### Partially Supported
- **Alpine Linux**: Docker mode only
- **Amazon Linux**: 2 with modifications
- **Other Unix-like**: May require dependency adjustments

## Network Requirements

### Data Downloads
- **Initial Setup**: ~3GB (reference genome + ChIP-seq data)
- **Null Models**: ~500MB (optional benchmark data)
- **Updates**: ~100MB (periodic updates)

### Proxy Configuration
The pipeline includes automatic proxy detection for enterprise environments:

```bash
# Automatic proxy detection
./check-proxy.sh

# Manual proxy configuration
export HTTP_PROXY="http://proxy.company.com:8080"
export HTTPS_PROXY="http://proxy.company.com:8080"
```

### Firewall Requirements
**Outbound connections needed**:
- HTTP/HTTPS (ports 80, 443): Data downloads
- Docker registry access: Container images
- GitHub access: Repository updates

## Performance Scaling

### Dataset Size Scaling

| Dataset Size | Memory Usage | Processing Time | Recommended Config |
|--------------|--------------|-----------------|-------------------|
| Demo (1K sequences) | 2GB | 5 minutes | Minimum setup |
| Standard (10K sequences) | 4-8GB | 15-30 minutes | Recommended setup |
| Large (35K sequences) | 8-16GB | 1-2 hours | High-performance setup |
| Research (100K+ sequences) | 16-32GB | 3-6 hours | Cluster/workstation |

### Parallelization Support
- **Multi-core processing**: Automatic detection and usage
- **Memory-efficient batching**: Large datasets processed in chunks
- **Disk I/O optimization**: Sequential and parallel file operations

## Storage Requirements

### Disk Space Breakdown
```
Initial data download:    ~3GB
  ├── hg38 reference:      2.5GB
  ├── CTCF ChIP-seq:       400MB
  └── Test datasets:       100MB

Working space:           ~2GB
  ├── Intermediate files:  1GB
  ├── Results storage:     500MB
  └── Temporary files:     500MB

Docker images:           ~1.5GB
  ├── Base R image:        800MB
  ├── Pipeline image:      500MB
  └── Dependencies:        200MB

Total recommended:       ~7GB free space
```

### I/O Performance
- **SSD recommended**: Significant performance improvement for large datasets
- **HDD acceptable**: For demo and small-scale usage
- **Network storage**: Supported but may impact performance

## Memory Management

### R Memory Configuration
```bash
# Set R memory limits
export R_MAX_VSIZE="16G"
export R_MAX_NUM_DLLS=500

# Docker memory limits
docker run --memory="8g" --memory-swap="12g" ctcf-predictor:latest
```

### Memory Usage Patterns
- **Peak usage**: During PWM construction and alignment
- **Sustained usage**: Throughout statistical analysis
- **Memory release**: Automatic garbage collection between phases

## Validation and Testing

### System Compatibility Test
```bash
# Quick system check
./scripts/check_system_requirements.sh

# Docker environment test
docker run --rm ctcf-predictor:latest Rscript --version

# Local environment test
Rscript -e "source('scripts/test_dependencies.R')"
```

### Performance Benchmarks
```bash
# Performance benchmark (5 minutes)
./run-in-docker.sh scripts/benchmark_system.R

# Expected outputs:
# - Processing speed: sequences/second
# - Memory efficiency: MB/thousand sequences
# - I/O performance: files/second
```

## Cloud and HPC Environments

### Cloud Platforms
- **AWS**: EC2 instances (r5.large or larger recommended)
- **Google Cloud**: Compute Engine (n2-standard-4 or larger)
- **Azure**: Virtual machines (Standard_D4s_v3 or larger)

### HPC Clusters
- **Slurm**: Job submission scripts provided
- **PBS/Torque**: Compatible with module loading
- **Docker/Singularity**: Container support for reproducibility

### Configuration Examples
```bash
# AWS EC2 user data script
#!/bin/bash
apt-get update
curl -fsSL https://get.docker.com | sh
git clone https://github.com/organization/ctcf-predictor.git
cd ctcf-predictor && ./smart-startup.sh
```

## Troubleshooting Common Issues

### Docker Issues
```bash
# Docker daemon not running
sudo systemctl start docker

# Permission denied
sudo usermod -aG docker $USER
newgrp docker

# Memory limit exceeded
docker system prune -f
```

### R Package Issues
```bash
# Package installation failures
Rscript -e "options(repos='https://cloud.r-project.org'); install.packages('Biostrings')"

# Bioconductor packages
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')"
```

### Network Issues
```bash
# Proxy detection failure
./check-proxy.sh --verbose
source ./set-proxy.sh

# DNS resolution issues
nslookup github.com
```

## Security Considerations

### Container Security
- Pipeline runs in isolated Docker containers
- No root privileges required for normal operation
- Network access limited to necessary connections

### Data Privacy
- All processing occurs locally
- No data transmitted to external services
- Genomic data remains on your infrastructure

---

**Next Steps**: Once your system meets these requirements, proceed to [Docker Setup](03-docker-setup.md) for containerized deployment or [User Guide](10-user-guide.md) for local installation.
