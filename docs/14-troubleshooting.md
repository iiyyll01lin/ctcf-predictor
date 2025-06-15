# Troubleshooting

> **🔧 Common issues, solutions, and debugging strategies**

## Quick Diagnostic Checklist

When issues arise, start with this systematic checklist:

### ✅ **System Health Check**
```bash
# 1. Docker status
docker --version && docker info

# 2. Container accessibility  
docker run --rm ctcf-predictor:latest echo "✅ Container OK"

# 3. R environment
docker run --rm ctcf-predictor:latest Rscript -e "cat('✅ R OK\n')"

# 4. R packages
docker run --rm ctcf-predictor:latest Rscript -e "library(Biostrings); cat('✅ Packages OK\n')"

# 5. File permissions
ls -la *.sh && echo "✅ Scripts executable"
```

### ✅ **Data Health Check**
```bash
# 1. Input files exist
ls -la data/ && echo "✅ Data directory OK"

# 2. Results directory writable
mkdir -p results && touch results/test.txt && rm results/test.txt && echo "✅ Results writable"

# 3. Network connectivity
./check-proxy.sh && echo "✅ Network OK"
```

## Common Issues and Solutions

### Issue 1: Docker Container Won't Start

**Symptoms**:
```
Error: Cannot connect to the Docker daemon
Error: docker: command not found
Permission denied while trying to connect to Docker
```

**Diagnosis**:
```bash
# Check Docker installation
which docker
docker --version

# Check Docker daemon status
sudo systemctl status docker  # Linux
# Docker Desktop status (Windows/Mac)

# Check user permissions
groups $USER | grep docker
```

**Solutions**:

**Docker not installed**:
```bash
# Linux (Ubuntu/Debian)
curl -fsSL https://get.docker.com | sh
sudo usermod -aG docker $USER
newgrp docker

# Windows: Install Docker Desktop
# Mac: Install Docker Desktop
```

**Docker daemon not running**:
```bash
# Linux
sudo systemctl start docker
sudo systemctl enable docker

# Windows/Mac: Start Docker Desktop application
```

**Permission issues**:
```bash
# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker

# Or run with sudo (not recommended)
sudo ./smart-startup.sh
```

### Issue 2: Memory Errors / Out of Memory

**Symptoms**:
```
Error: cannot allocate vector of size X GB
Container killed (exit code 137)
java.lang.OutOfMemoryError
```

**Diagnosis**:
```bash
# Check available memory
free -h

# Check Docker memory limits
docker system info | grep Memory

# Monitor memory usage during pipeline
docker stats
```

**Solutions**:

**Increase Docker memory limits**:
```bash
# Docker Desktop: Settings → Resources → Memory (set to 8-16GB)

# Command line limits
docker run --memory="8g" --memory-swap="12g" ctcf-predictor:latest
```

**Reduce pipeline memory usage**:
```bash
# Use smaller batch sizes
export CTCF_BATCH_SIZE="500"
export CTCF_MEMORY_LIMIT="4GB"

# Use efficient algorithms
./run-in-docker.sh scripts/efficient_aligned_pwm.R

# Process subsets instead of full dataset
./run-in-docker.sh scripts/build_subset_pwm.R 1000
```

**System-level solutions**:
```bash
# Close other applications
# Add swap space (Linux)
sudo fallocate -l 4G /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

### Issue 3: Network/Proxy Issues

**Symptoms**:
```
Curl: (7) Failed to connect to host
SSL certificate problem
Proxy authentication required
Download failed: timeout
```

**Diagnosis**:
```bash
# Test network connectivity
ping google.com
curl -I https://github.com

# Check proxy settings
echo $HTTP_PROXY $HTTPS_PROXY
./check-proxy.sh --verbose
```

**Solutions**:

**Configure proxy settings**:
```bash
# Automatic proxy detection
./check-proxy.sh
source ./set-proxy.sh

# Manual proxy configuration
export HTTP_PROXY="http://proxy.company.com:8080"
export HTTPS_PROXY="http://proxy.company.com:8080"
export NO_PROXY="localhost,127.0.0.1,.local"
```

**Docker proxy configuration**:
```bash
# Configure Docker daemon proxy
sudo mkdir -p /etc/systemd/system/docker.service.d
sudo tee /etc/systemd/system/docker.service.d/proxy.conf <<EOF
[Service]
Environment="HTTP_PROXY=$HTTP_PROXY"
Environment="HTTPS_PROXY=$HTTPS_PROXY"
EOF

sudo systemctl daemon-reload
sudo systemctl restart docker
```

**Alternative data download**:
```bash
# Use fallback download script
./download-with-fallback.sh

# Manual download with different tools
wget --no-check-certificate https://github.com/...
```

### Issue 4: Low PWM Quality (Information Content < 5 bits)

**Symptoms**:
```
PWM Information Content: 0.695 bits (POOR)
Conserved positions: 0
Quality grade: POOR - Requires significant improvement
```

**Diagnosis**:
```bash
# Check input data quality
./run-in-docker.sh scripts/analyze_sequence_quality.R data/training_sequences.fasta

# Check alignment quality
./run-in-docker.sh scripts/analyze_sequence_alignment.R
```

**Solutions**:

**Improve sequence alignment**:
```bash
# Try different alignment methods
./run-in-docker.sh scripts/advanced_alignment.R consensus
./run-in-docker.sh scripts/advanced_alignment.R length
./run-in-docker.sh scripts/advanced_alignment.R progressive
```

**Use quality-filtered subsets**:
```bash
# Build high-quality subset PWMs
./run-in-docker.sh scripts/build_subset_pwm.R 1000 2000 5000

# Expected improvement: 15-20 bits information content
```

**Check data source**:
```r
# Verify sequences contain CTCF binding sites
sequences <- readDNAStringSet("data/training_sequences.fasta")
cat("Number of sequences:", length(sequences), "\n")
cat("Average length:", mean(width(sequences)), "\n")
cat("GC content:", mean(letterFrequency(sequences, "GC", as.prob=TRUE)), "\n")

# Look for CTCF-like patterns manually
consensus_pattern <- "CCGC[GATC][GATC]GG[GATC]GGCAG"
matches <- vcountPattern(consensus_pattern, sequences, max.mismatch = 3)
cat("Sequences with CTCF-like patterns:", sum(matches > 0), "\n")
```

### Issue 5: R Package Installation Failures

**Symptoms**:
```
Error: package 'Biostrings' is not available
Warning: installation of package 'X' had non-zero exit status
Error in library(Biostrings): there is no package called 'Biostrings'
```

**Diagnosis**:
```bash
# Check R version
docker run --rm ctcf-predictor:latest Rscript -e "R.version.string"

# Check available packages
docker run --rm ctcf-predictor:latest Rscript -e "installed.packages()[,1]"

# Test Bioconductor
docker run --rm ctcf-predictor:latest Rscript -e "source('https://bioconductor.org/biocLite.R')"
```

**Solutions**:

**Rebuild Docker image**:
```bash
# Clean rebuild
docker system prune -f
docker build --no-cache -t ctcf-predictor:latest .
```

**Manual package installation**:
```bash
# Install missing packages
docker run -it --rm ctcf-predictor:latest bash
Rscript -e "
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install('Biostrings')
"
```

**Alternative package sources**:
```r
# Use different CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))
install.packages("ggplot2")

# Use development versions
devtools::install_github("bioconductor/Biostrings")
```

### Issue 6: File Permission Errors

**Symptoms**:
```
Permission denied
cannot create directory 'results'
Error in file.create(): cannot create file
bash: ./script.sh: Permission denied
```

**Diagnosis**:
```bash
# Check file permissions
ls -la *.sh
ls -la scripts/
ls -ld results/

# Check user/group ownership
id
```

**Solutions**:

**Fix script permissions**:
```bash
# Make scripts executable
chmod +x *.sh
chmod +x scripts/*.sh
find . -name "*.sh" -exec chmod +x {} \;
```

**Fix directory permissions**:
```bash
# Create and fix results directory
mkdir -p results
chmod 755 results

# Fix ownership if needed
sudo chown -R $USER:$USER .
```

**Docker volume permissions**:
```bash
# Run with user mapping
docker run --user $(id -u):$(id -g) -v $(pwd):/app ctcf-predictor:latest

# Or fix permissions after container run
sudo chown -R $USER:$USER results/
```

### Issue 7: Statistical Validation Failures

**Symptoms**:
```
P-value: 0.234 (NOT SIGNIFICANT)
Effect size: 0.12 (negligible)
No evidence above random expectation
```

**Diagnosis**:
```bash
# Check null model generation
ls -la results/null_models/
./run-in-docker.sh scripts/validate_null_models.R

# Analyze null model statistics
./run-in-docker.sh scripts/analyze_null_statistics.R
```

**Solutions**:

**Regenerate null models**:
```bash
# Clean regeneration
rm -rf results/null_models/
./run-in-docker.sh scripts/generate_null_models.R --replicates 100
```

**Check data quality**:
```bash
# Ensure using real CTCF binding sites
./run-in-docker.sh scripts/validate_binding_sites.R data/training_sequences.fasta

# Use experimentally validated sites
./run-in-docker.sh scripts/download_validated_sites.R
```

**Use quality subsets**:
```bash
# High-quality subsets should show significance
./run-in-docker.sh scripts/build_subset_pwm.R 1000
./run-in-docker.sh scripts/statistical_significance_test.R results/subset_pwm_size1000.rds
```

## Advanced Debugging

### Enable Debug Mode

**Verbose logging**:
```bash
# Enable detailed logging
export CTCF_DEBUG=true
export CTCF_VERBOSE=true
./run-in-docker.sh scripts/build_pwm_robust.R
```

**R debugging**:
```r
# Add to R scripts for debugging
options(error = traceback)
options(warn = 2)  # Convert warnings to errors

# Debug specific functions
debug(build_pwm_function)
browser()  # Interactive debugging
```

### Performance Profiling

**Memory profiling**:
```r
# Add to R scripts
library(pryr)
mem_used()  # Check memory usage

# Profile memory usage
Rprof("memory.prof", memory.profiling = TRUE)
# ... your code ...
Rprof(NULL)
summaryRprof("memory.prof")
```

**Time profiling**:
```bash
# Time entire pipeline
time ./test_pwm_improvements_with_null_analysis.sh

# Profile individual scripts
time ./run-in-docker.sh scripts/build_subset_pwm.R 1000
```

### Log Analysis

**Docker logs**:
```bash
# View container logs
docker logs $(docker ps -ql)

# Follow logs in real-time
docker logs -f $(docker ps -ql)

# Export logs for analysis
docker logs $(docker ps -ql) > pipeline.log 2>&1
```

**R session logs**:
```r
# Save R session for debugging
save.image("debug_session.RData")

# Capture warnings and errors
warnings()
traceback()
```

## Environment-Specific Issues

### Windows-Specific Issues

**WSL2 Integration**:
```powershell
# Enable WSL2 integration in Docker Desktop
# Settings → Resources → WSL Integration

# Check WSL2 status
wsl --list --verbose
wsl --set-default-version 2
```

**Path Issues**:
```powershell
# Use WSL paths in Docker commands
docker run -v /mnt/c/Users/username/ctcf-predictor:/app ctcf-predictor:latest
```

**Line Ending Issues**:
```bash
# Convert Windows line endings to Unix
dos2unix *.sh
dos2unix scripts/*.sh
```

### macOS-Specific Issues

**Apple Silicon (M1/M2)**:
```bash
# Build for ARM architecture
docker buildx build --platform linux/arm64 -t ctcf-predictor:arm64 .

# Run ARM container
docker run --platform linux/arm64 ctcf-predictor:arm64
```

**File Sharing Performance**:
```bash
# Use named volumes for better performance
docker volume create ctcf-data
docker run -v ctcf-data:/app/data ctcf-predictor:latest
```

### Linux Distribution-Specific Issues

**CentOS/RHEL**:
```bash
# Install Docker on CentOS
sudo yum install -y yum-utils
sudo yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
sudo yum install docker-ce
```

**Ubuntu/Debian Package Conflicts**:
```bash
# Remove conflicting packages
sudo apt-get remove docker docker-engine docker.io containerd runc
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

## Getting Help

### Information to Collect

When reporting issues, provide:

1. **System Information**:
   ```bash
   uname -a
   docker --version
   docker info
   ```

2. **Error Messages**:
   ```bash
   # Full error output
   ./run-in-docker.sh scripts/build_pwm.R 2>&1 | tee error.log
   ```

3. **Configuration**:
   ```bash
   # Current settings
   env | grep CTCF
   cat scripts/preprocess_config.json
   ```

4. **Data Information**:
   ```bash
   # Input data summary
   ls -la data/
   head -20 data/training_sequences.fasta
   ```

### Diagnostic Script

**Run comprehensive diagnostics**:
```bash
#!/bin/bash
# diagnostics.sh
echo "=== CTCF Pipeline Diagnostics ==="
echo "Date: $(date)"
echo "User: $(whoami)"
echo "System: $(uname -a)"
echo ""

echo "=== Docker Status ==="
docker --version
docker info | head -20
echo ""

echo "=== Container Test ==="
docker run --rm ctcf-predictor:latest echo "Container OK"
docker run --rm ctcf-predictor:latest Rscript -e "cat('R OK\n')"
echo ""

echo "=== File System ==="
ls -la *.sh | head -5
ls -la data/ | head -5
ls -la results/ | head -5
echo ""

echo "=== Network ==="
ping -c 3 google.com
echo ""

echo "=== Environment ==="
env | grep -E "(PROXY|CTCF|DOCKER)" | head -10
```

### Self-Help Resources

1. **Documentation**: Check relevant sections in this documentation
2. **Examples**: Look at working examples in `test_*.sh` scripts  
3. **Community**: Search GitHub issues for similar problems
4. **Logs**: Always check logs first before asking for help

---

Most issues can be resolved through systematic debugging and checking the basics first. When in doubt, start with the diagnostic checklist and work through the common solutions.

**Next Reading**: [API Reference](15-api-reference.md) for technical specifications and advanced usage patterns.
