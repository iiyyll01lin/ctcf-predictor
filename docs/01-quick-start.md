# Quick Start Guide

> **🚀 Get the CTCF PWM Testing Pipeline running in 5 minutes**

## Prerequisites

- **Docker** (recommended) or **R 4.0+** with required packages
- **8GB RAM** minimum, 16GB recommended
- **5GB free disk space** for full dataset
- **Stable internet connection** for data download

## Option 1: Docker Quick Start (Recommended)

### Step 1: Clone and Setup
```bash
git clone https://github.com/organization/ctcf-predictor.git
cd ctcf-predictor
```

### Step 2: Start Docker Environment
```bash
# Automatically detects proxy and starts containerized environment
./smart-startup.sh
```

### Step 3: Run Demo Pipeline
```bash
# Run the complete pipeline with demo data
./test_pipeline_chromosome_split.sh demo
```

**Expected Output**: Pipeline completes in ~10 minutes, generates PWM files in `results/` directory.

## Option 2: Local Installation Quick Start

### Step 1: Install R Dependencies
```bash
# Install required R packages
Rscript -e "install.packages(c('Biostrings', 'seqinr', 'ggplot2', 'pROC', 'jsonlite'))"
```

### Step 2: Download Data
```bash
# Download CTCF ChIP-seq data
bash scripts/download_data.sh
```

### Step 3: Run Basic Pipeline
```bash
# Run simple PWM building
Rscript scripts/simple_aligned_pwm.R
```

## Verify Installation

### Check Docker Setup
```bash
# Test Docker environment
docker run --rm ctcf-predictor:latest Rscript --version
```

### Check Local Setup
```bash
# Test R dependencies
Rscript -e "library(Biostrings); cat('✅ Dependencies OK\\n')"
```

## Quick Test Runs

### Demo Mode (Fast - 2 minutes)
```bash
# Test with small dataset
./run-in-docker.sh test_chromosome_split.R demo

# Expected: PWM with ~8-12 bits information content
```

### Production Mode (Full - 30 minutes)
```bash
# Run complete pipeline
./test_pwm_improvements_with_null_analysis.sh

# Expected: Multiple PWMs with quality reports
```

## Understanding the Output

### Key Files Generated
```
results/
├── simple_aligned_pwm.rds           # Basic PWM model
├── subset_pwm_size1000.rds          # High-quality subset PWM
├── robust_pwm_metadata.json         # Quality metrics
└── enhanced_pwm_comparison.html     # Visual comparison report
```

### Quality Indicators
- **Information Content > 8 bits**: Good quality PWM
- **Information Content > 15 bits**: Excellent quality PWM
- **Conserved Positions > 3**: Acceptable motif definition
- **Statistical Significance p < 0.05**: Validated against null models

## Next Steps

1. **Understand the Science**: Read [Biological Background](04-biological-background.md)
2. **Explore Results**: See [Results Analysis](12-results-analysis.md)
3. **Customize Pipeline**: Check [Configuration Guide](09-configuration.md)
4. **Add New Methods**: See [Extending the Pipeline](13-extending-pipeline.md)

## Common Quick Start Issues

### Docker Issues
```bash
# If Docker fails to start
docker system prune -f
./smart-startup.sh --rebuild
```

### Memory Issues
```bash
# Reduce memory usage
export CTCF_MEMORY_LIMIT="4GB"
./run-in-docker.sh test_chromosome_split.R demo
```

### Network Issues
```bash
# Manual proxy setup
export HTTP_PROXY="http://your-proxy:port"
./set-proxy.sh
```

## Performance Expectations

| Mode | Time | Memory | Disk Space | Output Quality |
|------|------|--------|------------|----------------|
| **Demo** | 2-5 min | 2GB | 100MB | Basic validation |
| **Standard** | 15-30 min | 8GB | 2GB | Production ready |
| **Full Analysis** | 1-2 hours | 16GB | 5GB | Research grade |

## Troubleshooting Quick Fixes

### "R package not found"
```bash
# Reinstall dependencies
./run-in-docker.sh scripts/install_dependencies.sh
```

### "Permission denied"
```bash
# Fix file permissions
chmod +x *.sh
sudo chmod +x scripts/*.sh
```

### "Network timeout"
```bash
# Configure proxy detection
./check-proxy.sh
source set-proxy.sh
```

## Success Verification

✅ **You've successfully completed quick start when you see:**
- PWM files in `results/` directory
- Information content > 8 bits reported
- No fatal errors in console output
- Docker container starts without issues (Docker mode)

🔬 **Ready for Science**: You now have a working CTCF PWM prediction pipeline!

---

**Next Recommended Reading**: [Biological Background](04-biological-background.md) to understand what CTCF is and why this pipeline matters for genomics research.
