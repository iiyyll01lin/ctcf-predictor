# Configuration Guide

> **⚙️ Parameters, Settings, and Customization**  
> Complete reference for configuring the CTCF PWM Testing Pipeline for different research scenarios and computational environments.

## 🎯 Configuration Philosophy

The pipeline uses a **hierarchical configuration system** that balances ease of use with powerful customization capabilities:

1. **Default Settings** - Optimized for most CTCF research scenarios
2. **Environment Variables** - System-wide configuration
3. **Configuration Files** - Project-specific settings
4. **Command-line Arguments** - Execution-specific overrides

## 🏗️ Configuration Hierarchy

### Priority Order (Highest to Lowest)

```
1. Command-line Arguments    --parameter value
2. Environment Variables     CTCF_PARAMETER=value
3. Configuration Files       config/pipeline.yml
4. Script Defaults          Built-in default values
```

## 📁 Configuration Files

### Main Configuration: `config/pipeline.yml`

**Location:** `config/pipeline.yml`

**Complete Configuration Template:**
```yaml
# CTCF PWM Testing Pipeline Configuration
# Version: 2.0

# ======================
# Data Processing
# ======================
data:
  input:
    format: "fasta"                    # Input format: fasta, bed, custom
    min_sequences: 100                 # Minimum sequences required
    max_sequences: 50000               # Maximum sequences to process
    
  preprocessing:
    min_length: 50                     # Minimum sequence length (bp)
    max_length: 500                    # Maximum sequence length (bp)
    quality_threshold: 0.8             # Quality score threshold (0-1)
    remove_duplicates: true            # Remove duplicate sequences
    filter_repetitive: true            # Filter repetitive elements
    gc_content_range: [0.3, 0.7]       # Acceptable GC content range
    complexity_threshold: 0.6          # Sequence complexity threshold
    
  validation:
    split_method: "chromosome"         # Validation method: chromosome, random, cv
    test_chromosomes: ["16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    train_ratio: 0.8                   # Training set ratio (for random split)
    cv_folds: 5                        # Cross-validation folds

# ======================
# Alignment Settings
# ======================
alignment:
  method: "integrated"                 # Alignment method: center, consensus, integrated
  center_window: 50                    # Window size for center alignment
  consensus_threshold: 0.8             # Consensus threshold (0-1)
  max_iterations: 10                   # Maximum alignment iterations
  convergence_threshold: 0.01          # Convergence threshold
  
  refinement:
    enable: true                       # Enable iterative refinement
    max_refinement_cycles: 3           # Maximum refinement cycles
    improvement_threshold: 0.1         # Minimum improvement required

# ======================
# PWM Construction
# ======================
pwm:
  pseudocount: 0.1                     # Pseudocount for probability calculation
  background_freq:                     # Background nucleotide frequencies
    A: 0.25
    C: 0.25
    G: 0.25
    T: 0.25
  
  information_content:
    min_total_ic: 8.0                  # Minimum total information content
    min_conserved_positions: 3         # Minimum conserved positions (>1 bit)
    conservation_threshold: 1.0         # Information content threshold for conservation
    
  quality_thresholds:
    excellent: 16.0                    # Excellent quality threshold (bits)
    good: 12.0                         # Good quality threshold (bits)
    acceptable: 8.0                    # Minimum acceptable threshold (bits)

# ======================
# Statistical Validation
# ======================
statistics:
  significance_level: 0.01             # Statistical significance level
  null_model_iterations: 1000          # Null model iterations
  bootstrap_iterations: 500            # Bootstrap iterations
  
  tests:
    - "wilcoxon"                       # Wilcoxon rank-sum test
    - "ks"                             # Kolmogorov-Smirnov test
    - "permutation"                    # Permutation test
    
  multiple_testing_correction: "bonferroni"  # bonferroni, fdr, holm

# ======================
# Performance Settings
# ======================
performance:
  threads: 4                           # Number of CPU threads
  memory_limit: "8G"                   # Memory limit
  batch_size: 1000                     # Batch size for processing
  use_parallel: true                   # Enable parallel processing
  temp_dir: "/tmp/ctcf"                # Temporary directory
  
  optimization:
    use_streaming: false               # Use streaming for large datasets
    cache_intermediate: true           # Cache intermediate results
    compress_temp_files: true          # Compress temporary files

# ======================
# Output Settings
# ======================
output:
  formats:                             # Output formats to generate
    - "meme"                           # MEME format
    - "jaspar"                         # JASPAR format
    - "transfac"                       # TRANSFAC format
    
  reports:
    generate_html: true                # Generate HTML reports
    generate_pdf: false                # Generate PDF reports
    include_plots: true                # Include visualization plots
    
  file_naming:
    timestamp: true                    # Include timestamp in filenames
    method_suffix: true                # Include method in filename
    quality_suffix: true               # Include quality grade in filename

# ======================
# Logging and Debugging
# ======================
logging:
  level: "INFO"                        # Log level: DEBUG, INFO, WARN, ERROR
  file: "logs/pipeline.log"            # Log file path
  console: true                        # Log to console
  rotate_logs: true                    # Rotate log files
  max_log_size: "100M"                 # Maximum log file size
  
  debug:
    save_intermediate: false           # Save intermediate results for debugging
    verbose_alignment: false           # Verbose alignment output
    profile_memory: false              # Profile memory usage

# ======================
# Docker Configuration
# ======================
docker:
  image: "ctcf-pipeline:latest"        # Docker image name
  memory_limit: "8g"                   # Docker memory limit
  cpu_limit: "4"                       # Docker CPU limit
  
  volumes:                             # Volume mappings
    data: "/data"
    results: "/results"
    config: "/config"
    
  environment:                         # Environment variables
    CTCF_THREADS: "4"
    CTCF_MEMORY: "8G"
```

### Specialized Configurations

#### High-Quality Configuration: `config/high_quality.yml`
```yaml
# Optimized for maximum quality, slower processing
data:
  preprocessing:
    quality_threshold: 0.95
    min_length: 100
    max_length: 200
    complexity_threshold: 0.8

alignment:
  method: "integrated"
  max_iterations: 20
  refinement:
    max_refinement_cycles: 5

pwm:
  min_total_ic: 12.0
  min_conserved_positions: 5
```

#### Fast Processing Configuration: `config/fast.yml`
```yaml
# Optimized for speed, moderate quality
data:
  preprocessing:
    quality_threshold: 0.7
    filter_repetitive: false

alignment:
  method: "center"
  max_iterations: 3
  refinement:
    enable: false

performance:
  use_parallel: true
  threads: 8
  batch_size: 2000
```

#### Memory-Efficient Configuration: `config/low_memory.yml`
```yaml
# Optimized for low memory usage
performance:
  memory_limit: "2G"
  batch_size: 100
  use_streaming: true
  cache_intermediate: false
  compress_temp_files: true

output:
  reports:
    include_plots: false
```

## 🌍 Environment Variables

### Core Environment Variables

**Data Paths:**
```bash
export CTCF_DATA_DIR="/path/to/data"           # Input data directory
export CTCF_OUTPUT_DIR="/path/to/results"      # Output directory
export CTCF_CONFIG_DIR="/path/to/config"       # Configuration directory
export CTCF_TEMP_DIR="/tmp/ctcf"               # Temporary files directory
export CTCF_LOG_DIR="/var/log/ctcf"            # Log files directory
```

**Processing Parameters:**
```bash
export CTCF_THREADS=4                          # Number of CPU threads
export CTCF_MEMORY="8G"                        # Memory limit
export CTCF_BATCH_SIZE=1000                    # Processing batch size
export CTCF_QUALITY_THRESHOLD=0.8              # Quality threshold
```

**Pipeline Behavior:**
```bash
export CTCF_ALIGNMENT_METHOD="integrated"      # Default alignment method
export CTCF_MIN_IC=8.0                        # Minimum information content
export CTCF_DEBUG_MODE="false"                # Enable debug mode
export CTCF_VERBOSE="false"                   # Verbose output
```

**Docker Settings:**
```bash
export CTCF_DOCKER_IMAGE="ctcf-pipeline:latest"  # Docker image
export CTCF_DOCKER_MEMORY="8g"                   # Docker memory limit
export CTCF_DOCKER_CPUS="4"                      # Docker CPU limit
```

### Environment Setup Scripts

#### `setup_environment.sh`
```bash
#!/bin/bash
# CTCF Pipeline Environment Setup Script
# Configures environment variables and directory structure

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== CTCF Pipeline Environment Setup ===${NC}"

# Detect operating system
OS=$(uname -s)
echo "Detected OS: $OS"

# Base directories
if [[ "$OS" == "MINGW"* ]] || [[ "$OS" == "CYGWIN"* ]] || [[ "$OS" == "MSYS"* ]]; then
    # Windows paths
    export CTCF_HOME="$(pwd)"
    export CTCF_DATA_DIR="$CTCF_HOME/data"
    export CTCF_OUTPUT_DIR="$CTCF_HOME/results"
    export CTCF_CONFIG_DIR="$CTCF_HOME/config"
    export CTCF_TEMP_DIR="/tmp/ctcf-$$"
    export CTCF_LOG_DIR="$CTCF_HOME/logs"
    export CTCF_SCRIPTS_DIR="$CTCF_HOME/scripts"
else
    # Unix-like paths
    export CTCF_HOME="/opt/ctcf-pipeline"
    export CTCF_DATA_DIR="$CTCF_HOME/data"
    export CTCF_OUTPUT_DIR="$CTCF_HOME/results"
    export CTCF_CONFIG_DIR="$CTCF_HOME/config"
    export CTCF_TEMP_DIR="/tmp/ctcf-$$"
    export CTCF_LOG_DIR="$CTCF_HOME/logs"
    export CTCF_SCRIPTS_DIR="$CTCF_HOME/scripts"
fi

# Auto-detect system resources
detect_system_resources() {
    echo -e "${YELLOW}Detecting system resources...${NC}"
    
    # CPU cores
    if command -v nproc &> /dev/null; then
        CPU_CORES=$(nproc)
    elif command -v sysctl &> /dev/null; then
        CPU_CORES=$(sysctl -n hw.ncpu)
    else
        CPU_CORES=4
    fi
    
    # Memory (in GB)
    if command -v free &> /dev/null; then
        TOTAL_MEMORY_GB=$(free -g | awk '/^Mem:/{print $2}')
        AVAILABLE_MEMORY_GB=$(echo "$TOTAL_MEMORY_GB * 0.8" | bc 2>/dev/null || echo "8")
    elif command -v vm_stat &> /dev/null; then
        # macOS
        TOTAL_MEMORY_GB=$(echo "$(vm_stat | grep "Pages free" | awk '{print $3}' | sed 's/\.//' ) * 4096 / 1024 / 1024 / 1024" | bc 2>/dev/null || echo "8")
        AVAILABLE_MEMORY_GB=$(echo "$TOTAL_MEMORY_GB * 0.8" | bc 2>/dev/null || echo "8")
    else
        TOTAL_MEMORY_GB=8
        AVAILABLE_MEMORY_GB=6
    fi
    
    echo "CPU cores: $CPU_CORES"
    echo "Total memory: ${TOTAL_MEMORY_GB}GB"
    echo "Available memory: ${AVAILABLE_MEMORY_GB}GB"
}

# Processing settings
detect_system_resources

export CTCF_THREADS=$CPU_CORES
export CTCF_MEMORY="${AVAILABLE_MEMORY_GB}G"
export CTCF_BATCH_SIZE=1000
export CTCF_QUALITY_THRESHOLD=0.8

# Pipeline behavior settings
export CTCF_ALIGNMENT_METHOD="integrated"
export CTCF_MIN_IC=8.0
export CTCF_DEBUG_MODE="false"
export CTCF_VERBOSE="false"

# Docker settings
export CTCF_DOCKER_IMAGE="ctcf-pipeline:latest"
export CTCF_DOCKER_MEMORY="${AVAILABLE_MEMORY_GB}g"
export CTCF_DOCKER_CPUS="$CPU_CORES"

# Create directory structure
echo -e "${YELLOW}Creating directory structure...${NC}"
create_directories() {
    local dirs=(
        "$CTCF_DATA_DIR"
        "$CTCF_OUTPUT_DIR"
        "$CTCF_CONFIG_DIR"
        "$CTCF_TEMP_DIR"
        "$CTCF_LOG_DIR"
        "$CTCF_SCRIPTS_DIR"
        "$CTCF_DATA_DIR/reference_genome"
        "$CTCF_OUTPUT_DIR/pwms"
        "$CTCF_OUTPUT_DIR/reports"
        "$CTCF_OUTPUT_DIR/plots"
        "$CTCF_OUTPUT_DIR/validation"
    )
    
    for dir in "${dirs[@]}"; do
        if [ ! -d "$dir" ]; then
            mkdir -p "$dir"
            echo "Created: $dir"
        else
            echo "Exists: $dir"
        fi
    done
}

create_directories

# Set permissions
echo -e "${YELLOW}Setting permissions...${NC}"
if [[ "$OS" != "MINGW"* ]] && [[ "$OS" != "CYGWIN"* ]] && [[ "$OS" != "MSYS"* ]]; then
    chmod 755 "$CTCF_TEMP_DIR" 2>/dev/null || true
    chmod 755 "$CTCF_LOG_DIR" 2>/dev/null || true
    chmod 755 "$CTCF_OUTPUT_DIR" 2>/dev/null || true
fi

# Create environment file
echo -e "${YELLOW}Creating environment file...${NC}"
cat > "$CTCF_HOME/.ctcf_env" << EOF
# CTCF Pipeline Environment Variables
# Source this file to load the environment: source .ctcf_env

# Base directories
export CTCF_HOME="$CTCF_HOME"
export CTCF_DATA_DIR="$CTCF_DATA_DIR"
export CTCF_OUTPUT_DIR="$CTCF_OUTPUT_DIR"
export CTCF_CONFIG_DIR="$CTCF_CONFIG_DIR"
export CTCF_TEMP_DIR="$CTCF_TEMP_DIR"
export CTCF_LOG_DIR="$CTCF_LOG_DIR"
export CTCF_SCRIPTS_DIR="$CTCF_SCRIPTS_DIR"

# Processing settings
export CTCF_THREADS=$CTCF_THREADS
export CTCF_MEMORY="$CTCF_MEMORY"
export CTCF_BATCH_SIZE=$CTCF_BATCH_SIZE
export CTCF_QUALITY_THRESHOLD=$CTCF_QUALITY_THRESHOLD

# Pipeline behavior
export CTCF_ALIGNMENT_METHOD="$CTCF_ALIGNMENT_METHOD"
export CTCF_MIN_IC=$CTCF_MIN_IC
export CTCF_DEBUG_MODE="$CTCF_DEBUG_MODE"
export CTCF_VERBOSE="$CTCF_VERBOSE"

# Docker settings
export CTCF_DOCKER_IMAGE="$CTCF_DOCKER_IMAGE"
export CTCF_DOCKER_MEMORY="$CTCF_DOCKER_MEMORY"
export CTCF_DOCKER_CPUS="$CTCF_DOCKER_CPUS"

# Add scripts to PATH
export PATH="\$CTCF_SCRIPTS_DIR:\$PATH"

echo "CTCF Pipeline environment loaded"
EOF

# Check for required tools
echo -e "${YELLOW}Checking for required tools...${NC}"
check_requirements() {
    local missing_tools=()
    
    # Check for R
    if ! command -v R &> /dev/null; then
        missing_tools+=("R")
    fi
    
    # Check for Docker (optional)
    if ! command -v docker &> /dev/null; then
        echo -e "${YELLOW}Warning: Docker not found (optional)${NC}"
    fi
    
    # Check for basic tools
    local tools=("awk" "sed" "grep" "sort" "uniq")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -gt 0 ]; then
        echo -e "${RED}Missing required tools: ${missing_tools[*]}${NC}"
        echo "Please install the missing tools before running the pipeline."
        return 1
    else
        echo -e "${GREEN}All required tools are available${NC}"
        return 0
    fi
}

check_requirements

# Create quick test script
echo -e "${YELLOW}Creating quick test script...${NC}"
cat > "$CTCF_HOME/test_environment.sh" << 'EOF'
#!/bin/bash
# Quick environment test

source .ctcf_env

echo "=== CTCF Pipeline Environment Test ==="
echo "CTCF_HOME: $CTCF_HOME"
echo "CTCF_DATA_DIR: $CTCF_DATA_DIR"
echo "CTCF_OUTPUT_DIR: $CTCF_OUTPUT_DIR"
echo "CTCF_THREADS: $CTCF_THREADS"
echo "CTCF_MEMORY: $CTCF_MEMORY"

echo ""
echo "Directory structure:"
find "$CTCF_HOME" -type d -name "ctcf-*" -o -name "data" -o -name "results" -o -name "config" -o -name "logs" 2>/dev/null | sort

echo ""
echo "Environment test completed successfully!"
EOF

chmod +x "$CTCF_HOME/test_environment.sh"

# Summary
echo -e "${GREEN}=== Environment Setup Complete ===${NC}"
echo ""
echo "Environment configuration:"
echo "  Home: $CTCF_HOME"
echo "  Data: $CTCF_DATA_DIR"
echo "  Output: $CTCF_OUTPUT_DIR"
echo "  Config: $CTCF_CONFIG_DIR"
echo "  Logs: $CTCF_LOG_DIR"
echo "  Threads: $CTCF_THREADS"
echo "  Memory: $CTCF_MEMORY"
echo ""
echo "To load the environment in future sessions:"
echo "  source $CTCF_HOME/.ctcf_env"
echo ""
echo "To test the environment:"
echo "  ./test_environment.sh"
echo ""
echo -e "${GREEN}Setup completed successfully!${NC}"
```
export CTCF_HOME="/opt/ctcf-pipeline"
export CTCF_DATA_DIR="$CTCF_HOME/data"
export CTCF_OUTPUT_DIR="$CTCF_HOME/results"
export CTCF_CONFIG_DIR="$CTCF_HOME/config"
export CTCF_TEMP_DIR="/tmp/ctcf-$$"
export CTCF_LOG_DIR="$CTCF_HOME/logs"

# Processing settings
export CTCF_THREADS=$(nproc)
export CTCF_MEMORY="$(free -g | awk '/^Mem:/{print int($2*0.8)}')G"
export CTCF_BATCH_SIZE=1000

# Create directories
mkdir -p "$CTCF_DATA_DIR" "$CTCF_OUTPUT_DIR" "$CTCF_CONFIG_DIR" \
         "$CTCF_TEMP_DIR" "$CTCF_LOG_DIR"

# Set permissions
chmod 755 "$CTCF_TEMP_DIR"
chmod 755 "$CTCF_LOG_DIR"

echo "CTCF Pipeline environment configured:"
echo "  Data: $CTCF_DATA_DIR"
echo "  Output: $CTCF_OUTPUT_DIR"
echo "  Threads: $CTCF_THREADS"
echo "  Memory: $CTCF_MEMORY"
```

## ⚙️ Command-Line Parameters

### Global Parameters (Available for Most Scripts)

**Input/Output:**
```bash
--input FILE              # Input file path
--output FILE             # Output file path
--config FILE             # Configuration file path
--format FORMAT           # Output format (meme, jaspar, transfac)
--log-file FILE           # Log file path
```

**Processing Control:**
```bash
--threads N               # Number of threads to use
--memory SIZE             # Memory limit (e.g., 8G)
--batch-size N            # Batch size for processing
--temp-dir DIR            # Temporary directory
--no-parallel             # Disable parallel processing
```

**Quality Control:**
```bash
--min-ic FLOAT            # Minimum information content
--quality-threshold FLOAT # Quality threshold (0-1)
--min-sequences N         # Minimum number of sequences
--max-sequences N         # Maximum number of sequences
```

**Behavior Modifiers:**
```bash
--verbose                 # Verbose output
--debug                   # Debug mode
--quiet                   # Suppress non-error output
--force                   # Force overwrite existing files
--dry-run                 # Show what would be done without executing
```

### Script-Specific Parameters

#### `build_pwm_robust.R`
```bash
--sequences FILE          # Input aligned sequences (required)
--output FILE             # Output PWM file (required)
--pseudocount FLOAT       # Pseudocount value (default: 0.1)
--background A,C,G,T      # Background frequencies (default: 0.25,0.25,0.25,0.25)
--min-conserved N         # Minimum conserved positions (default: 3)
--format FORMAT           # Output format: meme, jaspar, transfac (default: meme)
```

#### `advanced_alignment.R`
```bash
--sequences FILE          # Input sequences (required)
--method METHOD           # Alignment method: center, consensus, integrated (required)
--output FILE             # Output aligned sequences (required)
--window-size N           # Window size for center alignment (default: 50)
--consensus-threshold FLOAT # Consensus threshold (default: 0.8)
--max-iterations N        # Maximum iterations (default: 10)
--refinement              # Enable iterative refinement
```

#### `validate_pwm_quality.R`
```bash
--pwm FILE                # PWM file to validate (required)
--sequences FILE          # Original sequences (optional)
--threshold FLOAT         # Quality threshold (default: 8.0)
--report-format FORMAT    # Report format: html, pdf, text (default: html)
--include-plots           # Include quality plots
```

## 🎨 Customization Scenarios

### Scenario 1: High-Throughput Analysis

**Objective:** Process large datasets quickly with moderate quality

**Configuration:**
```yaml
# config/high_throughput.yml
data:
  max_sequences: 100000
  preprocessing:
    quality_threshold: 0.7
    filter_repetitive: false

alignment:
  method: "center"
  max_iterations: 3

performance:
  threads: 16
  batch_size: 5000
  use_parallel: true
  use_streaming: true

output:
  reports:
    include_plots: false
```

**Command:**
```bash
Rscript scripts/build_pwm_robust.R \
  --config config/high_throughput.yml \
  --sequences data/large_dataset.fa \
  --threads 16 \
  --batch-size 5000
```

### Scenario 2: Publication-Quality Analysis

**Objective:** Maximum quality for research publication

**Configuration:**
```yaml
# config/publication.yml
data:
  preprocessing:
    quality_threshold: 0.95
    complexity_threshold: 0.8
    gc_content_range: [0.4, 0.6]

alignment:
  method: "integrated"
  max_iterations: 20
  refinement:
    enable: true
    max_refinement_cycles: 5

pwm:
  min_total_ic: 15.0
  min_conserved_positions: 8

statistics:
  bootstrap_iterations: 2000
  null_model_iterations: 5000

output:
  reports:
    generate_html: true
    generate_pdf: true
    include_plots: true
```

### Scenario 3: Resource-Constrained Environment

**Objective:** Run on limited computational resources

**Configuration:**
```yaml
# config/minimal.yml
performance:
  threads: 2
  memory_limit: "2G"
  batch_size: 50
  use_streaming: true
  cache_intermediate: false

data:
  max_sequences: 5000

output:
  formats: ["meme"]
  reports:
    generate_html: false
    include_plots: false
```

### Scenario 4: Method Comparison Study

**Objective:** Compare different alignment methods systematically

**Script:** `scripts/compare_methods.sh`
```bash
#!/bin/bash
# Compare alignment methods

SEQUENCES="data/input_sequences.fa"
METHODS=("center" "consensus" "integrated")
CONFIG="config/comparison.yml"

for method in "${METHODS[@]}"; do
    echo "Processing method: $method"
    
    # Align sequences
    Rscript scripts/advanced_alignment.R \
        --sequences "$SEQUENCES" \
        --method "$method" \
        --output "data/aligned_${method}.fa" \
        --config "$CONFIG"
    
    # Build PWM
    Rscript scripts/build_pwm_robust.R \
        --sequences "data/aligned_${method}.fa" \
        --output "results/pwm_${method}.meme" \
        --config "$CONFIG"
    
    # Validate quality
    Rscript scripts/validate_pwm_quality.R \
        --pwm "results/pwm_${method}.meme" \
        --output "results/quality_${method}.html" \
        --config "$CONFIG"
done

# Compare results
Rscript scripts/enhanced_compare_pwms.R \
    --pwms results/pwm_*.meme \
    --output results/method_comparison.html \
    --config "$CONFIG"
```

## 🔧 Advanced Configuration

### Custom Quality Metrics

**Define Custom Quality Function:**
```r
# In config/custom_quality.R
custom_quality_assessment <- function(pwm, sequences, config) {
  # Custom quality metrics
  custom_score <- calculate_custom_metric(pwm)
  
  # Combine with standard metrics
  standard_metrics <- standard_quality_assessment(pwm, sequences, config)
  
  return(list(
    standard = standard_metrics,
    custom = custom_score,
    combined = weight_metrics(standard_metrics, custom_score, config$weights)
  ))
}
```

**Configuration:**
```yaml
# config/custom.yml
quality:
  custom_function: "config/custom_quality.R"
  weights:
    information_content: 0.4
    conservation: 0.3
    biological_relevance: 0.2
    custom_metric: 0.1
```

### Conditional Processing

**Configuration with Conditions:**
```yaml
# config/adaptive.yml
data:
  preprocessing:
    # Adaptive quality threshold based on dataset size
    quality_threshold: 
      - condition: "sequence_count < 1000"
        value: 0.7
      - condition: "sequence_count >= 1000 && sequence_count < 10000"
        value: 0.8
      - condition: "sequence_count >= 10000"
        value: 0.9

alignment:
  # Method selection based on available memory
  method:
    - condition: "available_memory > 8G"
      value: "integrated"
    - condition: "available_memory > 4G"
      value: "consensus"
    - condition: "available_memory <= 4G"
      value: "center"
```

## 🔍 Configuration Validation

### Validation Script: `scripts/validate_config.R`

```r
# Validate configuration files
validate_configuration <- function(config_file) {
  config <- load_config(config_file)
  
  # Required parameters check
  required_params <- c(
    "data.preprocessing.quality_threshold",
    "alignment.method",
    "pwm.min_total_ic",
    "performance.threads"
  )
  
  missing <- check_required_parameters(config, required_params)
  if (length(missing) > 0) {
    stop("Missing required parameters: ", paste(missing, collapse = ", "))
  }
  
  # Parameter range validation
  validate_ranges(config)
  
  # Dependency validation
  validate_dependencies(config)
  
  return(TRUE)
}
```

### Configuration Testing

**Test Configuration:**
```bash
# Validate configuration file
Rscript scripts/validate_config.R --config config/pipeline.yml

# Test run with configuration
Rscript scripts/build_pwm_robust.R \
  --config config/pipeline.yml \
  --dry-run \
  --verbose
```

---

## 📖 Next Reading

- **[User Guide](10-user-guide.md)** - Step-by-step usage instructions
- **[Testing & Validation](11-testing-validation.md)** - Quality assurance procedures
- **[Troubleshooting](14-troubleshooting.md)** - Common configuration issues

---

*This comprehensive configuration guide enables you to customize the CTCF PWM Testing Pipeline for any research scenario, computational environment, or quality requirement.*
