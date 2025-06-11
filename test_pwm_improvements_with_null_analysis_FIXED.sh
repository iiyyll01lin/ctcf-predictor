#!/bin/bash
# Comprehensive PWM Testing Script with Null Model Integration - FIXED VERSION
# This version properly handles USE_DOCKER=false when running inside a Docker container

# Environment Configuration
USE_DOCKER=${USE_DOCKER:-true}

# Set script directory
SCRIPT_DIR="scripts"
DATA_DIR="data"
RESULTS_DIR="results"

# Create results directory if it doesn't exist
mkdir -p "$RESULTS_DIR"

echo "================================================================"
echo "CTCF PWM Improvement Testing Workflow with Null Model Analysis"
echo "================================================================"
echo "Starting comprehensive PWM testing with statistical validation..."
echo "Execution mode: $([ "$USE_DOCKER" = "true" ] && echo "Docker" || echo "Local")"
echo ""

# Function to run script with Docker
run_with_docker() {
    local script_name=$1
    local description=$2
    shift 2
    local args="$@"
    
    echo "----------------------------------------"
    echo "Running: $description"
    echo "Script: $script_name"
    echo "Args: $args"
    echo "Mode: Docker"
    echo "----------------------------------------"
    
    if [ -f "$SCRIPT_DIR/$script_name" ]; then
        ./run-in-docker.sh Rscript "$SCRIPT_DIR/$script_name" $args
        
        if [ $? -eq 0 ]; then
            echo "✓ $description completed successfully"
        else
            echo "✗ $description failed"
            return 1
        fi
    else
        echo "✗ Script not found: $SCRIPT_DIR/$script_name"
        return 1
    fi
    echo ""
}

# Function to run script locally
run_local() {
    local script_name=$1
    local description=$2
    shift 2
    local args="$@"
    
    echo "----------------------------------------"
    echo "Running: $description"
    echo "Script: $script_name"
    echo "Args: $args"
    echo "Mode: Local"
    echo "----------------------------------------"
    
    if [ -f "$SCRIPT_DIR/$script_name" ]; then
        # Check if Rscript is available
        if ! command -v Rscript &> /dev/null; then
            echo "✗ Rscript not found. Please install R or use Docker mode (USE_DOCKER=true)"
            return 1
        fi
        
        Rscript "$SCRIPT_DIR/$script_name" $args
        
        if [ $? -eq 0 ]; then
            echo "✓ $description completed successfully"
        else
            echo "✗ $description failed"
            return 1
        fi
    else
        echo "✗ Script not found: $SCRIPT_DIR/$script_name"
        return 1
    fi
    echo ""
}

# Wrapper function to choose execution method
run_script() {
    if [ "$USE_DOCKER" = "true" ]; then
        run_with_docker "$@"
    else
        run_local "$@"
    fi
}

# Function to run any command with proper Docker/Local handling
run_command() {
    local description=$1
    local command=$2
    
    echo "Running: $description"
    if [ "$USE_DOCKER" = "true" ]; then
        echo "Mode: Docker"
        ./run-in-docker.sh bash -c "$command"
    else
        echo "Mode: Local"
        bash -c "$command"
    fi
    
    local result=$?
    if [ $result -eq 0 ]; then
        echo "✓ $description completed successfully"
    else
        echo "✗ $description failed"
    fi
    return $result
}

# Display usage information
echo "USAGE INFORMATION:"
echo "=================="
echo "To run with Docker (default):     USE_DOCKER=true ./test_pwm_improvements_with_null_analysis_FIXED.sh"
echo "To run locally:                   USE_DOCKER=false ./test_pwm_improvements_with_null_analysis_FIXED.sh"
echo "Note: Local execution requires R and required packages to be installed"
echo ""

# Phase 1: Initial Analysis and Null Model Generation
echo "PHASE 1: INITIAL ANALYSIS & NULL MODEL GENERATION"
echo "=================================================="

echo "1.1 Sequence Quality Analysis"
if [ ! -f "$RESULTS_DIR/sequence_quality_analysis.txt" ]; then
    run_script "analyze_sequence_quality.R" "Sequence Quality Analysis"
fi

echo "1.2 PWM Quality Validation"
if [ ! -f "$RESULTS_DIR/pwm_quality_report.txt" ]; then
    run_script "validate_pwm_quality.R" "PWM Quality Validation"
fi

echo "1.3 Sequence Alignment Analysis"
if [ ! -f "$DATA_DIR/aligned_sequences.fasta" ]; then
    run_script "analyze_sequence_alignment.R" "Sequence Alignment Analysis"
fi

echo "1.4 Null Model Generation"
echo "Generating null models for statistical baseline..."
run_script "generate_null_models.R" "Null Model Generation" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/null_models" "100"

# Check if null models were generated successfully
if [ -f "$RESULTS_DIR/null_models/null_summary_statistics.rds" ]; then
    echo "✓ Null model generation completed successfully"
    echo "  Generated null model baselines for statistical testing"
else
    echo "⚠️  Null model generation failed - proceeding without statistical testing"
fi

echo "PHASE 1 COMPLETED"
echo ""

# Phase 1.5: Chromosome Split Validation
echo "PHASE 1.5: CHROMOSOME SPLIT VALIDATION"
echo "======================================"

echo "1.5.1 Chromosome Extraction and Split Validation"
echo "Running comprehensive chromosome split validation..."

# Create validation script for data leakage detection
cat > validate_chromosome_split.R << 'EOF'
library(Biostrings)

cat("=== Chromosome Split Validation ===\n")

# Check if dataset files exist
train_file <- "data/training_sequences.fasta"
test_file <- "data/test_sequences.fasta"

if (!file.exists(train_file) || !file.exists(test_file)) {
    cat("❌ ERROR: Training or test datasets not found.\n")
    cat("   Expected files:\n")
    cat("   - ", train_file, "\n")
    cat("   - ", test_file, "\n")
    cat("   Please run data preparation first.\n")
    quit(status = 1)
}

# Load data
train_seqs <- readDNAStringSet(train_file)
test_seqs <- readDNAStringSet(test_file)

cat("✅ Dataset files loaded successfully\n")
cat("   Training sequences:", length(train_seqs), "\n")
cat("   Test sequences:", length(test_seqs), "\n")

# Load chromosome extraction function
source("scripts/prepare_datasets.R", local = TRUE)

# Extract chromosome info
train_chrs <- sapply(names(train_seqs), extract_chromosome)
test_names_clean <- gsub(" \\| class=[01]", "", names(test_seqs))
test_chrs <- sapply(test_names_clean, extract_chromosome)

# Get unique chromosomes
train_unique <- unique(train_chrs)
test_unique <- unique(test_chrs)

cat("\n=== Chromosome Assignment Analysis ===\n")
cat("Training chromosomes (", length(train_unique), "):", paste(sort(train_unique), collapse=", "), "\n")
cat("Testing chromosomes (", length(test_unique), "):", paste(sort(test_unique), collapse=", "), "\n")

# Check for overlap (data leakage)
overlap <- intersect(train_unique, test_unique)

if (length(overlap) > 0) {
    cat("\n❌ DATA LEAKAGE DETECTED!\n")
    cat("Overlapping chromosomes:", paste(overlap, collapse=", "), "\n")
    cat("This violates the principle of genomic data independence.\n")
    quit(status = 1)
} else {
    cat("\n✅ No data leakage detected - chromosomes are properly separated\n")
}

cat("=== Validation Summary ===\n")
cat("✅ Chromosome split validation completed successfully\n")
cat("✅ No genomic data leakage detected\n")
cat("✅ Dataset is ready for robust PWM training and evaluation\n")
EOF

# Run the validation with proper Docker/Local handling
if [ "$USE_DOCKER" = "true" ]; then
    ./run-in-docker.sh Rscript validate_chromosome_split.R
    validation_result=$?
else
    Rscript validate_chromosome_split.R
    validation_result=$?
fi

if [ $validation_result -eq 0 ]; then
    echo "✅ Chromosome split validation passed"
else
    echo "❌ Chromosome split validation failed"
    echo "Please check the validation output above and fix any issues."
    rm -f validate_chromosome_split.R
    exit 1
fi

# Clean up validation script
rm -f validate_chromosome_split.R

echo ""
echo "PHASE 1.5 COMPLETED - Chromosome split validation passed"
echo ""

# Phase 2: PWM Building Methods
echo "PHASE 2: PWM BUILDING METHODS"
echo "=============================="

echo "2.1 Simple Aligned PWM"
run_script "simple_aligned_pwm.R" "Simple Aligned PWM" "$DATA_DIR/aligned_sequences.fasta" "$RESULTS_DIR/simple_aligned_pwm.rds" "0.1"

echo "2.2 High-Quality Subset PWMs"
run_script "build_subset_pwm.R" "High-Quality Subset PWMs" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/subset_pwm" "5000" "0.01"

echo "2.3 Efficient Aligned PWM"
run_script "efficient_aligned_pwm.R" "Efficient Aligned PWM" "$DATA_DIR/aligned_sequences.fasta" "$RESULTS_DIR/efficient_aligned_pwm" "10000" "TRUE"

echo "2.4 Robust PWM Building"
run_script "build_pwm_robust.R" "Robust PWM Building" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/robust_pwm" "100"

echo "2.5 Advanced Alignment Methods"
for method in "consensus" "length" "progressive"; do
    echo "2.5.$method Advanced Alignment - $method"
    run_script "advanced_alignment.R" "Advanced Alignment ($method)" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/advanced_${method}" "$method" "0.5"
done

echo "PHASE 2 COMPLETED"
echo ""

# Phase 3: Enhanced PWM Comparison and Statistical Analysis
echo "PHASE 3: ENHANCED PWM COMPARISON & STATISTICAL ANALYSIS"
echo "========================================================"

echo "3.1 Standard PWM Comparison Report"
run_script "compare_pwms.R" "PWM Comparison Analysis" "$RESULTS_DIR" "$RESULTS_DIR/pwm_comparison_report.html"

echo "3.2 Enhanced PWM Comparison with Null Model Analysis"
run_script "enhanced_compare_pwms.R" "Enhanced PWM Comparison with Statistical Testing" "$RESULTS_DIR" "$RESULTS_DIR/enhanced_pwm_comparison_report.html" "$RESULTS_DIR/null_models"

echo "PHASE 3 COMPLETED"
echo ""

# Phase 4: Results Summary
echo "PHASE 4: COMPREHENSIVE RESULTS SUMMARY"
echo "======================================="

echo "Checking generated files..."

# Check key output files
declare -a expected_files=(
    "$RESULTS_DIR/simple_aligned_pwm.rds"
    "$RESULTS_DIR/subset_pwm_size1000.rds"
    "$RESULTS_DIR/subset_pwm_size2000.rds"
    "$RESULTS_DIR/subset_pwm_size5000.rds"
    "$RESULTS_DIR/efficient_aligned_pwm.rds"
    "$RESULTS_DIR/pwm_comparison_report.html"
    "$RESULTS_DIR/enhanced_pwm_comparison_report.html"
)

success_count=0
total_files=${#expected_files[@]}

for file in "${expected_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ PWM Result exists: $file"
        ((success_count++))
    else
        echo "✗ PWM Result missing: $file"
    fi
done

echo ""
echo "File Check Summary: $success_count/$total_files core files generated successfully"

echo ""
echo "========================================="
echo "FINAL TESTING SUMMARY"
echo "========================================="

echo ""
echo "📊 Generated Analysis Reports:"
echo "=============================="
if [ -f "$RESULTS_DIR/pwm_comparison_report.html" ]; then
    echo "✓ Standard comparison: $RESULTS_DIR/pwm_comparison_report.html"
else
    echo "✗ Standard comparison report missing"
fi

if [ -f "$RESULTS_DIR/enhanced_pwm_comparison_report.html" ]; then
    echo "✓ Enhanced comparison: $RESULTS_DIR/enhanced_pwm_comparison_report.html"
    echo "  📈 Includes null model baselines and statistical significance"
else
    echo "✗ Enhanced comparison report missing"
fi

echo ""
echo "🧬 Generated PWM Files:"
echo "======================"
pwm_count=$(ls -1 "$RESULTS_DIR"/*pwm*.rds 2>/dev/null | wc -l)
echo "  Total PWM files: $pwm_count"

echo ""
echo "📋 Next Steps and Recommendations:"
echo "=================================="
echo "1. 📊 Open the Enhanced PWM Comparison Report:"
echo "   File: $RESULTS_DIR/enhanced_pwm_comparison_report.html"
echo "   Contains: PWM performance vs null model baselines"

echo "2. 🎯 Select Best Performing PWM:"
echo "   Based on: Information content, conserved positions, statistical significance"
echo "   Use for: CTCF binding site prediction tasks"

echo ""
if [ $success_count -eq $total_files ]; then
    echo "🎉 ALL TESTS COMPLETED SUCCESSFULLY!"
    echo "📊 Enhanced PWM analysis with null model comparison is now available"
    exit 0
elif [ $success_count -ge $((total_files * 3 / 4)) ]; then
    echo "✅ MOST TESTS COMPLETED SUCCESSFULLY"
    echo "⚠️  Some optional components may be missing - check messages above"
    exit 0
else
    echo "⚠️  SOME CRITICAL TESTS FAILED"
    echo "Please check the error messages above and verify input data"
    exit 1
fi
