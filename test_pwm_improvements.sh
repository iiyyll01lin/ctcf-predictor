# Comprehensive PWM Testing Script
# This script runs all PWM improvement methods and provides complete testing workflow

# Set script directory
SCRIPT_DIR="scripts"
DATA_DIR="data"
RESULTS_DIR="results"

# Create results directory if it doesn't exist
mkdir -p "$RESULTS_DIR"

echo "========================================"
echo "CTCF PWM Improvement Testing Workflow"
echo "========================================"
echo "Starting comprehensive PWM testing..."
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

# Function to check if file exists
check_file() {
    local file=$1
    local description=$2
    
    if [ -f "$file" ]; then
        echo "✓ $description exists: $file"
        return 0
    else
        echo "✗ $description missing: $file"
        return 1
    fi
}

# Phase 1: Initial Analysis (if not already done)
echo "PHASE 1: INITIAL ANALYSIS"
echo "=========================="

if [ ! -f "$RESULTS_DIR/sequence_quality_analysis.txt" ]; then
    run_with_docker "analyze_sequence_quality.R" "Sequence Quality Analysis"
fi

if [ ! -f "$RESULTS_DIR/pwm_quality_report.txt" ]; then
    run_with_docker "validate_pwm_quality.R" "PWM Quality Validation"
fi

if [ ! -f "$DATA_DIR/aligned_sequences.fasta" ]; then
    run_with_docker "analyze_sequence_alignment.R" "Sequence Alignment Analysis"
fi

echo "PHASE 1 COMPLETED"
echo ""

# Phase 2: PWM Building Methods
echo "PHASE 2: PWM BUILDING METHODS"
echo "=============================="

# 2.1 Simple Aligned PWM
echo "2.1 Simple Aligned PWM"
run_with_docker "simple_aligned_pwm.R" "Simple Aligned PWM" "$DATA_DIR/aligned_sequences.fasta" "$RESULTS_DIR/simple_aligned_pwm.rds" "0.1"

# 2.2 High-Quality Subset PWMs
echo "2.2 High-Quality Subset PWMs"
run_with_docker "build_subset_pwm.R" "High-Quality Subset PWMs" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/subset_pwm" "5000" "0.01"

# 2.3 Efficient Aligned PWM
echo "2.3 Efficient Aligned PWM"
run_with_docker "efficient_aligned_pwm.R" "Efficient Aligned PWM" "$DATA_DIR/aligned_sequences.fasta" "$RESULTS_DIR/efficient_aligned_pwm" "10000" "TRUE"

# 2.4 Robust PWM Building
echo "2.4 Robust PWM Building"
run_with_docker "build_pwm_robust.R" "Robust PWM Building" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/robust_pwm" "100"

# 2.5 Advanced Alignment Methods
echo "2.5 Advanced Alignment Methods"

# Try different alignment methods
for method in "consensus" "length" "progressive"; do
    echo "2.5.$method Advanced Alignment - $method"
    run_with_docker "advanced_alignment.R" "Advanced Alignment ($method)" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/advanced_${method}" "$method" "0.5"
done

echo "PHASE 2 COMPLETED"
echo ""

# Phase 3: PWM Comparison and Analysis
echo "PHASE 3: PWM COMPARISON AND ANALYSIS"
echo "====================================="

echo "3.1 PWM Comparison Report"
run_with_docker "compare_pwms.R" "PWM Comparison Analysis" "$RESULTS_DIR" "$RESULTS_DIR/pwm_comparison_report.html"

echo "PHASE 3 COMPLETED"
echo ""

# Phase 4: Results Summary
echo "PHASE 4: RESULTS SUMMARY"
echo "========================="

echo "Checking generated files..."

# Check key output files
declare -a expected_files=(
    "$RESULTS_DIR/simple_aligned_pwm.rds"
    "$RESULTS_DIR/subset_pwm_size1000.rds"
    "$RESULTS_DIR/subset_pwm_size2000.rds"
    "$RESULTS_DIR/subset_pwm_size5000.rds"
    "$RESULTS_DIR/efficient_aligned_pwm.rds"
    "$RESULTS_DIR/efficient_aligned_pwm_report.txt"
    "$RESULTS_DIR/advanced_consensus_basic.rds"
    "$RESULTS_DIR/advanced_length_basic.rds"
    "$RESULTS_DIR/pwm_comparison_report.html"
)

success_count=0
total_files=${#expected_files[@]}

for file in "${expected_files[@]}"; do
    if check_file "$file" "PWM Result"; then
        ((success_count++))
    fi
done

echo ""
echo "File Check Summary: $success_count/$total_files files generated successfully"

# Generate final summary
echo ""
echo "========================================="
echo "FINAL TESTING SUMMARY"
echo "========================================="

if [ -f "$RESULTS_DIR/pwm_comparison_report.html" ]; then
    echo "✓ Complete PWM comparison report generated"
    echo "  View at: $RESULTS_DIR/pwm_comparison_report.html"
else
    echo "✗ PWM comparison report not generated"
fi

echo ""
echo "Generated PWM Files:"
ls -la "$RESULTS_DIR"/*.rds 2>/dev/null | wc -l | xargs echo "  PWM files:"
ls -la "$RESULTS_DIR"/*pwm*.rds 2>/dev/null || echo "  No PWM files found"

echo ""
echo "Analysis Reports:"
ls -la "$RESULTS_DIR"/*.txt "$RESULTS_DIR"/*.html 2>/dev/null || echo "  No report files found"

echo ""
echo "Quality Metrics Summary:"
echo "========================"

# Extract key metrics if comparison report exists
if [ -f "$RESULTS_DIR/pwm_comparison_report_data.rds" ]; then
    echo "✓ Detailed comparison data available"
    echo "  File: $RESULTS_DIR/pwm_comparison_report_data.rds"
else
    echo "✗ Detailed comparison data not available"
fi

echo ""
echo "Next Steps:"
echo "==========="
echo "1. Open the HTML comparison report: $RESULTS_DIR/pwm_comparison_report.html"
echo "2. Review the detailed analysis for best-performing PWM"
echo "3. Use the best PWM for CTCF binding site prediction"
echo "4. Consider further optimization based on results"

echo ""
if [ $success_count -eq $total_files ]; then
    echo "🎉 ALL TESTS COMPLETED SUCCESSFULLY!"
    exit 0
else
    echo "⚠️  SOME TESTS FAILED OR FILES MISSING"
    echo "Please check the error messages above"
    exit 1
fi
