# Comprehensive PWM Testing Script with Null Model Integration
# 這個腳本運行所有PWM改進方法並提供完整的測試工作流程，現在包含null model分析

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

# Phase 1: Initial Analysis and Null Model Generation
echo "PHASE 1: INITIAL ANALYSIS & NULL MODEL GENERATION"
echo "=================================================="

echo "1.1 Sequence Quality Analysis"
if [ ! -f "$RESULTS_DIR/sequence_quality_analysis.txt" ]; then
    run_with_docker "analyze_sequence_quality.R" "Sequence Quality Analysis"
fi

echo "1.2 PWM Quality Validation"
if [ ! -f "$RESULTS_DIR/pwm_quality_report.txt" ]; then
    run_with_docker "validate_pwm_quality.R" "PWM Quality Validation"
fi

echo "1.3 Sequence Alignment Analysis"
if [ ! -f "$DATA_DIR/aligned_sequences.fasta" ]; then
    run_with_docker "analyze_sequence_alignment.R" "Sequence Alignment Analysis"
fi

echo "1.4 Null Model Generation"
echo "Generating null models for statistical baseline..."
run_with_docker "generate_null_models.R" "Null Model Generation" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/null_models" "100"

# Check if null models were generated successfully
if [ -f "$RESULTS_DIR/null_models/null_summary_statistics.rds" ]; then
    echo "✓ Null model generation completed successfully"
    echo "  Generated null model baselines for statistical testing"
else
    echo "⚠️  Null model generation failed - proceeding without statistical testing"
fi

echo "PHASE 1 COMPLETED"
echo ""

# Phase 2: PWM Building Methods
echo "PHASE 2: PWM BUILDING METHODS"
echo "=============================="

echo "2.1 Simple Aligned PWM"
run_with_docker "simple_aligned_pwm.R" "Simple Aligned PWM" "$DATA_DIR/aligned_sequences.fasta" "$RESULTS_DIR/simple_aligned_pwm.rds" "0.1"

echo "2.2 High-Quality Subset PWMs"
run_with_docker "build_subset_pwm.R" "High-Quality Subset PWMs" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/subset_pwm" "5000" "0.01"

echo "2.3 Efficient Aligned PWM"
run_with_docker "efficient_aligned_pwm.R" "Efficient Aligned PWM" "$DATA_DIR/aligned_sequences.fasta" "$RESULTS_DIR/efficient_aligned_pwm" "10000" "TRUE"

echo "2.4 Robust PWM Building"
run_with_docker "build_pwm_robust.R" "Robust PWM Building" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/robust_pwm" "100"

echo "2.5 Advanced Alignment Methods"
# Try different alignment methods
for method in "consensus" "length" "progressive"; do
    echo "2.5.$method Advanced Alignment - $method"
    run_with_docker "advanced_alignment.R" "Advanced Alignment ($method)" "$DATA_DIR/training_sequences.fasta" "$RESULTS_DIR/advanced_${method}" "$method" "0.5"
done

echo "PHASE 2 COMPLETED"
echo ""

# Phase 3: Enhanced PWM Comparison and Statistical Analysis
echo "PHASE 3: ENHANCED PWM COMPARISON & STATISTICAL ANALYSIS"
echo "========================================================"

echo "3.1 Standard PWM Comparison Report"
run_with_docker "compare_pwms.R" "PWM Comparison Analysis" "$RESULTS_DIR" "$RESULTS_DIR/pwm_comparison_report.html"

echo "3.2 Enhanced PWM Comparison with Null Model Analysis"
run_with_docker "enhanced_compare_pwms.R" "Enhanced PWM Comparison with Statistical Testing" "$RESULTS_DIR" "$RESULTS_DIR/enhanced_pwm_comparison_report.html" "$RESULTS_DIR/null_models"

echo "3.3 Detailed Statistical Significance Testing"
if [ -f "$RESULTS_DIR/null_models/null_summary_statistics.rds" ]; then
    run_with_docker "statistical_significance_test.R" "Statistical Significance Analysis" "$RESULTS_DIR" "$RESULTS_DIR/null_models" "$RESULTS_DIR/statistical_significance_report.html"
    echo "✓ Comprehensive statistical analysis completed"
else
    echo "⚠️  Skipping detailed statistical testing - null models not available"
fi

echo "PHASE 3 COMPLETED"
echo ""

# Phase 4: Results Summary and Validation
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
    "$RESULTS_DIR/efficient_aligned_pwm_report.txt"
    "$RESULTS_DIR/advanced_consensus_basic.rds"
    "$RESULTS_DIR/advanced_length_basic.rds"
    "$RESULTS_DIR/pwm_comparison_report.html"
    "$RESULTS_DIR/enhanced_pwm_comparison_report.html"
)

success_count=0
total_files=${#expected_files[@]}

for file in "${expected_files[@]}"; do
    if check_file "$file" "PWM Result"; then
        ((success_count++))
    fi
done

echo ""
echo "File Check Summary: $success_count/$total_files core files generated successfully"

# Check null model and statistical files
echo ""
echo "Statistical Analysis Files:"
if [ -f "$RESULTS_DIR/null_models/null_summary_statistics.rds" ]; then
    echo "✓ Null model statistics: $RESULTS_DIR/null_models/null_summary_statistics.rds"
    null_files=$(ls "$RESULTS_DIR/null_models/"*.fasta 2>/dev/null | wc -l)
    echo "✓ Null model sequences: $null_files FASTA files"
else
    echo "✗ Null model statistics not found"
fi

if [ -f "$RESULTS_DIR/statistical_significance_report.html" ]; then
    echo "✓ Statistical significance report: $RESULTS_DIR/statistical_significance_report.html"
else
    echo "✗ Statistical significance report not generated"
fi

if [ -f "$RESULTS_DIR/enhanced_pwm_comparison_report.html" ]; then
    echo "✓ Enhanced comparison report: $RESULTS_DIR/enhanced_pwm_comparison_report.html"
else
    echo "✗ Enhanced comparison report not generated"
fi

echo ""
echo "========================================="
echo "FINAL TESTING SUMMARY"
echo "========================================="

echo ""
echo "📊 Generated Analysis Reports:"
echo "=============================="
echo "1. Standard PWM Comparison:"
if [ -f "$RESULTS_DIR/pwm_comparison_report.html" ]; then
    echo "   ✓ $RESULTS_DIR/pwm_comparison_report.html"
else
    echo "   ✗ Standard comparison report missing"
fi

echo "2. Enhanced PWM Comparison with Statistical Analysis:"
if [ -f "$RESULTS_DIR/enhanced_pwm_comparison_report.html" ]; then
    echo "   ✓ $RESULTS_DIR/enhanced_pwm_comparison_report.html"
    echo "   📈 Includes null model baselines and statistical significance"
else
    echo "   ✗ Enhanced comparison report missing"
fi

echo "3. Detailed Statistical Significance Analysis:"
if [ -f "$RESULTS_DIR/statistical_significance_report.html" ]; then
    echo "   ✓ $RESULTS_DIR/statistical_significance_report.html"
    echo "   📊 Includes p-values, effect sizes, and interpretations"
else
    echo "   ✗ Statistical significance report missing"
fi

echo ""
echo "🧬 Generated PWM Files:"
echo "======================"
pwm_count=$(ls -1 "$RESULTS_DIR"/*pwm*.rds 2>/dev/null | wc -l)
echo "  Total PWM files: $pwm_count"
ls -la "$RESULTS_DIR"/*pwm*.rds 2>/dev/null | head -10

if [ $pwm_count -gt 10 ]; then
    echo "  ... and $((pwm_count - 10)) more PWM files"
fi

echo ""
echo "📈 Statistical Validation Components:"
echo "====================================="
if [ -f "$RESULTS_DIR/null_models/null_summary_statistics.rds" ]; then
    null_types=$(Rscript -e "
        null_data <- readRDS('$RESULTS_DIR/null_models/null_summary_statistics.rds')
        cat('Null model types:', length(null_data), '\\n')
        cat('Types:', paste(names(null_data), collapse=', '), '\\n')
    " 2>/dev/null)
    echo "  ✓ Null model baselines generated"
    echo "  $null_types"
    
    # Count null model replicates
    null_files=$(ls "$RESULTS_DIR/null_models/"*.fasta 2>/dev/null | wc -l)
    echo "  ✓ Null model replicates: $null_files sequence sets"
else
    echo "  ✗ Null model baselines not available"
    echo "  ⚠️  Statistical significance testing was not performed"
fi

echo ""
echo "🎯 Key Quality Metrics Summary:"
echo "=============================="

# Extract key metrics if enhanced comparison data exists
if [ -f "$RESULTS_DIR/enhanced_pwm_comparison_report_data.rds" ]; then
    echo "✓ Enhanced comparison data available"
    
    # Try to extract and display key findings
    Rscript -e "
        tryCatch({
            data <- readRDS('$RESULTS_DIR/enhanced_pwm_comparison_report_data.rds')
            metrics <- data\$comparison_results\$metrics_df
            
            cat('Best PWM Performance:\\n')
            best_total <- metrics[which.max(metrics\$total_info), ]
            cat('  Highest total information:', best_total\$name, '-', round(best_total\$total_info, 3), 'bits\\n')
            
            best_conserved <- metrics[which.max(metrics\$conserved_1bit), ]
            cat('  Most conserved positions:', best_conserved\$name, '-', best_conserved\$conserved_1bit, 'positions\\n')
            
            # Statistical significance summary
            if (!is.null(data\$significance_results)) {
                sig_count <- 0
                for (pwm in names(data\$significance_results)) {
                    pwm_sig <- data\$significance_results[[pwm]]
                    if (!is.null(pwm_sig)) {
                        for (null_type in names(pwm_sig)) {
                            sig_data <- pwm_sig[[null_type]]
                            if (any(sig_data\$total_info\$is_significant,
                                    sig_data\$conserved_positions\$is_significant,
                                    sig_data\$avg_info\$is_significant)) {
                                sig_count <- sig_count + 1
                                break
                            }
                        }
                    }
                }
                cat('\\nStatistical Significance:\\n')
                cat('  PWMs with significant improvement:', sig_count, 'out of', length(data\$significance_results), '\\n')
            }
        }, error = function(e) {
            cat('Error extracting metrics summary\\n')
        })
    " 2>/dev/null
else
    echo "✗ Enhanced comparison data not available"
fi

echo ""
echo "📋 Next Steps and Recommendations:"
echo "=================================="
echo "1. 📊 Open the Enhanced PWM Comparison Report:"
echo "   File: $RESULTS_DIR/enhanced_pwm_comparison_report.html"
echo "   Contains: PWM performance vs null model baselines"

echo "2. 📈 Review Statistical Significance Analysis:"
echo "   File: $RESULTS_DIR/statistical_significance_report.html"
echo "   Contains: Detailed p-values and effect sizes"

echo "3. 🎯 Select Best Performing PWM:"
echo "   Based on: Information content, conserved positions, and statistical significance"
echo "   Use for: CTCF binding site prediction tasks"

echo "4. 🔬 Consider Further Optimization:"
echo "   If needed: Adjust parameters based on statistical analysis results"
echo "   Focus on: PWMs showing significant improvement over null models"

echo "5. 📝 Document Findings:"
echo "   Archive: PWM files and analysis reports for reproducibility"
echo "   Share: Statistical validation results with research team"

echo ""
if [ $success_count -eq $total_files ] && [ -f "$RESULTS_DIR/enhanced_pwm_comparison_report.html" ]; then
    echo "🎉 ALL TESTS COMPLETED SUCCESSFULLY WITH STATISTICAL VALIDATION!"
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
