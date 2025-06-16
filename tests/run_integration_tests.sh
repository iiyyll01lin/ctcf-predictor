#!/bin/bash

#' Integration Testing Script for CTCF PWM Testing Pipeline
#' 
#' This script runs comprehensive integration tests to validate the entire
#' pipeline from input data to final results.

set -e  # Exit on any error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TEST_DATA_DIR="$PROJECT_ROOT/tests/test_data"
TEST_RESULTS_DIR="$PROJECT_ROOT/tests/test_results"
REFERENCE_DIR="$PROJECT_ROOT/tests/reference_outputs"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test configuration
QUICK_TEST_SEQUENCES=50
STANDARD_TEST_SEQUENCES=200
COMPREHENSIVE_TEST_SEQUENCES=1000

# Logging
LOG_FILE="$TEST_RESULTS_DIR/integration_test.log"

# Initialize logging
log() {
    local level="$1"
    shift
    local message="$*"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    echo "[$timestamp] [$level] $message" | tee -a "$LOG_FILE"
    
    case "$level" in
        "ERROR")
            echo -e "${RED}✗ $message${NC}" >&2
            ;;
        "SUCCESS")
            echo -e "${GREEN}✓ $message${NC}"
            ;;
        "WARNING")
            echo -e "${YELLOW}⚠ $message${NC}"
            ;;
        "INFO")
            echo -e "$message"
            ;;
    esac
}

# Setup test environment
setup_test_environment() {
    log "INFO" "Setting up test environment..."
    
    # Create test directories
    mkdir -p "$TEST_DATA_DIR"
    mkdir -p "$TEST_RESULTS_DIR"
    mkdir -p "$REFERENCE_DIR"
    
    # Initialize log file
    echo "Integration Test Log - $(date)" > "$LOG_FILE"
    
    log "SUCCESS" "Test environment setup complete"
}

# Generate test data
generate_test_data() {
    local n_sequences="$1"
    local output_file="$2"
    
    log "INFO" "Generating $n_sequences test sequences..."
    
    # Check if R is available
    if ! command -v Rscript &> /dev/null; then
        log "ERROR" "Rscript not found. Please install R."
        return 1
    fi
    
    # Generate synthetic CTCF-like sequences
    cat > "$TEST_DATA_DIR/generate_test_sequences.R" << 'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript generate_test_sequences.R <n_sequences> <output_file>")
}

n_sequences <- as.numeric(args[1])
output_file <- args[2]

# Generate synthetic CTCF-like sequences
set.seed(42)  # For reproducible results

# CTCF consensus pattern
ctcf_core <- "CCGCGGGGGCGCT"
sequence_length <- 150

sequences <- character(n_sequences)

for (i in 1:n_sequences) {
    # Generate random background
    seq <- paste(sample(c("A", "C", "G", "T"), sequence_length, 
                       replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25)), 
                collapse = "")
    
    # Insert CTCF motif in 70% of sequences
    if (runif(1) < 0.7) {
        pos <- sample(1:(sequence_length - nchar(ctcf_core) + 1), 1)
        
        # Add some variation to the motif (1-2 mutations in 20% of cases)
        motif <- ctcf_core
        if (runif(1) < 0.2) {
            mut_pos <- sample(1:nchar(motif), sample(1:2, 1))
            for (mp in mut_pos) {
                substr(motif, mp, mp) <- sample(c("A", "C", "G", "T"), 1)
            }
        }
        
        substr(seq, pos, pos + nchar(motif) - 1) <- motif
    }
    
    sequences[i] <- paste0(">seq_", i, "\n", seq)
}

# Write FASTA file
writeLines(sequences, output_file)
cat("Generated", n_sequences, "sequences in", output_file, "\n")
EOF
    
    Rscript "$TEST_DATA_DIR/generate_test_sequences.R" "$n_sequences" "$output_file"
    
    if [[ -f "$output_file" ]]; then
        log "SUCCESS" "Test data generated: $output_file"
        return 0
    else
        log "ERROR" "Failed to generate test data"
        return 1
    fi
}

# Test data validation
test_data_validation() {
    log "INFO" "Testing data validation..."
    
    local test_file="$TEST_DATA_DIR/test_sequences.fa"
    local validation_script="$PROJECT_ROOT/scripts/validate_input_data.R"
    
    if [[ ! -f "$validation_script" ]]; then
        log "WARNING" "Data validation script not found: $validation_script"
        return 0
    fi
    
    if Rscript "$validation_script" --input "$test_file" --output "$TEST_RESULTS_DIR/data_validation.html" &>> "$LOG_FILE"; then
        log "SUCCESS" "Data validation test passed"
        return 0
    else
        log "ERROR" "Data validation test failed"
        return 1
    fi
}

# Test PWM construction
test_pwm_construction() {
    local method="$1"
    local input_file="$2"
    local output_file="$3"
    
    log "INFO" "Testing PWM construction with $method method..."
    
    # Check if main PWM construction scripts exist
    local main_script=""
    if [[ -f "$PROJECT_ROOT/scripts/build_pwm_robust.R" ]]; then
        main_script="$PROJECT_ROOT/scripts/build_pwm_robust.R"
    elif [[ -f "$PROJECT_ROOT/scripts/simple_aligned_pwm.R" ]]; then
        main_script="$PROJECT_ROOT/scripts/simple_aligned_pwm.R"
    else
        log "WARNING" "PWM construction script not found"
        return 0
    fi
    
    # Run PWM construction
    if Rscript "$main_script" --input "$input_file" --method "$method" --output "$output_file" &>> "$LOG_FILE"; then
        log "SUCCESS" "PWM construction test ($method) passed"
        return 0
    else
        log "ERROR" "PWM construction test ($method) failed"
        return 1
    fi
}

# Test bootstrap validation
test_bootstrap_validation() {
    log "INFO" "Testing bootstrap validation..."
    
    local bootstrap_script="$PROJECT_ROOT/scripts/bootstrap_validation.R"
    local input_file="$TEST_DATA_DIR/test_sequences.fa"
    
    if [[ ! -f "$bootstrap_script" ]]; then
        log "WARNING" "Bootstrap validation script not found: $bootstrap_script"
        return 0
    fi
    
    if Rscript "$bootstrap_script" --input "$input_file" --bootstrap-samples 100 --output "$TEST_RESULTS_DIR/bootstrap_validation.html" &>> "$LOG_FILE"; then
        log "SUCCESS" "Bootstrap validation test passed"
        return 0
    else
        log "ERROR" "Bootstrap validation test failed"
        return 1
    fi
}

# Test biological validation
test_biological_validation() {
    log "INFO" "Testing biological validation..."
    
    local bio_validation_script="$PROJECT_ROOT/scripts/biological_validation.R"
    local input_file="$TEST_DATA_DIR/test_sequences.fa"
    
    if [[ ! -f "$bio_validation_script" ]]; then
        log "WARNING" "Biological validation script not found: $bio_validation_script"
        return 0
    fi
    
    if Rscript "$bio_validation_script" --input "$input_file" --output "$TEST_RESULTS_DIR/biological_validation.html" &>> "$LOG_FILE"; then
        log "SUCCESS" "Biological validation test passed"
        return 0
    else
        log "ERROR" "Biological validation test failed"
        return 1
    fi
}

# Test performance benchmarking
test_performance_benchmark() {
    log "INFO" "Testing performance benchmarking..."
    
    local benchmark_script="$PROJECT_ROOT/scripts/benchmark_performance.R"
    
    if [[ ! -f "$benchmark_script" ]]; then
        log "WARNING" "Performance benchmark script not found: $benchmark_script"
        return 0
    fi
    
    if Rscript "$benchmark_script" --dataset-sizes "50,100" --methods "center" --repetitions 3 --output "$TEST_RESULTS_DIR/performance_benchmark.html" &>> "$LOG_FILE"; then
        log "SUCCESS" "Performance benchmark test passed"
        return 0
    else
        log "ERROR" "Performance benchmark test failed"
        return 1
    fi
}

# Test full pipeline integration
test_full_pipeline() {
    log "INFO" "Testing full pipeline integration..."
    
    local input_file="$TEST_DATA_DIR/test_sequences.fa"
    local config_file="$PROJECT_ROOT/config/pipeline.yml"
    local pipeline_script=""
    
    # Find main pipeline script
    if [[ -f "$PROJECT_ROOT/scripts/run_pipeline.R" ]]; then
        pipeline_script="$PROJECT_ROOT/scripts/run_pipeline.R"
    elif [[ -f "$PROJECT_ROOT/scripts/main_pipeline.R" ]]; then
        pipeline_script="$PROJECT_ROOT/scripts/main_pipeline.R"
    else
        log "WARNING" "Main pipeline script not found"
        return 0
    fi
    
    if Rscript "$pipeline_script" --input "$input_file" --config "$config_file" --output "$TEST_RESULTS_DIR/full_pipeline_result.rds" &>> "$LOG_FILE"; then
        log "SUCCESS" "Full pipeline integration test passed"
        return 0
    else
        log "ERROR" "Full pipeline integration test failed"
        return 1
    fi
}

# Run syntax checks
test_syntax_checks() {
    log "INFO" "Running syntax checks on all R scripts..."
    
    local r_scripts=($(find "$PROJECT_ROOT/scripts" -name "*.R" -type f))
    local syntax_errors=0
    
    for script in "${r_scripts[@]}"; do
        if ! Rscript -e "source('$script')" &> /dev/null; then
            log "ERROR" "Syntax error in $(basename "$script")"
            ((syntax_errors++))
        fi
    done
    
    if [[ $syntax_errors -eq 0 ]]; then
        log "SUCCESS" "All R scripts passed syntax check"
        return 0
    else
        log "ERROR" "$syntax_errors scripts have syntax errors"
        return 1
    fi
}

# Main integration test runner
run_integration_tests() {
    local test_level="$1"
    
    log "INFO" "Starting integration tests (level: $test_level)..."
    
    # Determine test parameters based on level
    local n_sequences
    case "$test_level" in
        "quick")
            n_sequences=$QUICK_TEST_SEQUENCES
            ;;
        "standard")
            n_sequences=$STANDARD_TEST_SEQUENCES
            ;;
        "comprehensive")
            n_sequences=$COMPREHENSIVE_TEST_SEQUENCES
            ;;
        *)
            log "ERROR" "Invalid test level: $test_level. Use: quick, standard, comprehensive"
            return 1
            ;;
    esac
    
    local test_file="$TEST_DATA_DIR/test_sequences.fa"
    local total_tests=0
    local passed_tests=0
    local failed_tests=0
    
    # Test 1: Generate test data
    ((total_tests++))
    if generate_test_data "$n_sequences" "$test_file"; then
        ((passed_tests++))
    else
        ((failed_tests++))
    fi
    
    # Test 2: Syntax checks
    ((total_tests++))
    if test_syntax_checks; then
        ((passed_tests++))
    else
        ((failed_tests++))
    fi
    
    # Test 3: Data validation
    ((total_tests++))
    if test_data_validation; then
        ((passed_tests++))
    else
        ((failed_tests++))
    fi
    
    # Test 4: PWM construction methods
    for method in center consensus integrated; do
        ((total_tests++))
        if test_pwm_construction "$method" "$test_file" "$TEST_RESULTS_DIR/pwm_${method}.rds"; then
            ((passed_tests++))
        else
            ((failed_tests++))
        fi
    done
    
    # Test 5: Bootstrap validation
    ((total_tests++))
    if test_bootstrap_validation; then
        ((passed_tests++))
    else
        ((failed_tests++))
    fi
    
    # Test 6: Biological validation
    ((total_tests++))
    if test_biological_validation; then
        ((passed_tests++))
    else
        ((failed_tests++))
    fi
    
    # Test 7: Performance benchmarking (only for standard and comprehensive)
    if [[ "$test_level" != "quick" ]]; then
        ((total_tests++))
        if test_performance_benchmark; then
            ((passed_tests++))
        else
            ((failed_tests++))
        fi
    fi
    
    # Test 8: Full pipeline integration (only for comprehensive)
    if [[ "$test_level" == "comprehensive" ]]; then
        ((total_tests++))
        if test_full_pipeline; then
            ((passed_tests++))
        else
            ((failed_tests++))
        fi
    fi
    
    # Generate summary
    log "INFO" "============== INTEGRATION TEST SUMMARY =============="
    log "INFO" "Test Level: $test_level"
    log "INFO" "Total Tests: $total_tests"
    log "INFO" "Passed: $passed_tests"
    log "INFO" "Failed: $failed_tests"
    log "INFO" "Success Rate: $(bc -l <<< "scale=1; $passed_tests * 100 / $total_tests")%"
    log "INFO" "=================================================="
    
    if [[ $failed_tests -eq 0 ]]; then
        log "SUCCESS" "All integration tests passed!"
        return 0
    else
        log "ERROR" "$failed_tests integration tests failed!"
        return 1
    fi
}

# Show usage
show_usage() {
    echo "Usage: $0 [quick|standard|comprehensive]"
    echo ""
    echo "Test levels:"
    echo "  quick        - Basic tests with small dataset (50 sequences)"
    echo "  standard     - Standard tests with medium dataset (200 sequences)"
    echo "  comprehensive - Full tests with large dataset (1000 sequences)"
    echo ""
    echo "Examples:"
    echo "  $0 quick"
    echo "  $0 standard"
    echo "  $0 comprehensive"
}

# Main execution
main() {
    local test_level="${1:-standard}"
    
    case "$test_level" in
        "-h"|"--help"|"help")
            show_usage
            exit 0
            ;;
        "quick"|"standard"|"comprehensive")
            ;;
        *)
            log "ERROR" "Invalid test level: $test_level"
            show_usage
            exit 1
            ;;
    esac
    
    # Setup environment
    setup_test_environment
    
    # Run tests
    if run_integration_tests "$test_level"; then
        exit 0
    else
        exit 1
    fi
}

# Run main if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
