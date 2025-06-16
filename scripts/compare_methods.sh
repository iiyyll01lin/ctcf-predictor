#!/bin/bash
# Method Comparison Script
# Compare different alignment methods systematically

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default parameters
SEQUENCES=""
METHODS=("center" "consensus" "integrated")
CONFIG="config/comparison.yml"
OUTPUT_DIR="results/method_comparison"
THREADS=4
VERBOSE=false

# Function to print usage
print_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Compare different alignment methods systematically for CTCF PWM construction.

OPTIONS:
    -s, --sequences FILE    Input sequences file (required)
    -m, --methods LIST      Comma-separated list of methods to compare
                           Available: center,consensus,integrated
                           Default: center,consensus,integrated
    -c, --config FILE       Configuration file (default: config/comparison.yml)
    -o, --output DIR        Output directory (default: results/method_comparison)
    -t, --threads N         Number of threads (default: 4)
    -v, --verbose           Verbose output
    -h, --help              Show this help message

EXAMPLES:
    $0 -s data/ctcf_sequences.fa
    $0 -s data/input.fa -m "center,integrated" -t 8
    $0 -s data/input.fa -c config/custom_comparison.yml -o results/my_comparison

EOF
}

# Function to log messages
log_message() {
    local level=$1
    shift
    local message="$*"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo -e "${GREEN}[INFO]${NC} ${timestamp} - $message"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} ${timestamp} - $message"
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} ${timestamp} - $message"
            ;;
        "DEBUG")
            if [ "$VERBOSE" = true ]; then
                echo -e "${BLUE}[DEBUG]${NC} ${timestamp} - $message"
            fi
            ;;
    esac
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -s|--sequences)
                SEQUENCES="$2"
                shift 2
                ;;
            -m|--methods)
                IFS=',' read -ra METHODS <<< "$2"
                shift 2
                ;;
            -c|--config)
                CONFIG="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            -h|--help)
                print_usage
                exit 0
                ;;
            *)
                echo "Unknown option: $1"
                print_usage
                exit 1
                ;;
        esac
    done
}

# Validate inputs
validate_inputs() {
    log_message "INFO" "Validating inputs..."
    
    # Check if sequences file exists
    if [ -z "$SEQUENCES" ]; then
        log_message "ERROR" "Input sequences file is required"
        exit 1
    fi
    
    if [ ! -f "$SEQUENCES" ]; then
        log_message "ERROR" "Sequences file not found: $SEQUENCES"
        exit 1
    fi
    
    # Check if config file exists
    if [ ! -f "$CONFIG" ]; then
        log_message "ERROR" "Configuration file not found: $CONFIG"
        exit 1
    fi
    
    # Validate methods
    valid_methods=("center" "consensus" "integrated")
    for method in "${METHODS[@]}"; do
        if [[ ! " ${valid_methods[@]} " =~ " ${method} " ]]; then
            log_message "ERROR" "Invalid method: $method. Valid methods: ${valid_methods[*]}"
            exit 1
        fi
    done
    
    # Check if R is available
    if ! command -v Rscript &> /dev/null; then
        log_message "ERROR" "Rscript not found. Please install R."
        exit 1
    fi
    
    # Check if required R scripts exist
    required_scripts=("scripts/advanced_alignment.R" "scripts/build_pwm_robust.R" "scripts/validate_pwm_quality.R")
    for script in "${required_scripts[@]}"; do
        if [ ! -f "$script" ]; then
            log_message "ERROR" "Required script not found: $script"
            exit 1
        fi
    done
    
    log_message "INFO" "Input validation completed successfully"
}

# Create output directory structure
setup_output_directory() {
    log_message "INFO" "Setting up output directory: $OUTPUT_DIR"
    
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR/alignments"
    mkdir -p "$OUTPUT_DIR/pwms"
    mkdir -p "$OUTPUT_DIR/quality_reports"
    mkdir -p "$OUTPUT_DIR/plots"
    mkdir -p "$OUTPUT_DIR/logs"
    
    # Create comparison log file
    COMPARISON_LOG="$OUTPUT_DIR/logs/comparison.log"
    echo "CTCF Method Comparison Log" > "$COMPARISON_LOG"
    echo "Started: $(date)" >> "$COMPARISON_LOG"
    echo "Input sequences: $SEQUENCES" >> "$COMPARISON_LOG"
    echo "Methods: ${METHODS[*]}" >> "$COMPARISON_LOG"
    echo "Configuration: $CONFIG" >> "$COMPARISON_LOG"
    echo "----------------------------------------" >> "$COMPARISON_LOG"
}

# Process single method
process_method() {
    local method=$1
    local start_time=$(date +%s)
    
    log_message "INFO" "Processing method: $method"
    
    # Define output files
    local aligned_file="$OUTPUT_DIR/alignments/aligned_${method}.fa"
    local pwm_file="$OUTPUT_DIR/pwms/pwm_${method}.meme"
    local quality_file="$OUTPUT_DIR/quality_reports/quality_${method}.html"
    local method_log="$OUTPUT_DIR/logs/${method}.log"
    
    # Step 1: Align sequences
    log_message "DEBUG" "Aligning sequences with method: $method"
    if ! Rscript scripts/advanced_alignment.R \
        --sequences "$SEQUENCES" \
        --method "$method" \
        --output "$aligned_file" \
        --config "$CONFIG" \
        --threads "$THREADS" \
        --log-file "$method_log" 2>&1; then
        
        log_message "ERROR" "Alignment failed for method: $method"
        echo "FAILED: Alignment" >> "$COMPARISON_LOG"
        return 1
    fi
    
    # Step 2: Build PWM
    log_message "DEBUG" "Building PWM for method: $method"
    if ! Rscript scripts/build_pwm_robust.R \
        --sequences "$aligned_file" \
        --output "$pwm_file" \
        --config "$CONFIG" \
        --threads "$THREADS" \
        --log-file "$method_log" 2>&1; then
        
        log_message "ERROR" "PWM construction failed for method: $method"
        echo "FAILED: PWM construction" >> "$COMPARISON_LOG"
        return 1
    fi
    
    # Step 3: Validate quality
    log_message "DEBUG" "Validating PWM quality for method: $method"
    if ! Rscript scripts/validate_pwm_quality.R \
        --pwm "$pwm_file" \
        --sequences "$aligned_file" \
        --output "$quality_file" \
        --config "$CONFIG" \
        --include-plots \
        --log-file "$method_log" 2>&1; then
        
        log_message "ERROR" "Quality validation failed for method: $method"
        echo "FAILED: Quality validation" >> "$COMPARISON_LOG"
        return 1
    fi
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log_message "INFO" "Method $method completed successfully in ${duration}s"
    echo "SUCCESS: $method (${duration}s)" >> "$COMPARISON_LOG"
    
    return 0
}

# Process all methods
process_all_methods() {
    log_message "INFO" "Processing ${#METHODS[@]} methods: ${METHODS[*]}"
    
    local failed_methods=()
    local successful_methods=()
    
    for method in "${METHODS[@]}"; do
        if process_method "$method"; then
            successful_methods+=("$method")
        else
            failed_methods+=("$method")
        fi
    done
    
    # Report results
    if [ ${#successful_methods[@]} -gt 0 ]; then
        log_message "INFO" "Successfully processed methods: ${successful_methods[*]}"
    fi
    
    if [ ${#failed_methods[@]} -gt 0 ]; then
        log_message "WARN" "Failed methods: ${failed_methods[*]}"
    fi
    
    # Update comparison log
    echo "----------------------------------------" >> "$COMPARISON_LOG"
    echo "Successful methods: ${successful_methods[*]}" >> "$COMPARISON_LOG"
    echo "Failed methods: ${failed_methods[*]}" >> "$COMPARISON_LOG"
    
    return ${#failed_methods[@]}
}

# Generate comparison report
generate_comparison_report() {
    log_message "INFO" "Generating comprehensive comparison report"
    
    # Find all PWM files
    local pwm_files=()
    for method in "${METHODS[@]}"; do
        local pwm_file="$OUTPUT_DIR/pwms/pwm_${method}.meme"
        if [ -f "$pwm_file" ]; then
            pwm_files+=("$pwm_file")
        fi
    done
    
    if [ ${#pwm_files[@]} -eq 0 ]; then
        log_message "ERROR" "No PWM files found for comparison"
        return 1
    fi
    
    # Generate comparison report
    local comparison_report="$OUTPUT_DIR/method_comparison_report.html"
    log_message "DEBUG" "Creating comparison report: $comparison_report"
    
    if ! Rscript scripts/enhanced_compare_pwms.R \
        --pwms "${pwm_files[@]}" \
        --output "$comparison_report" \
        --config "$CONFIG" \
        --include-plots \
        --format "html" 2>&1; then
        
        log_message "ERROR" "Failed to generate comparison report"
        return 1
    fi
    
    log_message "INFO" "Comparison report generated: $comparison_report"
    return 0
}

# Generate summary statistics
generate_summary() {
    log_message "INFO" "Generating comparison summary"
    
    local summary_file="$OUTPUT_DIR/comparison_summary.txt"
    
    cat > "$summary_file" << EOF
CTCF Method Comparison Summary
==============================

Input file: $SEQUENCES
Methods compared: ${METHODS[*]}
Configuration: $CONFIG
Output directory: $OUTPUT_DIR
Threads used: $THREADS

Results:
--------
EOF
    
    # Add method-specific results
    for method in "${METHODS[@]}"; do
        local pwm_file="$OUTPUT_DIR/pwms/pwm_${method}.meme"
        local quality_file="$OUTPUT_DIR/quality_reports/quality_${method}.html"
        
        echo "" >> "$summary_file"
        echo "Method: $method" >> "$summary_file"
        
        if [ -f "$pwm_file" ]; then
            echo "  PWM file: $(basename "$pwm_file")" >> "$summary_file"
            echo "  Status: SUCCESS" >> "$summary_file"
            
            # Extract quality metrics if available
            if [ -f "$quality_file" ]; then
                echo "  Quality report: $(basename "$quality_file")" >> "$summary_file"
            fi
        else
            echo "  Status: FAILED" >> "$summary_file"
        fi
    done
    
    echo "" >> "$summary_file"
    echo "Generated: $(date)" >> "$summary_file"
    
    log_message "INFO" "Summary saved to: $summary_file"
}

# Main function
main() {
    echo -e "${GREEN}=== CTCF Method Comparison Pipeline ===${NC}"
    echo ""
    
    # Parse arguments
    parse_arguments "$@"
    
    # Validate inputs
    validate_inputs
    
    # Setup output directory
    setup_output_directory
    
    # Process all methods
    log_message "INFO" "Starting method comparison pipeline"
    local start_time=$(date +%s)
    
    if ! process_all_methods; then
        log_message "WARN" "Some methods failed, but continuing with comparison"
    fi
    
    # Generate comparison report
    if ! generate_comparison_report; then
        log_message "ERROR" "Failed to generate comparison report"
        exit 1
    fi
    
    # Generate summary
    generate_summary
    
    local end_time=$(date +%s)
    local total_duration=$((end_time - start_time))
    
    # Final report
    echo ""
    echo -e "${GREEN}=== Comparison Complete ===${NC}"
    echo "Total time: ${total_duration}s"
    echo "Output directory: $OUTPUT_DIR"
    echo "Comparison report: $OUTPUT_DIR/method_comparison_report.html"
    echo "Summary: $OUTPUT_DIR/comparison_summary.txt"
    echo ""
    
    log_message "INFO" "Method comparison pipeline completed successfully"
}

# Run main function with all arguments
main "$@"
