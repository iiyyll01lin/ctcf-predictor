#!/bin/bash

# Run Parallel Alignment Methods
# Executes multiple alignment methods in parallel and compares results
# Author: CTCF PWM Testing Pipeline Team

set -euo pipefail

# Default parameters
INPUT_FILE=""
OUTPUT_DIR="results/parallel_alignment"
METHODS="center,consensus,integrated"
CONFIG_FILE=""
THREADS=4
VERBOSE=false
CLEANUP=true

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local color=$1
    local message=$2
    echo -e "${color}[$(date +'%Y-%m-%d %H:%M:%S')] ${message}${NC}"
}

# Function to show usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Run multiple alignment methods in parallel for PWM testing pipeline.

OPTIONS:
    -i, --input FILE        Input FASTA file with sequences (required)
    -o, --output DIR        Output directory [default: results/parallel_alignment]
    -m, --methods LIST      Comma-separated alignment methods [default: center,consensus,integrated]
    -c, --config FILE       Configuration file
    -t, --threads NUM       Number of parallel threads [default: 4]
    -v, --verbose           Enable verbose output
    --no-cleanup            Don't remove temporary files
    -h, --help              Show this help message

ALIGNMENT METHODS:
    center       - Center-based alignment
    consensus    - Consensus-driven alignment
    integrated   - Hybrid approach combining methods

EXAMPLES:
    $0 -i data/sequences.fasta
    $0 -i data/sequences.fasta -o results/my_alignment -m center,consensus
    $0 -i data/sequences.fasta -t 8 -v

EOF
}

# Function to check if required tools are available
check_dependencies() {
    local missing_deps=()
    
    if ! command -v Rscript &> /dev/null; then
        missing_deps+=("Rscript")
    fi
    
    if ! command -v parallel &> /dev/null; then
        print_status "$YELLOW" "GNU parallel not found. Will use background processes instead."
    fi
    
    if [ ${#missing_deps[@]} -ne 0 ]; then
        print_status "$RED" "Missing dependencies: ${missing_deps[*]}"
        exit 1
    fi
}

# Function to validate input parameters
validate_inputs() {
    if [ -z "$INPUT_FILE" ]; then
        print_status "$RED" "Input file is required"
        usage
        exit 1
    fi
    
    if [ ! -f "$INPUT_FILE" ]; then
        print_status "$RED" "Input file does not exist: $INPUT_FILE"
        exit 1
    fi
    
    if [ ! -r "$INPUT_FILE" ]; then
        print_status "$RED" "Input file is not readable: $INPUT_FILE"
        exit 1
    fi
    
    # Check if input file is not empty
    if [ ! -s "$INPUT_FILE" ]; then
        print_status "$RED" "Input file is empty: $INPUT_FILE"
        exit 1
    fi
    
    # Validate threads parameter
    if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [ "$THREADS" -lt 1 ]; then
        print_status "$RED" "Invalid threads parameter: $THREADS"
        exit 1
    fi
}

# Function to setup output directory
setup_output_dir() {
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR/alignments"
    mkdir -p "$OUTPUT_DIR/logs"
    mkdir -p "$OUTPUT_DIR/reports"
    
    print_status "$BLUE" "Output directory: $OUTPUT_DIR"
}

# Function to run alignment method
run_alignment_method() {
    local method=$1
    local input_file=$2
    local output_dir=$3
    local config_file=$4
    local verbose=$5
    
    local method_output_dir="$output_dir/alignments/$method"
    local log_file="$output_dir/logs/${method}.log"
    local error_file="$output_dir/logs/${method}.error"
    
    mkdir -p "$method_output_dir"
    
    print_status "$BLUE" "Starting alignment method: $method"
    
    # Determine R script based on method
    local script_file=""
    case $method in
        "center")
            script_file="scripts/align_center.R"
            ;;
        "consensus")
            script_file="scripts/align_consensus.R"
            ;;
        "integrated")
            script_file="scripts/align_integrated.R"
            ;;
        *)
            print_status "$RED" "Unknown alignment method: $method"
            return 1
            ;;
    esac
    
    # Check if script exists
    if [ ! -f "$script_file" ]; then
        print_status "$YELLOW" "Script not found: $script_file. Using alignment fallback."
        script_file="scripts/alignment_fallback.R"
        
        if [ ! -f "$script_file" ]; then
            print_status "$RED" "Fallback script not found: $script_file"
            return 1
        fi
    fi
    
    # Build R command
    local r_cmd="Rscript $script_file -i \"$input_file\" -o \"$method_output_dir/aligned_sequences.fasta\""
    
    if [ -n "$config_file" ]; then
        r_cmd="$r_cmd -c \"$config_file\""
    fi
    
    if [ "$method" = "center" ] || [ "$method" = "consensus" ] || [ "$method" = "integrated" ]; then
        r_cmd="$r_cmd -m \"$method\""
    fi
    
    if [ "$verbose" = true ]; then
        r_cmd="$r_cmd -v"
    fi
    
    # Execute command
    local start_time=$(date +%s)
    
    if eval "$r_cmd" > "$log_file" 2> "$error_file"; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        
        print_status "$GREEN" "Completed $method alignment in ${duration}s"
        
        # Save method metadata
        cat > "$method_output_dir/metadata.json" << EOF
{
  "method": "$method",
  "input_file": "$input_file",
  "start_time": "$start_time",
  "end_time": "$end_time",
  "duration_seconds": $duration,
  "status": "success",
  "log_file": "$log_file"
}
EOF
        
        return 0
    else
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        
        print_status "$RED" "Failed $method alignment after ${duration}s"
        
        # Save failure metadata
        cat > "$method_output_dir/metadata.json" << EOF
{
  "method": "$method",
  "input_file": "$input_file",
  "start_time": "$start_time",
  "end_time": "$end_time",
  "duration_seconds": $duration,
  "status": "failed",
  "log_file": "$log_file",
  "error_file": "$error_file"
}
EOF
        
        return 1
    fi
}

# Function to run methods in parallel using background processes
run_methods_parallel_bg() {
    local methods_array=($1)
    local input_file=$2
    local output_dir=$3
    local config_file=$4
    local verbose=$5
    
    local pids=()
    local method_status=()
    
    # Start background jobs
    for method in "${methods_array[@]}"; do
        run_alignment_method "$method" "$input_file" "$output_dir" "$config_file" "$verbose" &
        local pid=$!
        pids+=($pid)
        method_status+=("$method:$pid")
        
        # Limit concurrent processes
        if [ ${#pids[@]} -ge $THREADS ]; then
            wait ${pids[0]}
            pids=("${pids[@]:1}")  # Remove first element
        fi
    done
    
    # Wait for remaining processes
    for pid in "${pids[@]}"; do
        wait $pid
    done
    
    print_status "$BLUE" "All alignment methods completed"
}

# Function to run methods using GNU parallel
run_methods_parallel_gnu() {
    local methods_array=($1)
    local input_file=$2
    local output_dir=$3
    local config_file=$4
    local verbose=$5
    
    # Create function for parallel execution
    export -f run_alignment_method print_status
    export INPUT_FILE="$input_file"
    export OUTPUT_DIR="$output_dir"
    export CONFIG_FILE="$config_file"
    export VERBOSE="$verbose"
    export RED GREEN YELLOW BLUE NC
    
    # Run in parallel
    printf '%s\n' "${methods_array[@]}" | \
    parallel -j $THREADS run_alignment_method {} "$input_file" "$output_dir" "$config_file" "$verbose"
    
    print_status "$BLUE" "All alignment methods completed (GNU parallel)"
}

# Function to generate comparison report
generate_comparison_report() {
    local output_dir=$1
    local methods_array=($2)
    
    print_status "$BLUE" "Generating comparison report..."
    
    local report_file="$output_dir/reports/alignment_comparison.json"
    local summary_file="$output_dir/reports/alignment_summary.txt"
    
    # Start JSON report
    echo "{" > "$report_file"
    echo "  \"timestamp\": \"$(date -Iseconds)\"," >> "$report_file"
    echo "  \"input_file\": \"$INPUT_FILE\"," >> "$report_file"
    echo "  \"methods\": [" >> "$report_file"
    
    # Start summary report
    cat > "$summary_file" << EOF
CTCF PWM Pipeline - Parallel Alignment Results
==============================================

Input File: $INPUT_FILE
Timestamp: $(date)
Methods Tested: ${methods_array[*]}

Results Summary:
---------------
EOF
    
    local first=true
    local successful_methods=()
    local failed_methods=()
    
    for method in "${methods_array[@]}"; do
        local metadata_file="$output_dir/alignments/$method/metadata.json"
        
        if [ ! "$first" = true ]; then
            echo "," >> "$report_file"
        fi
        first=false
        
        if [ -f "$metadata_file" ]; then
            cat "$metadata_file" >> "$report_file"
            
            # Extract status for summary
            local status=$(grep '"status"' "$metadata_file" | cut -d'"' -f4)
            local duration=$(grep '"duration_seconds"' "$metadata_file" | cut -d':' -f2 | tr -d ' ,')
            
            if [ "$status" = "success" ]; then
                successful_methods+=("$method")
                printf "✓ %-12s Success (${duration}s)\n" "$method" >> "$summary_file"
            else
                failed_methods+=("$method")
                printf "✗ %-12s Failed  (${duration}s)\n" "$method" >> "$summary_file"
            fi
        else
            echo "    {\"method\": \"$method\", \"status\": \"unknown\"}" >> "$report_file"
            failed_methods+=("$method")
            printf "? %-12s Unknown\n" "$method" >> "$summary_file"
        fi
    done
    
    # Close JSON report
    echo "" >> "$report_file"
    echo "  ]" >> "$report_file"
    echo "}" >> "$report_file"
    
    # Complete summary report
    cat >> "$summary_file" << EOF

Summary:
--------
Successful: ${#successful_methods[@]} (${successful_methods[*]})
Failed:     ${#failed_methods[@]} (${failed_methods[*]})

Output Files:
------------
Alignments: $output_dir/alignments/
Logs:       $output_dir/logs/
Reports:    $output_dir/reports/

EOF
    
    print_status "$GREEN" "Comparison report generated: $report_file"
    print_status "$GREEN" "Summary report generated: $summary_file"
    
    # Display summary
    cat "$summary_file"
}

# Function to cleanup temporary files
cleanup_temp_files() {
    if [ "$CLEANUP" = true ]; then
        print_status "$BLUE" "Cleaning up temporary files..."
        
        # Remove empty log files
        find "$OUTPUT_DIR/logs" -name "*.log" -size 0 -delete 2>/dev/null || true
        find "$OUTPUT_DIR/logs" -name "*.error" -size 0 -delete 2>/dev/null || true
        
        print_status "$BLUE" "Cleanup completed"
    fi
}

# Main execution function
main() {
    print_status "$BLUE" "Starting parallel alignment pipeline..."
    
    # Check dependencies
    check_dependencies
    
    # Validate inputs
    validate_inputs
    
    # Setup output directory
    setup_output_dir
    
    # Parse methods
    IFS=',' read -ra METHODS_ARRAY <<< "$METHODS"
    
    print_status "$BLUE" "Running alignment methods: ${METHODS_ARRAY[*]}"
    print_status "$BLUE" "Using $THREADS parallel threads"
    
    # Record start time
    local overall_start_time=$(date +%s)
    
    # Run alignment methods in parallel
    if command -v parallel &> /dev/null; then
        run_methods_parallel_gnu "${METHODS_ARRAY[*]}" "$INPUT_FILE" "$OUTPUT_DIR" "$CONFIG_FILE" "$VERBOSE"
    else
        run_methods_parallel_bg "${METHODS_ARRAY[*]}" "$INPUT_FILE" "$OUTPUT_DIR" "$CONFIG_FILE" "$VERBOSE"
    fi
    
    # Record end time
    local overall_end_time=$(date +%s)
    local overall_duration=$((overall_end_time - overall_start_time))
    
    # Generate comparison report
    generate_comparison_report "$OUTPUT_DIR" "${METHODS_ARRAY[*]}"
    
    # Cleanup
    cleanup_temp_files
    
    print_status "$GREEN" "Parallel alignment completed in ${overall_duration}s"
    print_status "$GREEN" "Results available in: $OUTPUT_DIR"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -m|--methods)
            METHODS="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
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
        --no-cleanup)
            CLEANUP=false
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            print_status "$RED" "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Run main function
main
