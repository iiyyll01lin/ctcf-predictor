#!/bin/bash

# This script downloads a CTCF ChIP-seq peak file (BED format)
# for the K562 cell line from the ENCODE project, reference genome
# (hg38 full or chr21), and then uses bedtools to extract sequences from the peaks.
#
# By default, this script will download the full hg38 reference genome for real analysis.
# You can use the -d or --demo flag to download only chr21 for demonstration purposes.
# The -f or --force flag can be used to force re-download of files even if they exist.
# Usage:
#   ./download_data.sh         # Download full genome (recommended for real analysis)
#   ./download_data.sh -d      # Download only chr21 (demonstration mode)
#   ./download_data.sh --demo  # Download only chr21 (demonstration mode)
#   ./download_data.sh -f      # Force re-download of all files
#   ./download_data.sh -d -f   # Force re-download in demo mode

# --- Function Definitions --- 
# These must be at the top before they're called

# Function to check file integrity
check_fasta_integrity() {
    local fasta_file="$1"
    # Check if file exists and is not empty
    if [ -s "$fasta_file" ]; then
        # For FASTA files, first line should start with ">"
        if head -n 1 "$fasta_file" | grep -q "^>"; then
            return 0 # Valid FASTA file
        else
            echo "ERROR: Invalid FASTA format in $fasta_file"
            return 1 # Invalid FASTA file
        fi
    else
        echo "ERROR: File $fasta_file does not exist or is empty"
        return 1
    fi
}

# Function to check file with MD5
verify_md5() {
    local file="$1"
    local expected_md5="$2"
    local actual_md5=""
    
    if [ -f "$file" ]; then
        if command -v md5sum > /dev/null; then
            actual_md5=$(md5sum "$file" | cut -d' ' -f1)
        elif command -v md5 > /dev/null; then
            # For MacOS
            actual_md5=$(md5 -q "$file")
        else
            echo "WARNING: md5sum/md5 not available. Skipping checksum verification."
            return 0
        fi
        
        if [ "$actual_md5" = "$expected_md5" ]; then
            echo "MD5 checksum verified for $file"
            return 0
        else
            echo "WARNING: MD5 checksum mismatch for $file"
            echo "  Expected: $expected_md5"
            echo "  Actual:   $actual_md5"
            return 1
        fi
    else
        echo "ERROR: File $file does not exist"
        return 1
    fi
}

# --- Configuration ---
# Process command line arguments
DEMO_MODE=false
FORCE=false

# Process all command-line arguments
for arg in "$@"; do
    case $arg in
        -d|--demo)
            DEMO_MODE=true
            ;;
        -f|--force)
            FORCE=true
            ;;
    esac
done

# Print active modes
if [ "$DEMO_MODE" == true ]; then
    echo "Running in demo mode (chr21 only)"
fi
if [ "$FORCE" == true ]; then
    echo "Force mode enabled - will re-download files even if they exist"
fi

# --- MD5 Checksum Constants ---
# These are the actual MD5 checksums from the files
if [ "$DEMO_MODE" == true ]; then
    # Demo mode - chr21 only
    EXPECTED_REF_GENOME_MD5="184df2bd9b812b6e6b6da16c6021369e" # Actual MD5 for chr21
else
    # Full mode
    EXPECTED_REF_GENOME_MD5="1c9dcaddfa41027f17cd8f7a82c7293b" # Actual MD5 for full hg38
fi
# ENCODE peak file expected checksum
EXPECTED_ENCODE_MD5="3f36655eb9a421a542ff031fc7b5e9aa" # Actual MD5 for ENCODE peaks

# Absolute path of this script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Store results inside the repo (<repo>/data)
OUTPUT_DIR="$SCRIPT_DIR/../data"
mkdir -p "$OUTPUT_DIR"
# ENCODE file details
ENCODE_URL="https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz"
ENCODE_FILENAME="K562_CTCF_peaks.bed.gz"
ENCODE_BED_PATH="$OUTPUT_DIR/${ENCODE_FILENAME%.gz}"
# Reference Genome details
REF_GENOME_DIR="$OUTPUT_DIR/reference_genome"

# Set reference genome details based on mode
if [ "$DEMO_MODE" == true ]; then
    # Demo mode - download only chr21
    REF_GENOME_FILENAME="hg38.chr21.fa.gz"
    REF_GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
    REF_GENOME_FASTA_PATH="$REF_GENOME_DIR/${REF_GENOME_FILENAME%.gz}"
else
    # Full analysis mode - download complete genome
    REF_GENOME_FILENAME="hg38.fa.gz"
    REF_GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    REF_GENOME_FASTA_PATH="$REF_GENOME_DIR/${REF_GENOME_FILENAME%.gz}"
fi

# Output Fasta file from bedtools
EXTRACTED_FASTA_PATH="$OUTPUT_DIR/extracted_sequences.fasta"

# --- Dependency Checks ---
echo "Checking dependencies..."
if ! command -v gunzip > /dev/null; then
    echo "Error: gunzip is required but not found. Please install it." 
    exit 1
fi
if ! command -v wget > /dev/null && ! command -v curl > /dev/null; then
    echo "Error: wget or curl is required for downloading but neither is found. Please install one."
    exit 1
fi
if ! command -v bedtools > /dev/null; then
    echo "Error: bedtools is required for sequence extraction but not found. Please install it (e.g., 'conda install bedtools' or 'sudo apt-get install bedtools')."
    exit 1
fi
echo "Dependencies found."

# --- Download ENCODE Peaks ---
# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

ENCODE_GZ_PATH="$OUTPUT_DIR/$ENCODE_FILENAME"
if [ ! -f "$ENCODE_BED_PATH" ] && [ ! -f "$ENCODE_GZ_PATH" ] || [ "$FORCE" == true ]; then
    # If force flag is set and files exist, remove them
    if [ "$FORCE" == true ]; then
        [ -f "$ENCODE_BED_PATH" ] && rm "$ENCODE_BED_PATH" && echo "Removed existing file: $ENCODE_BED_PATH"
        [ -f "$ENCODE_GZ_PATH" ] && rm "$ENCODE_GZ_PATH" && echo "Removed existing file: $ENCODE_GZ_PATH"
    fi
    
    echo -e "\nDownloading example CTCF ChIP-seq peaks for K562..."
    if command -v wget > /dev/null; then
        wget -O "$ENCODE_GZ_PATH" "$ENCODE_URL"
    elif command -v curl > /dev/null; then
        curl -L -o "$ENCODE_GZ_PATH" "$ENCODE_URL"
    fi
    # Check if download was successful
    if [ $? -ne 0 ]; then
        echo "Error: ENCODE Peak download failed."
        exit 1
    fi
    
    # Verify the downloaded file
    if ! verify_md5 "$ENCODE_GZ_PATH" "$EXPECTED_ENCODE_MD5"; then
        echo "WARNING: ENCODE peak file MD5 checksum does not match expected value."
        echo "The file may be corrupted or has been updated upstream."
        # Optional: uncomment to force exit on checksum failure
        # echo "Exiting due to checksum verification failure."
        # exit 1
    fi
else
    echo -e "\nENCODE peak file already exists (or was downloaded previously). Skipping download."
    
    # Optional: verify existing files if not forcing download
    if [ -f "$ENCODE_GZ_PATH" ]; then
        verify_md5 "$ENCODE_GZ_PATH" "$EXPECTED_ENCODE_MD5"
    elif [ -f "$ENCODE_BED_PATH" ]; then
        echo "Note: Cannot verify MD5 for uncompressed BED file. Consider using --force to re-download."
    fi
fi

# Decompress ENCODE peaks if necessary
if ([ -f "$ENCODE_GZ_PATH" ] && [ ! -f "$ENCODE_BED_PATH" ]) || ([ "$FORCE" == true ] && [ -f "$ENCODE_GZ_PATH" ]); then
    echo "Decompressing $ENCODE_FILENAME..."
    # If force and BED exists, remove it first
    if [ "$FORCE" == true ] && [ -f "$ENCODE_BED_PATH" ]; then
        rm "$ENCODE_BED_PATH" && echo "Removed existing file: $ENCODE_BED_PATH"
    fi
    
    gunzip -f "$ENCODE_GZ_PATH"
    if [ $? -ne 0 ]; then
        echo "Error: ENCODE Peak decompression failed."
        exit 1
    else
        echo "ENCODE Peak decompression complete: $ENCODE_BED_PATH"
        # Verify BED file is not empty
        if [ -s "$ENCODE_BED_PATH" ]; then
            echo "ENCODE BED file verified (non-empty)."
            # Optional: Add more specific BED format validation here
        else
            echo "WARNING: ENCODE BED file is empty after decompression."
            echo "Consider downloading it again with --force."
            exit 1
        fi
    fi
elif [ -f "$ENCODE_BED_PATH" ]; then
    echo "ENCODE Peak file already uncompressed: $ENCODE_BED_PATH"
    # Verify BED file is not empty
    if [ ! -s "$ENCODE_BED_PATH" ]; then
        echo "WARNING: Existing ENCODE BED file appears to be empty."
        echo "Consider using --force to re-download and decompress it."
    fi
fi

# --- Download Reference Genome (Example: hg38 chr21) ---
mkdir -p "$REF_GENOME_DIR"
REF_GENOME_GZ_PATH="$REF_GENOME_DIR/$REF_GENOME_FILENAME"

if [ ! -f "$REF_GENOME_FASTA_PATH" ] && [ ! -f "$REF_GENOME_GZ_PATH" ] || [ "$FORCE" == true ]; then
    # If force flag is set and files exist, remove them
    if [ "$FORCE" == true ]; then
        [ -f "$REF_GENOME_FASTA_PATH" ] && rm "$REF_GENOME_FASTA_PATH" && echo "Removed existing file: $REF_GENOME_FASTA_PATH"
        [ -f "$REF_GENOME_GZ_PATH" ] && rm "$REF_GENOME_GZ_PATH" && echo "Removed existing file: $REF_GENOME_GZ_PATH"
    fi
    
    echo -e "\nDownloading reference genome (${REF_GENOME_FILENAME%.gz})..."
    if command -v wget > /dev/null; then
        wget -O "$REF_GENOME_GZ_PATH" "$REF_GENOME_URL"
    elif command -v curl > /dev/null; then
        curl -L -o "$REF_GENOME_GZ_PATH" "$REF_GENOME_URL"
    fi
    # Check download success
    if [ $? -ne 0 ]; then
        echo "Error: Reference genome download failed. Sequence extraction cannot proceed without it."
        exit 1
    fi
    
    # Verify the downloaded file
    if ! verify_md5 "$REF_GENOME_GZ_PATH" "$EXPECTED_REF_GENOME_MD5"; then
        echo "WARNING: Reference genome file MD5 checksum does not match expected value."
        echo "The file may be corrupted or has been updated upstream."
        # Optional: uncomment to force exit on checksum failure
        # echo "Exiting due to checksum verification failure."
        # exit 1
    fi
else
    echo -e "\nReference genome file already exists (or was downloaded previously). Skipping download."
    
    # Optional: verify existing files if not forcing download
    if [ -f "$REF_GENOME_GZ_PATH" ]; then
        verify_md5 "$REF_GENOME_GZ_PATH" "$EXPECTED_REF_GENOME_MD5"
    elif [ -f "$REF_GENOME_FASTA_PATH" ]; then
        # Check if valid FASTA
        if check_fasta_integrity "$REF_GENOME_FASTA_PATH"; then
            echo "Reference genome FASTA file format verified."
        else
            echo "WARNING: Reference genome FASTA file appears to be invalid."
            echo "Consider using --force to re-download it."
        fi
    fi
fi

# Decompress reference genome if necessary
if ([ -f "$REF_GENOME_GZ_PATH" ] && [ ! -f "$REF_GENOME_FASTA_PATH" ]) || ([ "$FORCE" == true ] && [ -f "$REF_GENOME_GZ_PATH" ]); then
    echo "Decompressing $REF_GENOME_FILENAME..."
    # If force and FASTA exists, remove it first
    if [ "$FORCE" == true ] && [ -f "$REF_GENOME_FASTA_PATH" ]; then
        rm "$REF_GENOME_FASTA_PATH" && echo "Removed existing file: $REF_GENOME_FASTA_PATH"
    fi
    
    gunzip -f "$REF_GENOME_GZ_PATH"
    if [ $? -ne 0 ]; then
        echo "Error: Reference genome decompression failed."
        exit 1
    else
        echo "Reference genome decompression complete: $REF_GENOME_FASTA_PATH"
        # Verify FASTA file integrity
        if check_fasta_integrity "$REF_GENOME_FASTA_PATH"; then
            echo "Reference genome FASTA file format verified."
        else
            echo "WARNING: Reference genome FASTA file appears to be invalid after decompression."
            echo "Consider downloading it again with --force."
            exit 1
        fi
    fi
elif [ -f "$REF_GENOME_FASTA_PATH" ]; then
    echo "Reference genome file already uncompressed: $REF_GENOME_FASTA_PATH"
    # Verify FASTA file integrity if not already done
    if ! check_fasta_integrity "$REF_GENOME_FASTA_PATH"; then
        echo "WARNING: Existing reference genome FASTA file appears to be invalid."
        echo "Consider using --force to re-download and decompress it."
    fi
fi

# --- Extract Sequences using bedtools getfasta ---

# Check if both necessary files exist before attempting extraction
if [ -f "$ENCODE_BED_PATH" ] && [ -f "$REF_GENOME_FASTA_PATH" ]; then
    echo -e "\nExtracting sequences using bedtools getfasta..."
    echo "Input BED: $ENCODE_BED_PATH"
    echo "Reference FASTA: $REF_GENOME_FASTA_PATH"
    echo "Output FASTA: $EXTRACTED_FASTA_PATH"
    
    if [ "$DEMO_MODE" != true ]; then
        echo "Using full genome - extraction may take longer than with demo data"
    fi
    
    # getfasta command:
    # to extract DNA sequences corresponding to CTCF binding sites
    # Creates extracted_sequences.fasta containing actual genomic sequences from the ENCODE peak locations
    # Reports statistics on the number of sequences extracted (44,217 sequences in the full genome case)

    # Run bedtools getfasta
    # -fi: input FASTA file
    # -bed: input BED file
    # -fo: output FASTA file
    bedtools getfasta -fi "$REF_GENOME_FASTA_PATH" -bed "$ENCODE_BED_PATH" -fo "$EXTRACTED_FASTA_PATH"
    
    # Check if bedtools command was successful
    if [ $? -ne 0 ]; then
        echo "Error: bedtools getfasta command failed."
    else
        echo "Sequence extraction complete. Output saved to: $EXTRACTED_FASTA_PATH"
        echo "Note: This FASTA file can now be used as input for training/testing (after potential splitting and length filtering)."
        
        # Report statistics on the extracted sequences
        SEQ_COUNT=$(grep -c "^>" "$EXTRACTED_FASTA_PATH")
        echo "Extracted $SEQ_COUNT sequences from the reference genome."
    fi
else
    echo -e "\nSkipping sequence extraction: Required input files not found."
    if [ ! -f "$ENCODE_BED_PATH" ]; then echo "Missing: $ENCODE_BED_PATH"; fi
    if [ ! -f "$REF_GENOME_FASTA_PATH" ]; then echo "Missing: $REF_GENOME_FASTA_PATH"; fi
fi

# Provide next steps information
echo -e "\nScript finished."
if [ -f "$EXTRACTED_FASTA_PATH" ]; then
    echo -e "\nNext steps:"
    echo "1. [Optional] Preprocess the extracted sequences:"
    echo "   Rscript scripts/preprocess_sequences.R $EXTRACTED_FASTA_PATH data/preprocessed_sequences.fasta"
    echo "2. Prepare training and test datasets:"
    echo "   Rscript scripts/prepare_datasets.R"
    echo "3. Build the PWM model:"
    echo "   Rscript scripts/build_pwm.R"
    echo "4. Evaluate the model and find the optimal threshold"
fi

# Print rerun options
echo -e "\nRerun options:"
echo "- For full genome: ./scripts/download_data.sh"
echo "- For demo mode (chr21 only): ./scripts/download_data.sh -d"
echo "- To force re-download: ./scripts/download_data.sh -f"
echo "- To force re-download in demo mode: ./scripts/download_data.sh -d -f"
