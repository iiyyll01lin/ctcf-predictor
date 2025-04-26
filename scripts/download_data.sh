#!/bin/bash

# This script downloads an example CTCF ChIP-seq peak file (BED format)
# for the K562 cell line from the ENCODE project, an example reference genome
# (hg38 chr21), and then uses bedtools to extract sequences from the peaks.

# --- Configuration ---
# Output directory
OUTPUT_DIR="../data"
# ENCODE file details
ENCODE_URL="https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz"
ENCODE_FILENAME="K562_CTCF_peaks.bed.gz"
ENCODE_BED_PATH="$OUTPUT_DIR/${ENCODE_FILENAME%.gz}"
# Reference Genome details
REF_GENOME_DIR="$OUTPUT_DIR/reference_genome"
REF_GENOME_FILENAME="hg38.chr21.fa.gz"
REF_GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
REF_GENOME_FASTA_PATH="$REF_GENOME_DIR/${REF_GENOME_FILENAME%.gz}"
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
if [ ! -f "$ENCODE_BED_PATH" ] && [ ! -f "$ENCODE_GZ_PATH" ]; then
    echo "\nDownloading example CTCF ChIP-seq peaks for K562..."
    # ... existing download logic using wget/curl for ENCODE_GZ_PATH ...
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
else
    echo "\nENCODE peak file already exists (or was downloaded previously). Skipping download."
fi

# Decompress ENCODE peaks if necessary
if [ -f "$ENCODE_GZ_PATH" ] && [ ! -f "$ENCODE_BED_PATH" ]; then
    echo "Decompressing $ENCODE_FILENAME..."
    gunzip -f "$ENCODE_GZ_PATH"
    if [ $? -ne 0 ]; then
        echo "Error: ENCODE Peak decompression failed."
        exit 1
    else
        echo "ENCODE Peak decompression complete: $ENCODE_BED_PATH"
    fi
elif [ -f "$ENCODE_BED_PATH" ]; then
     echo "ENCODE Peak file already uncompressed: $ENCODE_BED_PATH"
fi

# --- Download Reference Genome (Example: hg38 chr21) ---
mkdir -p "$REF_GENOME_DIR"
REF_GENOME_GZ_PATH="$REF_GENOME_DIR/$REF_GENOME_FILENAME"

if [ ! -f "$REF_GENOME_FASTA_PATH" ] && [ ! -f "$REF_GENOME_GZ_PATH" ]; then
    echo "\nDownloading example reference genome (hg38 chr21)..."
    # ... existing download logic using wget/curl for REF_GENOME_GZ_PATH ...
    if command -v wget > /dev/null; then
        wget -O "$REF_GENOME_GZ_PATH" "$REF_GENOME_URL"
    elif command -v curl > /dev/null; then
        curl -L -o "$REF_GENOME_GZ_PATH" "$REF_GENOME_URL"
    fi
    # Check download success
    if [ $? -ne 0 ]; then
        echo "Warning: Reference genome download failed. Sequence extraction cannot proceed without it."
        # Decide whether to exit or continue
        # exit 1 
    fi
else
    echo "\nReference genome file already exists (or was downloaded previously). Skipping download."
fi

# Decompress reference genome if necessary
if [ -f "$REF_GENOME_GZ_PATH" ] && [ ! -f "$REF_GENOME_FASTA_PATH" ]; then
    echo "Decompressing $REF_GENOME_FILENAME..."
    gunzip -f "$REF_GENOME_GZ_PATH"
    if [ $? -ne 0 ]; then
        echo "Error: Reference genome decompression failed."
    else
        echo "Reference genome decompression complete: $REF_GENOME_FASTA_PATH"
    fi
elif [ -f "$REF_GENOME_FASTA_PATH" ]; then
    echo "Reference genome file already uncompressed: $REF_GENOME_FASTA_PATH"
fi

# --- Extract Sequences using bedtools getfasta ---

# Check if both necessary files exist before attempting extraction
if [ -f "$ENCODE_BED_PATH" ] && [ -f "$REF_GENOME_FASTA_PATH" ]; then
    echo "\nExtracting sequences using bedtools getfasta..."
    echo "Input BED: $ENCODE_BED_PATH"
    echo "Reference FASTA: $REF_GENOME_FASTA_PATH"
    echo "Output FASTA: $EXTRACTED_FASTA_PATH"
    
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
    fi
else
    echo "\nSkipping sequence extraction: Required input files not found."
    if [ ! -f "$ENCODE_BED_PATH" ]; then echo "Missing: $ENCODE_BED_PATH"; fi
    if [ ! -f "$REF_GENOME_FASTA_PATH" ]; then echo "Missing: $REF_GENOME_FASTA_PATH"; fi
fi

echo "\nScript finished."
