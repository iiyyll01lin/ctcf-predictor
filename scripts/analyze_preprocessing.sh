#!/bin/bash

echo "=== CTCF Sequence Preprocessing Analysis ==="
echo "Analyzing sequence characteristics to optimize preprocessing parameters"
echo ""

# Basic file analysis
echo "1. INPUT FILE ANALYSIS:"
echo "   File: data/extracted_sequences.fasta"
echo "   Total lines: $(wc -l < data/extracted_sequences.fasta)"
echo "   Total sequences: $(grep -c '^>' data/extracted_sequences.fasta)"
echo ""

# Extract some example lengths from coordinates
echo "2. LENGTH ANALYSIS (from first 10 sequences):"
head -20 data/extracted_sequences.fasta | grep '^>' | head -10 | while read line; do
    # Extract coordinates and calculate length
    coords=$(echo "$line" | cut -d':' -f2)
    if [[ "$coords" == *"-"* ]]; then
        start=$(echo "$coords" | cut -d'-' -f1)
        end=$(echo "$coords" | cut -d'-' -f2)
        length=$((end - start + 1))
        echo "   $line -> ${length}bp"
    fi
done
echo ""

echo "3. CURRENT PREPROCESSING RESULTS:"
if [ -f "data/preprocessed_sequences.fasta" ]; then
    input_count=$(grep -c '^>' data/extracted_sequences.fasta)
    output_count=$(grep -c '^>' data/preprocessed_sequences.fasta)
    retention_rate=$(echo "scale=2; $output_count * 100 / $input_count" | bc -l)
    echo "   Input: $input_count sequences"
    echo "   Output: $output_count sequences"
    echo "   Retention rate: ${retention_rate}%"
else
    echo "   No preprocessed file found yet"
fi
echo ""

echo "4. OPTIMIZED PREPROCESSING RESULTS:"
if [ -f "data/preprocessed_sequences_optimized.fasta" ]; then
    input_count=$(grep -c '^>' data/extracted_sequences.fasta)
    output_count=$(grep -c '^>' data/preprocessed_sequences_optimized.fasta)
    retention_rate=$(echo "scale=2; $output_count * 100 / $input_count" | bc -l)
    echo "   Input: $input_count sequences"
    echo "   Output: $output_count sequences"
    echo "   Retention rate: ${retention_rate}%"
else
    echo "   Optimized preprocessing not run yet"
fi
echo ""

echo "5. PARAMETER COMPARISON:"
echo "   Current config (preprocess_config.json):"
echo "     max_length: $(grep -o '"max_length": [0-9]*' scripts/preprocess_config.json | cut -d' ' -f2)"
echo "     max_n_percent: $(grep -o '"max_n_percent": [0-9]*' scripts/preprocess_config.json | cut -d' ' -f2)"
echo "     low_complexity_filter: $(grep -o '"low_complexity_filter": [a-z]*' scripts/preprocess_config.json | cut -d' ' -f2)"
echo "     mask_repeats: $(grep -o '"mask_repeats": [a-z]*' scripts/preprocess_config.json | cut -d' ' -f2)"
echo ""

echo "6. RECOMMENDATIONS:"
echo "   - CTCF binding sites are typically 200-300bp long"
echo "   - Current max_length=100 filters out ~99% of sequences"
echo "   - Recommended max_length=300 for initial training"
echo "   - Consider relaxing complexity filters for more data"
echo "   - Target retention rate: 10-30% for sufficient training data"
echo ""

echo "=== Analysis Complete ==="
