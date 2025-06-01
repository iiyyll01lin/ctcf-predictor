# Sequence Analysis and Parameter Optimization Guide

This guide explains how to use the sequence analysis script to optimize your preprocessing parameters for the CTCF binding site prediction pipeline.

## Step 1: Run the Sequence Analysis

Execute the analysis script to understand your sequence characteristics:

```bash
# Run basic analysis (output to console)
./run-in-docker.sh Rscript scripts/analyze_sequences.R data/extracted_sequences.fasta

# Run analysis and save detailed report to file
./run-in-docker.sh Rscript scripts/analyze_sequences.R data/extracted_sequences.fasta results/sequence_analysis_report.txt
```

## Step 2: Interpret the Analysis Results

The analysis provides 6 key sections:

### 1. Sequence Length Analysis
- **What it shows**: Distribution of sequence lengths in your dataset
- **Key metrics**: Average, median, min/max lengths
- **What to look for**: 
  - If most sequences are >100bp, increase `max_length`
  - Check what percentage would be retained with different length limits

### 2. N Base Analysis  
- **What it shows**: Percentage of undefined nucleotides (N) in sequences
- **Key metrics**: How many sequences contain N bases and their percentages
- **What to look for**:
  - If many sequences have >15% N bases, increase `max_n_percent`
  - High N content may indicate lower quality sequences

### 3. Sequence Complexity Analysis
- **What it shows**: Entropy measurements indicating sequence complexity
- **Key metrics**: Distribution of entropy values (low entropy = repetitive/simple)
- **What to look for**:
  - If many sequences have low entropy, consider disabling complexity filtering
  - Entropy <1.5 indicates very repetitive sequences

### 4. Repeat Sequence Analysis
- **What it shows**: Presence of homopolymer and dinucleotide repeats
- **Key metrics**: Percentage of sequences containing repeats
- **What to look for**:
  - High repeat content may justify repeat masking
  - Consider if repeats are biologically relevant for CTCF binding

### 5. Preprocessing Parameter Optimization
- **What it shows**: Simulation of different parameter combinations
- **Key metrics**: Number and percentage of sequences retained under each scenario
- **What to look for**:
  - Scenarios that retain 1000+ sequences for reliable model training
  - Balance between data quality and quantity

### 6. Recommendations
- **What it shows**: Specific parameter suggestions based on your data
- **Key output**: Recommended config.json updates

## Step 3: Optimize Your Parameters

Based on the analysis results, follow these guidelines:

### Target Number of Sequences
- **Minimum**: 200-500 sequences (basic PWM training)
- **Good**: 1000-2000 sequences (reliable training)
- **Ideal**: 3000+ sequences (robust model)

### Parameter Adjustment Strategy

1. **Start with length filtering only**:
   ```json
   {
     "length_filter": true,
     "min_length": 11,
     "max_length": 200,  // Adjust based on analysis
     "low_complexity_filter": false,
     "mask_repeats": false,
     "max_n_percent": 50  // Very permissive initially
   }
   ```

2. **Gradually add restrictions**:
   - If you have >5000 sequences, tighten N base filtering
   - If you have >3000 sequences, consider enabling complexity filtering
   - Only enable repeat masking if biological justification exists

3. **Monitor the trade-offs**:
   - More filtering = higher quality but fewer sequences
   - Less filtering = more sequences but potentially lower quality

## Step 4: Test Your Optimized Parameters

1. **Create optimized config** (already provided as `preprocess_config_optimized.json`):
   ```bash
   # Test with optimized parameters
   ./run-in-docker.sh Rscript scripts/preprocess_sequences.R data/extracted_sequences.fasta data/preprocessed_sequences_optimized.fasta scripts/preprocess_config_optimized.json
   ```

2. **Compare results**:
   ```bash
   # Check how many sequences were retained
   ./run-in-docker.sh Rscript -e "
   library(Biostrings)
   orig <- readDNAStringSet('data/extracted_sequences.fasta')
   proc <- readDNAStringSet('data/preprocessed_sequences_optimized.fasta')
   cat('Original sequences:', length(orig), '\n')
   cat('Processed sequences:', length(proc), '\n')
   cat('Retention rate:', round(length(proc)/length(orig)*100, 2), '%\n')
   "
   ```

## Step 5: Iterative Optimization

If the results aren't satisfactory:

1. **Too few sequences retained**:
   - Increase `max_length` (try 300, 500, or even 1000)
   - Increase `max_n_percent` (try 30%, 40%)
   - Disable `low_complexity_filter`
   - Disable `mask_repeats`

2. **Too many low-quality sequences**:
   - Decrease `max_n_percent`
   - Enable `low_complexity_filter` with `entropy_threshold: 1.0`
   - Enable `mask_repeats`

3. **Re-run analysis after changes**:
   ```bash
   # Analyze the preprocessed sequences
   ./run-in-docker.sh Rscript scripts/analyze_sequences.R data/preprocessed_sequences_optimized.fasta results/processed_analysis.txt
   ```

## Step 6: Validate Model Performance

After finding optimal parameters:

1. **Proceed with the pipeline**:
   ```bash
   # Continue with dataset preparation
   ./run-in-docker.sh Rscript scripts/prepare_datasets.R
   
   # Build PWM model
   ./run-in-docker.sh Rscript scripts/build_pwm.R
   
   # Evaluate model performance
   ./run-in-docker.sh Rscript scripts/evaluate_models.R
   ```

2. **Monitor model quality**:
   - Check PWM matrix for reasonable base frequencies
   - Evaluate ROC curves and AUC values
   - If performance is poor, consider adjusting parameters

## Common Parameter Combinations

### Conservative (High Quality)
```json
{
  "max_length": 150,
  "max_n_percent": 10,
  "low_complexity_filter": true,
  "entropy_threshold": 1.5,
  "mask_repeats": true
}
```

### Balanced (Good Trade-off)
```json
{
  "max_length": 200,
  "max_n_percent": 20,
  "low_complexity_filter": false,
  "mask_repeats": false
}
```

### Permissive (Maximum Retention)
```json
{
  "max_length": 500,
  "max_n_percent": 30,
  "low_complexity_filter": false,
  "mask_repeats": false
}
```

## Troubleshooting

### Issue: Analysis script runs too slowly
**Solution**: The script automatically samples large datasets. For very large files (>10K sequences), it uses a sample of 1000 sequences for complexity analysis.

### Issue: Very few sequences pass any filter
**Solution**: Your dataset might have unusual characteristics. Try the "Length_Only" scenario or investigate the source data quality.

### Issue: All sequences are similar length
**Solution**: Your sequences might be pre-processed. Check if they're already standardized to a specific length.

### Issue: High N content across all sequences
**Solution**: This might indicate masking from previous processing steps. Consider if these N's represent repetitive elements that should be preserved.

Remember: The goal is to balance data quality with quantity sufficient for reliable PWM model training. Start permissive and gradually add restrictions based on your specific dataset characteristics.
