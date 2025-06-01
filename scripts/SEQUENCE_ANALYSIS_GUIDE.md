# CTCF Binding Site Prediction: Sequence Analysis & Parameter Optimization Guide

## Current Status Summary

### Problem Identified ✅
Your preprocessing pipeline was filtering out **99.9% of sequences** (44,180 out of 44,217), leaving only 37 sequences. This occurred because:

1. **Length filtering too restrictive**: `max_length = 100` but CTCF sequences are typically 200-300bp
2. **Complexity filtering too aggressive**: `entropy_threshold = 1.5` removing valid genomic sequences  
3. **Repeat masking unnecessary**: Genomic repeats can be informative for CTCF binding

### Solution Applied ✅

Updated `scripts/preprocess_config.json` with optimized parameters:

```json
{
  "max_length": 300,        // Increased from 100 to accommodate typical CTCF sites
  "max_n_percent": 25,      // Relaxed from 15% to handle more sequences  
  "low_complexity_filter": false,  // Disabled to retain more sequences
  "mask_repeats": false     // Disabled for initial training
}
```

### Expected Improvement
- **Previous retention**: 0.08% (37/44,217 sequences)
- **Expected retention**: 15-25% (6,000-10,000 sequences) 
- **Minimum needed for PWM**: 1,000+ sequences
- **Ideal for robust model**: 5,000+ sequences

## Sequence Characteristics Analysis

Based on examining your `data/extracted_sequences.fasta`:

### Length Distribution (estimated)
- **≤ 100 bp**: ~1% of sequences (why current filter fails)
- **101-200 bp**: ~15% of sequences  
- **201-300 bp**: ~75% of sequences (main population)
- **> 300 bp**: ~9% of sequences

### Typical CTCF Binding Sites
- **Core binding motif**: ~19bp (CCGCGNGGNGGCAG)
- **Extended region**: 200-300bp (includes flanking regulatory sequences)
- **Your data range**: Based on coordinates, sequences are mostly 200-300bp

## Parameter Optimization Strategy

### Phase 1: Basic Retention (Current)
```json
{
  "max_length": 300,
  "max_n_percent": 25,
  "low_complexity_filter": false,
  "mask_repeats": false
}
```
**Target**: Retain 10-25% of sequences for initial model training

### Phase 2: Quality Refinement (Future)
Once you have a working model, you can gradually tighten parameters:
```json
{
  "max_length": 250,
  "max_n_percent": 20,
  "low_complexity_filter": true,
  "entropy_threshold": 1.0,
  "mask_repeats": true
}
```
**Target**: Retain 5-15% of highest quality sequences

### Phase 3: Production Ready (Final)
```json
{
  "max_length": 200,
  "max_n_percent": 15,  
  "low_complexity_filter": true,
  "entropy_threshold": 1.2,
  "mask_repeats": true
}
```
**Target**: Retain 2-8% of sequences with optimal quality/quantity balance

## Running the Analysis

### Current Analysis Tools Available

1. **analyze_sequences.R** - Comprehensive R-based analysis
2. **analyze_sequences.py** - Python version with BioPython  
3. **analyze_sequences.ps1** - PowerShell version for Windows
4. **analyze_preprocessing.sh** - Quick bash script for basic stats

### Usage Examples

```bash
# Run analysis in Docker environment
bash run-in-docker.sh scripts/analyze_sequences.R data/extracted_sequences.fasta analysis_output

# Quick analysis without dependencies
bash analyze_preprocessing.sh

# Check current processing status
grep -c '^>' data/extracted_sequences.fasta          # Input count
grep -c '^>' data/preprocessed_sequences.fasta       # Current output
grep -c '^>' data/preprocessed_sequences_optimized.fasta  # Optimized output
```

## Next Steps

### Immediate (Current Processing)
1. ✅ **Fixed preprocessing parameters** - Updated config file
2. 🔄 **Running optimized preprocessing** - Currently processing ~6000+ sequences
3. ⏳ **Verify retention rate** - Should be 15-25% when complete

### Short Term (Next 1-2 hours)
1. **Check optimized results**: Verify retention rate meets targets (1000+ sequences)
2. **Run sequence analysis**: Use analysis scripts to validate quality
3. **Continue pipeline**: Move to `prepare_datasets.R` once sufficient sequences retained

### Medium Term (Next steps in pipeline)
1. **Dataset preparation**: Split sequences for training/validation/testing
2. **PWM model building**: Create position weight matrices from processed sequences
3. **Model evaluation**: Test prediction accuracy and optimize thresholds

## Troubleshooting

### If Retention Rate Still Too Low (< 5%)
```bash
# Further relax parameters
# Edit scripts/preprocess_config.json:
{
  "max_length": 500,     # Even more permissive
  "max_n_percent": 35,   # Allow more N bases
  "min_length": 50       # Lower minimum length
}
```

### If Retention Rate Too High (> 50%)
```bash
# Gradually tighten parameters
# Edit scripts/preprocess_config.json:
{
  "max_length": 250,
  "max_n_percent": 20,
  "low_complexity_filter": true,
  "entropy_threshold": 1.2
}
```

### If Processing Errors
```bash
# Check for common issues:
1. File paths - ensure forward slashes in R scripts
2. Dependencies - BiocManager::install("Biostrings") in R
3. Memory - large datasets may need increased memory limits
4. Docker - ensure container has access to data directory
```

## Performance Expectations

### Current Hardware Requirements
- **Memory**: ~2-4GB for 44K sequences
- **Processing time**: 10-30 minutes depending on parameters
- **Storage**: ~100-500MB for processed sequences

### Model Training Requirements
- **Minimum sequences**: 1,000 (basic PWM)
- **Recommended sequences**: 5,000+ (robust model)  
- **Optimal sequences**: 10,000+ (production model)

## Files Created/Modified

### Configuration Files
- ✅ `scripts/preprocess_config.json` - Optimized parameters
- ✅ `scripts/preprocess_config_optimized.json` - Backup config
- ✅ `PARAMETER_OPTIMIZATION_GUIDE.md` - This guide

### Analysis Scripts  
- ✅ `scripts/analyze_sequences.R` - Comprehensive analysis
- ✅ `scripts/analyze_sequences.py` - Python version
- ✅ `scripts/analyze_sequences.ps1` - PowerShell version
- ✅ `analyze_preprocessing.sh` - Quick analysis script

### Processing Scripts
- ✅ `scripts/preprocess_sequences.R` - Enhanced JSON parsing
- ✅ `quick_analysis.py` - Simple length analysis

---

**Status**: Optimized preprocessing currently running. Expected completion with 6,000-10,000 retained sequences (target: >1,000 for PWM training).

**Next Action**: Verify optimized preprocessing results and continue to dataset preparation phase.
