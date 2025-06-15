# Scripts Reference

> **📜 Complete Guide to All R Scripts and Utilities**  
> Detailed documentation of every script in the CTCF PWM Testing Pipeline, including parameters, usage, and integration points.

## 📁 Script Organization

The pipeline's scripts are organized by functionality in the `scripts/` directory:

```
scripts/
├── Data Preprocessing
│   ├── preprocess_sequences.R
│   ├── prepare_datasets.R
│   └── analyze_sequences.R
├── Alignment Methods
│   ├── advanced_alignment.R
│   ├── analyze_sequence_alignment.R
│   └── simple_aligned_pwm.R
├── PWM Construction
│   ├── build_pwm_robust.R
│   ├── build_aligned_pwm.R
│   ├── build_subset_pwm.R
│   └── efficient_aligned_pwm.R
├── Analysis & Validation
│   ├── validate_pwm_quality.R
│   ├── compare_pwms.R
│   ├── enhanced_compare_pwms.R
│   └── statistical_significance_test.R
├── Model Evaluation
│   ├── evaluate_models.R
│   ├── evaluate_models_with_cv.R
│   └── generate_null_models.R
├── Utilities
│   ├── predict_ctcf.R
│   ├── optimize_threshold.R
│   └── analyze_sequence_quality.R
└── Configuration & Guides
    ├── PARAMETER_OPTIMIZATION_GUIDE.md
    └── SEQUENCE_ANALYSIS_GUIDE.md
```

## 🔧 Core Pipeline Scripts

### Data Preprocessing

#### `preprocess_sequences.R`

**Purpose:** Initial sequence preprocessing and quality control

**Usage:**
```bash
Rscript scripts/preprocess_sequences.R \
  --input data/raw_sequences.fa \
  --output data/processed_sequences.fa \
  --config scripts/preprocess_config.json
```

**Key Parameters:**
- `--input`: Input FASTA file path
- `--output`: Output file path
- `--config`: Configuration file (JSON format)
- `--min-length`: Minimum sequence length (default: 50)
- `--max-length`: Maximum sequence length (default: 500)
- `--quality-threshold`: Quality score threshold (default: 0.8)

**Key Functions:**
```r
# Main preprocessing function
preprocess_sequences <- function(sequences, config) {
  # Filter by length
  filtered <- filter_by_length(sequences, config$min_length, config$max_length)
  
  # Quality scoring
  scored <- calculate_quality_scores(filtered)
  
  # Remove low-quality sequences
  high_quality <- filter_by_quality(scored, config$quality_threshold)
  
  return(high_quality)
}
```

**Output Files:**
- Processed sequence FASTA file
- Quality report (JSON)
- Processing log

#### `prepare_datasets.R`

**Purpose:** Prepare datasets for different validation approaches

**Usage:**
```bash
Rscript scripts/prepare_datasets.R \
  --input data/processed_sequences.fa \
  --method chromosome_split \
  --output-dir data/datasets/
```

**Key Functions:**
- `prepare_chromosome_split()` - Split by chromosome for validation
- `prepare_random_split()` - Random train/test split
- `prepare_cross_validation()` - K-fold cross-validation setup

**Dataset Types:**
- **Chromosome Split**: Chromosomes 1-15 for training, 16-22+X+Y for testing
- **Random Split**: 80/20 random split
- **Cross-Validation**: 5-fold or 10-fold CV

#### `analyze_sequences.R`

**Purpose:** Comprehensive sequence analysis and statistics

**Usage:**
```bash
Rscript scripts/analyze_sequences.R \
  --input data/processed_sequences.fa \
  --output results/sequence_analysis.html
```

**Analysis Components:**
- Length distribution analysis
- Nucleotide composition
- Sequence complexity scores
- Motif content assessment
- Quality distribution

### Alignment Methods

#### `advanced_alignment.R`

**Purpose:** Advanced alignment algorithms for sequence optimization

**Key Methods:**
```r
# Multiple alignment strategies
align_sequences <- function(sequences, method = "integrated") {
  switch(method,
    "center" = center_based_alignment(sequences),
    "consensus" = consensus_based_alignment(sequences),
    "integrated" = integrated_alignment(sequences),
    "iterative" = iterative_alignment(sequences)
  )
}
```

**Alignment Strategies:**
- **Center-based**: Align by sequence center
- **Consensus-based**: Align to consensus pattern
- **Integrated**: Hybrid approach combining methods
- **Iterative**: Progressive refinement alignment

#### `analyze_sequence_alignment.R`

**Purpose:** Evaluate alignment quality and optimization

**Quality Metrics:**
- Information content improvement
- Alignment stability
- Pattern recognition accuracy
- Conservation enhancement

**Usage:**
```bash
Rscript scripts/analyze_sequence_alignment.R \
  --sequences data/sequences.fa \
  --method integrated \
  --output results/alignment_analysis.pdf
```

### PWM Construction

#### `build_pwm_robust.R` ⭐

**Purpose:** Main PWM construction with robust quality control

**Usage:**
```bash
Rscript scripts/build_pwm_robust.R \
  --sequences data/aligned_sequences.fa \
  --output results/ctcf_pwm.meme \
  --format meme \
  --pseudocount 0.1
```

**Key Parameters:**
- `--sequences`: Input aligned sequences
- `--output`: Output PWM file
- `--format`: Output format (meme, jaspar, transfac)
- `--pseudocount`: Pseudocount for probability calculation
- `--min-ic`: Minimum information content threshold
- `--background`: Background nucleotide frequencies

**Core Functions:**
```r
# Main PWM construction
build_robust_pwm <- function(sequences, config) {
  # Calculate frequency matrix
  freq_matrix <- calculate_frequencies(sequences)
  
  # Add pseudocounts
  pseudo_matrix <- add_pseudocounts(freq_matrix, config$pseudocount)
  
  # Convert to probabilities
  prob_matrix <- normalize_probabilities(pseudo_matrix)
  
  # Calculate information content
  ic_vector <- calculate_information_content(prob_matrix)
  
  # Quality assessment
  quality <- assess_pwm_quality(ic_vector, config)
  
  return(list(
    pwm = prob_matrix,
    information_content = ic_vector,
    quality_metrics = quality
  ))
}
```

**Quality Thresholds:**
- Total IC > 8.0 bits (minimum acceptable)
- Total IC > 12.0 bits (good quality)
- Total IC > 16.0 bits (excellent quality)
- Conserved positions > 3 (minimum)

#### `build_aligned_pwm.R`

**Purpose:** PWM construction specifically for aligned sequences

**Optimizations:**
- Alignment-aware processing
- Position-specific quality control
- Enhanced information content calculation

#### `efficient_aligned_pwm.R`

**Purpose:** Memory and performance optimized PWM construction

**Features:**
- Streaming sequence processing
- Reduced memory footprint
- Parallel processing support
- Large dataset handling

### Analysis & Validation

#### `validate_pwm_quality.R`

**Purpose:** Comprehensive PWM quality validation

**Validation Tests:**
```r
# Quality validation suite
validate_quality <- function(pwm, sequences) {
  tests <- list(
    information_content = test_information_content(pwm),
    conservation = test_conservation_pattern(pwm),
    biological_relevance = test_biological_pattern(pwm),
    statistical_significance = test_statistical_significance(pwm, sequences)
  )
  
  return(compile_validation_report(tests))
}
```

**Quality Metrics:**
- Total information content
- Conserved positions count
- Pattern recognition accuracy
- Statistical significance

#### `compare_pwms.R`

**Purpose:** Compare multiple PWMs across different methods

**Comparison Metrics:**
- Information content comparison
- Motif similarity scores
- Performance benchmarking
- Statistical comparisons

**Usage:**
```bash
Rscript scripts/compare_pwms.R \
  --pwm1 results/method1_pwm.meme \
  --pwm2 results/method2_pwm.meme \
  --output results/comparison_report.html
```

#### `statistical_significance_test.R`

**Purpose:** Statistical validation of PWM quality

**Statistical Tests:**
- Wilcoxon rank-sum test
- Kolmogorov-Smirnov test
- Chi-square goodness of fit
- Permutation testing

### Model Evaluation

#### `evaluate_models.R`

**Purpose:** Comprehensive model evaluation and benchmarking

**Evaluation Metrics:**
- ROC curve analysis
- Precision-recall curves
- F1 scores
- Matthews correlation coefficient

**Usage:**
```bash
Rscript scripts/evaluate_models.R \
  --pwm results/ctcf_pwm.meme \
  --test-sequences data/test_sequences.fa \
  --output results/evaluation_report.html
```

#### `evaluate_models_with_cv.R`

**Purpose:** Cross-validation based model evaluation

**CV Methods:**
- K-fold cross-validation
- Leave-one-chromosome-out
- Stratified sampling
- Bootstrap validation

#### `generate_null_models.R`

**Purpose:** Generate null models for statistical comparison

**Null Model Types:**
- Shuffled sequences
- Random nucleotide models
- Uniform distribution models
- Background frequency models

### Utility Scripts

#### `predict_ctcf.R`

**Purpose:** Use trained PWM to predict CTCF binding sites

**Usage:**
```bash
Rscript scripts/predict_ctcf.R \
  --pwm results/ctcf_pwm.meme \
  --genome data/genome.fa \
  --threshold 0.8 \
  --output results/predictions.bed
```

**Prediction Parameters:**
- `--threshold`: Score threshold for predictions
- `--window-size`: Sliding window size
- `--step-size`: Step size for scanning
- `--output-format`: BED, GFF, or custom format

#### `optimize_threshold.R`

**Purpose:** Optimize score thresholds for predictions

**Optimization Methods:**
- ROC curve optimization
- F1 score maximization
- Precision-recall optimization
- Custom objective functions

#### `analyze_sequence_quality.R`

**Purpose:** Quality control and assessment of input sequences

**Quality Metrics:**
- Sequence complexity
- Repetitive element content
- GC content distribution
- Length statistics

## 🔧 Configuration Files

### `preprocess_config.json`

**Standard Configuration:**
```json
{
  "min_length": 50,
  "max_length": 500,
  "quality_threshold": 0.8,
  "remove_duplicates": true,
  "filter_repetitive": true,
  "gc_content_range": [0.3, 0.7]
}
```

### `preprocess_config_optimized.json`

**High-Quality Configuration:**
```json
{
  "min_length": 100,
  "max_length": 300,
  "quality_threshold": 0.9,
  "remove_duplicates": true,
  "filter_repetitive": true,
  "gc_content_range": [0.4, 0.6],
  "complexity_threshold": 0.7
}
```

### `preprocess_config_relaxed.json`

**Permissive Configuration:**
```json
{
  "min_length": 30,
  "max_length": 1000,
  "quality_threshold": 0.6,
  "remove_duplicates": false,
  "filter_repetitive": false,
  "gc_content_range": [0.2, 0.8]
}
```

## 🚀 Usage Patterns

### Standard Pipeline Execution

```bash
# 1. Preprocess sequences
Rscript scripts/preprocess_sequences.R \
  --input data/raw_sequences.fa \
  --output data/processed_sequences.fa

# 2. Prepare datasets
Rscript scripts/prepare_datasets.R \
  --input data/processed_sequences.fa \
  --method chromosome_split

# 3. Advanced alignment
Rscript scripts/advanced_alignment.R \
  --input data/processed_sequences.fa \
  --method integrated \
  --output data/aligned_sequences.fa

# 4. Build robust PWM
Rscript scripts/build_pwm_robust.R \
  --sequences data/aligned_sequences.fa \
  --output results/ctcf_pwm.meme

# 5. Validate quality
Rscript scripts/validate_pwm_quality.R \
  --pwm results/ctcf_pwm.meme \
  --sequences data/aligned_sequences.fa
```

### Quick Analysis Pipeline

```bash
# Quick quality check
Rscript scripts/analyze_sequences.R --input data/sequences.fa

# Quick PWM build and validation
Rscript scripts/build_pwm_robust.R --sequences data/sequences.fa
Rscript scripts/validate_pwm_quality.R --pwm results/ctcf_pwm.meme
```

### Comparison Analysis

```bash
# Compare different alignment methods
for method in center consensus integrated; do
  Rscript scripts/advanced_alignment.R \
    --method $method \
    --output data/aligned_${method}.fa
  
  Rscript scripts/build_pwm_robust.R \
    --sequences data/aligned_${method}.fa \
    --output results/pwm_${method}.meme
done

# Compare resulting PWMs
Rscript scripts/enhanced_compare_pwms.R \
  --pwms results/pwm_*.meme \
  --output results/method_comparison.html
```

## 🔍 Error Handling

### Common Error Messages

**File Not Found:**
```
Error: Input file 'data/sequences.fa' not found
Solution: Check file path and ensure file exists
```

**Memory Issues:**
```
Error: Cannot allocate vector of size X MB
Solution: Use efficient_aligned_pwm.R or reduce batch size
```

**Quality Threshold Failures:**
```
Warning: PWM quality below threshold (IC=4.2, required=8.0)
Solution: Improve alignment or increase dataset size
```

### Debug Mode

Most scripts support debug mode:
```bash
Rscript scripts/build_pwm_robust.R --debug --verbose
```

## 📊 Performance Optimization

### Memory Usage

**High Memory Scripts:**
- `build_pwm_robust.R` - Use `--batch-size` parameter
- `advanced_alignment.R` - Use `--streaming` mode
- `evaluate_models.R` - Use `--memory-limit` parameter

**Low Memory Alternatives:**
- `efficient_aligned_pwm.R` instead of `build_aligned_pwm.R`
- `simple_aligned_pwm.R` for basic analysis

### Parallel Processing

**Multi-core Support:**
```bash
# Use multiple cores
export CTCF_THREADS=4
Rscript scripts/build_pwm_robust.R --threads 4
```

**Batch Processing:**
```bash
# Process large datasets in batches
Rscript scripts/build_pwm_robust.R --batch-size 1000
```

---

## 📖 Next Reading

- **[Configuration Guide](09-configuration.md)** - Complete parameter reference
- **[User Guide](10-user-guide.md)** - Step-by-step usage instructions
- **[Testing & Validation](11-testing-validation.md)** - Quality assurance procedures

---

*This comprehensive script reference ensures you can effectively use every component of the CTCF PWM Testing Pipeline for your specific research needs.*
