# PWM Quality Validation Script

# --- Dependencies ---
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Package 'Biostrings' is needed. Please install via BiocManager.")
}
library(Biostrings)

# --- Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
pwm_file <- if (length(args) >= 1) args[1] else "results/generated_pwm.rds"

# --- Main Analysis ---
cat("=== PWM Quality Validation Report ===\n")
cat("PWM file:", pwm_file, "\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Validate input file
if (!file.exists(pwm_file)) {
  stop("PWM file not found: ", pwm_file)
}

# Load PWM
cat("Loading PWM...\n")
pwm <- readRDS(pwm_file)

if (!is.matrix(pwm)) {
  stop("PWM object is not a matrix")
}

if (nrow(pwm) != 4) {
  stop("PWM should have 4 rows (A, C, G, T)")
}

cat("PWM dimensions:", nrow(pwm), "x", ncol(pwm), "(bases x positions)\n")
cat("PWM length:", ncol(pwm), "bp\n\n")

# 1. Basic PWM Validation
cat("1. BASIC PWM VALIDATION\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# Check if probabilities sum to 1
col_sums <- colSums(pwm)
cat("Column sums (should be ~1.0):\n")
cat("Min:", round(min(col_sums), 4), "\n")
cat("Max:", round(max(col_sums), 4), "\n")
cat("Mean:", round(mean(col_sums), 4), "\n")

if (all(abs(col_sums - 1.0) < 0.01)) {
  cat("✅ All columns sum to ~1.0 - VALID probability matrix\n")
} else {
  cat("❌ Some columns don't sum to 1.0 - INVALID probability matrix\n")
}

# Check for zero probabilities
zero_count <- sum(pwm == 0)
if (zero_count > 0) {
  cat("⚠️  Warning:", zero_count, "zero probabilities found\n")
  cat("   This may cause issues in log-likelihood scoring\n")
} else {
  cat("✅ No zero probabilities - good for scoring\n")
}

# 2. Information Content Analysis
cat("\n2. INFORMATION CONTENT ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# Calculate information content for each position
info_content <- apply(pwm, 2, function(x) {
  # Replace zeros with small value to avoid log(0)
  x[x == 0] <- 1e-10
  # Information content = sum(p * log2(p/0.25))
  sum(x * log2(x/0.25))
})

cat("Information content statistics (bits):\n")
cat("Total information content:", round(sum(info_content), 3), "bits\n")
cat("Average per position:", round(mean(info_content), 3), "bits\n")
cat("Maximum position:", round(max(info_content), 3), "bits\n")
cat("Minimum position:", round(min(info_content), 3), "bits\n")
cat("Standard deviation:", round(sd(info_content), 3), "bits\n\n")

# Classify information content quality
total_info <- sum(info_content)
avg_info <- mean(info_content)

cat("PWM Quality Assessment:\n")
if (total_info < 5) {
  cat("❌ POOR: Very low total information content (<5 bits)\n")
  cat("   Recommendation: Need more or higher quality sequences\n")
} else if (total_info < 10) {
  cat("⚠️  FAIR: Low total information content (5-10 bits)\n")
  cat("   Recommendation: Consider improving training data\n")
} else if (total_info < 20) {
  cat("✅ GOOD: Moderate total information content (10-20 bits)\n")
} else {
  cat("✅ EXCELLENT: High total information content (>20 bits)\n")
}

if (avg_info < 0.5) {
  cat("❌ POOR: Very low average information per position\n")
} else if (avg_info < 1.0) {
  cat("⚠️  FAIR: Low average information per position\n")
} else {
  cat("✅ GOOD: Adequate average information per position\n")
}

# 3. Motif Core Analysis
cat("\n3. MOTIF CORE ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# Find high-information positions (conserved core)
high_info_threshold <- 1.0
high_info_positions <- which(info_content > high_info_threshold)

cat("High-information positions (>", high_info_threshold, "bits):", 
    length(high_info_positions), "positions\n")

if (length(high_info_positions) > 0) {
  cat("Core positions:", paste(high_info_positions, collapse = ", "), "\n")
  
  # Extract core motif sequence
  core_motif <- ""
  for (pos in high_info_positions) {
    dominant_base <- rownames(pwm)[which.max(pwm[, pos])]
    core_motif <- paste0(core_motif, dominant_base)
  }
  cat("Core motif sequence:", core_motif, "\n")
  
  # Calculate core conservation
  core_conservation <- mean(apply(pwm[, high_info_positions, drop=FALSE], 2, max))
  cat("Average conservation in core:", round(core_conservation * 100, 1), "%\n")
  
} else {
  cat("❌ No strongly conserved core found\n")
  cat("   Recommendation: Check if sequences are properly aligned\n")
}

# 4. Position-by-Position Analysis
cat("\n4. TOP 10 MOST INFORMATIVE POSITIONS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# Sort positions by information content
top_positions <- order(info_content, decreasing = TRUE)[1:min(10, ncol(pwm))]

cat("Position | Info(bits) | Dominant Base | Probability | Consensus\n")
cat(paste(rep("-", 60), collapse=""), "\n")

for (i in 1:length(top_positions)) {
  pos <- top_positions[i]
  dominant_idx <- which.max(pwm[, pos])
  dominant_base <- rownames(pwm)[dominant_idx]
  dominant_prob <- pwm[dominant_idx, pos]
  info_bits <- info_content[pos]
  
  # Create consensus representation
  probs <- pwm[, pos]
  consensus <- paste0(rownames(pwm), ":", round(probs * 100, 0), "%", collapse = " ")
  
  cat(sprintf("%8d | %10.3f | %13s | %11.3f | %s\n", 
              pos, info_bits, dominant_base, dominant_prob, 
              paste0(dominant_base, "(", round(dominant_prob*100, 0), "%)")))
}

# 5. Positional Bias Analysis
cat("\n5. POSITIONAL BIAS ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# Check if information is concentrated in specific regions
first_third <- 1:floor(ncol(pwm)/3)
middle_third <- (floor(ncol(pwm)/3)+1):(2*floor(ncol(pwm)/3))
last_third <- (2*floor(ncol(pwm)/3)+1):ncol(pwm)

info_first <- mean(info_content[first_third])
info_middle <- mean(info_content[middle_third])
info_last <- mean(info_content[last_third])

cat("Average information content by region:\n")
cat("First third (positions 1-", max(first_third), "):", round(info_first, 3), "bits\n")
cat("Middle third (positions", min(middle_third), "-", max(middle_third), "):", round(info_middle, 3), "bits\n")
cat("Last third (positions", min(last_third), "-", max(last_third), "):", round(info_last, 3), "bits\n")

# Identify the most informative region
max_region_info <- max(info_first, info_middle, info_last)
if (max_region_info == info_first) {
  cat("Most informative region: First third\n")
} else if (max_region_info == info_middle) {
  cat("Most informative region: Middle third\n")
} else {
  cat("Most informative region: Last third\n")
}

# 6. Overall Recommendations
cat("\n6. OVERALL RECOMMENDATIONS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

cat("Based on PWM analysis:\n\n")

# Information content recommendations
if (total_info >= 15 && length(high_info_positions) >= 5) {
  cat("✅ PWM QUALITY: Excellent - ready for prediction\n")
} else if (total_info >= 8 && length(high_info_positions) >= 3) {
  cat("✅ PWM QUALITY: Good - suitable for prediction with caution\n")
} else if (total_info >= 5) {
  cat("⚠️  PWM QUALITY: Fair - may work but consider improvements\n")
  cat("   Suggestions:\n")
  cat("   - Increase training data size\n")
  cat("   - Improve sequence alignment\n")
  cat("   - Check sequence quality\n")
} else {
  cat("❌ PWM QUALITY: Poor - not recommended for prediction\n")
  cat("   Critical issues to address:\n")
  cat("   - Very low information content\n")
  cat("   - Insufficient sequence conservation\n")
  cat("   - Need better training data\n")
}

# Core motif recommendations
if (length(high_info_positions) >= 5) {
  cat("✅ MOTIF CORE: Well-defined core region identified\n")
} else if (length(high_info_positions) >= 2) {
  cat("⚠️  MOTIF CORE: Weak core region - check alignment\n")
} else {
  cat("❌ MOTIF CORE: No clear core - major alignment issues\n")
}

cat("\n=== PWM Validation Complete ===\n")
