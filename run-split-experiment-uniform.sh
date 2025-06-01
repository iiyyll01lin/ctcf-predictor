#!/usr/bin/env bash
# run-split-experiment-uniform.sh
# -------------------------------------------------------------
# 比較四種資料切分策略對 Uniform Null PWM (全位置均一 0.25) 模型預測效能的影響：
#   1. random split
#   2. chrom split（染色體劃分）
#   3. stratified k-fold CV
#   4. groupk（染色體為群組的 k-fold CV）
# -------------------------------------------------------------

set -euo pipefail
LOG_DIR=logs/uniform
mkdir -p "$LOG_DIR"

# 1. Random Split ------------------------------------------------------------
echo "[1] Random split (Uniform Null PWM)"
# 1-1. 切分原始資料 (正+負)，此步驟與真實 PWM 相同
Rscript scripts/split_strategy.R \
        --strategy random \
        --input data/complete_labeled_82.fasta \
        --out_dir data/splits/random \
        | tee "$LOG_DIR/split_random.log"

# 1-2. 產生 Uniform Null PWM (長度 82)
Rscript scripts/build_null_pwm_uniform.R \
        --length 82 \
        --output results/null_random_uniform.rds \
        | tee "$LOG_DIR/build_null_random_uniform.log"

# 1-3. 在同一個 Random test 集上評估 Uniform Null
Rscript scripts/evaluate_models.R \
        --pwm   results/null_random_uniform.rds \
        --test  data/splits/random/test_sequences.fasta \
        --out   results/metrics_random_uniform.csv \
        | tee "$LOG_DIR/eval_random_uniform.log"


# 2. Chromosome Split -------------------------------------------------------
echo "[2] Chromosome-based split (Uniform Null PWM)"
# 2-1. 切分原始資料 (正+負)
Rscript scripts/split_strategy.R \
        --strategy chrom \
        --input data/complete_labeled_82.fasta \
        --train_chr chr1,chr2,chr3,chr4,chr5,chr6 \
        --val_chr   chr7 \
        --test_chr  chr8,chr9 \
        --out_dir data/splits/chrom \
        | tee "$LOG_DIR/split_chrom.log"

# 2-2. 產生 Uniform Null PWM (長度 82)
Rscript scripts/build_null_pwm_uniform.R \
        --length 82 \
        --output results/null_chrom_uniform.rds \
        | tee "$LOG_DIR/build_null_chrom_uniform.log"

# 2-3. 在 Chrom test 集上評估 Uniform Null
Rscript scripts/evaluate_models.R \
        --pwm   results/null_chrom_uniform.rds \
        --test  data/splits/chrom/test_sequences.fasta \
        --out   results/metrics_chrom_uniform.csv \
        | tee "$LOG_DIR/eval_chrom_uniform.log"


# 3. Stratified CV（k-fold CV）------------------------------------------------
echo "[3] Stratified k-fold CV (Uniform Null PWM)"

# 3-1. 產生 Uniform Null PWM (長度 82)
Rscript scripts/build_null_pwm_uniform.R \
        --length 82 \
        --output results/null_cv_uniform.rds \
        | tee "$LOG_DIR/build_null_cv_uniform.log"

# 3-2. 用 CV 評估 Uniform Null
Rscript scripts/evaluate_models_with_cv.R \
        --fasta      data/complete_labeled_82.fasta \
        --pwm_dir    results \
        --pwm_files  null_cv_uniform.rds \
        --k          5 \
        --repeats    3 \
        | tee "$LOG_DIR/eval_cv_uniform.log"

mv results/metrics_cv.csv results/metrics_cv_uniform.csv

# 4. GroupKFold（染色體為群組的 k-fold CV）-------------------------------
echo "[4] GroupKFold (chromosome-based CV, Uniform Null PWM)"

# 4-1. 先對原始資料 (正+負) 做 GroupKFold 拆分
Rscript scripts/split_strategy.R \
        --strategy groupk \
        --input data/complete_labeled_82.fasta \
        --out_dir data/splits/groupk/standard \
        --k 5 \
        | tee "$LOG_DIR/split_groupk_standard.log"

# 4-2. Uniform Null 下的 GroupKFold：呼叫專用的 evaluate_models_with_groupk_uniform.R
#      （此腳本內部會針對每個 fold 都直接用 uniform PWM 計算 AUC）
Rscript scripts/evaluate_models_with_groupk_uniformPWM.R \
        --fold_dir data/splits/groupk/standard \
        --out_dir  results \
        | tee "$LOG_DIR/eval_groupk_uniform.log"


# Visualize ------------------------------------------------------------------
echo "[✔] Visualizing Uniform Null results"
Rscript scripts/visualize_metrics.R \
        --metrics_dir results \
        --out_dir    results/plots_uniform \
        | tee "$LOG_DIR/plot_uniform.log"

echo "[✔] All Uniform Null experiments done."
