# run-split-experiment-final.sh
# -------------------------------------------------------------
# 比較四種資料切分策略對 PWM 模型預測效能的影響：
#   1. random split
#   2. chrom split（染色體劃分）
#   3. stratified k-fold CV（隨機交叉驗證）
#   4. groupk（以染色體為群組的 k-fold CV）
# -------------------------------------------------------------

set -euo pipefail
LOG_DIR=logs
mkdir -p "$LOG_DIR"

# 1. Random Split ------------------------------------------------------------
echo "[1] Random split"
Rscript scripts/split_strategy.R \
        --strategy random \
        --input data/complete_labeled_82.fasta \
        --out_dir data/splits/random \
        | tee "$LOG_DIR/split_random.log"

Rscript scripts/build_pwm.R \
        --input  data/splits/random/train_sequences.fasta \
        --output results/pwm_random.rds \
        | tee "$LOG_DIR/build_pwm_random.log"

Rscript scripts/evaluate_models.R \
        --pwm results/pwm_random.rds \
        --test data/splits/random/test_sequences.fasta \
        --out results/metrics_random.csv \
        | tee "$LOG_DIR/eval_random.log"

# 2. Chromosome Split -------------------------------------------------------
echo "[2] Chromosome-based split"
Rscript scripts/split_strategy.R \
        --strategy chrom \
        --input data/complete_labeled_82.fasta \
        --train_chr chr1,chr2,chr3,chr4,chr5,chr6 \
        --val_chr chr7 \
        --test_chr chr8,chr9 \
        --out_dir data/splits/chrom \
        | tee "$LOG_DIR/split_chrom.log"

Rscript scripts/build_pwm.R \
        --input  data/splits/chrom/train_sequences.fasta \
        --output results/pwm_chrom.rds \
        | tee "$LOG_DIR/build_pwm_chrom.log"

Rscript scripts/evaluate_models.R \
        --pwm results/pwm_chrom.rds \
        --test data/splits/chrom/test_sequences.fasta \
        --out results/metrics_chrom.csv \
        | tee "$LOG_DIR/eval_chrom.log"

# 3. Stratified CV（舊稱 groupk，實為一般 CV） -----------------------------
echo "[3] Stratified CV"
Rscript scripts/build_pwm.R \
        --input  data/complete_labeled_82.fasta \
        --output results/pwm_cv.rds \
        | tee "$LOG_DIR/build_pwm_cv.log"

Rscript scripts/evaluate_models_with_cv.R \
        --fasta   data/complete_labeled_82.fasta \
        --pwm_dir results \
        --pwm_files   pwm_cv.rds \
        --k       5 \
        --repeats 3 \
        | tee "$LOG_DIR/eval_cv.log"

mv results/metrics_cv.csv results/metrics_cv_real.csv

# 4. GroupKFold（染色體分組 CV） -------------------------------------------
echo "[4] GroupKFold (chromosome-based CV)"
Rscript scripts/split_strategy.R \
        --strategy groupk \
        --input data/complete_labeled_82.fasta \
        --out_dir data/splits/groupk \
        --k 5 \
        | tee "$LOG_DIR/split_groupk.log"

Rscript scripts/evaluate_models_with_groupk.R \
        --fold_dir data/splits/groupk \
        --out_dir results \
        | tee "$LOG_DIR/eval_groupk.log"

# Visualize ------------------------------------------------------------------
echo "[✔] Visualizing results"
Rscript scripts/visualize_metrics.R \
        --metrics_dir results \
        --out_dir results \
        | tee "$LOG_DIR/plot_metrics.log"

echo "[✔] All done."
