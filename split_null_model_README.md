```markdown
# CTCF 結合位點預測管線說明

## 專案目標

本專案旨在建立一條完整的自動化分析管線，用以預測 DNA 序列中 CTCF（CCCTC‐binding factor）轉錄因子之潛在結合位點。CTCF 在染色質組織與基因表現調控上扮演重要角色。本管線從資料下載、前處理、標記正負例、資料切分、模型訓練、評估，到最終比較不同切分策略與 null model（Uniform PWM）的表現，逐步完成整個流程。

---

## 專案目錄結構

```

/
├── data/
│   ├── reference\_genome/
│   │   ├── hg38.chr21.fa             # Demo 用 chr21 參考序列
│   │   └── hg38.fa                   # 完整 hg38 參考序列
│   ├── K562\_CTCF\_peaks.bed           # 下載之 ChIP‐seq 峰值 BED 檔
│   ├── extracted\_sequences.fasta     # 由 bedtools getfasta 擷取之正向序列
│   ├── preprocessed\_sequences.fasta  # 經過 preprocess\_sequences.R 處理後的序列
│   └── splits/                       # split\_strategy.R 生成之切分結果
│       ├── random/
│       │   ├── train\_sequences.fasta
│       │   └── test\_sequences.fasta
│       ├── chrom/
│       │   ├── train\_sequences.fasta
│       │   ├── val\_sequences.fasta   # 若有指定 val\_chr
│       │   └── test\_sequences.fasta
│       └── groupk/
│           ├── fold1/train\_sequences.fasta
│           ├── fold1/test\_sequences.fasta
│           ├── fold2/…
│           └── foldK/…
├── results/
│   ├── pwm\_random.rds                # build\_pwm.R 產出之一般 PWM (random)
│   ├── pwm\_chrom.rds                 # build\_pwm.R 產出之一般 PWM (chrom)
│   ├── pwm\_cv.rds                    # build\_pwm.R 產出之一般 PWM (CV)
│   ├── uniform\_pwm.rds               # build\_null\_pwm\_uniform.R 產出之 Uniform PWM
│   ├── metrics\_random.csv            # evaluate\_models.R (random) 產出之 AUC 檔
│   ├── metrics\_chrom.csv             # evaluate\_models.R (chrom) 產出之 AUC 檔
│   ├── metrics\_cv\_real.csv           # evaluate\_models\_with\_cv.R (CV) 產出之 AUC 檔（已改名）
│   ├── metrics\_groupk.csv            # evaluate\_models\_with\_groupk.R 產出之 AUC 檔
│   ├── summary\_auc.csv               # visualize\_metrics.R 產出之各策略統計摘要
│   ├── auc\_boxplot.png               # visualize\_metrics.R 產出之 Boxplot 圖
│   ├── plots\_uniform/                # run-split-experiment-uniform.sh 執行後產生
│   │   ├── random\_uniform\_auc.csv
│   │   ├── chrom\_uniform\_auc.csv
│   │   ├── cv\_uniform\_auc.csv
│   │   ├── groupk\_uniform\_auc.csv
│   │   ├── summary\_auc.csv
│   │   └── auc\_boxplt.png
│   └── predictions.tsv               # predict\_ctcf.R 之範例預測輸出
├── scripts/
│   ├── download\_data.sh              # 下載範例 CTCF 峰值資料與參考基因組
│   ├── preprocess\_sequences.R        # 序列前處理腳本 (長度/Ｎ-base/低複雜度/重複遮罩)
│   ├── preprocess\_config.json        # 前處理參數範例
│   ├── preprocess\_sequences\_with\_stats.R # 前處理並計算統計資訊（備用）
│   ├── prepare\_datasets\_full.R       # 以預處理後序列產生完整標記資料 (正負例)
│   ├── prepare\_datasets\_full\_with\_bggcFasta.R # 額外版本：結合 bggcFasta 做標記
│   ├── split\_strategy.R              # 多種切分策略 (random / chrom / groupk)
│   ├── build\_pwm.R                   # 由訓練序列計算一般 PWM
│   ├── build\_null\_pwm\_uniform.R      # 由均勻機率 (0.25) 建構 Uniform PWM (Null Model)
│   ├── evaluate\_models.R             # 一般 PWM + 單一切分測試集之 ROC‐AUC 評估
│   ├── evaluate\_models\_with\_cv.R     # Stratified K‐fold CV 評估 (一般 PWM)
│   ├── evaluate\_models\_with\_groupk.R # GroupKFold (染色體分組 CV) + 一般 PWM 評估
│   ├── evaluate\_models\_with\_groupk\_uniformPWM.R # GroupKFold + Uniform PWM (Null Model) 評估
│   ├── optimize\_threshold.R          # 為單一模型優化分類門檻
│   ├── predict\_ctcf.R                # 使用 PWM (一般或 Uniform) 進行滑動窗掃描並輸出預測
│   ├── visualize\_metrics.R           # 繪製多種切分策略下之 AUC Boxplot
│   ├── run-split-experiment-final.sh  # Pipeline：比較不同切分策略 + 一般 PWM
│   ├── run-split-experiment-uniform.sh # Pipeline：比較不同切分策略 + Uniform PWM
│   └── visualize\_metrics.r           # 已載入於 metrics\_visualization/ 下
└── README.md                         # 本檔案

````

---

## 功能說明

以下分章節說明各個步驟與腳本的功能，並以範例指令示範執行流程。

---

### 1. 資料下載與擷取

1. **scripts/download_data.sh**  
   - 功能：下載示範用的 CTCF ChIP‐seq 峰值（BED 格式）及參考基因組（chr21 demo 或完整 hg38）。  
   - 內部使用 `bedtools getfasta` 擷取 BED 區間對應的 FASTA 序列，輸出至 `data/extracted_sequences.fasta`。  
   - 執行範例：  
     ```bash
     chmod +x scripts/download_data.sh

     # 僅下載 chr21 Demo
     ./scripts/download_data.sh -d

     # 完整 hg38 (需較長時間，建議於真實環境使用)
     ./scripts/download_data.sh

     # 如果要強制重新下載，可加 -f
     ./scripts/download_data.sh -f
     ```
   - 參數說明：  
     - `-d`：僅下載 chr21 demo 資料、快速完成  
     - `-f`：強制覆蓋原先已下載之檔案  
   - 輸出：  
     - `data/extracted_sequences.fasta`

---

### 2. 序列前處理（Sequence Preprocessing）

2. **scripts/preprocess_sequences.R**  
   - 功能：針對 `data/extracted_sequences.fasta` 進行多種前處理，包括：  
     - 長度過濾（`min_length`、`max_length` 或指定 `target_length`）  
     - N-base 處理（遮罩或去除）  
     - 低複雜度序列過濾（基於熵值計算）  
     - 重複序列/同核苷酸（homopolymer）遮罩或去除  
   - 可透過 `scripts/preprocess_config.json` 以 JSON 方式自訂參數，也可改為命令列參數。  
   - 範例執行：  
     ```bash
     Rscript scripts/preprocess_sequences.R \
       data/extracted_sequences.fasta \
       data/preprocessed_sequences.fasta \
       scripts/preprocess_config.json
     ```
   - 輸出：  
     - `data/preprocessed_sequences.fasta`（前處理後序列）  
     - 序列標頭中會註記每個步驟的處理狀態，方便後續除錯與篩選。  

3. **scripts/preprocess_sequences_with_stats.R**  
   - 功能：與 `preprocess_sequences.R` 類似，但同時會計算並輸出各序列在各步驟之統計資訊（例如長度分布、N‐base 比例、熵值）。  
   - 可視需求選擇使用，並不一定每次都要跑。  

---

### 3. 標記正負例（Prepare Full Labeled Dataset）

4. **scripts/prepare_datasets_full.R**  
   - 功能：以 `data/preprocessed_sequences.fasta` 為基礎，為每條序列標記「正例 (class=1)」。  
   - 建立對應數量之負例，支援以下方式：  
     1. 隨機基因組序列擷取  
     2. 保留二核苷酸分布的洗牌（dinucleotide shuffle）  
     3. 完全隨機洗牌（shuffle）  
   - 最後將正例與負例打散、合併，輸出完整標記資料集（fasta header 中含 `class=0/1`）。  
   - 範例執行：  
     ```bash
     Rscript scripts/prepare_datasets_full.R \
       --input_fasta data/preprocessed_sequences.fasta \
       --output_fasta data/complete_labeled_82.fasta \
       --neg_method dinuc_shuf \
       --neg_ratio 1.0
     ```
   - 參數說明：  
     - `--input_fasta`：前處理後之正向序列  
     - `--output_fasta`：輸出之完整標記資料集  
     - `--neg_method`：負例生成方法（`random`、`shuffle`、`dinuc_shuf`）  
     - `--neg_ratio`：負例與正例比例（1.0 表示 1:1）  
   - 輸出：  
     - `data/complete_labeled_82.fasta`（每條序列的 header 都會有 `class=0` 或 `class=1`）

5. **scripts/prepare_datasets_full_with_bggcFasta.R**  
   - 功能：類似 `prepare_datasets_full.R`，但保留額外格式或與 bggcFasta 整合的版本，視需求使用。  

---

### 4. 資料切分（Split Strategies）

6. **scripts/split_strategy.R**  
   - 功能：針對完整標記的 `data/complete_labeled_82.fasta` 產生三種切分策略所需之訓練/測試檔案（FASTA）。  
   - 支援三種切分方式：  
     1. **random**：隨機將一定比例（預設 80%）歸為訓練集，其餘為測試集。  
     2. **chrom**：依 `--train_chr`、`--val_chr`、`--test_chr` 所指定之染色體名稱，將同一染色體上的序列全部分配到相同集合。  
     3. **groupk**：以染色體名稱作為 group，將所有染色體隨機打亂後拆成 K 個 group，每個 group 依序當測試集，其餘為訓練集，共產生 K 組 fold。  
   - 參數說明：  
     - `--strategy {random|chrom|groupk}`：選擇切分策略  
     - `--input`：輸入 `data/complete_labeled_82.fasta`  
     - `--out_dir`：指定輸出資料夾路徑，會自動建立對應子資料夾與檔案  
     - `--train_ratio`：在 `random` 模式下，訓練集比例（預設 0.8）  
     - `--seed`：亂數種子，確保結果可重現  
     - `--train_chr`、`--val_chr`、`--test_chr`：在 `chrom` 模式下，指定各染色體清單（以逗號分隔）  
     - `--k`：在 `groupk` 模式下，K‐fold 數量（預設 5）  
   - 範例執行：  
     ```bash
     # 1. 隨機切分 (80% train, 20% test)
     Rscript scripts/split_strategy.R \
       --strategy random \
       --input data/complete_labeled_82.fasta \
       --out_dir data/splits/random \
       --train_ratio 0.8 \
       --seed 42

     # 2. 染色體切分：chr1-6 為 train，chr7 為 val，chr8-9 為 test
     Rscript scripts/split_strategy.R \
       --strategy chrom \
       --input data/complete_labeled_82.fasta \
       --out_dir data/splits/chrom \
       --train_chr chr1,chr2,chr3,chr4,chr5,chr6 \
       --val_chr chr7 \
       --test_chr chr8,chr9

     # 3. Group K-Fold (以染色體為 group, K=5)
     Rscript scripts/split_strategy.R \
       --strategy groupk \
       --input data/complete_labeled_82.fasta \
       --out_dir data/splits/groupk \
       --k 5 \
       --seed 123
     ```
   - 輸出：  
     - `data/splits/random/train_sequences.fasta`  
     - `data/splits/random/test_sequences.fasta`  
     - `data/splits/chrom/train_sequences.fasta`、`val_sequences.fasta`、`test_sequences.fasta`  
     - `data/splits/groupk/fold1/train_sequences.fasta`、`fold1/test_sequences.fasta`、…、`foldK/`

---

### 5. 建立 PWM 模型與 Uniform Null Model

7. **scripts/build_pwm.R**  
   - 功能：讀取某個訓練集（例如 `data/splits/random/train_sequences.fasta`），計算每個位置上 A/C/G/T 的出現頻率，加上 pseudocount（如 +1），並將其轉換成 log₂‐likelihood ratio 矩陣（即 PWM）。  
   - 輸出：  
     - RDS 檔：`results/pwm_*.rds`  
   - 範例執行：  
     ```bash
     Rscript scripts/build_pwm.R \
       --input data/splits/random/train_sequences.fasta \
       --output results/pwm_random.rds
     ```

8. **scripts/build_null_pwm_uniform.R**  
   - 功能：基於「均勻機率 (0.25, 0.25, 0.25, 0.25)」對應到每個序列位置，直接生成一個 Uniform PWM，也就是 null model。  
   - 使用者須指定 `--length`（目標序列長度），或者可自行計算訓練集中所有序列的長度並傳入。  
   - 輸出：  
     - `results/uniform_pwm.rds`（RDS 格式）  
   - 範例執行：  
     ```bash
     Rscript scripts/build_null_pwm_uniform.R \
       --length 82 \
       --output results/uniform_pwm.rds
     ```

---

### 6. 模型評估（ROC‐AUC 計算）

9. **scripts/evaluate_models.R**  
   - 功能：對「已建立的 PWM」進行單一切分測試集（例：`data/splits/random/test_sequences.fasta`）的 ROC‐AUC 評估。  
   - 參數：  
     - `--pwm`：指定 RDS 檔案  
     - `--test`：指定測試集 FASTA 檔案  
     - `--out`：指定輸出 CSV 檔案  
   - 範例執行：  
     ```bash
     Rscript scripts/evaluate_models.R \
       --pwm results/pwm_random.rds \
       --test data/splits/random/test_sequences.fasta \
       --out results/metrics_random.csv
     ```

10. **scripts/evaluate_models_with_cv.R**  
    - 功能：對「完整標記資料集」(`data/complete_labeled_82.fasta`) 以 Stratified K‐Fold CV (預設 K=5, 可重複多次) 進行交叉驗證評估。  
    - 參數：  
      - `--fasta`：指定完整標記 FASTA 檔案  
      - `--pwm_dir`：指定存放 RDS 檔（如 `pwm_cv.rds`）的資料夾  
      - `--pwm_files`：指定要評估的 RDS 檔名稱（例如 `pwm_cv.rds`）  
      - `--k`：CV fold 數（預設 5）  
      - `--repeats`：重複次數  
    - 範例執行：  
      ```bash
      Rscript scripts/evaluate_models_with_cv.R \
        --fasta data/complete_labeled_82.fasta \
        --pwm_dir results \
        --pwm_files pwm_cv.rds \
        --k 5 \
        --repeats 3
      ```
    - 執行後預設會在 `results/metrics_cv.csv` 產生 AUC 結果。若需要改名，可加一句：  
      ```bash
      mv results/metrics_cv.csv results/metrics_cv_real.csv
      ```

11. **scripts/evaluate_models_with_groupk.R**  
    - 功能：當您使用 `split_strategy.R --strategy groupk` 生成了多個 fold (`data/splits/groupk/fold{i}`) 後，對每個 fold 的 train/test，分別建立一般 PWM，並計算 ROC‐AUC。  
    - 參數：  
      - `--fold_dir`：指定 GroupKFold 切分後的資料夾（例如 `data/splits/groupk`）  
      - `--out_dir`：指定要把結果輸到哪裡（如 `results`）  
    - 範例執行：  
      ```bash
      Rscript scripts/evaluate_models_with_groupk.R \
        --fold_dir data/splits/groupk \
        --out_dir results
      ```
    - 執行後會在 `results/metrics_groupk.csv` 產生 AUC 結果。

12. **scripts/evaluate_models_with_groupk_uniformPWM.R**  
    - 功能：同上，但針對 Uniform PWM (Null Model) 做 GroupKFold 評估。  
    - 參數：  
      - `--fold_dir`：指定 GroupKFold 切分後的資料夾（例如 `data/splits/groupk/standard`）  
      - `--out_dir`：指定要把結果輸到哪裡（如 `results`）  
    - 範例執行：  
      ```bash
      Rscript scripts/evaluate_models_with_groupk_uniformPWM.R \
        --fold_dir data/splits/groupk/standard \
        --out_dir results
      ```
    - 執行後會在 `results/metrics_groupk_uniform.csv` 產生 AUC 結果。

13. **scripts/optimize_threshold.R**  
    - 功能：針對單一 PWM 模型與指定測試集，計算各種分數門檻下的敏感度、特異度、F₁ 等指標，並找出最佳門檻 (Youden’s J、Sensitivity=Specificity、F1 最佳等)。  
    - 範例執行：  
      ```bash
      Rscript scripts/optimize_threshold.R \
        --test_fasta data/splits/random/test_sequences.fasta \
        --pwm_rds results/pwm_random.rds \
        --method youden \
        --output_json results/threshold_optimization.json
      ```

---

### 7. 預測與結果匯出

14. **scripts/predict_ctcf.R**  
    - 功能：讀取指定的 PWM (一般或 Uniform)，對某段或整個 `data/extracted_sequences.fasta` 進行滑動視窗掃描，計算每一個可能的結合位點分數，依使用者給定的閾值 (threshold) 篩選正例，輸出預測結果 TSV。  
    - 範例執行：  
      ```bash
      Rscript scripts/predict_ctcf.R \
        --input_fasta data/extracted_sequences.fasta \
        --pwm_rds results/pwm_random.rds \
        --threshold 5.0 \
        --output_tsv results/predictions.tsv
      ```

---

### 8. 四種切分策略下的 AUC 比較（Visualize Metrics）

本節說明如何使用 `scripts/visualize_metrics.R`，將先前跑好的各種切分策略與模型評量結果 CSV 合併、計算統計值（平均、標準差），並繪製 AUC 的 Boxplot 圖。

#### 腳本功能

1. **自動掃描**  
   - 以 `--metrics_dir` 指定一個資料夾，程式會在其中尋找所有檔名符合 `--pattern`（預設為 `^metrics_.*\.csv$`）的 CSV 檔。  
   - 在我們的 pipeline 中，`run-split-experiment-final.sh` 與 `run-split-experiment-uniform.sh` 最後會把 `metrics_random.csv`、`metrics_chrom.csv`、`metrics_cv_real.csv`、`metrics_groupk.csv`、以及 Uniform 版本的 `metrics_random_uniform.csv`、`metrics_chrom_uniform.csv`、`metrics_cv_uniform.csv`、`metrics_groupk_uniform.csv` 都放在 `results/` 底下。

2. **計算統計值**  
   - 針對讀入的每個檔案（假設都具有 `AUC` 欄位），程式會先把檔名（去掉前綴和副檔名）當作 `strategy`，再把所有檔案依照 `strategy` 合併。  
   - 接著計算每個 `strategy` 底下的樣本數 (`n`)、平均 AUC (`mean_metric`)、標準差 (`sd_metric`)。

3. **輸出結果**  
   - 產生統計摘要 CSV：`results/summary_auc.csv`  
   - 繪製 Boxplot 圖：`results/auc_boxplot.png`

#### 使用範例

在執行完所有 pipeline 後，統一執行：

```bash
Rscript scripts/visualize_metrics.R \
  --metrics_dir results \
  --out_dir results
````

* `--metrics_dir results`
  指定放置所有 `metrics_*.csv` 檔案的資料夾。
* `--out_dir results`
  指定輸出結果的資料夾，不需預先建立，程式會自動建立。
* 其餘參數（`--pattern`、`--metric_col`）皆使用預設值：

  * `--pattern "^metrics_.*\\.csv$"`
  * `--metric_col "AUC"`

執行完成後，`results/` 目錄底下會包含：

```
results/
├── metrics_random.csv
├── metrics_chrom.csv
├── metrics_cv_real.csv
├── metrics_groupk.csv
├── metrics_random_uniform.csv
├── metrics_chrom_uniform.csv
├── metrics_cv_uniform.csv
├── metrics_groupk_uniform.csv
├── summary_auc.csv
└── auc_boxplot.png
```

* **summary\_auc.csv**（格式範例）

  | strategy        | n  | mean\_metric | sd\_metric |
  | --------------- | -- | ------------ | ---------- |
  | random          | 50 | 0.912        | 0.017      |
  | chrom           | 1  | NaN          | NaN        |
  | cv              | 15 | 0.905        | 0.014      |
  | groupk          | 5  | 0.899        | 0.020      |
  | random\_uniform | 50 | 0.500        | 0.000      |
  | chrom\_uniform  | 1  | NaN          | NaN        |
  | cv\_uniform     | 15 | 0.500        | 0.000      |
  | groupk\_uniform | 5  | 0.500        | 0.000      |

  * `strategy`：由檔名 `metrics_{strategy}.csv` 自動抽取（例如 `random`、`chrom`、`cv`、`groupk`，以及後綴 `_uniform` 的四種 Uniform Null 策略）。
  * `n`：該策略底下的讀入樣本筆數（即每個檔案裡多少行有 `AUC`）。
  * `mean_metric`、`sd_metric`：該策略底下所有 AUC 的平均值與標準差，若 `n=1` 或所有 AUC 為 `NA`，則顯示 `NaN`。

* **auc\_boxplot.png**

  * 橫軸：`strategy`
  * 縱軸：`AUC`
  * 每個 Boxplot 同時顯示 jitter（透明小點），以便觀察各筆 AUC 值。

---

### 9. 一鍵式實驗管線 (Pipeline)

下面提供兩個一鍵式腳本，分別用來比較：「四種切分策略 + 一般 PWM 模型」以及「四種切分策略 + Uniform PWM (Null Model)」。

#### 9.1 比較「四種切分策略 + 一般 PWM 模型」

**腳本檔案：** `scripts/run-split-experiment-final.sh`

```bash
#!/usr/bin/env bash
# run-split-experiment-final.sh
# -------------------------------------------------------------
# 比較四種資料切分策略對 PWM 模型預測效能的影響：
#   1. random split
#   2. chrom split（染色體劃分）
#   3. stratified k-fold CV
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
        --input data/splits/random/train_sequences.fasta \
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
        --input data/splits/chrom/train_sequences.fasta \
        --output results/pwm_chrom.rds \
        | tee "$LOG_DIR/build_pwm_chrom.log"

Rscript scripts/evaluate_models.R \
        --pwm results/pwm_chrom.rds \
        --test data/splits/chrom/test_sequences.fasta \
        --out results/metrics_chrom.csv \
        | tee "$LOG_DIR/eval_chrom.log"


# 3. Stratified CV ----------------------------------------------------------
echo "[3] Stratified CV"
Rscript scripts/build_pwm.R \
        --input data/complete_labeled_82.fasta \
        --output results/pwm_cv.rds \
        | tee "$LOG_DIR/build_pwm_cv.log"

Rscript scripts/evaluate_models_with_cv.R \
        --fasta data/complete_labeled_82.fasta \
        --pwm_dir results \
        --pwm_files pwm_cv.rds \
        --k 5 \
        --repeats 3 \
        | tee "$LOG_DIR/eval_cv.log"

mv results/metrics_cv.csv results/metrics_cv_real.csv


# 4. GroupKFold (染色體分組 CV) ---------------------------------------------
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
```

##### 使用方式

1. **確保已安裝相依套件（R 套件、bedtools 等）**
2. 設定腳本可執行：

   ```bash
   chmod +x scripts/run-split-experiment-final.sh
   ```
3. 執行：

   ```bash
   ./scripts/run-split-experiment-final.sh
   ```
4. 執行結束後，`results/` 目錄下會出現：

   ```
   results/
   ├── pwm_random.rds
   ├── pwm_chrom.rds
   ├── pwm_cv.rds
   ├── metrics_random.csv
   ├── metrics_chrom.csv
   ├── metrics_cv_real.csv
   ├── metrics_groupk.csv
   ├── summary_auc.csv
   └── auc_boxplot.png
   ```

   * `metrics_groupk.csv` 來自 `evaluate_models_with_groupk.R`
   * `summary_auc.csv`、`auc_boxplot.png` 由 `visualize_metrics.R` 自動產生，整合所有 `metrics_*.csv` 統計並作圖。

---

#### 9.2 比較「四種切分策略 + Uniform PWM (Null Model)」

**腳本檔案：** `scripts/run-split-experiment-uniform.sh`

```bash
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
        --pwm results/null_random_uniform.rds \
        --test data/splits/random/test_sequences.fasta \
        --out results/metrics_random_uniform.csv \
        | tee "$LOG_DIR/eval_random_uniform.log"


# 2. Chromosome Split --------------------------------------------------------
echo "[2] Chromosome-based split (Uniform Null PWM)"
# 2-1. 切分原始資料 (正+負)
Rscript scripts/split_strategy.R \
        --strategy chrom \
        --input data/complete_labeled_82.fasta \
        --train_chr chr1,chr2,chr3,chr4,chr5,chr6 \
        --val_chr chr7 \
        --test_chr chr8,chr9 \
        --out_dir data/splits/chrom \
        | tee "$LOG_DIR/split_chrom.log"

# 2-2. 產生 Uniform Null PWM (長度 82)
Rscript scripts/build_null_pwm_uniform.R \
        --length 82 \
        --output results/null_chrom_uniform.rds \
        | tee "$LOG_DIR/build_null_chrom_uniform.log"

# 2-3. 在 Chrom test 集上評估 Uniform Null
Rscript scripts/evaluate_models.R \
        --pwm results/null_chrom_uniform.rds \
        --test data/splits/chrom/test_sequences.fasta \
        --out results/metrics_chrom_uniform.csv \
        | tee "$LOG_DIR/eval_chrom_uniform.log"


# 3. Stratified CV (k-fold CV) -----------------------------------------------
echo "[3] Stratified k-fold CV (Uniform Null PWM)"

# 3-1. 產生 Uniform Null PWM (長度 82)
Rscript scripts/build_null_pwm_uniform.R \
        --length 82 \
        --output results/null_cv_uniform.rds \
        | tee "$LOG_DIR/build_null_cv_uniform.log"

# 3-2. 用 CV 評估 Uniform Null
Rscript scripts/evaluate_models_with_cv.R \
        --fasta data/complete_labeled_82.fasta \
        --pwm_dir results \
        --pwm_files null_cv_uniform.rds \
        --k 5 \
        --repeats 3 \
        | tee "$LOG_DIR/eval_cv_uniform.log"

mv results/metrics_cv.csv results/metrics_cv_uniform.csv


# 4. GroupKFold (染色體為群組的 k-fold CV) ------------------------------------
echo "[4] GroupKFold (chromosome-based CV, Uniform Null PWM)"

# 4-1. 先對原始資料 (正+負) 做 GroupKFold 拆分
Rscript scripts/split_strategy.R \
        --strategy groupk \
        --input data/complete_labeled_82.fasta \
        --out_dir data/splits/groupk/standard \
        --k 5 \
        | tee "$LOG_DIR/split_groupk_standard.log"

# 4-2. Uniform Null 下的 GroupKFold：呼叫專用的 evaluate_models_with_groupk_uniformPWM.R
#      （此腳本會對每個 fold 都直接用 uniform PWM 計算 AUC）
Rscript scripts/evaluate_models_with_groupk_uniformPWM.R \
        --fold_dir data/splits/groupk/standard \
        --out_dir results \
        | tee "$LOG_DIR/eval_groupk_uniform.log"


# Visualize ------------------------------------------------------------------
echo "[✔] Visualizing Uniform Null results"
Rscript scripts/visualize_metrics.R \
        --metrics_dir results \
        --out_dir results/plots_uniform \
        | tee "$LOG_DIR/plot_uniform.log"

echo "[✔] All Uniform Null experiments done."
```

##### 使用方式

1. **請先完成完整標記資料（`data/complete_labeled_82.fasta`）之建置**

   * 也就是把 `prepare_datasets_full.R` 執行完畢，得到 `data/complete_labeled_82.fasta`。
2. **先執行一般 PWM Pipeline**

   ```bash
   chmod +x scripts/run-split-experiment-final.sh
   ./scripts/run-split-experiment-final.sh
   ```
3. **再執行 Uniform Null Pipeline**

   ```bash
   chmod +x scripts/run-split-experiment-uniform.sh
   ./scripts/run-split-experiment-uniform.sh
   ```
4. **結果會依序生成**

   * 在 `results/` 底下原先的：

     ```
     pwm_random.rds
     pwm_chrom.rds
     pwm_cv.rds
     metrics_random.csv
     metrics_chrom.csv
     metrics_cv_real.csv
     metrics_groupk.csv
     summary_auc.csv
     auc_boxplot.png
     ```
   * 再新增：

     ```
     null_random_uniform.rds
     metrics_random_uniform.csv
     null_chrom_uniform.rds
     metrics_chrom_uniform.csv
     null_cv_uniform.rds
     metrics_cv_uniform.csv
     metrics_groupk_uniform.csv
     plots_uniform/
       ├── random_uniform_auc.csv
       ├── chrom_uniform_auc.csv
       ├── cv_uniform_auc.csv
       ├── groupk_uniform_auc.csv
       ├── summary_auc.csv
       └── auc_boxplt.png
     ```

---

## 注意事項

1. **相依套件與環境**

   * R 套件：

     * `BiocManager`, `Biostrings`（處理 FASTA）、
     * `pROC`（計算 ROC‐AUC）、`optparse`（解析命令列參數）、
     * `readr`, `dplyr`, `ggplot2`（繪圖與資料處理）。
   * 系統工具：

     * `bedtools`（擷取 FASTA）、`wget` 或 `curl`、`gunzip` 等。
   * 在不同系統上安裝 R 套件時，請先執行：

     ```r
     if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
     BiocManager::install("Biostrings")
     install.packages(c("pROC","optparse","readr","dplyr","ggplot2"), repos="https://cloud.r-project.org")
     ```

2. **目標序列長度 (Length) 一致性**

   * 若要將 Uniform PWM 與一般 PWM 做直接比較，必須確保 Uniform PWM 的 `--length` 與一般 PWM 訓練所用序列長度一致。例如若 `train_sequences.fasta` 中的所有序列長度均為 82，則 `--length 82`。

3. **檔案路徑與命名**

   * 請確認 `data/`、`scripts/`、`results/` 等資料夾結構與本 README 中一致，若有修改須同步更新參數與路徑。
   * 特別注意「完整標記資料」檔名：`data/complete_labeled_82.fasta`。

4. **Split 策略差異**

   * `random`：隨機切分最簡單，但可能造成同染色體或同區域序列既在 train 也在 test，影響生物學意義。
   * `chrom`：完全依染色體名稱分配，具有實際生物學邏輯。
   * `stratified CV`：對整體標記資料做分層交叉驗證，可衡量模型在不同 subset 上的穩定性。
   * `groupk (chrom 版 K‐fold)`：將整染色體作為一個 group，避免資料跨染色體洩漏，適合基因組範疇分析。

5. **Uniform PWM (Null Model)**

   * 由均勻機率設定的 Uniform PWM 當作基線 (baseline)，其預測 AUC 通常約為 0.5。若一般 PWM 與 Uniform PWM 差距不大，代表此資料集難以區分正負例或模型需更複雜的特徵。

6. **測試集僅含單一類別時導致 AUC 為 NA**

   * 在「Chromosome split」策略下，若所指定的 `--test_chr`（例如只用 chr8, chr9）這些染色體對應的序列全屬於同一類別（全部正例或全部負例），則在該測試集上計算 ROC‐AUC 時會出現 `NA`，因為沒有辦法分辨正負類。
   * 同理在 Uniform Null PWM 的 Chromosome split 比對時，若該測試集也只含一種類別，計算 AUC 依舊會是 `NA`。
   * 這是資料本身的「天生限制」，並不影響整體實驗結論，因為其他策略 (random、CV、groupk) 都能提供合理的分數比較。

7. **結果解讀與視覺化**

   * `results/summary_auc.csv` 中，若某策略對應的 `n=1` 或所有 AUC 為 `NA`，則會顯示 `mean_metric = NaN`、`sd_metric = NaN`。
   * 這正是因為該測試集只包含單一類別，導致無法計算 AUC。
   * `results/auc_boxplot.png` 依然會顯示其他策略的分布，並且只要至少一個策略有意義的分數，就能完整比較各種策略之間的差異。

---

## 開發與改版建議

1. **更多切分方法**

   * 未來可考慮「基於 GC 含量匹配之背景剪裁 (GC‐matched negative sampling)」、「功能區域切分 (functional partition)」等，探討不同生物學前提下之評估差異。

2. **負例生成強化**

   * 增加「依據已知 CTCF‐depleted 區域擷取非結合序列」或「不同 species 間比對」等，更具生物學意義的負例策略。

3. **結果報表自動化**

   * 可使用 RMarkdown、Shiny App 或 Jupyter Notebook 產出動態報表，自動顯示各種比較圖表、AUC 值、門檻最佳化結果，便於分享與重現。

4. **參數化管理**

   * 將所有腳本所需參數集中管理（例如放在一個 `config.yaml` 或 JSON），並在 pipeline 中統一讀取，方便團隊成員快速部署與共享。

5. **擴充模型**

   * 除了 PWM，後續可進一步加入 Deep Learning、SVM、Random Forest 等模型，並比較其與 PWM/Uniform PWM 在同樣資料切分策略下的差異。

---

**最後更新：2025-06-01**


```
```
