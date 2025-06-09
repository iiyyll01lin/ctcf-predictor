#!/usr/bin/env Rscript
# build_null_pwm_uniform.R — 產生 Uniform PWM (4×L 機率皆 0.25)，輸出 RDS
# 用法：Rscript build_null_pwm_uniform.R --length <L> --output <out_rds>

suppressPackageStartupMessages({
  # 不需要額外套件，純 R base 就能做
})

# 解析參數
args <- commandArgs(trailingOnly = TRUE)
# 預設
L <- NA
output_rds <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--length") {
    L <- as.integer(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--output") {
    output_rds <- args[i + 1]
    i <- i + 2
  } else {
    stop("Unknown option: ", args[i])
  }
}

if (is.na(L) || is.null(output_rds)) {
  cat("用法：Rscript build_null_pwm_uniform.R --length <L> --output <out_rds>\n")
  stop("缺少參數：必須指定 --length 以及 --output。")
}

cat("[build_null_pwm_uniform] motif length =", L, "\n")
# 產生一個 4×L 的矩陣，行名為 A,C,G,T，欄名 Pos1~PosL，所有元素皆 0.25
bases <- c("A", "C", "G", "T")
pwm <- matrix(0.25, nrow = 4, ncol = L,
              dimnames = list(bases, paste0("Pos", seq_len(L))))
cat("[build_null_pwm_uniform] PWM 維度：", dim(pwm)[1], "×", dim(pwm)[2], "\n")

cat("[build_null_pwm_uniform] 儲存成：", output_rds, "\n")
saveRDS(pwm, file = output_rds)

cat("[build_null_pwm_uniform] 完成。\n")
