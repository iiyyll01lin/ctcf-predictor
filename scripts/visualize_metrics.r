# --- 1. auto-install & load ----------------------------------------------------
pkgs <- c("readr", "dplyr", "ggplot2", "optparse")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))

# --- 2. parse CLI options ------------------------------------------------------
option_list <- list(
  make_option("--metrics_dir", type="character", default="results"),
#   make_option("--pattern",     type="character", default="^(metrics_.*|cv_evaluation_results)\\.csv$",
#             help = "檔名 regex pattern"),
  make_option("--pattern", type="character", default="^metrics_.*\\.csv$"),
  make_option("--out_dir",     type="character", default="results"),
  make_option("--metric_col",  type="character", default="AUC")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$metrics_dir))
  stop("metrics_dir not found: ", opt$metrics_dir)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. load & merge -----------------------------------------------------------
files <- list.files(opt$metrics_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) == 0)
  stop("找不到符合 pattern 的 metrics CSV 檔案")

df_list <- lapply(files, function(f) {
  tbl <- readr::read_csv(f, show_col_types = FALSE)

  # 如果沒有 AUC，就跳過
  if (!opt$metric_col %in% names(tbl)) {
    message(" – Skip ", basename(f), " (no '", opt$metric_col, "' column)")
    return(NULL)
  }

  tbl %>% dplyr::mutate(strategy = gsub("^metrics_(.*)\\.csv$", "\\1", basename(f)))
})

# df <- bind_rows(lapply(files, function(f) {
#   tbl <- read_csv(f, show_col_types = FALSE)
#   if (!opt$metric_col %in% names(tbl))
#     stop("在 ", basename(f), " 找不到欄位 '", opt$metric_col, "'")
#   tbl %>% mutate(strategy = gsub("^metrics_(.*)\\.csv$", "\\1", basename(f)))
# }))
df <- dplyr::bind_rows(df_list)   # ← 直接把剛剛的 list 合併

# --- 4. summary ----------------------------------------------------------------
summary_df <- df %>%
  group_by(strategy) %>%
  summarise(n = n(),
            mean_metric = mean(.data[[opt$metric_col]]),
            sd_metric   = sd(.data[[opt$metric_col]]),
            .groups = "drop")

metric_name <- tolower(opt$metric_col)
summary_path <- file.path(opt$out_dir, paste0("summary_", metric_name, ".csv"))
write_csv(summary_df, summary_path)

# --- 5. plot -------------------------------------------------------------------
plot_path <- file.path(opt$out_dir, paste0(metric_name, "_boxplot.png"))
p <- ggplot(df, aes(x = strategy, y = .data[[opt$metric_col]])) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  labs(title = paste(opt$metric_col, "comparison of split strategies"),
       x = "Split strategy", y = opt$metric_col) +
  theme_minimal(base_size = 13)
ggsave(plot_path, p, width = 6, height = 4, dpi = 300)

cat("✔ Finished!\n  • Summary :", summary_path,
    "\n  • Plot    :", plot_path, "\n")
