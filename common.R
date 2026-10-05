# Shared setup for the heart disease scripts.
.libPaths(c("C:/Users/b7993/diamonds_eda/rlib", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

BASE <- "C:/Users/b7993/heart_disease_eda"
CSV <- "C:/Users/b7993/.cache/kagglehub/datasets/oktayrdeki/heart-disease/versions/1/heart_disease.csv"

SEP <- strrep("=", 70)

# the target, and the graded categoricals in real order (not alphabetical)
TARGET  <- "Heart Disease Status"
ORDERS  <- list(
  `Exercise Habits`     = c("Low", "Medium", "High"),
  `Stress Level`        = c("Low", "Medium", "High"),
  `Sugar Consumption`   = c("Low", "Medium", "High"),
  `Alcohol Consumption` = c("None", "Low", "Moderate", "High")
)
# these are 0/1 Yes-No risk flags, already numeric in spirit
FLAGS <- c("Smoking", "Family Heart Disease", "Diabetes",
           "High Blood Pressure", "Low HDL Cholesterol", "High LDL Cholesterol")

theme_set(theme_minimal(base_size = 12) +
            theme(plot.title = element_text(face = "bold", size = 13)))

# save a ggplot into <stage>/, numbered in write order
png <- function(stage, name, p, w = 8, h = 5) {
  d <- file.path(BASE, stage)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(d, sprintf("%02d_%s.png", length(list.files(d, pattern = "[.]png$")) + 1, name))
  ggsave(f, p, width = w, height = h, dpi = 110, bg = "white")
  cat("  saved:", f, "\n")
  invisible(f)
}

# VIF via auxiliary R2: regress each design-matrix column on all the others.
vif_table <- function(formula, data) {
  X <- as.data.frame(model.matrix(terms(formula), data)[, -1, drop = FALSE])
  X <- X[, sapply(X, function(cc) length(unique(cc)) > 1), drop = FALSE]
  vapply(names(X), function(col) {
    r2 <- summary(lm(X[[col]] ~ ., data = X[setdiff(names(X), col)]))$r.squared
    if (!is.finite(r2) || r2 >= 1) Inf else 1 / (1 - r2)
  }, numeric(1))
}

# classification metrics on a 0/1 response (accuracy, precision, recall, F1, AUC)
clf_metrics <- function(y, prob, cut = 0.5) {
  pred <- as.integer(prob >= cut)
  tp <- sum(pred == 1 & y == 1); tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0); fn <- sum(pred == 0 & y == 1)
  r <- c(accuracy  = (tp + tn) / length(y),
         precision = if (tp + fp) tp / (tp + fp) else NA,
         recall    = if (tp + fn) tp / (tp + fn) else NA,
         f1        = NA, prevalence = mean(y), predicted_rate = mean(pred),
         tp = tp, tn = tn, fp = fp, fn = fn)
  r["f1"] <- if (is.na(r["precision"]) || is.na(r["recall"]) || r["precision"] + r["recall"] == 0)
    NA else 2 * r["precision"] * r["recall"] / (r["precision"] + r["recall"])
  r
}

# disease rate for every level of a factor, including NA levels.
# tapply() drops empty levels and silently misaligns the result, so index it.
# binary 0/1 flags are numeric (correlations need that) - grade_factor() turns
# them back into a labelled factor for the category charts.
grade_factor <- function(x) {
  if (is.factor(x)) return(x)
  factor(ifelse(is.na(x), NA, ifelse(x == 1, "Yes", "No")), levels = c("No", "Yes"))
}

rate_table <- function(data, var, target) {
  v <- grade_factor(data[[var]])
  data.frame(feature = var,
             grade    = levels(v),
             n        = as.integer(table(v)[levels(v)]),
             rate     = as.numeric(tapply(data[[target]], v, function(z) mean(z, na.rm = TRUE))))
}

read_data <- function() {
  cat("dataset:", CSV, "\n")
  d <- readr::read_csv(CSV, show_col_types = FALSE, na = c("", "NA"))
  for (nm in names(ORDERS)) d[[nm]] <- factor(d[[nm]], levels = ORDERS[[nm]])
  # Yes/No flags and the target become 0/1 so lm() can be used directly
  for (nm in c(FLAGS, TARGET))
    d[[nm]] <- as.integer(factor(d[[nm]], levels = c("No", "Yes"))) - 1L
  d <- d %>% mutate(across(where(is.character), ~ factor(.x)))
  d
}

# 80/20 row split, fixed seed, shared by both model scripts
split_rows <- function(n) {
  set.seed(42)
  tr <- sample(n, 0.8 * n)
  list(train = tr, test = setdiff(seq_len(n), tr))
}
