# Stage 2: SIMPLE linear regression on the binary target.
#   heart_disease ~ age            (single continuous predictor)
# Also reports the single best categorical predictor, for contrast.
# Console: Rscript 02_linear_regression.R > linear_regression_report.txt
source("common.R")
df <- read_data()

# lm() drops rows with any NA; do it once here so every model sees the same data
feat <- setdiff(names(df), TARGET)
df <- df[stats::complete.cases(df[, c(TARGET, feat)]), ]
cat("rows after dropping NA:", nrow(df), "of 10000\n")

sp <- split_rows(nrow(df))
tr <- df[sp$train, ]
te <- df[sp$test, ]
cat("split:", length(sp$train), "train /", length(sp$test), "test\n")
base_rate <- mean(te[[TARGET]])
cat(sprintf("test prevalence (always predict No): %.4f\n", base_rate))

cat("\n", SEP, "\n1. MODEL: heart_disease = b0 + b1 * age\n", SEP, "\n")
m <- lm(`Heart Disease Status` ~ Age, data = tr)
print(summary(m))
cat("\ncoefficients:\n")
print(coef(summary(m)))
cat(sprintf("\nequation: P(disease) = %.4f + %.6f * age\n", coef(m)[1], coef(m)[2]))
cat("at age 70 the model says", round(predict(m, data.frame(Age = 70)), 4),
    "- above prevalence, but so is every fitted value\n")

p <- predict(m, te)
mt <- clf_metrics(te[[TARGET]], p)
cat("\n", SEP, "\n2. TEST-SET METRICS\n", SEP, "\n")
print(round(mt, 4))
cat(sprintf("\nRMSE  = %.4f\n", sqrt(mean((te[[TARGET]] - p)^2))))
cat(sprintf("AUC   = %.4f  (0.5 = coin flip)\n",
            as.numeric(pROC_auc <- {
              o <- order(p); r <- rank(p)
              n1 <- sum(te[[TARGET]] == 1); n0 <- sum(te[[TARGET]] == 0)
              (sum(r[te[[TARGET]] == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
            })))
cat(sprintf("train prevalence = %.4f\n", mean(tr[[TARGET]])))

cat("\n", SEP, "\n3. DOES ANY SINGLE PREDICTOR DO BETTER?\n", SEP, "\n")
singles <- data.frame()
for (v in feat) {
  if (is.numeric(df[[v]]) && !is.null(df[[v]])) {
    f <- as.formula(sprintf("`%s` ~ `%s`", TARGET, v))
    mi <- lm(f, data = df)
    pi <- predict(mi, df)
    singles <- rbind(singles, data.frame(feature = v, slope = coef(mi)[2],
                                        auc = {
                                          r <- rank(pi); n1 <- sum(df[[TARGET]] == 1); n0 <- sum(df[[TARGET]] == 0)
                                          (sum(r[df[[TARGET]] == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
                                        }))
  }
}
singles <- singles[order(-abs(singles$auc - 0.5)), ]
print(head(singles, 12), row.names = FALSE, digits = 5)
cat("\nbest single-feature AUC is", round(max(abs(singles$auc - 0.5)) + 0.5, 4),
    "- the distance from 0.5 is the whole story here\n")

cat("\n", SEP, "\n4. SINGLE CATEGORICAL PREDICTOR (Smoking)\n", SEP, "\n")
ms <- lm(`Heart Disease Status` ~ Smoking, data = tr)
print(summary(ms))
cat("\nmeans by smoking status:\n")
print(round(tapply(tr[[TARGET]], grade_factor(tr$Smoking), mean, na.rm = TRUE), 4))
cat("\ndisease rate by every Yes/No flag:\n")
for (v in FLAGS) {
  t <- round(tapply(df[[TARGET]], grade_factor(df[[v]]), mean, na.rm = TRUE), 4)
  cat(sprintf("%-24s %s\n", v, paste(names(t), t, sep = "=", collapse = "  ")))
}

cat("\n", SEP, "\n5. PREDICTED PROBABILITY BY AGE\n", SEP, "\n")
grid <- data.frame(Age = c(18, 30, 40, 50, 60, 70, 80))
grid$predicted_prob <- round(predict(m, grid), 4)
grid$lower_ci <- round(predict(m, grid, interval = "confidence")[, 1], 4)
grid$upper_ci <- round(predict(m, grid, interval = "confidence")[, 3], 4)
print(grid, row.names = FALSE)

cat("\n", SEP, "\n6. RESIDUAL DIAGNOSTICS\n", SEP, "\n")
cat(sprintf("residual mean = %.5f (should be ~0)\n", mean(te[[TARGET]] - p)))
cat(sprintf("residual sd   = %.5f (sqrt(p(1-p)) at base rate = %.5f)\n",
            sd(te[[TARGET]] - p), sqrt(base_rate * (1 - base_rate))))
cat(sprintf("fitted range  = %.4f to %.4f\n", min(p), max(p)))
cat(sprintf("p-values of the single predictor: %s\n",
            paste(sprintf("%.3f", coef(summary(m))[, 4]), collapse = ", ")))
cat("a p above 0.05 means the feature is not distinguishable from noise\n")

# ---------------------------------------------------------------- visuals
cat("\n", SEP, "\n7. PLOTS\n", SEP, "\n")

# 1 scatter of disease vs age
png("images_r/linear_regression", "disease_by_age",
    ggplot(tr, aes(Age, `Heart Disease Status`)) +
      geom_jitter(height = .06, alpha = .12, color = "#4C72B0", size = .8) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#C44E52",
                  fill = "#C44E52", alpha = .15) +
      labs(title = sprintf("Simple LR: P(disease) = %.4f + %.6f x age   p = %.3f",
                           coef(m)[1], coef(m)[2], coef(summary(m))[2, 4]),
           x = "age", y = "heart disease (0/1)"), 8, 5)

# 2 single-feature AUC ladder
top12 <- head(singles, 12)
png("images_r/linear_regression", "single_feature_auc",
    ggplot(top12, aes(auc, reorder(feature, auc), fill = auc > 0.5)) +
      geom_col() +
      geom_vline(xintercept = 0.5, color = "#C44E52", linewidth = .9) +
      geom_text(aes(label = sprintf("%.4f", auc)), hjust = -0.15, size = 3.2) +
      scale_fill_manual(values = c("TRUE" = "#4C72B0", "FALSE" = "#C44E52"), guide = "none") +
      scale_x_continuous(limits = c(0.45, 0.56)) +
      labs(title = "AUC of every single predictor (red line = 0.5 coin flip)",
           x = "in-sample AUC", y = NULL) +
      theme_minimal(base_size = 11), 8, 5.5)

# 3 predicted probability curve with CI
png("images_r/linear_regression", "predicted_probability_by_age",
    ggplot(grid, aes(Age, predicted_prob)) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "#4C72B0", alpha = .2) +
      geom_line(color = "#C44E52", linewidth = 1) +
      geom_hline(yintercept = base_rate, linetype = "dashed", color = "grey40") +
      labs(title = "Predicted probability stays pinned near the 0.20 base rate",
           x = "age", y = "predicted probability") +
      theme_minimal(base_size = 12), 8, 5)

# 4 disease rate by age decile
agedf <- tr %>% mutate(dec = ntile(Age, 10)) %>%
  group_by(dec) %>%
  summarise(age_mid = round(median(Age), 1), rate = mean(`Heart Disease Status`), n = n())
png("images_r/linear_regression", "rate_by_age_decile",
    ggplot(agedf, aes(age_mid, rate)) +
      geom_line(color = "#C44E52", linewidth = .9) +
      geom_point(color = "#4C72B0", size = 2.5) +
      geom_hline(yintercept = base_rate, linetype = "dashed", color = "grey40") +
      labs(title = "Observed rate per age decile: noise around 0.20",
           x = "median age in decile", y = "observed disease rate") +
      theme_minimal(base_size = 12), 8, 5)

# 5 rate by every Yes/No flag
fl <- bind_rows(lapply(FLAGS, function(v) rate_table(df, v, TARGET)))
png("images_r/linear_regression", "rate_by_risk_flag",
    ggplot(fl, aes(grade, rate, fill = grade)) +
      geom_col(position = "dodge") +
      facet_wrap(~ feature, ncol = 3) +
      geom_hline(yintercept = base_rate, linetype = "dashed", color = "#C44E52") +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      labs(title = "Disease rate by risk flag - every split lands on the base rate",
           x = NULL, y = "rate") +
      theme_minimal(base_size = 10), 9, 5)

# 6 rate by graded category
gr <- bind_rows(lapply(names(ORDERS), function(v) rate_table(df, v, TARGET)))
png("images_r/linear_regression", "rate_by_graded_category",
    ggplot(gr, aes(grade, rate, fill = grade)) +
      geom_col() +
      facet_wrap(~ feature, ncol = 2) +
      geom_hline(yintercept = base_rate, linetype = "dashed", color = "#C44E52") +
      scale_fill_brewer(palette = "Set3", guide = "none") +
      labs(title = "Disease rate by Low/Medium/High grade - flat everywhere",
           x = NULL, y = "rate") +
      theme_minimal(base_size = 10), 9, 6)

# 7 predicted vs actual, binned
binned <- te %>% mutate(bin = ntile(p, 10)) %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(p), obs_rate = mean(`Heart Disease Status`), n = n())
png("images_r/linear_regression", "predicted_vs_observed_bins",
    ggplot(binned, aes(mean_pred, obs_rate)) +
      geom_abline(color = "#C44E52", linewidth = .8) +
      geom_line(color = "grey50", linewidth = .5) +
      geom_point(size = 3, color = "#4C72B0") +
      geom_text(aes(label = n), vjust = -1, size = 3) +
      labs(title = "Model prediction vs observed rate, in 10 bins - points sit on the base rate",
           x = "mean predicted probability", y = "observed rate") +
      theme_minimal(base_size = 12), 8, 5)

# 8 residual spread
png("images_r/linear_regression", "residual_spread",
    ggplot(data.frame(r = te[[TARGET]] - p, fit = p), aes(fit, r)) +
      geom_hline(yintercept = 0, color = "#C44E52") +
      geom_point(alpha = .2, size = .8, color = "#4C72B0") +
      labs(title = sprintf("Residuals: sd %.4f vs Bernoulli sd %.4f at base rate",
                           sd(te[[TARGET]] - p), sqrt(base_rate * (1 - base_rate))),
           x = "fitted probability", y = "residual") +
      theme_minimal(base_size = 12), 8, 5)

# 9 coefficient with CI
png("images_r/linear_regression", "coefficient_ci",
    ggplot(data.frame(term = "age", est = coef(m)[2], se = coef(summary(m))[2, 2],
                      lo = coef(m)[2] - 1.96 * coef(summary(m))[2, 2],
                      hi = coef(m)[2] + 1.96 * coef(summary(m))[2, 2]),
           aes(est, term, xmin = lo, xmax = hi)) +
      geom_errorbar(orientation = "y", width = .1, color = "#4C72B0", linewidth = 1.2) +
      geom_point(color = "#C44E52", size = 3) +
      geom_vline(xintercept = 0, color = "grey40", linewidth = .8) +
      labs(title = "Age coefficient, 95% CI - it straddles zero",
           x = "coefficient", y = NULL) +
      theme_minimal(base_size = 12), 8, 4)

stopifnot(length(list.files(file.path(BASE, "images_r/linear_regression"), pattern = "[.]png$")) >= 9,
          mt["accuracy"] > 0.75 && mt["accuracy"] < 0.87)
cat("\nSimple linear regression done. accuracy =", round(mt["accuracy"], 4), "\n")
