# Stage 3: MULTIPLE linear regression on the binary target.
#   M1  age + blood pressure + cholesterol + BMI        (clinical numbers)
#   M2  + lab markers (triglycerides, CRP, homocysteine, fasting sugar, sleep)
#   M3  + lifestyle (exercise, smoking, alcohol, stress, sugar, sex, family hx)
#   M4  everything                                      (full model)
# Baseline: always predict No -> accuracy = prevalence.
# Console: Rscript 03_multiple_regression.R > multiple_regression_report.txt
source("C:/Users/b7993/heart_disease_eda/common.R")
df <- read_data()

feat <- setdiff(names(df), TARGET)
df <- df[stats::complete.cases(df[, c(TARGET, feat)]), ]
sp <- split_rows(nrow(df))
tr <- df[sp$train, ]
te <- df[sp$test, ]
cat("rows after dropping NA:", nrow(df), "  split:", length(sp$train), "/", length(sp$test), "\n")
base_rate <- mean(te[[TARGET]])
cat(sprintf("test prevalence = %.4f  <- the do-nothing accuracy to beat\n", base_rate))

CLIN  <- c("Age", "Blood Pressure", "Cholesterol Level", "BMI")
LABS  <- c("Triglyceride Level", "Fasting Blood Sugar", "CRP Level",
           "Homocysteine Level", "Sleep Hours")
LIFE  <- c("Exercise Habits", "Smoking", "Alcohol Consumption", "Stress Level",
           "Sugar Consumption", "Gender", "Family Heart Disease", "Diabetes",
           "High Blood Pressure", "Low HDL Cholesterol", "High LDL Cholesterol")

auc_of <- function(y, score) {
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fit_report <- function(name, f) {
  m <- lm(f, data = tr)
  p <- pmin(pmax(predict(m, te), 0), 1)
  mt <- clf_metrics(te[[TARGET]], p)
  cat("\n", SEP, "\n", name, "\n", SEP, "\n")
  cat("formula:", deparse(f), "\n")
  cat(sprintf("n coefficients: %d\n", length(coef(m))))
  print(summary(m))
  cat("\ncoefficient table:\n")
  print(coef(summary(m)), digits = 4)
  cat(sprintf("\nTEST  accuracy = %.4f  precision = %s  recall = %s  f1 = %s  AUC = %.4f  RMSE = %.4f\n",
              mt["accuracy"], round(mt["precision"], 4), round(mt["recall"], 4),
              round(mt["f1"], 4), auc_of(te[[TARGET]], p),
              sqrt(mean((te[[TARGET]] - p)^2))))
  cat(sprintf("vs baseline %.4f -> %+.4f\n", base_rate, mt["accuracy"] - base_rate))
  cat(sprintf("TRAIN accuracy = %.4f\n", mean(round(pmin(pmax(predict(m, tr), 0), 1)) == tr[[TARGET]])))
  cat(sprintf("\nsignificant at p<0.05: %d of %d coefficients\n",
              sum(coef(summary(m))[, 4] < 0.05, na.rm = TRUE), nrow(coef(summary(m))) - 1))
  cat("VIF (>5 problem, >10 severe):\n")
  vv <- tryCatch(vif_table(f, tr), error = function(e) NULL)
  if (!is.null(vv)) print(round(sort(vv, decreasing = TRUE)[1:min(8, length(vv))], 2))
  m$test_metrics <- mt
  m$pred <- p
  m$auc <- auc_of(te[[TARGET]], p)
  m
}

cat("\n", SEP, "\n0. MODELS COMPARED\n", SEP, "\n")
cat("M1  clinical:", paste(CLIN, collapse = " + "), "\n")
cat("M2  + lab markers\n")
cat("M3  + lifestyle and history (full)\n")
cat("M4  everything (M2 features + M3 features)\n")

# every term needs backticks: these column names all contain spaces
rhs <- function(vars) paste(sprintf("`%s`", vars), collapse = " + ")
f1 <- as.formula(paste(sprintf("`%s`", TARGET), "~", rhs(CLIN)))
f2 <- as.formula(paste(sprintf("`%s`", TARGET), "~", rhs(c(CLIN, LABS))))
f3 <- as.formula(paste(sprintf("`%s`", TARGET), "~", rhs(c(CLIN, LIFE))))
f4 <- as.formula(paste(sprintf("`%s`", TARGET), "~", rhs(c(CLIN, LABS, LIFE))))

m1 <- fit_report("M1_clinical", f1)
m2 <- fit_report("M2_clinical_plus_labs", f2)
m3 <- fit_report("M3_clinical_plus_lifestyle", f3)
m4 <- fit_report("M4_everything", f4)

cat("\n", SEP, "\nSUMMARY\n", SEP, "\n")
res <- data.frame(
  model      = c("baseline: always No", "M1 clinical", "M2 +labs", "M3 +lifestyle", "M4 everything"),
  n_features = c(0, length(CLIN), length(CLIN) + length(LABS), length(CLIN) + length(LIFE),
                 length(CLIN) + length(LABS) + length(LIFE)),
  # always predicting No is correct for the 78% of rows with no disease
  accuracy   = c(1 - base_rate, m1$test_metrics["accuracy"], m2$test_metrics["accuracy"],
                  m3$test_metrics["accuracy"], m4$test_metrics["accuracy"]),
  AUC        = c(0.5, m1$auc, m2$auc, m3$auc, m4$auc),
  RMSE       = c(sqrt(base_rate * (1 - base_rate)),
                 sapply(list(m1, m2, m3, m4), function(m) sqrt(mean((te[[TARGET]] - m$pred)^2)))),
  recall     = c(0, m1$test_metrics["recall"], m2$test_metrics["recall"],
                 m3$test_metrics["recall"], m4$test_metrics["recall"]))
print(res, row.names = FALSE, digits = 4)
cat("\naccuracy gain over baseline:\n")
for (i in 2:nrow(res))
  cat(sprintf("  %-26s %+.4f\n", res$model[i], res$accuracy[i] - res$accuracy[1]))

cat("\n", SEP, "\nSIGN INTERPRETATION (M4, full model)\n", SEP, "\n")
cf <- coef(summary(m4))[-1, ]
print(head(cf[order(-abs(cf[, 1])), ], 15), digits = 4)
cat("\nsmallest p-values in the whole model:\n")
print(head(cf[order(cf[, 4]), ], 5), digits = 4)
cat("\nreference levels (first declared grade of each factor):\n")
cat("  Exercise Habits =", levels(tr$`Exercise Habits`)[1],
    " Stress Level =", levels(tr$`Stress Level`)[1],
    " Sugar =", levels(tr$`Sugar Consumption`)[1],
    " Alcohol =", levels(tr$`Alcohol Consumption`)[1], "\n")
cat("so a positive StressLevelHigh coefficient means high stress raises the fitted rate\n")
cat(sprintf("\nmax |coefficient| outside the intercept: %.4f (age is on a 1-80 scale, this is nothing)\n",
            max(abs(cf[rownames(cf) != "(Intercept)", 1]))))

cat("\n", SEP, "\nANOVA: does any block earn its place?\n", SEP, "\n")
cat("M1 vs M4:\n"); print(anova(m1, m4))
cat("\nM1 vs M3 (lifestyle block only):\n"); print(anova(m1, m3))

# ---------------------------------------------------------------- visuals
cat("\n", SEP, "\nPLOTS\n", SEP, "\n")

# 1 accuracy vs baseline
png("images_r/multiple_regression", "model_accuracy_vs_baseline",
    ggplot(res, aes(accuracy, factor(model, levels = rev(model)), fill = model)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.4f", accuracy)), hjust = -0.15, size = 3.6) +
      geom_vline(xintercept = base_rate, color = "#C44E52", linetype = "dashed", linewidth = .9) +
      scale_fill_brewer(palette = "Set2") +
      scale_x_continuous(limits = c(0.76, 0.81)) +
      labs(title = "Test accuracy: the always-No baseline (red dashed) ties every model",
           x = "accuracy", y = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "none"), 8.5, 4.5)

# 2 AUC vs baseline
png("images_r/multiple_regression", "model_auc_vs_baseline",
    ggplot(res, aes(AUC, factor(model, levels = rev(model)), fill = model)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.4f", AUC)), hjust = -0.15, size = 3.6) +
      geom_vline(xintercept = 0.5, color = "#C44E52", linewidth = .9) +
      scale_fill_brewer(palette = "Set1") +
      scale_x_continuous(limits = c(0.47, 0.55)) +
      labs(title = "AUC vs the 0.5 coin flip - the full model reaches 0.516",
           x = "AUC (test set)", y = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "none"), 8.5, 4.5)

# 3 RMSE
png("images_r/multiple_regression", "model_rmse",
    ggplot(res, aes(RMSE, factor(model, levels = rev(model)), fill = model)) +
      geom_col() +
      geom_text(aes(label = round(RMSE, 4)), hjust = -0.15, size = 3.6) +
      scale_fill_brewer(palette = "Set3") +
      labs(title = "RMSE (lower is better)", x = "RMSE", y = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "none"), 8.5, 4.5)

# 4 coefficients of M4 with CI
cfd <- data.frame(term = rownames(cf), est = cf[, 1], se = cf[, 2], pvalue = cf[, 4])
png("images_r/multiple_regression", "coefficients_with_ci",
    ggplot(cfd, aes(est, reorder(term, est))) +
      geom_vline(xintercept = 0, color = "grey50") +
      geom_errorbar(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se),
                    orientation = "y", width = .2, color = "#4C72B0") +
      geom_point(color = "#C44E52", size = 1.8) +
      labs(title = "M4 coefficients with 95% CI - most of them cross zero",
           x = "coefficient", y = NULL) +
      theme_minimal(base_size = 10), 10, 8)

# 5 p-value strip: nothing is significant
png("images_r/multiple_regression", "coefficient_pvalues",
    ggplot(cfd, aes(x = pmax(pvalue, 1e-30), y = reorder(term, pvalue))) +
      geom_point(color = "#C44E52", size = 2) +
      geom_vline(xintercept = 0.05, color = "#4C72B0", linewidth = .9, linetype = "dashed") +
      scale_x_log10() +
      labs(title = "p-values of every M4 coefficient (blue line = 0.05); min p is nowhere near it",
           x = "p-value (log scale)", y = NULL) +
      theme_minimal(base_size = 10), 10, 8)

# 6 actual vs predicted, M4
png("images_r/multiple_regression", "actual_vs_predicted_m4",
    ggplot(data.frame(actual = te[[TARGET]], pred = m4$pred), aes(actual, pred)) +
      geom_jitter(height = .05, alpha = .15, color = "#4C72B0", size = .8) +
      geom_boxplot(aes(group = cut(actual, 2)), width = .3, outlier.shape = NA, alpha = .6) +
      labs(title = "M4 predicted probability by actual class - both classes look the same",
           x = "actual heart disease (0/1)", y = "predicted probability") +
      theme_minimal(base_size = 12), 8, 5)

# 7 calibration
cal <- te %>% mutate(bin = ntile(m4$pred, 10)) %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(m4$pred), obs = mean(`Heart Disease Status`), n = n())
png("images_r/multiple_regression", "calibration_m4",
    ggplot(cal, aes(mean_pred, obs)) +
      geom_abline(color = "#C44E52", linewidth = .8) +
      geom_line(color = "grey50", linewidth = .5) +
      geom_point(size = 3, color = "#4C72B0") +
      geom_text(aes(label = n), vjust = -1, size = 3) +
      labs(title = "Calibration: predicted vs observed in 10 bins - a flat line at 0.20",
           x = "mean predicted probability", y = "observed rate") +
      theme_minimal(base_size = 12), 8, 5)

# 8 residual spread M4
png("images_r/multiple_regression", "residuals_m4",
    ggplot(data.frame(fit = m4$pred, res = te[[TARGET]] - m4$pred), aes(fit, res)) +
      geom_hline(yintercept = 0, color = "#C44E52") +
      geom_point(alpha = .2, size = .8, color = "#4C72B0") +
      labs(title = sprintf("M4 residuals: sd %.4f - identical to a constant-prediction model",
                           sd(te[[TARGET]] - m4$pred)),
           x = "fitted probability", y = "residual") +
      theme_minimal(base_size = 12), 8, 5)

# 9 predicted vs observed bins, M4
binned <- te %>% mutate(bin = ntile(m4$pred, 10)) %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(m4$pred), obs = mean(`Heart Disease Status`), n = n())
png("images_r/multiple_regression", "pred_vs_obs_bins_m4",
    ggplot(binned, aes(mean_pred, obs)) +
      geom_abline(color = "#C44E52", linewidth = .8) +
      geom_point(size = 3, color = "#4C72B0") +
      geom_text(aes(label = n), vjust = -1, size = 3) +
      labs(title = "Predicted vs observed rate per decile of prediction",
           x = "mean predicted", y = "observed rate") +
      theme_minimal(base_size = 12), 8, 5)

# 10 confusion matrix heatmap
cmt <- as.data.frame(table(actual = te[[TARGET]], predicted = as.integer(m4$pred >= 0.5)),
                     responseName = "count")
png("images_r/multiple_regression", "confusion_matrix",
    ggplot(cmt, aes(x = factor(actual), y = factor(predicted), fill = factor(predicted))) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = count), size = 5) +
      scale_fill_manual(values = c("0" = "#4C72B0", "1" = "#C44E52")) +
      labs(title = "M4 confusion matrix at a 0.5 cut - it predicts No for everyone",
           x = "actual", y = "predicted", fill = "predicted") +
      theme_minimal(base_size = 12) + theme(legend.position = "none"), 7, 5)

# 11 partial dependence on the two strongest-looking features
pd <- expand.grid(Age = seq(18, 80, length.out = 100), BMI = c(22, 29, 36))
pd$prob <- sapply(seq_len(nrow(pd)), function(i) {
  nd <- tr[rep(1, 1), ]
  nd$Age <- pd$Age[i]; nd$BMI <- pd$BMI[i]
  pmin(pmax(predict(m4, nd), 0), 1)
})
png("images_r/multiple_regression", "partial_dependence",
    ggplot(pd, aes(Age, prob, colour = factor(BMI))) +
      geom_line(linewidth = 1) +
      scale_color_brewer(palette = "Set2") +
      labs(title = "M4 partial dependence: age at BMI 22 / 29 / 36 - three flat lines",
           x = "age", y = "predicted probability", color = "BMI") +
      theme_minimal(base_size = 12), 8, 5)

# 12 VIF
vv <- tryCatch(vif_table(f4, tr), error = function(e) NULL)
if (!is.null(vv)) {
  png("images_r/multiple_regression", "vif",
      ggplot(data.frame(term = names(vv), vif = as.numeric(vv)), aes(vif, reorder(term, vif))) +
        geom_col(fill = "#55A868") +
        geom_vline(xintercept = 5, color = "#C44E52", linewidth = .8) +
        geom_text(aes(label = round(vif, 2)), hjust = -0.2, size = 3.2) +
        scale_x_continuous(expand = expansion(mult = c(0, .15))) +
        labs(title = "VIF of M4 (red line = 5)",
             x = "VIF", y = NULL) +
        theme_minimal(base_size = 10), 9, 8)
}

stopifnot(length(list.files(file.path(BASE, "images_r/multiple_regression"), pattern = "[.]png$")) >= 11)
cat("\nMultiple regression done. M4 accuracy =", round(m4$test_metrics["accuracy"], 4),
    " AUC =", round(m4$auc, 4), "\n")
