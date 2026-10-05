# Stage 1: EDA on oktayrdeki/heart-disease.
# Console: Rscript 01_eda.R > eda_report.txt
source("C:/Users/b7993/heart_disease_eda/common.R")
df <- read_data()

cat(SEP, "\n1. STRUCTURE\n", SEP, "\n")
cat("rows:", nrow(df), " cols:", ncol(df), "\n")
cat("memory:", format(object.size(df), units = "MB"), "\n")
cat("duplicate rows:", sum(duplicated(df)), "\n\n")
cat("column types:\n"); print(sapply(df, function(x) class(x)[1]))

cat("\n", SEP, "\n2. MISSING DATA\n", SEP, "\n")
nn <- colSums(is.na(df))
print(sort(nn, decreasing = TRUE))
cat("\nrows with at least one NA:", sum(!complete.cases(df)),
    sprintf("(%.1f%%)\n", 100 * mean(!complete.cases(df))))
cat("columns fully missing:", paste(names(nn)[nn == nrow(df)], collapse = ", "), "\n")

cat("\n", SEP, "\n3. FIRST 5 ROWS\n", SEP, "\n")
print(head(as.data.frame(df), 5))

cat("\n", SEP, "\n4. SUMMARY (numeric)\n", SEP, "\n")
print(summary(as.data.frame(df[, sapply(df, is.numeric)])))

cat("\n", SEP, "\n5. SUMMARY (categorical)\n", SEP, "\n")
for (v in names(ORDERS)) {
  cat("\n--", v, "\n"); print(as.data.frame(table(df[[v]], useNA = "ifany")))
}
cat("\n-- Gender\n"); print(table(df$Gender, useNA = "ifany"))

cat("\n", SEP, "\n6. TARGET BALANCE\n", SEP, "\n")
cat(sprintf("heart disease: %d yes / %d no  -> prevalence %.4f\n",
            sum(df[[TARGET]]), sum(df[[TARGET]] == 0), mean(df[[TARGET]], na.rm = TRUE)))
cat("a model predicting 'No' always would score accuracy 0.8000 - that is the floor to beat\n")

cat("\n", SEP, "\n7. FEATURE SIGNAL vs TARGET\n", SEP, "\n")
num <- as.data.frame(df %>% select(where(is.numeric)))
cr <- sort(sapply(num, function(v) cor(v, df[[TARGET]], use = "complete.obs")), decreasing = TRUE)
print(round(cr, 4))
cat("\nfull matrix:\n"); print(round(cor(num), 3))

cat("\n", SEP, "\n8. DISEASE RATE BY CATEGORY\n", SEP, "\n")
grp <- c("Gender", "Smoking", "Diabetes", "Family Heart Disease", "High Blood Pressure",
         "Low HDL Cholesterol", "High LDL Cholesterol", names(ORDERS))
for (v in grp) {
  t <- round(tapply(df[[TARGET]], df[[v]], mean, na.rm = TRUE), 4)
  cat(sprintf("%-22s %s\n", v, paste(names(t), t, sep = "=", collapse = "  ")))
}
cat("\nby age decade:\n")
print(as.data.frame(df %>% mutate(dec = cut(Age, c(17, 30, 40, 50, 60, 70, 90))) %>%
                       group_by(dec) %>%
                       summarise(n = n(), rate = round(mean(!!sym(TARGET), na.rm = TRUE), 4)), row.names = FALSE))

cat("\n", SEP, "\n9. OUTLIERS (IQR rule, 1.5x)\n", SEP, "\n")
fences <- list()
for (v in names(num)) {
  q <- quantile(num[[v]], c(.25, .75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  k <- sum(num[[v]] < q[1] - 1.5 * iqr | num[[v]] > q[2] + 1.5 * iqr, na.rm = TRUE)
  if (k > 0) {
    cat(sprintf("%-22s outliers=%-6d lo=%9.2f hi=%9.2f\n", v, k, q[1] - 1.5 * iqr, q[2] + 1.5 * iqr))
    fences[[v]] <- c(q[1] - 1.5 * iqr, q[2] + 1.5 * iqr)
  }
}

cat("\n", SEP, "\n10. HOW FAR DOES THE TARGET LIVE IN THE FILE ORDER?\n", SEP, "\n")
blk <- (seq_len(nrow(df)) - 1) %/% 1000
print(round(tapply(df[[TARGET]], blk, mean, na.rm = TRUE), 3))
cat("if this is 0 then 0 then 1 the label is positional, not clinical - see the report\n")

# ---------------------------------------------------------------- visuals
cat("\n", SEP, "\n11. PLOTS\n", SEP, "\n")

# 1 class balance
bal <- as.data.frame(table(df[[TARGET]]), responseName = "n")
bal$grp <- factor(c("No", "Yes")[match(as.character(bal[[1]]), c("0", "1"))],
                  levels = c("No", "Yes"))
png("images_r/eda", "class_balance",
    ggplot(bal, aes(grp, n, fill = grp)) +
      geom_col(width = .6) +
      geom_text(aes(label = n), vjust = -0.4, size = 4) +
      labs(title = "80/20 class imbalance - accuracy 0.80 is the do-nothing baseline",
           x = "heart disease", y = "count") +
      theme_minimal(base_size = 12) + theme(legend.position = "none"), 7, 5)

# 2 correlation with target
cp <- data.frame(feature = names(cr), r = as.numeric(cr))
cp <- cp[!cp$feature %in% c(TARGET, "Heart Disease Status"), ]
png("images_r/eda", "corr_with_target",
    ggplot(cp, aes(r, reorder(feature, r), fill = r)) +
      geom_col() + geom_vline(xintercept = 0, color = "grey40") +
      geom_text(aes(label = sprintf("%+.3f", r)), hjust = ifelse(cp$r > 0, -0.2, 1.2), size = 3.5) +
      scale_fill_gradient(low = "#C44E52", high = "#4C72B0", limits = c(-0.05, 0.05)) +
      labs(x = "Pearson correlation with heart disease", y = NULL,
           title = "Every |r| is under 0.02 - no feature carries signal on its own"), 8, 5.5)

# 3 correlation heatmap
cm <- cor(num)
lg <- data.frame(x = rep(colnames(cm), times = ncol(cm)),
                 y = rep(colnames(cm), each = nrow(cm)),
                 value = as.vector(cm))
png("images_r/eda", "corr_heatmap",
    ggplot(lg, aes(x, y, fill = value)) +
      geom_tile(color = "white", linewidth = .3) +
      scale_fill_gradient2(low = "#C44E52", mid = "white", high = "#4C72B0", limits = c(-1, 1)) +
      labs(x = NULL, y = NULL, fill = "r", title = "Correlation matrix") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)), 9, 7.5)

# 4 missing data
png("images_r/eda", "missing_data",
    ggplot(data.frame(col = names(nn), n = as.numeric(nn)), aes(reorder(col, n), n)) +
      geom_col(fill = "#DD8452") +
      geom_text(aes(label = n), hjust = -0.15, size = 3.2) +
      scale_y_continuous(expand = expansion(mult = c(0, .12))) +
      labs(title = "Missing values per column (500 rows incomplete, none fully empty)",
           x = NULL, y = "NA count") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)), 9, 6)

# 5 target rate by categorical feature (the important EDA chart)
cat_rate <- bind_rows(lapply(setdiff(grp, "Gender"), function(v)
  rate_table(df, v, TARGET)))
png("images_r/eda", "disease_rate_by_category",
    ggplot(cat_rate, aes(grade, rate, fill = grade)) +
      geom_col() +
      facet_wrap(~ feature, scales = "free_x", ncol = 4) +
      geom_hline(yintercept = mean(df[[TARGET]], na.rm = TRUE), color = "#C44E52",
                 linetype = "dashed", linewidth = .6) +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      labs(title = "Disease rate per category - every bar sits on the 0.20 line",
           x = NULL, y = "heart disease rate") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)), 11, 7)

# 6 rate vs age
agedf <- df %>% mutate(dec = cut(Age, c(17, 30, 40, 50, 60, 70, 90))) %>%
  group_by(dec) %>% summarise(n = n(), rate = mean(!!sym(TARGET), na.rm = TRUE))
png("images_r/eda", "rate_by_age_decade",
    ggplot(agedf, aes(dec, rate)) +
      geom_col(fill = "#4C72B0") +
      geom_text(aes(label = sprintf("%.3f\nn=%d", rate, n)), vjust = -0.4, size = 3) +
      geom_hline(yintercept = mean(df[[TARGET]], na.rm = TRUE), color = "#C44E52", linetype = "dashed") +
      labs(title = "Disease rate by age decade - flat", x = "age band", y = "rate") +
      theme_minimal(base_size = 12), 8, 5)

# 7 numeric distributions
png("images_r/eda", "numeric_distributions",
    ggplot(pivot_longer(num %>% select(-!!sym(TARGET)), names_to = "var", values_to = "v", cols = everything()) %>%
             filter(!is.na(v)), aes(v)) +
      geom_histogram(bins = 30, fill = "#4C72B0", color = NA) +
      facet_wrap(~ var, scales = "free", ncol = 5) +
      labs(title = "Every continuous feature is a uniform-ish block of noise",
           x = NULL, y = "count") +
      theme_minimal(base_size = 9), 12, 6)

# 8 outlier overlay on the strongest-looking feature
bmi_fit <- lm(`Heart Disease Status` ~ BMI, df)
bmi_slope <- coef(summary(bmi_fit))[2, 1]
png("images_r/eda", "outliers_bmi",
    ggplot(df, aes(BMI, `Heart Disease Status`)) +
      geom_point(alpha = .2, color = "#4C72B0", size = .8) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#C44E52",
                  fill = "#C44E52", alpha = .15) +
      labs(title = sprintf("BMI vs disease: slope %+.5f - flat", bmi_slope),
           x = "BMI", y = "heart disease (0/1)"), 8, 5)

# 9 label vs row position
png("images_r/eda", "label_by_row_position",
    ggplot(data.frame(i = seq_len(nrow(df)), y = df[[TARGET]]), aes(i, y)) +
      geom_point(alpha = .05, size = .5, color = "#4C72B0") +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#C44E52", linewidth = .9) +
      labs(title = "Label vs row index: 0 for the first ~8000 rows, then 1 - the label is positional",
           x = "row number in the file", y = "heart disease (0/1)"), 8, 5)

# 10 categorical counts
png("images_r/eda", "category_counts",
    ggplot(bind_rows(lapply(c(grp, "Gender"), function(v) {
      data.frame(var = v, grade = levels(grade_factor(df[[v]])),
                 n = as.integer(table(grade_factor(df[[v]]))[levels(grade_factor(df[[v]]))]))
    })), aes(grade, n, fill = var)) +
      geom_col() +
      facet_wrap(~ var, scales = "free_x", ncol = 5) +
      labs(title = "Categorical frequencies - all near-uniform", x = NULL, y = "count") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"), 12, 7)

# 11 duplicate check
png("images_r/eda", "duplicate_rows",
    ggplot(data.frame(k = "duplicate rows", v = sum(duplicated(df))), aes(k, v, label = v)) +
      geom_col(fill = "#55A868", width = .4) +
      geom_text(vjust = -0.5, size = 4) +
      scale_y_continuous(expand = expansion(mult = c(0, .15)), limits = c(0, 1)) +
      labs(title = "Exact duplicate rows: 0", x = NULL, y = "count") +
      theme_minimal(base_size = 12), 6, 4)

stopifnot(length(list.files(file.path(BASE, "images_r/eda"), pattern = "[.]png$")) >= 11)
cat("\nEDA done.\n")
