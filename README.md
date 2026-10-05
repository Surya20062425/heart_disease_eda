# Heart Disease EDA + Regression (R)

Exploratory Data Analysis and regression modeling on the [Kaggle Heart Disease dataset](https://www.kaggle.com/datasets/oktayrdeki/heart-disease) (10,000 rows, 21 features).

## Dataset
- **Source**: `oktayrdeki/heart-disease` via `kagglehub`
- **Target**: `Heart Disease Status` (Yes/No, ~20% positive)
- **Features**: 9 numeric, 12 categorical
- **Missing**: ~29% rows have at least one NA; `Alcohol Consumption` has 25.8% missing

## Project Structure
```
heart_disease_eda/
├── common.R                      # Shared: data loading, split, metrics, png()
├── 01_eda.R           -> eda_report.txt              + images_r/eda/
├── 02_linear_regression.R -> linear_regression_report.txt + images_r/linear_regression/
├── 03_multiple_regression.R -> multiple_regression_report.txt + images_r/multiple_regression/
├── heart_disease.csv             # Local copy of dataset
├── .gitignore
└── README.md
```

## Quick Start
```bash
# Prerequisites: R installed (winget install --id RProject.R)

# Run all stages
Rscript 01_eda.R
Rscript 02_linear_regression.R
Rscript 03_multiple_regression.R

# Reports are in *_report.txt, charts in images_r/
```

## Key Findings (EDA)
- No numeric feature has |r| > 0.1 with target — **no linear signal**
- Categorical target rates range 18.7–21.4% — small differences
- All numerics near-uniform (skew ≈ 0, kurt ≈ -1.2) — synthetic appearance
- No multicollinearity (max |r| < 0.7)

## Regression Results
| Model | Predictors | R² / Pseudo-R² | Notes |
|-------|------------|----------------|-------|
| Simple LR | Single best numeric | ~0.000 | No signal |
| Multiple LR (numeric) | 9 numeric | ~0.004 | Negligible |
| Multiple LR + one-hot | All 21 features | ~0.02 | Marginal |

See `*_report.txt` for full diagnostics: VIF, anova block tests, per-segment fit.

## Requirements
- R ≥ 4.0
- Packages: `data.table`, `ggplot2`, `car` (for VIF) — installed on first run

## License
MIT