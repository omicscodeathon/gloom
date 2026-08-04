
## `models/`

This folder stores trained machine-learning models.

### Purpose

Training can take time. This folder keeps the trained models so you can reuse them.

### Typical contents

- `.joblib` model files
- `best_model.joblib` / `best_model_calibrated.joblib`
- `best_model_uncalibrated.joblib`
- `best_model_name.txt`
- `cv_results.csv`

### When to use this folder

Use this folder when:

- you want to score new data,
- you want to reload a trained model,
- you want to compare models or deploy one later.

## `/models/` Directory

Contains trained machine learning models, evaluation results, and preprocessing objects used in the analysis pipeline.

| File Name                            | Description                                                                                                     | Usage & Key Differences                                                                                                            |
| ------------------------------------ | --------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `best_model.joblib`                | Current deployment-ready best model selected in this folder.                                                    | Identical to `best_model_calibrated.joblib` and `model_xgboost.joblib`; a calibrated `xgboost.XGBClassifier` for final inference. |
| `best_model_calibrated.joblib`     | Explicit calibrated copy of the selected best model.                                                            | Same artifact as `best_model.joblib`; kept as a clear calibrated export.                                                           |
| `best_model_name.txt`              | Text file indicating the selected best model name.                                                              | Current value: `xgboost`.                                                                                                          |
| `best_model_uncalibrated.joblib`   | Separate uncalibrated model artifact saved for comparison/debugging.                                            | Loads as `GradientBoostingClassifier`; not the same artifact as the calibrated best model.                                         |
| `cv_results.csv`                   | Five-fold cross-validation summary with AUROC, AUPRC, F1, and MCC for each candidate model.                    | Useful for tradeoff analysis: `xgboost` leads mean AUPRC, `gradient_boosting` mean AUROC, `extra_trees` mean F1, `random_forest` mean MCC. |
| `model_extra_trees.joblib`         | Calibrated Extra Trees classifier.                                                                              | Strong tree-ensemble option; highest mean F1 in `cv_results.csv`.                                                                  |
| `model_gradient_boosting.joblib`   | Calibrated Gradient Boosting classifier.                                                                        | Best mean AUROC in `cv_results.csv`.                                                                                                |
| `model_logistic_regression.joblib` | Calibrated Logistic Regression baseline model.                                                                  | Most interpretable baseline, but lowest overall scores among the saved candidates.                                                  |
| `model_random_forest.joblib`       | Calibrated Random Forest classifier.                                                                            | Best mean MCC in `cv_results.csv`.                                                                                                  |
| `model_xgboost.joblib`             | Calibrated XGBoost classifier.                                                                                  | Current best-model artifact; identical to `best_model.joblib`.                                                                     |
| `robust_scaler.joblib`             | `RobustScaler` preprocessing object using median and IQR scaling.                                               | Load it when the inference pipeline expects the same robust scaling applied during training.                                        |

---

### Usage Tips

- To load a model for prediction or analysis, use joblib's `load()` function in Python:

  ```python
  import joblib
  model = joblib.load('models/best_model.joblib')
  ```
- The current best-model marker in `best_model_name.txt` is `xgboost`.
- The `cv_results.csv` provides detailed cross-validation metrics for comparing ranking criteria, not just one overall winner.

### Usage Summary

- These saved models include a baseline (`logistic_regression`) and several calibrated tree-based ensembles (`gradient_boosting`, `extra_trees`, `random_forest`, `xgboost`).
- The current deployment choice is `best_model.joblib`, which resolves to the calibrated XGBoost artifact.
- The validation table shows that "best" depends on the metric: `xgboost` wins mean AUPRC, `gradient_boosting` wins mean AUROC, `extra_trees` wins mean F1, and `random_forest` wins mean MCC.
- Use the models to generate predictions on new data, compare performance, or interpret feature importance based on your analysis goals.
