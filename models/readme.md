
## `models/`

This folder stores trained machine-learning models.

### Purpose

Training can take time. This folder keeps the trained models so you can reuse them.

### Typical contents

- `.joblib` model files
- `best_model.joblib`
- `best_model_name.txt`
- cross-validation results

### When to use this folder

Use this folder when:

- you want to score new data,
- you want to reload a trained model,
- you want to compare models or deploy one later.

## `/models/` Directory

Contains trained machine learning models, evaluation results, and preprocessing objects used in the analysis pipeline.

| File Name                            | Description                                                                                                     | Usage & Key Differences                                                                                                            |
| ------------------------------------ | --------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `best_model.joblib`                | The overall best-performing model selected based on validation metrics.                                         | Used for final predictions on unseen data. Combines the strengths of multiple algorithms optimized during training.                |
| `best_model_name.txt`              | Text file indicating the name or description of the best model.                                                 | Clarifies which model is considered the best for deployment or further testing.                                                    |
| `cv_results.csv`                   | Cross-validation results including metrics such as AUROC, accuracy, precision, recall, and hyperparameters.     | Used for performance comparison, hyperparameter tuning insights, and model validation.                                             |
| `model_extra_trees.joblib`         | Extra Trees classifier, an ensemble method similar to Random Forest but with more randomness for decorrelation. | Fast training, robust to overfitting, good for high-dimensional data. Suitable when feature importance interpretability is needed. |
| `model_gradient_boosting.joblib`   | Gradient Boosting classifier, builds models sequentially to correct previous errors.                            | Often achieves high accuracy, but may require careful tuning to prevent overfitting. Good for complex patterns.                    |
| `model_logistic_regression.joblib` | Logistic Regression, a linear model for binary classification.                                                  | Simple, interpretable, fast, suitable for baseline or when relationships are linear. Less effective with complex patterns.         |
| `model_random_forest.joblib`       | Random Forest classifier, an ensemble of decision trees built with bootstrap sampling.                          | Generally robust, handles feature interactions well, less prone to overfitting, provides feature importance.                       |
| `model_svm.joblib`                 | Support Vector Machine, finds optimal separating hyperplane with kernel options.                                | Effective in high-dimensional spaces, good with clear margins, but computationally intensive with large datasets.                  |
| `robust_scaler.joblib`             | Scaler that scales features using median and IQR, reducing the impact of outliers.                              | Used to normalize features before model training and prediction to improve model performance.                                      |

---

### Usage Tips

- To load a model for prediction or analysis, use joblib's `load()` function in Python:

  ```python
  import joblib
  model = joblib.load('models/best_model.joblib')
  ```
- The `cv_results.csv` provides detailed cross-validation metrics useful for model evaluation and comparison.

### Usage Summary

- These models are trained and saved for various purposes: some are baseline models (`logistic_regression`), others are more complex ensemble methods (`gradient_boosting`, `extra_trees`, `random_forest`), and SVM for high-dimensional data.
- The `best_model.joblib` is the optimal model selected after validation, suitable for deployment.
- Use the models to generate predictions on new data, compare performance, or interpret feature importance based on your analysis goals.
