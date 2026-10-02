# Benchmark results

40 samples; 5 classes; 3 x 4-fold nested cross-validation.

| Method | Mean balanced accuracy |
|---|---:|
| A only | 0.683 |
| B only | 1.000 |
| Concatenated | 0.967 |
| DIABLO | 0.983 |
| DIABLO: A available | 0.533 |
| DIABLO: B available | 0.983 |
| Inner-selected single | 1.000 |
| Panel 10/block/component | 0.983 |
| Panel 5/block/component | 0.983 |

DIABLO minus inner-selected single-omics baseline: -1.7 percentage points.

Repeated estimates share samples and training sets; their spread is descriptive, not a confidence interval.
Panel budgets apply per block per component. Actual panel size is the union across components.
Missing-block scenarios reuse the complete-data-trained DIABLO model without retraining.
Selection frequency measures model reuse across overlapping training sets, not biomarker validation.
Small training classes (<8 observations in this demo) can make glmnet estimates unstable.
These results describe this dataset and search grid; no method is assumed to win.
