# SWORD reproducibility code

Code accompanying *"Oblique Random Forests for Regression via Weighted Support Vector
Machine"* (SWORD).

**Authors**: Andrea Carta, Luca Frigau
Department of Economics and Business Sciences, University of Cagliari, Cagliari, Italy

Paper submitted to *Advances in Data Analysis and Classification* (ADAC).

## Repository structure

```
R/                Core SWORD algorithm (SVM.ROT_Functions_ADAC.R)
data/             Canonical input data: simulated datasets, 248 real benchmark
                  datasets, dataset metadata
Simulation/       Simulation study: scripts that reproduce every simulation-based
                  table/figure, plus the (slower) scripts that generate the
                  underlying results from scratch
Benchmark/        Benchmark on the 248 real datasets: same pattern, for SWORD and
                  the 10 competing methods
results_cache/    Light cached results used as input by the scripts above
```

Every `fig#_*.R` / `tab#_*.R` script reproduces one result from the paper directly
from the shipped caches — the number matches the figure/table number in the paper
(an `A` prefix means it's numbered in the Appendix, e.g. `tabA1`). Every `run_*.R`
script regenerates the underlying results from scratch.

| Script | Reproduces |
|---|---|
| `Simulation/fig1_variable_importance.R` | Fig. 1 — OIW-VI vs RF variable importance |
| `Simulation/fig2_tabA2_compute_time.R` | Fig. 2 + Table A2 — tree-building time (+ scaling exponents) |
| `Simulation/tab1_wilcox_nrmse.R` | Table 1 — NRMSE Wilcoxon tests |
| `Simulation/tab2_stability_bootstrap.R` | Table 2 — OIW-VI bootstrap stability |
| `Simulation/tab3_tabA1_wilcox_corstrength.R` | Table 3 + Table A1 — cor/strength Wilcoxon tests |
| `Benchmark/fig3_fig4_benchmark_rmse_rga.R` | Fig. 3 (RMSE) + Fig. 4 (RGA) — 248-dataset benchmark |
| `Benchmark/fig5_benchmark_rgr.R` | Fig. 5 — RGR robustness |
| `Benchmark/analisi_RGE_explainability/fig6_rge_explainability.R` | Fig. 6 — RGE explainability |
| `Benchmark/tab4_summary_realds.R` | Table 4 — benchmark dataset summary |

Figures A1/A2 (algorithm flowcharts) are not data-driven and have no reproduction
script.

## Example

```r
source("R/SVM.ROT_Functions_ADAC.R")
library(data.tree)  # required by mean_variable_importance_SVM_ROT()

set.seed(1)
n <- 200; p <- 5
X <- as.data.frame(matrix(rnorm(n * p), n, p))
y <- 2 * X$V1 - X$V2 + rnorm(n, sd = 0.5)   # V1, V2 informative; V3-V5 noise

mod <- SVM.ROT.RF.OOB(
  Covariates = X, y = y,
  nmin = 5, cp = 0, n_perc = 1, n_topCor = 2, threshold_COR = 1,
  m = 100, rf_var = p, rand_ntopcor = TRUE, relation = "Pearson",
  Weight_Scheme = "scale", type_of_svm = "C-classification",
  cost_C = 1, cost_nu = 0.5, seed_BS = 1, OOB = TRUE, parallel = FALSE
)

mod$RMSE                                    # out-of-bag RMSE
vi <- mean_variable_importance_SVM_ROT(mod)
vi / sum(vi)                                # OIW-VI: V1, V2 dominate, as expected
```

## Environment

R 4.5.2. Package versions with which we have verified that the code reproduces the
results in the paper exactly:

| Package | Version | | Package | Version |
|---|---|---|---|---|
| data.table | 1.18.2.1 | | RRF | 1.9.4.1 |
| WeightSVM | 1.7-16 | | aorsf | 0.1.6 |
| infotheo | 1.2.0.1 | | partykit | 1.2-27 |
| corpcor | 1.6.10 | | extraTrees | 1.0.5 |
| MASS | 7.3-65 | | rJava | 1.0-16 |
| Matrix | 1.7-5 | | e1071 | 1.7-17 |
| data.tree | 1.2.0 | | caret | 7.0-1 |
| stringr | 1.6.0 | | fastDummies | 1.7.5 |
| Metrics | 0.1.4 | | dplyr | 1.2.0 |
| randomForest | 4.7-1.2 | | tidyr | 1.3.2 |
| ODRF | 0.0.5 | | ggplot2 | 4.0.2 |
| foreach | 1.5.2 | | ggpubr | 0.6.3 |
| doParallel | 1.0.17 | | cowplot | 1.2.0 |
| parallel | 4.5.2 (base) | | patchwork | 1.3.2 |
| rpart | 4.1.27 | | DiagrammeR | 1.0.11 |
