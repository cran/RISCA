RISCA: an R Package for Causal Inference and Prediction in Cohort-Based Analyses
================

## What is the ‘RISCA’ package?

Description: Numerous functions for cohort-based analyses, either for prediction or causal inference. For causal inference, it includes Inverse Probability Weighting and G-computation for marginal estimation of an exposure effect when confounders are expected. We deal with binary outcomes, times-to-events, competing events, and multi-state data. For multistate data, semi-Markov model with interval censoring may be considered, and we propose the possibility to consider the excess of mortality related to the disease compared to reference lifetime tables. For predictive studies, we propose a set of functions to estimate time-dependent receiver operating characteristic (ROC) curves with the possible consideration of right-censoring times-to-events or the presence of confounders. Finally, several functions are available to assess time-dependent ROC curves or survival curves from aggregated data.

## Installation

Install the latest release from CRAN:

``` r
install.packages("RISCA")
```

Install the development version from GitHub:

``` r
remotes::install_github("foucher-y/RISCA")
```
