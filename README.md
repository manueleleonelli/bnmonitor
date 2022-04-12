
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bnmonitor

`bnmonitor` is a package for sensitivity analysis and robustness in
Bayesian networks (BNs).

## Installation

The package `bnmonitor` can be installed from CRAN using the command

``` r
install.packages("bnmonitor")
```

and loaded in R with

``` r
library(bnmonitor)
```

Note that `bnmonitor` requires the package `gRain` which, while on CRAN,
depends on packages that are on Bioconductor both directly and through
the `gRbase` package, which depends on `RBGL`:

``` r
install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz", "RBGL"))
install.packages("gRain")
```

## Overview

`bnmonitor` provides a suite of function to investigate either a
data-learnt or an expert elicited BN. Its functions can be classified
into three main areas:

-   *Robustness in discrete BNs*: checking how well a BN represents a
    dataset;
-   *Sensitivity in discrete BNs*: assessing the effect of changes in
    the discrete BN’s probabilities;
-   *Sensitivity in continuous BNs*: assessing the effect of changes in
    the continuous BN’s probabilities, either in the standard or
    model-preserving framework

<!-- The prequential diagnostics examine the forecasts that flow from a model in sequence. -->
<!-- Each monitor given below indicates the probability of a particular observation based on the previous observations and the model structure.  -->
<!-- In the prequential mindset, we compute a probability of each subsequent observation based on all previous data points.  -->
<!-- These observations are then scored, and in this package we use the logarithmic score function. -->
<!-- The observations are then standardized to give a z-score statistic.  -->
<!-- Following the recommendation of Cowell (2007), scores indicate a poor fit where |z| > 1.96  -->
<!-- We demonstrate the efficacy of the prequential monitors with the Asia data set from the bnlearn package. Details of the variables (nodes) can be found in the documentation for bnlearn. -->

Refer to the articles section for guidance on each of these areas.

## Papers where bnmonitor is used

-   Görgen, Christiane, and Manuele Leonelli. “Model-Preserving
    Sensitivity Analysis for Families of Gaussian Distributions.” J.
    Mach. Learn. Res. 21 (2020): 84-1.

-   Leonelli, Manuele, and Eva Riccomagno. “A geometric characterisation
    of sensitivity analysis in monomial models.” arXiv preprint
    arXiv:1901.02058 (2018).

-   Leonelli, Manuele, Ramsiya Ramanathan, and Rachel L. Wilkerson.
    “Sensitivity and robustness analysis in Bayesian networks with the
    bnmonitor R package.” arXiv preprint arXiv:2107.11785 (2021).
