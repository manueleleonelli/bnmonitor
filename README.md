
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bnmonitor

`bnmonitor` is a package for sensitivity analysis and robustness in
Bayesian networks (BNs). If you use the package in your work please
consider citing as

``` r
citation("bnmonitor")
#> To cite package 'bnmonitor' in publications use:
#> 
#>   Leonelli M, Ramanathan R, Wilkerson RL (2023). "Sensitivity and
#>   robustness analysis in Bayesian networks with the bnmonitor R
#>   package." _Knowledge-Based Systems_, *278*, 110882.
#>   doi:10.1016/j.knosys.2023.110882
#>   <https://doi.org/10.1016/j.knosys.2023.110882>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Sensitivity and robustness analysis in {Bayesian} networks with the bnmonitor R package},
#>     author = {Manuele Leonelli and Ramsiya Ramanathan and Rachel L. Wilkerson},
#>     journal = {Knowledge-Based Systems},
#>     year = {2023},
#>     volume = {278},
#>     pages = {110882},
#>     doi = {10.1016/j.knosys.2023.110882},
#>   }
```

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
into the following main areas:

- **Parametric sensitivity analysis**: Investigate the effect of changes
  in some of the parameter values in a Bayesian network and quantify the
  difference between the original and perturbed Bayesian networks using
  dissimilarity measures (both for discrete and Gaussian BNs).

- **Robustness to data**: Verify how well a Bayesian network fits a
  specific dataset that was used either for learning or for testing
  (only for discrete BNs).

- **Node influence**: Quantify how much the nodes of a Bayesian network
  influence an output node of interest (only for discrete BNs).

- **Edge strength**: Assess the strength of the edges of a Bayesian
  network (only for discrete BNs).

- **Other investigations**: Including the diameter of the conditional
  probability tables, measures of asymmetric independence, and level
  amalgamation.

<!-- The prequential diagnostics examine the forecasts that flow from a model in sequence. -->
<!-- Each monitor given below indicates the probability of a particular observation based on the previous observations and the model structure.  -->
<!-- In the prequential mindset, we compute a probability of each subsequent observation based on all previous data points.  -->
<!-- These observations are then scored, and in this package we use the logarithmic score function. -->
<!-- The observations are then standardized to give a z-score statistic.  -->
<!-- Following the recommendation of Cowell (2007), scores indicate a poor fit where |z| > 1.96  -->
<!-- We demonstrate the efficacy of the prequential monitors with the Asia data set from the bnlearn package. Details of the variables (nodes) can be found in the documentation for bnlearn. -->

Refer to the articles section for case studies showcasing the use of the
`bnmonitor` functions.

## Papers where bnmonitor is used

- Görgen, C., & Leonelli, M. (2020). Model-preserving sensitivity
  analysis for families of Gaussian distributions. Journal of Machine
  Learning Research, 21(84), 1-32.

- Leonelli, M., & Riccomagno, E. (2022). A geometric characterization of
  sensitivity analysis in monomial models. International Journal of
  Approximate Reasoning, 151, 64-84.

- Leonelli, M., Ramanathan, R., & Wilkerson, R. L. (2023). Sensitivity
  and robustness analysis in Bayesian networks with the bnmonitor R
  package. Knowledge-Based Systems, 278, 110882.

- Leonelli, M., Smith, J. Q., & Wright, S. K. (2024). The diameter of a
  stochastic matrix: A new measure for sensitivity analysis in Bayesian
  networks. arXiv preprint arXiv:2407.04667.
