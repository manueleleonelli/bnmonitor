
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bnmonitor

`bnmonitor` is a package for sensitivity analysis and robustness in
Bayesian networks.

## Installation

The package `bnmonitor` can be installed from GitHub using the command

``` r
# install.packages("devtools")
devtools::install_github("manueleleonelli/bnmonitor")
```

and loaded in R with

``` r
library(bnmonitor)
```

Note that `bnmonitor` requires the package `gRain` which, while on CRAN,
depends on packages that are on Bioconductor both directly and through
the `gRbase` package, which depends on `RBGL`:

``` r
BiocManager::install()
BiocManager::install(c("graph", "Rgraphviz", "RBGL"))
install.packages("gRain")
```

## Overview

The prequential diagnostics examine the forecasts that flow from a model
in sequence. Each monitor given below indicates the probability of a
particular observation based on the previous observations and the model
structure. In the prequential mindset, we compute a probability of each
subsequent observation based on all previous data points. These
observations are then scored, and in this package we use the logarithmic
score function. The observations are then standardized to give a z-score
statistic. Following the recommendation of Cowell (2007), scores
indicate a poor fit where |z| \> 1.96

We demonstrate the efficacy of the prequential monitors with the Asia
data set from the bnlearn package. Details of the variables (nodes) can
be found in the documentation for bnlearn.

``` r
library(bnlearn)
data(asia)
summary(asia)
#>    A          S          T          L          B          E          X       
#>  no :4958   no :2485   no :4956   no :4670   no :2451   no :4630   no :4431  
#>  yes:  42   yes:2515   yes:  44   yes: 330   yes:2549   yes: 370   yes: 569  
#>    D       
#>  no :2650  
#>  yes:2350
```

Lauritzen defines two candidate models we replicate
below.

``` r
asia_bn <- bnlearn::model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
bnlearn::graphviz.plot(asia_bn) 
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="50%" />

``` r
asia_bn_alt <-  bnlearn::model2network("[A][S][T|A][L|S][B|S][E|T:L][X|E][D|B:E:S]")
bnlearn::graphviz.plot(asia_bn_alt) 
```

<img src="man/figures/README-unnamed-chunk-6-2.png" width="50%" />

### Global monitor

The global monitor is equivalent to the Bayes Factor and offers a good
assessment of the model as a whole. It is primarily useful for
differentiating between candidate models as the global monitor pinpoints
the nodes with different contributions for two candidate models. We can
assess the log likelihood contribution from each of the nodes and
determine which nodes contribute the most to the global monitor.
global.monitor.graph provides a quick visual for the entire network. The
darker the color, the more substantial the contribution of that node to
the log likelihood.

The global monitor assesses the probability of observing the data given
the model and as such does not particularly fit the prequential view.
However, it is a quick and useful first diagnostic assessment.

``` r
glob1 <- global_monitor(asia_bn, asia, alpha = 3)
glob2 <- global_monitor(asia_bn_alt, asia, alpha = 3)
glob1
#>   Vertex      Score
#> 1      A  249.74568
#> 2      B 3020.72273
#> 3      D 2144.16059
#> 4      E   26.31263
#> 5      L 1101.55381
#> 6      S 3469.43743
#> 7      T  260.38452
#> 8      X  851.22579
glob2
#>   Vertex      Score
#> 1      A  249.74568
#> 2      B 3020.72273
#> 3      D 2153.19484
#> 4      E   26.31263
#> 5      L 1101.55381
#> 6      S 3469.43743
#> 7      T  260.38452
#> 8      X  851.22579
```

In the alternative model, Dysnopea contributes more to the
log-likelihood.

``` r
plot(glob1)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />
