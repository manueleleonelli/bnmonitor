
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

## Gaussian Bayesian Networks

The functionalities of `bnmonitor` for sensitivity analysis in Gaussian
Bayesian networks are illustrated next using the `mathmarks` dataset
bundled within the package.

``` r
data(mathmarks)
head(mathmarks)
#>   mechanics vectors algebra analysis statistics
#> 1        77      82      67       67         81
#> 2        63      78      80       70         81
#> 3        75      73      71       66         81
#> 4        55      72      63       70         68
#> 5        63      63      65       70         63
#> 6        53      61      72       64         73
```

The data includes the grades (out of 100) of students in five maths
exams: mechanics, vectors, algebra, analysis and statistics.

The structure of a Bayesian network for this data is first learnt using
the package `bnlearn` and the maximum likelihood estimate of its
parameters is computed and stored in `bnfit`.

``` r
library(bnlearn)
bn <- hc(mathmarks)
plot(bn)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="50%" />

``` r
bnfit <-bn.fit(bn,mathmarks)
```

To start the sensitivity analysis for the parameters of the learnt
Bayesian network, one first need to transform `bnfit` to objects of
class `GBN` (for standard sensitivity analysis) and `CI` (for
model-preserving sensitivity). This can be done using the functions
`bn2gbn` and `bn2ci` respectively.

``` r
gbn <- bn2gbn(bnfit)
ci <-  bn2ci(bnfit)
c(class(gbn),class(ci))
#> [1] "GBN" "CI"
```

#### Perturbation of the mean vector

A varied GBN after a perturbation of an entry of the mean vector can be
obtained with the function `mean_var`, which can only be applied to an
object of class `GBN`. Below, we vary the fifth entry of the mean vector
(statistics), \(\mu_5\) by an additive factor \(\delta = 10\).

``` r
rbind( t(gbn$mean),t(mean_var(gbn,entry = 5, delta = 10)$mean))
#>          [,1]     [,2]     [,3]     [,4]     [,5]
#> [1,] 38.95455 50.59091 50.60227 46.68182 42.30682
#> [2,] 38.95455 50.59091 50.60227 46.68182 52.30682
```

The overall effect of such variations can be assessed in terms of
dissimilarity measures: the Kullback-Leibler divergence (`KL`) and
Jeffrey’s divergence (`Jeffreys`). For instance, let’s see what’s the
effect of variations in the mean of the statistics exam.

``` r
mean_var5 <- KL(gbn, "mean", entry=5, delta = seq(-10,10,0.1))
mean_var5$plot
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="35%" />

More interestingly, one can check the different effect of variations of
different paramenters.
<img src="man/figures/README-unnamed-chunk-10-1.png" width="45%" /><img src="man/figures/README-unnamed-chunk-10-2.png" width="45%" />

Therefore, misspecifications of the mean of the algebra exam would have
the biggest effect on the distribution of the Gaussian Bayesian network.

#### Perturbation of the covariance matrix

Care must be taken when performing perturbations of the covariance
matrix, for two reasons: (1) the perturbed matrix may not be positive
semidefinite; (2) the perturbed matrix may not respect the conditional
indepedences of the underlying Bayesian network.
