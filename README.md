
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bnmonitor

`bnmonitor` is a package for sensitivity analysis and robustness in
Bayesian networks (BNs).

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

`bnmonitor` provides a suite of function to investigate either a
data-learnt or an expert elicited BN. Its functions can be classified
into three main areas:

  - *Robustness in discrete BNs*: checking how well a BN represents a
    dataset;
  - *Sensitivity in discrete BNs*: assessing the effect of changes in
    the discrete BN’s probabilities;
  - *Sensitivity in continuous BNs*: assessing the effect of changes in
    the continuous BN’s probabilities, either in the standard or
    model-preserving
framework

<!-- The prequential diagnostics examine the forecasts that flow from a model in sequence. -->

<!-- Each monitor given below indicates the probability of a particular observation based on the previous observations and the model structure.  -->

<!-- In the prequential mindset, we compute a probability of each subsequent observation based on all previous data points.  -->

<!-- These observations are then scored, and in this package we use the logarithmic score function. -->

<!-- The observations are then standardized to give a z-score statistic.  -->

<!-- Following the recommendation of Cowell (2007), scores indicate a poor fit where |z| > 1.96  -->

<!-- We demonstrate the efficacy of the prequential monitors with the Asia data set from the bnlearn package. Details of the variables (nodes) can be found in the documentation for bnlearn. -->

## Robustness

Consider the `asia` dataset available from `bnlearn`. Details of the
variables can be found in the `bnlearn` documentation.

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

Consider two candidate BN models (Lauritzen and Spiegelhalter, 1988),
which only differ in the edge from `S` to
`D`.

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

<!-- The global monitor is equivalent to the Bayes Factor and offers a good assessment of the model as a whole.  -->

<!-- It is primarily useful for differentiating between candidate models as the global monitor pinpoints the nodes with different contributions for two candidate models. -->

<!-- We can assess the log likelihood contribution from each of the nodes and determine which nodes contribute the most to the global monitor. -->

<!-- global.monitor.graph provides a quick visual for the entire network.  -->

<!-- The darker the color, the more substantial the contribution of that node to the log likelihood. -->

<!-- The global monitor assesses the probability of observing the data given the model and as such does not particularly fit the prequential view.  -->

<!-- However, it is a quick and useful first diagnostic assessment. -->

A first useful diagnostic is the `global_monitor`, reporting the
contribution of each vertex to the log-likelihood of the model.

``` r
glob_asia <- global_monitor(asia_bn, asia, alpha = 3)
glob_asia_alt <- global_monitor(asia_bn_alt, asia, alpha = 3)
glob_asia
#>   Vertex      Score
#> 1      A  249.74568
#> 2      B 3020.72273
#> 3      D 2144.16059
#> 4      E   26.31263
#> 5      L 1101.55381
#> 6      S 3469.43743
#> 7      T  260.38452
#> 8      X  851.22579
glob_asia_alt
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

Global monitors can be plotted giving a quick view of the decomposition
of the log-likelihood. The darker the color, the more substantial the
contribution of a vertex.

``` r
plot(glob_asia)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Node monitor

There are two variants of node monitors.

  - The marginal node monitor computes the probability of the \(i\)th
    observation in the data set in turn after passing the evidence of
    the \(i-1\)th cases in the data set.

  - The conditional node monitor computes the probability of the \(i\)th
    observation in the data set after passing evidence of the \(i-1\)th
    cases in the data set, and the evidence for all nodes in the \(i\)th
    case except the node of interest.

As a quick survey of the nodes, the `node_monitor` command computes the
marginal and conditional monitors for the final observation in the data
set.

``` r
node_asia <- node_monitor(asia_bn, asia)
node_asia
#>   node marg.z.score cond.z.score
#> 1    A   -0.1029759   -0.1029862
#> 2    S -114.8350581 -118.7946397
#> 3    T   -0.1054182   -0.1054288
#> 4    L   -0.3013303   -0.3013649
#> 5    B  -34.6833789  -35.0407002
#> 6    E   -0.3226028   -0.3226405
#> 7    X   -0.4234203   -0.4234785
#> 8    D  -10.4048397  -10.3678671
```

The scores indicate a poor fit of the probability distributions
specified for the Smoking, Bronchitis, and Dysnopea nodes, since these
are larger than 1.96 in absolute value. Plots can also be created to
give a visual counterpart of the node
monitors.

``` r
plot(node_asia, which = "marginal")
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

<!-- The modeller may decide which of these distributions is of the most interest.  -->

<!-- For the purposes of this vignette, we assume that we would like to check the forecasts for Dysnopea.  -->

As an illustration we investigate further the fit of the variable
Dysnopea.

The sequential marginal monitor `seq_marg_monitor` gives us a closer
look at which particular forecasts in the data set might cause this poor
fit. We examine the sequential monitor for both candidate models.

``` r
seq_asia <- seq_marg_monitor(asia_bn, asia, "D")
seq_asia_alt <- seq_marg_monitor(asia_bn_alt, asia, "D")
seq_asia
#> Marginal Node Monitor for D 
#>  Minimum      -2.182895 
#>  Maximum      0.446613
seq_asia_alt
#> Marginal Node Monitor for D 
#>  Minimum      -2.408245 
#>  Maximum      0.1543135
```

``` r
library(gridExtra)
grid.arrange(plot(seq_asia),plot(seq_asia_alt),ncol=2)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="50%" />

Both sequential marginal node monitors fall within the recommendation of
|z| \> 1.96, indicating that both models are appropriate. However, for
later forecasts in the data set, Dysnopea tends to predict higher than
average observations. The forecasts flowing from the alternative model
are more accurate.

``` r
seq_asia <- seq_cond_monitor(asia_bn, asia, "D")
seq_asia_alt <- seq_cond_monitor(asia_bn_alt, asia, "D")
seq_asia
#> Conditional Node Monitor for D 
#>  Minimum      -1.658465 
#>  Maximum      1.862883
seq_asia_alt
#> Conditional Node Monitor for D 
#>  Minimum      -1.534081 
#>  Maximum      2.118901
```

``` r
grid.arrange(plot(seq_asia),plot(seq_asia_alt),ncol=2)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="50%" />
