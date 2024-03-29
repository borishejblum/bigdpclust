
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `bigdpclust`

<!-- badges: start -->

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/bigdpclust)](https://cran.r-project.org/package=bigdpclust)
[![Travis-CI Build
Status](https://travis-ci.org/borishejblum/bigdpclust.svg?branch=master)](https://travis-ci.org/borishejblum/bigdpclust)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/borishejblum/bigdpclust?branch=master&svg=true)](https://ci.appveyor.com/project/borishejblum/bigdpclust)
[![Downloads](https://cranlogs.r-pkg.org/badges/bigdpclust?color=blue)](https://www.r-pkg.org/pkg/bigdpclust)
<!-- badges: end -->

`bigdpclust` performs clustering of tall data using a Bayesian
nonparametric Gaussian Dirichlet process mixture model.

## Installation

You can install the development version of `bigdpclust` from
[GitHub](https://github.com/bigdpclust) with:

``` r
#install.packages("devtools")
devtools::install_github("borishejblum/bigdpclust")
```

`bigdpclust` depends on the weightedobs branch from the `NPflow`
package, which can be installed through the following command:

``` r
devtools::install_github(repo = "borishejblum/NPflow", ref = "weightedobs")
```

``` r
library(ggplot2)
library(bigdpclust)

n1 <- 100000
n2 <- 100
mydata <- rbind(cbind(rnorm(n1), rnorm(n = n1)),
                cbind(rnorm(n2, m=10), rnorm(n = n2, m=10)))
plot(mydata)
```

<img src="man/figures/README-unnamed-chunk-3-1.jpeg" width="100%" />

``` r

res <- bigdpclust(mydata, nclumps=100, 
                  Nmcmc = 1000, plotevery = 2000, burnin = 500)
table(res$cluster[1:n1])
#> 
#>      2 
#> 100000
table(res$cluster[n1 + 1:n2])
#> 
#>   1 
#> 100
```

– Boris Hejblum & Paul Kirk
