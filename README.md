
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
library(NPflow)
#rm(list=ls())

#Number of datapoints
n <- 500

#seed
#set.seed(2019)

# Sample data
m <- matrix(nrow=2, ncol=4, c(-1, 1, 1.5, 2, 2, -2, -1.5, -2))
p <- c(0.2, 0.1, 0.4, 0.3) # cluster frequencies

sdev <- array(dim=c(2,2,4))
sdev[, ,1] <- matrix(nrow=2, ncol=2, c(0.3, 0, 0, 0.3))
sdev[, ,2] <- matrix(nrow=2, ncol=2, c(0.1, 0, 0, 0.3))
sdev[, ,3] <- matrix(nrow=2, ncol=2, c(0.3, 0.15, 0.15, 0.3))
sdev[, ,4] <- .3*diag(2)
c <- rep(0,n)
z <- matrix(0, nrow=2, ncol=n)
#cat(k, "/", n, " observations simulated\n", sep="")
for(k in 1:n){
    c[k] = which(rmultinom(n=1, size=1, prob=p)!=0)
    z[,k] <- m[, c[k]] + sdev[, , c[k]]%*%matrix(rnorm(2, mean = 0, sd = 1), nrow=2, ncol=1)
    #cat(k, "/", n, " observations simulated\n", sep="")
}

d<-2
# Set parameters of G0
hyperG0 <- list()
hyperG0[["mu"]] <- rep(0,d)
hyperG0[["kappa"]] <- 0.001
hyperG0[["nu"]] <- d+2
hyperG0[["lambda"]] <- diag(d)/10

# hyperprior on the Scale parameter of DPM
a <- 0.0001
b <- 0.0001

# Number of iterations
N <- 30

# do some plots
doPlot <- TRUE
nbclust_init <- 30

z2 <- z
c2 <- c
z2 <- z2[, c!=3]
c2 <- c2[c!=3]
z2 <- cbind(z2, rowMeans(z[, c==3]))
my_obs_weights <- rep(1, length(c2))
c2 <- c(c2, 3)
table(c)
#> c
#>   1   2   3   4 
#> 104  47 209 140
table(c2)
#> c2
#>   1   2   3   4 
#> 104  47   1 140
my_obs_weights <- c(my_obs_weights, sum(c==3))


MCMCsample <- DPMGibbsN(z, hyperG0, a, b, N=1000, doPlot, nbclust_init, plotevery=1000,
                        gg.add=list(theme_bw(),
                                    guides(shape=guide_legend(override.aes = list(fill="grey45")))),
                        diagVar=FALSE)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

``` r
MCMCsample_w <- DPMGibbsN(z2, hyperG0, a, b, N=1000, doPlot, nbclust_init, plotevery=1000,
                          gg.add=list(theme_bw(),
                                      guides(shape=guide_legend(override.aes = list(fill="grey45")))),
                          diagVar=FALSE, obs_weights = my_obs_weights)
```

<img src="man/figures/README-unnamed-chunk-3-3.png" width="100%" /><img src="man/figures/README-unnamed-chunk-3-4.png" width="100%" />

``` r

plot_ConvDPM(MCMCsample_w, from=2)

s <- summary(MCMCsample_w, burnin = 200, thin=2, posterior_approx=FALSE,
             lossFn = "MBinderN")
s
#> summaryDPMMclust object with 400 observations sampled from the posterior:
#> ------------------------------------------------------------------------
#> Burnin: 200 MCMC iterations discarded
#> 
#> Point estimate of the partition with a gaussian mixture:
#>       1       5      20      25 
#> 0.0034%   0.16%   0.36%   0.48%
plot(s)
```

<img src="man/figures/README-unnamed-chunk-3-5.png" width="100%" /><img src="man/figures/README-unnamed-chunk-3-6.png" width="100%" />

– Boris Hejblum & Paul Kirk
