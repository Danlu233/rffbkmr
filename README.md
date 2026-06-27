
# rffbkmr

<!-- badges: start -->
<!-- badges: end -->

The goal of rffbkmr is to use supervised random Fourier features to 
    approximate the Gaussian process to recude computation times.

## Installation

You can install the development version of rffbkmr like so:

``` r
# install.packages("devtools")
devtools::install_github("Danlu233/rffbkmr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(rffbkmr)

# load sample data
data(simudat)

# fit model with 10 basis functions and 1000 iterations
fit <- rffbkmr(y = simudat$y,
    X = simudat[,1:5],
    Z = simudat[,6:10],
    Z.eq = as.formula("~ z1 + as.factor(z2) + z3 + z4 + as.factor(z5)"),
    K = 1000,
    J = 10,
    verbose = FALSE)

# exposure-response function for marginal effect
est_all <- allHnew(fit = fit)

# univariate exposure-response function for x1
est_uni <- univarHnew(fit = fit, whichX = "x1")

# bivariate exposure-response function for x1 and x2
est_bi <- bivarHnew(fit = fit, whichX = c("x1", "x2"))
```

