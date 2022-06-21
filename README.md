# adjustr <a href="https://corymccartan.github.io/adjustr/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/CoryMcCartan/adjustr/workflows/R-CMD-check/badge.svg)](https://github.com/CoryMcCartan/adjustr/actions)
[![Codecov test coverage](https://codecov.io/gh/CoryMcCartan/adjustr/branch/master/graph/badge.svg)](https://codecov.io/gh/CoryMcCartan/adjustr?branch=master)
<!-- badges: end -->

Sensitivity analysis is a critical component of a good modeling workflow. Yet
as the number and power of Bayesian computational tools has increased, the
options for sensitivity analysis have remained largely the same: compute
importance sampling weights manually, or fit a large number of similar models,
dramatically increasing computation time. Neither option is satisfactory for
most applied modeling.

**adjustr** is an R package which aims to make sensitivity analysis faster
and easier, and works with Bayesian models fitted with [Stan](https://mc-stan.org). 
Users provide a series of alternate sampling specifications, and the package
uses Pareto-smoothed importance sampling to estimate the posterior under each
specification. The package also provides functions to summarize and plot how
posterior quantities quantities change across specifications.

The package provides simple interface that makes it as easy as possible
for researchers to try out various adjustments to their Stan models, without
needing to write any specific Stan code or even recompile or rerun their model.

The package works by parsing Stan model code, so everything works best if the
model was written by the user. Models made using **brms** may in principle be
used as well. Models made using **rstanarm** are constructed using more 
complex model templates, and cannot be used.

## Getting Started

The basic __adjustr__ workflow is as follows:

1. Use [`make_spec`](https://corymccartan.github.io/adjustr/reference/make_spec.html) 
to specify the set of alternative model specifications you'd like to fit.

2. Use [`adjust_weights`](https://corymccartan.github.io/adjustr/reference/adjust_weights.html) 
to calculate importance sampling weights which approximate the posterior of each
alternative specification.

3. Use [`summarize`](https://corymccartan.github.io/adjustr/reference/summarize.adjustr_weighted.html) 
and [`spec_plot`](https://corymccartan.github.io/adjustr/reference/spec_plot.html) 
to examine posterior quantities of interest for each alternative specification,
in order to assess the sensitivity of the underlying model.

To illustrate, the package lets us do the following:
```r
extract_samp_stmts(eightschools_m)
#> Sampling statements for model 2c8d1d8a30137533422c438f23b83428:
#>   parameter   eta ~ std_normal()
#>   data        y ~ normal(theta, sigma)

make_spec(eta ~ student_t(0, 1, df), df=1:10) %>%
    adjust_weights(eightschools_m) %>%
    summarize(wasserstein(mu)) 
#> # A tibble: 11 x 5
#>       df .samp                     .weights      .pareto_k `wasserstein(mu)`
#>    <int> <chr>                     <list>            <dbl>             <dbl>
#>  1     1 eta ~ student_t(df, 0, 1) <dbl [1,000]>     1.02              0.928
#>  2     2 eta ~ student_t(df, 0, 1) <dbl [1,000]>     1.03              0.736
#>  3     3 eta ~ student_t(df, 0, 1) <dbl [1,000]>     0.915             0.534
#>  4     4 eta ~ student_t(df, 0, 1) <dbl [1,000]>     0.856             0.411
#>  5     5 eta ~ student_t(df, 0, 1) <dbl [1,000]>     0.826             0.341
#>  6     6 eta ~ student_t(df, 0, 1) <dbl [1,000]>     0.803             0.275
#>  7     7 eta ~ student_t(df, 0, 1) <dbl [1,000]>     0.782             0.234
#>  8     8 eta ~ student_t(df, 0, 1) <dbl [1,000]>     0.753             0.195
#>  9     9 eta ~ student_t(df, 0, 1) <dbl [1,000]>     0.736             0.166
#> 10    10 eta ~ student_t(df, 0, 1) <dbl [1,000]>     0.721             0.151
#> 11    NA <original model>          <dbl [1,000]>  -Inf                 0
```

The tutorial [vignette](https://corymccartan.github.io/adjustr/articles/eight-schools.html) 
walks through a full sensitivity analysis for this 8-schools example.
Smaller examples are also included in the package 
[documentation](https://corymccartan.github.io/adjustr/reference/index.html). 

## Installation

Install the latest version from **GitHub**:

```r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("corymccartan/adjustr@*release")
```

## References

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2015). 
Pareto smoothed importance sampling. 
_[arXiv preprint arXiv:1507.02646](https://arxiv.org/abs/1507.02646)_.
