# adjustr

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/CoryMcCartan/adjustr/branch/master/graph/badge.svg)](https://codecov.io/gh/CoryMcCartan/adjustr?branch=master)
[![Travis build status](https://travis-ci.org/CoryMcCartan/adjustr.svg?branch=master)](https://travis-ci.org/CoryMcCartan/adjustr)
<!-- badges: end -->

**adjustr** is an R package which provides functions to help assess the
sensitivity of a Bayesian model (fitted with [Stan](https://mc-stan.org)) to the
specification of its likelihood and priors. Users provide a series of alternate
sampling specifications, and the package uses Pareto-smoothed importance
sampling to estimate posterior quantities of interest under each specification.
The package also provides functions to summarize and plot how these quantities
change across specifications.

The package aims to provide simple interface that makes it as easy as possible
for modellers to try out various adjustments to their Stan models, without
needing to write any specific Stan code or even recompile or rerun their model.

The package works by parsing Stan model code, so everything works best if the
model was written by the user. The complexity of **rstanarm** and **brms**
models makes it harder to make interpretable adjustments, but in principle 
these are possible too.

## Getting Started

The tutorial [vignettes](https://corymccartan.github.io/adjustr/articles/index.html) 
walk through a full sensitivity analysis for the classic 8-schools example.
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
