# adjustr

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.org/CoryMcCartan/adjustr.svg?branch=master)](https://travis-ci.org/CoryMcCartan/adjustr)
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
for modellers to try out various adjustments to their Stan models, without
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

The tutorial [vignette](https://corymccartan.github.io/adjustr/articles/eight-schools.html) 
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
