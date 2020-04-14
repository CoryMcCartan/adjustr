#' adjustr: Stan Model Adjustments and Sensitivity Analyses using Importance
#' Sampling
#'
#'
#' Functions to help assess the sensitivity of a Bayesian model to the
#' specification of its likelihood and priors, estimated using the rstan
#' package. Users provide a series of alternate sampling specifications, and the
#' package uses Pareto-smoothed importance sampling to estimate posterior
#' quantities of interest under each specification.
#'
#' See the list of key functions and the example below.
#' Full package documentation available at \link{https://corymccartan.github.io/adjustr/}.
#'
#' @section Key Functions:
#' \itemize{
#'     \item \code{\link{adjustment_weights}}
#' }
#'
#' @import rlang
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr map
#' @import dplyr
#' @import ggplot2
#'
#' @docType package
#' @name adjustr
NULL
#> NULL