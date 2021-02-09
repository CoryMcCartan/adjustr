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
#' Full package documentation available at \url{https://corymccartan.github.io/adjustr/}.
#'
#' @section Key Functions:
#' \itemize{
#'     \item \code{\link{make_spec}}
#'     \item \code{\link{adjust_weights}}
#'     \item \code{\link{summarize}}
#'     \item \code{\link{plot}}
#' }
#'
#' @importFrom methods is
#' @import rlang
#' @importFrom purrr map_chr map map2
#' @import dplyr
#' @import ggplot2
#'
#' @docType package
#' @name adjustr
NULL

# to store shared package objects
pkg_env = new_environment()

.onLoad = function(libname, pkgname) {  # nocov start
    # create the Stan parser
    #tryCatch(get_parser(), error = function(e) {})

    utils::globalVariables(c("name", "pos", "value", ".y", ".y_ol", ".y_ou",
                             ".y_il", ".y_iu", ".y_med", "quantile", "median"))
}  # nocov end
#> NULL