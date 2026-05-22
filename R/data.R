#' Pre-fitted Eight Schools Model
#'
#' A small \code{\link[rstan]{stanfit}} object from the classic
#' "eight schools" hierarchical model (Rubin, 1981), included for use in
#' examples and tests.  The model was fit with 2 chains of 1000 iterations
#' (500 warmup) for a compact object size.
#'
#' @format A \code{\link[rstan]{stanfit}} object.
#' @examples
#' extract_samp_stmts(eightschools_m)
"eightschools_m"
