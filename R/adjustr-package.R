#' @keywords internal
"_PACKAGE"

#' @import rlang
#' @importFrom purrr map_chr map map2
#' @import dplyr
NULL

# internal; to store shared package objects
pkg_env = new_environment()

.onLoad = function(libname, pkgname) {  # nocov start
    # create the Stan parser
    #tryCatch(get_parser(), error = function(e) {})

    utils::globalVariables(c("name", "pos", "value", ".y", ".y_ol", ".y_ou",
                             ".y_il", ".y_iu", ".y_med", "quantile", "median"))

    # Grab even more distributions from `extraDistr` if available
    distrs_onload()
}  # nocov end
#> NULL