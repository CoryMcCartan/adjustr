#' Set Up Model Adjustment Specifications
#'
#' Takes a set of new sampling statements, which can be parametrized by other
#' arguments, data frames, or lists, and creates an \code{adjustr_spec} object
#' suitable for use in \code{\link{adjust_weights}}.
#'
#' @param ... Model specification. Each argument can either be a formula,
#'   a named vector, data frames, or lists.
#'
#'   Formula arguments provide new sampling statements to replace their
#'   counterparts in the original Stan model. All such formulas must be of the
#'   form \code{variable ~ distribution(parameters)}, where \code{variable} and
#'   \code{parameters} are Stan data variables or parameters, or are provided by
#'   other arguments to this function (see below), and where \code{distribution}
#'   matches one of the univariate
#'   \href{https://mc-stan.org/docs/2_22/functions-reference/conventions-for-probability-functions.html}{Stan distributions}.
#'   Arithmetic expressions of parameters are also allowed, but care must be
#'   taken with multivariate parameter arguments.  Since specifications are
#'   passed as formulas, R's arithmetic operators are used, not Stan's. As a
#'   result, matrix and elementwise multiplcation in Stan sampling statments may
#'   not be interpreted correctly. Moving these computations out of sampling
#'   statements and into a local variables will ensure correct results.
#'
#'   For named vector arguments, each entry of the vector will be substituted
#'   into the corresponding parameter in the sampling statements. For data
#'   frame, each entry in each column will be substituted into the corresponding
#'   parameter in the sampling statements.
#'
#'   List arguments are coerced to data frame. They can either be lists of named
#'   vectors, or lists of lists of single-element named vector.
#'
#'   The lengths of all parameter arguments must be consistent.  Named vectors
#'   can have length 1 or must have length equal to the number of rows in all
#'   data frame arguments and the length of list arguments.
#'
#' @return An object of class \code{adjustr_spec}, which is essentially a list
#' with two elements: \code{samp}, which is a list of sampling formulas, and
#' \code{params}, which is a list of lists of parameters. Core
#' \link[=filter.adjustr_spec]{dplyr verbs} which don't involve grouping
#' (\code{\link[dplyr]{filter}}, \code{\link[dplyr]{arrange}},
#' \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}},
#' \code{\link[dplyr]{rename}}, and \code{\link[dplyr]{slice}}) are
#' supported and operate on the underlying table of specification parameters.
#'
#' @examples
#' make_spec(eta ~ cauchy(0, 1))
#'
#' make_spec(eta ~ student_t(df, 0, 1), df=1:10)
#'
#' params = tidyr::crossing(df=1:10, infl=c(1, 1.5, 2))
#' make_spec(eta ~ student_t(df, 0, 1),
#'           y ~ normal(theta, infl*sigma),
#'           params)
#'
#' @export
make_spec = function(...) {
    args = dots_list(..., .check_assign=T)

    spec_samp = purrr::keep(args, is_formula)
    if (length(spec_samp) == 0) warning("No sampling statements provided.")
    spec_params = purrr::imap(args, function(value, name) {
        if (!is.numeric(name) && name != "") { # named arguments are preserved as is
            list2(!!name := value)
        } else if (is(value, "data.frame")) {
            as.list(value)
        } else if (is_list(value)) {
            if (is_list(value[[1]])) {
                tryCatch(as.list(do.call(bind_rows, value)), error = function(e) {
                    stop("List-of-list arguments must be coercible to data frames.")
                })
            } else if (is_vector(value[[1]])) {
                tryCatch(do.call(bind_cols, value), error = function(e) {
                    stop("List-of-vector arguments must be coercible to data frames.")
                })
                value
            } else {
                stop("List arguments must be lists of lists or lists of ",
                    "vectors, and coercible to data frames.")
            }
        } else if (is_formula(value)) {
            NULL
        } else {
            stop("Arguments must be formulas, named vectors, data frames, ",
                "or lists. Use ?make_spec to see documentations.")
        }
    }) %>%
        purrr::compact() %>% # remove NULLS
        purrr::flatten() %>%
        as_tibble

    if (any(is.na(spec_params)))
        stop("NAs found. Check input parameters and format.")

    spec_obj = structure(list(
        samp = spec_samp,
        params = if (nrow(spec_params) > 0)
                purrr::transpose(as.list(spec_params))
            else
                list(list())
    ),  class="adjustr_spec")
    spec_obj
}

# GENERIC FUNCTIONS for `adjustr_spec`
is.adjustr_spec = function(x) inherits(x, "adjustr_spec")
#' @export
print.adjustr_spec = function(x, ...) {
    cat("Sampling specifications:\n")
    purrr::walk(x$samp, print)
    if (length(x$params[[1]]) > 0) {
        cat("\nSpecification parameters:\n")
        df = as.data.frame(do.call(rbind, x$params))
        print(df, row.names=F, max=15*ncol(df))
    }
}
#' @export
length.adjustr_spec = function(x) length(x$params)

#' Convert an \code{adjustr_spec} Object Into a Data Frame
#'
#' Returns the data frame of specification parameters, with added columns of
#' the form \code{.samp_1}, \code{.samp_2}, ... for each sampling statement
#' (or just \code{.samp} if there is only one sampling statement).
#'
#' @param x the \code{adjustr_spec} object
#' @param ... additional arguments to underlying method
#'
#' @export
as.data.frame.adjustr_spec = function(x, ...) {
    if (length(x$params[[1]]) == 0) {
        params_df = as.data.frame(matrix(nrow=1, ncol=0))
    } else {
        params_df = do.call(bind_rows, x$params) %>%
            as_tibble %>%
            as.data.frame
    }

    n_samp = length(x$samp)
    if (n_samp == 1) {
        params_df$.samp = format(x$samp[[1]])
    } else {
        for (i in 1:n_samp) {
            colname = paste0(".samp_", i)
            params_df[colname] = format(x$samp[[i]])
        }
    }

    params_df
}


#' \code{dplyr} Methods for \code{adjustr_spec} Objects
#'
#' Core \code{\link[dplyr]{dplyr}} verbs which don't involve grouping
#'  (\code{\link[dplyr]{filter}}, \code{\link[dplyr]{arrange}},
#'   \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}},
#'   \code{\link[dplyr]{rename}}, and \code{\link[dplyr]{slice}}) are
#' implemented and operate on the underlying table of specification parameters.
#'
#' @param .data the \code{adjustr_spec} object
#' @param ... additional arguments to underlying method
#' @param .preserve as in \code{filter} and \code{slice}
#' @name dplyr.adjustr_spec
#'
#' @examples \dontrun{
#' spec = make_spec(eta ~ student_t(df, 0, 1), df=1:10)
#'
#' arrange(spec, desc(df))
#' slice(spec, 4:7)
#' filter(spec, df == 2)
#' }
NULL
# dplyr generics
dplyr_handler = function(dplyr_func, x, ...) {
    if (length(x$params[[1]]) == 0) return(x)
    x$params = do.call(bind_rows, x$params) %>%
        as_tibble %>%
        dplyr_func(...) %>%
        as.list %>%
        purrr::transpose()
    x
}

# no @export because R CMD CHECK didn't like it
#' @rdname dplyr.adjustr_spec
filter.adjustr_spec = function(.data, ..., .preserve=FALSE) {
    dplyr_handler(dplyr::filter, .data, ..., .preserve=.preserve)
}
#' @rdname dplyr.adjustr_spec
#' @export
arrange.adjustr_spec = function(.data, ...) {
    dplyr_handler(dplyr::arrange, .data, ...)
}
#' @rdname dplyr.adjustr_spec
#' @export
rename.adjustr_spec = function(.data, ...) {
    dplyr_handler(dplyr::rename, .data, ...)
}
#' @rdname dplyr.adjustr_spec
#' @export
select.adjustr_spec = function(.data, ...) {
    dplyr_handler(dplyr::select, .data, ...)
}
#' @rdname dplyr.adjustr_spec
#' @export
slice.adjustr_spec = function(.data, ..., .preserve=FALSE) {
    dplyr_handler(dplyr::slice, .data, ..., .preserve=.preserve)
}

