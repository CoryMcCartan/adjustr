#' Create a model specification object
#'
#' Takes a set of new sampling statements, which can be parametrized by other
#' arguments, data frames, or lists, and creates an \code{adjustr_spec} object
#' suitable for use in \code{\link{adjustment_weights}}.
#'
#' @param ... Model specification. Each argument can either be a formula,
#'   a named vector, data frames,
#'   or lists.
#'
#'
#' data frame, or list of lists, containing parameters which will
#'   be substituted into the sampling statements contained in \code{spec_stamp}.
#'   If a data frame, columns will be substituted. If a list of lists, the named
#'   entries of each sublist will be substituted.
#' spec_samp A formula, or list of formulas, of new sampling statements
#'   to replace their counterparts in the original Stan model. Sampling
#'   distributions can be parametrized by constants, Stan parameters and data,
#'   or by variables from the \code{specs} object.
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
        params = purrr::transpose(as.list(spec_params))
    ),  class="adjustr_spec")
    spec_obj
}

# GENERIC FUNCTIONS for `adjustr_spec`
is.adjustr_spec = function(x) inherits(x, "adjustr_spec")
print.adjustr_spec = function(x, ...) {
    cat("Sampling specifications:\n")
    purrr::walk(x$samp, print)
    if (!is_empty(x$params)) {
        cat("\nSpecification parameters:\n")
        do.call(rbind, x$params) %>%
            as.data.frame %>%
            print(row.names=F, max=15*ncol(.))
    }
}
length.adjustr_spec = function(x) length(x$params)
#' Convert an \code{adjustr_spec} object into a data frame
#'
#' Returns the data frame of specification parameters, with added columns of
#' the form \code{.samp_1}, \code{.samp_2}, ... for each sampling statement
#' (or just \code{.samp} if there is only one sampling statement).
as.data.frame.adjustr_spec = function(x, ...) {
    if (length(x$params) == 0) {
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
# dplyr generics
dplyr_handler = function(dplyr_func, x, ...) {
    x$params = do.call(bind_rows, x$params) %>%
        as_tibble %>%
        dplyr_func(...) %>%
        as.list %>%
        purrr::transpose()
    x
}
#' \code{dplyr} methods for \code{adjustr_spec} objects
#'
#' Core \code{\link[dplyr]{dplyr}} verbs which don't involve grouping
#'  (\code{\link[dplyr]{filter}}, \code{\link[dplyr]{arrange}},
#'   \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}},
#'   \code{\link[dplyr]{rename}}, and \code{\link[dplyr]{slice}}) are
#' implemented and operate on the underlying table of specification parameters.
#'
#' @param .data the \code{adjustr_spec} object
#' @param ... additional arguments to underlying
#'
#' @rdname filter.adjust_spec
filter.adjustr_spec = function(.data, ..., .preserve=F) {
    dplyr_handler(dplyr::filter, .data, ..., .preserve=.preserve)
}
#' @rdname filter.adjust_spec
arrange.adjustr_spec = function(.data, ...) {
    dplyr_handler(dplyr::arrange, .data, ...)
}
#' @rdname filter.adjust_spec
rename.adjustr_spec = function(.data, ...) {
    dplyr_handler(dplyr::rename, .data, ...)
}
#' @rdname filter.adjust_spec
select.adjustr_spec = function(.data, ...) {
    dplyr_handler(dplyr::select, .data, ...)
}
#' @rdname filter.adjust_spec
slice.adjustr_spec = function(.data, ..., .preserve=F) {
    dplyr_handler(dplyr::slice, .data, ..., .preserve=.preserve)
}

