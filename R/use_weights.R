#' Get Importance Resampling Indices From Weights
#'
#' Takes a vector of weights, or data frame or list containing sets of weights,
#' and resamples indices for use in later computation.
#'
#' @param x A vector of weights, a list of weight vectors, or a data frame of
#'   type \code{adjustr_weighted} containing a \code{.weights} list-column
#'   of weights.
#' @param frac A real number giving the fraction of draws to resample; the
#'   default, 1, resamples all draws. Smaller values should be used when
#'   \code{replace=FALSE}.
#' @param replace Whether sampling should be with replacement. When weights
#'   are extreme it may make sense to use \code{replace=FALSE}, but accuracy
#'   is not guaranteed in these cases.
#'
#' @return A vector, list, or data frame, depending of the type of \code{x},
#' containing the sampled indices. If any weights are \code{NA}, the indices
#' will also be \code{NA}.
#'
#' @export
get_resampling_idxs = function(x, frac=1, replace=T) {
    if (frac < 0) stop("`frac` parameter must be nonnegative")
    get_idxs = function(w) {
        if (all(is.na(w))) return(NA_integer_)
        sample.int(length(w), size=round(frac*length(w)), replace=replace, prob=w)
    }

    if (is(x, "list")) {
        map(x, get_idxs)
    } else if (is(x, "adjustr_weighted")) {
        x$.idxs = map(x$.weights, get_idxs)
        x
    } else {
        get_idxs(x)
    }
}

#' Summarize Posterior Distributions Under Alternative Model Specifications
#'
#' Uses weights computed in \code{\link{adjust_weights}} to compute posterior
#' summary statistics. These statistics can be compared against their reference
#' values to quantify the sensitivity of the model to aspects of its
#' specification.
#'
#' @param .data An \code{adjustr_weighted} object.
#' @param ... Name-value pairs of expressions. The name of each argument will be
#'   the name of a new variable, and the value will be computed for the
#'   posterior distribution of eight alternative specification. For example,
#'   a value of \code{mean(theta)} will compute the posterior mean of
#'   \code{theta} for each alternative specification.
#' @param .resampling Wether to compute summary statistics by first resampling
#'  the data according to the weights. Defaults to \code{FALSE}, but will be
#'  used for any summary statistic that is not \code{mean}, \code{var} or
#'  \code{sd}.
#' @param .model_data Stan model data, if not provided in the earlier call to
#'   \code{\link{adjust_weights}}.
#'
#' @return An \code{adjustr_weighted} object, wth the new columns specified in
#' \code{...} added.
#'
#' @rdname summarize.adjustr_weighted
#' @export
summarise.adjustr_weighted = function(.data, ..., .resampling=F, .model_data=NULL) {
    stopifnot(is.adjustr_weighted(.data)) # just in case called manually
    args = enexprs(...)

    broadcast = function(x) {
        dims = c(dim(as.array(x)), iter)
        x = array(rep(x, iter), dim=dims)
        aperm(x, c(length(dims), 2:length(dims) - 1))
    }
    iter = attr(.data, "iter")
    if (!is_null(.model_data)) attr(.data, "data") = .model_data
    data = append(attr(.data, "draws"), map(attr(.data, "data"), broadcast))

    n_args = length(args)
    for (i in seq_along(args)) {
        name = names(args)[i]
        if (name == "") name = expr_name(args[[i]])

        call = args[[i]]
        if (!.resampling && exists(call_name(call), funs_env)) {
            fun = call_fn(call, funs_env)
        } else {
            fun = function(x, ...) apply(x, 2, call_fn(call), ...)
            .resampling = T
        }

        expr = expr_deparse(call_args(call)[[1]])
        expr = stringr::str_replace_all(expr, "\\[(\\d)", "[,\\1")
        expr = stringr::str_replace_all(expr, "(?<![a-zA-Z0-9._])mean\\(", "rowMeans(")
        expr = stringr::str_replace_all(expr, "(?<![a-zA-Z0-9._])sum\\(", "rowSum(")
        computed = as.array(eval_tidy(parse_expr(expr), data))
        if (length(dim(computed)) == 1) dim(computed) = c(dim(computed), 1)

        if (!.resampling) {
            new_col = map(.data$.weights, ~ fun(computed, .))
        } else {
            idxs = map(.data$.weights, ~ sample.int(iter, iter, replace=T, prob=.))
            new_col = map(idxs, function(idx) {
                comp = as.array(computed[idx,])
                if (length(dim(comp)) == 1) dim(comp) = c(dim(comp), 1)
                exec(fun, comp, !!!map(call_args(call)[-1], eval_tidy))
            })
        }
        if (length(new_col[[1]]) == 1 && is.numeric(new_col[[1]]))
            new_col = as.numeric(new_col)
        .data[[name]] = new_col
    }

    .data
}
#' @rdname summarize.adjustr_weighted
#' @export
summarize.adjustr_weighted = summarise.adjustr_weighted

# Weighted summary functions that work on arrays
wtd_array_mean = function(arr, wgt) colSums(as.array(arr)*wgt) / sum(wgt)
wtd_array_var = function(arr, wgt)  wtd_array_mean((arr - wtd_array_mean(arr, wgt))^2, wgt)
wtd_array_sd = function(arr, wgt) sqrt(wtd_array_var(arr, wgt))

funs_env = new_environment(list(
    mean = wtd_array_mean,
    var = wtd_array_var,
    sd = wtd_array_sd
))
