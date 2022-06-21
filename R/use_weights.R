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
#' @examples \dontrun{
#' spec = make_spec(eta ~ student_t(df, 0, 1), df=1:10)
#' adjusted = adjust_weights(spec, eightschools_m)
#'
#' get_resampling_idxs(adjusted)
#' get_resampling_idxs(adjusted, frac=0.1, replace=FALSE)
#' }
#'
#' @export
get_resampling_idxs = function(x, frac=1, replace=TRUE) {
    if (frac < 0) stop("`frac` parameter must be nonnegative")
    get_idxs = function(w) {
        if (all(is.na(w))) return(NA_integer_)
        sample.int(length(w), size=round(frac*length(w)), replace=replace, prob=w)
    }

    if (inherits(x, "list")) {
        map(x, get_idxs)
    } else if (inherits(x, "adjustr_weighted")) {
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
#'
#'   Also supported is the custom function \code{wasserstein}, which computes
#'   the Wasserstein-p distance between the posterior distribution of the
#'   provided expression under the new model and under the original model, with
#'   \code{p=1} the default. Lower the \code{spacing} parameter from the
#'   default of 0.005 to compute a finer (but slower) approximation.
#'
#'   The arguments in \code{...} are automatically quoted and evaluated in the
#'   context of \code{.data}. They support unquoting and splicing.
#' @param .resampling Whether to compute summary statistics by first resampling
#'  the data according to the weights. Defaults to \code{FALSE}, but will be
#'  used for any summary statistic that is not \code{mean}, \code{var} or
#'  \code{sd}.
#' @param .model_data Stan model data, if not provided in the earlier call to
#'   \code{\link{adjust_weights}}.
#'
#' @return An \code{adjustr_weighted} object, with the new columns specified in
#' \code{...} added.
#'
#' @seealso \code{\link{adjust_weights}}, \code{\link{spec_plot}}
#'
#' @examples \dontrun{
#' model_data = list(
#'     J = 8,
#'     y = c(28, 8, -3, 7, -1, 1, 18, 12),
#'     sigma = c(15, 10, 16, 11, 9, 11, 10, 18)
#' )
#'
#' spec = make_spec(eta ~ student_t(df, 0, 1), df=1:10)
#' adjusted = adjust_weights(spec, eightschools_m)
#'
#' summarize(adjusted, mean(mu), var(mu))
#' summarize(adjusted, wasserstein(mu, p=2))
#' summarize(adjusted, diff_1 = mean(y[1] - theta[1]), .model_data=model_data)
#' summarize(adjusted, quantile(tau, probs=c(0.05, 0.5, 0.95)))
#' }
#'
#' @rdname summarize.adjustr_weighted
#' @export
summarise.adjustr_weighted = function(.data, ..., .resampling=FALSE, .model_data=NULL) {
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
        if (!is_call(call)) {
            stop("Expressions must summarize posterior draws; `",
                 expr_text(call), "` has a different value for each draw.\n",
                 "  Try `mean(", expr_text(call), ")` or `sd(", expr_text(call), ")`.")
        }
        if (!.resampling && exists(call_name(call), funs_env)) {
            fun = funs_env[[call_name(call)]]
        } else {
            fun = function(x, ...) apply(x, 2, call_fn(call), ...)
            .resampling = T
        }

        expr = expr_deparse(call_args(call)[[1]])
        expr = stringr::str_replace_all(expr, "\\[(\\d)", "[,\\1")
        expr = stringr::str_replace_all(expr, "(?<![a-zA-Z0-9._])mean\\(", "rowMeans(")
        expr = stringr::str_replace_all(expr, "(?<![a-zA-Z0-9._])sum\\(", "rowSums(")
        computed = as.array(eval_tidy(parse_expr(expr), data))
        if (length(dim(computed)) == 1) dim(computed) = c(dim(computed), 1)

        if (!.resampling) {
            new_col = map(.data$.weights, function(w) {
                exec(fun, computed, w, !!!map(call_args(call)[-1], eval_tidy))
            })
        } else {
            n_idx = max(min(5*iter, 20e3), iter)
            idxs = map(.data$.weights, ~ sample.int(iter, n_idx, replace=TRUE, prob=.))
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

# Weighted ECDF
weighted.ecdf = function(samp, wgt=rep(1, length(samp))) {
    or = order(samp)
    y = cumsum(wgt[or])/sum(wgt)
    f = stats::stepfun(samp[or], c(0, y))
    class(f) = c("weighted.ecdf", "ecdf", class(f))
    attr(f, "call") = sys.call()
    f
}
quantile.weighted.ecdf = function(f, q) {
    x = environment(f)$x
    y = environment(f)$y
    purrr::map_dbl(q, function(q) {
        if (q == 0) return(x[1])
        if (q == 1) return(utils::tail(x, 1))
        idx = which(y > q)[1]
        if (idx == 1) return(x[1])
        stats::approx(y[idx-0:1], x[idx-0:1], q)$y
    })
}

weighted.wasserstein = function(samp, wgt, p=1, spacing=0.005) {
    f = weighted.ecdf(samp, wgt)
    q = seq(0, 1, spacing)
    W = mean(abs(stats::quantile(samp, q, names=FALSE, type=4) - quantile.weighted.ecdf(f, q))^p)
    if (W < .Machine$double.eps) 0 else W^(1/p)
}

# Weighted summary functions that work on arrays
wtd_array_mean = function(arr, wgt) colSums(as.array(arr)*wgt) / sum(wgt)
wtd_array_var = function(arr, wgt)  wtd_array_mean((arr - wtd_array_mean(arr, wgt))^2, wgt)
wtd_array_sd = function(arr, wgt) sqrt(wtd_array_var(arr, wgt))
wtd_array_quantile = function(arr, wgt, probs=c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    apply(arr, 2, function(x) quantile.weighted.ecdf(weighted.ecdf(x, wgt), probs))
}
wtd_array_median = function(arr, wgt) wtd_array_quantile(arr, wgt, 0.5)
wtd_array_wasserstein = function(arr, wgt, ...) {
    apply(arr, 2, function(x) weighted.wasserstein(x, wgt, ...))
}

funs_env = new_environment(list(
    mean = wtd_array_mean,
    var = wtd_array_var,
    sd = wtd_array_sd,
    quantile = wtd_array_quantile,
    median = wtd_array_median,
    wasserstein = wtd_array_wasserstein
))


#' Plot Posterior Quantities of Interest Under Alternative Model Specifications
#'
#' Uses weights computed in \code{\link{adjust_weights}} to plot posterior
#' quantities of interest versus specification parameters
#'
#' @param x An \code{adjustr_weighted} object.
#' @param by The x-axis variable, which is usually one of the specification
#'   parameters. Can be set to \code{1} if there is only one specification.
#'   Automatically quoted and evaluated in the context of \code{x}.
#' @param post The posterior quantity of interest, to be computed for each
#'   resampled draw of each specification. Should evaluate to a single number
#'   for each draw. Automatically quoted and evaluated in the context of \code{x}.
#' @param only_mean Whether to only plot the posterior mean. May be more stable.
#' @param ci_level The inner credible interval to plot. Central
#'   100*ci_level% intervals are computed from the quantiles of the resampled
#'   posterior draws.
#' @param outer_level The outer credible interval to plot.
#' @param ... Ignored.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object which can be further
#'   customized with the \strong{ggplot2} package.
#'
#' @seealso \code{\link{adjust_weights}}, \code{\link{summarize.adjustr_weighted}}
#'
#' @examples \dontrun{
#' spec = make_spec(eta ~ student_t(df, 0, scale),
#'                  df=1:10, scale=seq(2, 1, -1/9))
#' adjusted = adjust_weights(spec, eightschools_m)
#'
#' spec_plot(adjusted, df, theta[1])
#' spec_plot(adjusted, df, mu, only_mean=TRUE)
#' spec_plot(adjusted, scale, tau)
#' }
#'
#' @export
spec_plot = function(x, by, post, only_mean=FALSE, ci_level=0.8,
                                 outer_level=0.95, ...) {
    if (!requireNamespace("ggplot2", quietly=TRUE)) { # nocov start
        stop("Package `ggplot2` must be installed to plot posterior quantities of interest.")
    } # nocov end
    if (ci_level > outer_level) stop("`ci_level` should be less than `outer_level`")

    post = enexpr(post)
    orig_row = filter(x, if_any(starts_with(".samp"), ~ . == "<original model>"))
    if (!only_mean) {
        outer = (1 - outer_level) / 2
        inner = (1 - ci_level) / 2
        q_probs = c(outer, inner, 0.5, 1-inner, 1-outer)
        sum_arg = quo(quantile(!!post, probs = !!q_probs))

        filter(x, if_any(starts_with(".samp"), ~ . != "<original model>")) %>%
            summarise.adjustr_weighted(.y = !!sum_arg) %>%
            rowwise() %>%
            mutate(.y_ol = .y[1],
                   .y_il = .y[2],
                   .y_med = .y[3],
                   .y_iu = .y[4],
                   .y_ou = .y[5]) %>%
        ggplot2::ggplot(ggplot2::aes({{ by }}, .y_med)) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin=.y_ol, ymax=.y_ou), alpha=0.4) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin=.y_il, ymax=.y_iu), alpha=0.5) +
            { if (nrow(orig_row) == 1)
                ggplot2::geom_hline(yintercept=summarise.adjustr_weighted(orig_row, .y = median(!!post))$`.y`,
                           lty="dashed")
            } +
            ggplot2::geom_line() +
            ggplot2::geom_point(size=3) +
            ggplot2::theme_minimal() +
            ggplot2::labs(y = expr_name(post))
    } else {
        filter(x, if_any(starts_with(".samp"), ~ . != "<original model>")) %>%
            summarise.adjustr_weighted(.y = mean(!!post)) %>%
        ggplot2::ggplot(ggplot2::aes({{ by }}, .y)) +
            { if (nrow(orig_row) == 1)
                ggplot2::geom_hline(yintercept=summarise.adjustr_weighted(orig_row, .y = mean(!!post))$`.y`,
                           lty="dashed")
            } +
            ggplot2::geom_line() +
            ggplot2::geom_point(size=3) +
            ggplot2::theme_minimal() +
            ggplot2::labs(y = expr_name(post))
    }
}