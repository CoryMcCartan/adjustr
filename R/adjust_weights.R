#' Compute Pareto-smoothed Importance Weights for Alternative Model
#' Specifications
#'
#' Given a set of new sampling statements, which can be parametrized by a
#' data frame or list, compute Pareto-smoothed importance weights and attach
#' them to the specification object, for further calculation and plotting.
#'
#' @param spec An object of class \code{adjustr_spec}, probably produced by
#'   \code{\link{make_spec}}, containing the new sampling sampling statements
#'   to replace their counterparts in the original Stan model, and the data,
#'   if any, by which these sampling statements are parametrized.
#' @param object A \code{\link[rstan]{stanfit}} model object.
#' @param data The data that was used to fit the model in \code{object}.
#'   Required only if one of the new sampling specifications involves Stan data
#'   variables.
#' @param keep_bad When \code{FALSE} (the strongly recommended default),
#'   alternate specifications which deviate too much from the original
#'   posterior, and which as a result cannot be reliably estimated using
#'   importance sampling (i.e., if the Pareto shape parameter is larger than
#'   0.7), have their weights discarded.
#' @param incl_orig When \code{TRUE}, include a row for the original
#'   model specification, with all weights equal. Can facilitate comaprison
#'   and plotting later.
#'
#' @return A tibble, produced by converting the provided \code{specs} to a
#'   tibble (see \code{\link{as.data.frame.adjustr_spec}}), and adding columns
#'   \code{.weights}, containing vectors of weights for each draw, and
#'   \code{.pareto_k}, containing the diagnostic Pareto shape parameters. Values
#'   greater than 0.7 indicate that importance sampling is not reliable.
#'   If \code{incl_orig} is \code{TRUE}, a row is added for the original model
#'   specification. Weights can be extracted with the
#'   \code{\link{pull.adjustr_weighted}} method. The returned object also
#'   includes the model sample draws, in the \code{draws} attribute.
#'
#' @examples \dontrun{
#' model_data = list(
#'     J = 8,
#'     y = c(28, 8, -3, 7, -1, 1, 18, 12),
#'     sigma = c(15, 10, 16, 11, 9, 11, 10, 18)
#' )
#'
#' spec = make_spec(eta ~ student_t(df, 0, 1), df=1:10)
#' adjust_weights(spec, eightschools_m)
#' adjust_weights(spec, eightschools_m, keep_bad=TRUE)
#'
#' spec = make_spec(y ~ student_t(df, theta, sigma), df=1:10)
#' adjust_weights(spec, eightschools_m, data=model_data)
#' # will throw an error because `y` and `sigma` aren't provided
#' adjust_weights(spec, eightschools_m)
#' }
#'
#' @export
adjust_weights = function(spec, object, data=NULL, keep_bad=FALSE, incl_orig=TRUE) {
    # CHECK ARGUMENTS
    object = get_fit_obj(object)
    model_code = object@stanmodel@model_code
    stopifnot(is.adjustr_spec(spec))

    parsed = parse_model(model_code)

    # if no model data provided, we can only resample distributions of parameters
    if (is.null(data)) {
        samp_vars = map_chr(parsed$samp, ~ as.character(f_lhs(.)))
        prior_vars = parsed$vars[samp_vars] != "data"
        parsed$samp = parsed$samp[prior_vars]
        data = list()
    }

    matched_samp = match_sampling_stmts(spec$samp, parsed$samp)
    original_lp = calc_original_lp(object, matched_samp, parsed$vars, data)
    specs_lp = calc_specs_lp(object, spec$samp, parsed$vars, data, spec$params)

    # compute weights
    wgts = map(specs_lp, function(spec_lp) {
        lratio = spec_lp - original_lp
        dim(lratio) = c(dim(lratio), 1)
        r_eff = loo::relative_eff(as.array(exp(-lratio)))
        psis_wgt = suppressWarnings(loo::psis(lratio, r_eff=r_eff))
        pareto_k = loo::pareto_k_values(psis_wgt)
        if (all(psis_wgt$log_weights == psis_wgt$log_weights[1])) {
            warning("New specification equal to old specification.", call.=FALSE)
            pareto_k = -Inf
        }

        list(
            weights = loo::weights.importance_sampling(psis_wgt, log=FALSE),
            pareto_k = pareto_k
        )
    })

    adjust_obj = as_tibble(spec)
    class(adjust_obj) = c("adjustr_weighted", class(adjust_obj))
    adjust_obj$.weights = map(wgts, ~ as.numeric(.$weights))
    adjust_obj$.pareto_k = purrr::map_dbl(wgts, ~ .$pareto_k)
    if (!keep_bad)
        adjust_obj$.weights[adjust_obj$.pareto_k > 0.7] = list(NA_real_)
    attr(adjust_obj, "draws") = rstan::extract(object)
    attr(adjust_obj, "data") = data
    attr(adjust_obj, "iter") = object@sim$chains * (object@sim$iter - object@sim$warmup)
    if (incl_orig) {
        adjust_obj = bind_rows(adjust_obj, tibble(
            .weights=list(rep(1, attr(adjust_obj, "iter"))),
            .pareto_k = -Inf))
        samp_cols = stringr::str_detect(names(adjust_obj), "\\.samp")
        adjust_obj[nrow(adjust_obj), samp_cols] = "<original model>"
    }

    adjust_obj
}


# Generic methods
is.adjustr_weighted = function(x) inherits(x, "adjustr_weighted")
#' Extract Weights From an \code{adjustr_weighted} Object
#'
#' This function modifies the default behavior of \code{dplyr::pull} to extract
#' the \code{.weights} column.
#'
#' @param .data A table of data
#' @param var A variable, as in \code{\link[dplyr]{pull}}. The default returns
#'   the \code{.weights} column, and if there is only one row, it returns the
#'   first element of that column
#'
#' @export
pull.adjustr_weighted = function(.data, var=".weights") {
    var = tidyselect::vars_pull(names(.data), !!enquo(var))
    if (nrow(.data) == 1 && var == ".weights") {
        .data$.weights[[1]]
    } else {
        .data[[var]]
    }
}

#' Extract Model Sampling Statements From a Stan Model.
#'
#' Prints a list of sampling statements extracted from the \code{model} block of
#' a Stan program, with each labelled "parameter" or "data" depending on the
#' type of variable being sampled.
#'
#' @param object A \code{\link[rstan]{stanfit}} model object.
#'
#' @return Invisbly returns a list of sampling formulas.
#'
#' @examples \dontrun{
#' extract_samp_stmts(eightschools_m)
#' #> Sampling statements for model 2c8d1d8a30137533422c438f23b83428:
#' #>   parameter   eta ~ std_normal()
#' #>   data        y ~ normal(theta, sigma)
#' }
#' @export
extract_samp_stmts = function(object) {
    model_code = get_fit_obj(object)@stanmodel@model_code

    parsed = parse_model(model_code)

    samp_vars = map_chr(parsed$samp, ~ as.character(f_lhs(.)))
    type = map_chr(samp_vars, function(var) {
        if (stringr::str_ends(parsed$vars[var], "data")) "data" else "parameter"
    })
    print_order = order(type, samp_vars, decreasing=c(T, F))

    cat(paste0("Sampling statements for model ", object@model_name, ":\n"))
    purrr::walk(print_order, ~ cat(sprintf("  %-9s   %s\n", type[.], as.character(parsed$samp[.]))))
    invisible(parsed$samp)
}

# Check that the model object is correct, and extract its Stan code
get_fit_obj = function(object) {
    if (is(object, "stanfit")) {
        object
    } else if (is(object, "stanreg")) {
        object$stanfit
    } else if (is(object, "brmsfit")) {
        object$fit
    } else {
        stop("`object` must be of class `stanfit`, `stanreg`, or `brmsfit`.")
    }
}
