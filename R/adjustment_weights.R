#' Compute Pareto-smoothed importance weights for alternative model
#' specifications
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
#'
#' @return If no \code{specs} were provided, a data frame with columns
#'   \code{.samp_stmt}, containing the sampling statements, and \code{.weights},
#'   containing weight objects of class \code{\link[loo]{psis}}. If \code{specs}
#'   is a data frame or list, then a \code{.weights} entry is added for each row
#'   or sublist.
#'
#'
#' @export
adjustment_weights = function(spec, object, data=NULL) {
    # CHECK ARGUMENTS
    model_code = get_model_code(object)
    #spec_samp =

    parsed_model = parse_model(model_code)
    parsed_vars = get_variables(parsed_model)
    parsed_samp = get_sampling_stmts(parsed_model)

    # if no model data provided, we can only resample distributions of parameters
    if (is.null(data)) {
        samp_vars = map_chr(parsed_samp, ~ as.character(f_lhs(.)))
        prior_vars = parsed_vars[samp_vars] != "data"
        parsed_samp = parsed_samp[prior_vars]
        data = list()
    }

    matched_samp = match_sampling_stmts(spec_samp, parsed_samp)
    original_lp = calc_original_lp(object, matched_samp, parsed_vars, data)
    specs_lp = calc_specs_lp(object, spec_samp, parsed_vars, data, specs)
    # compute weights
    map(specs_lp, function(spec_lp) {
        lratio = spec_lp - original_lp
        dim(lratio) = c(dim(lratio), 1)
        r_eff = loo::relative_eff(as.array(exp(-lratio)))
        psis_wgt = suppressWarnings(loo::psis(lratio, r_eff=r_eff))
        print(loo::pareto_k_values(psis_wgt))
        loo::weights.importance_sampling(psis_wgt, log=F)
    })
}

#' Extract model sampling statements from a Stan model.
#'
#' Prints a list of sampling statements extracted from the \code{model} block of
#' a Stan program, with each labelled "parameter" or "data" depending on the
#' type of variable being sampled.
#'
#' @param object A \code{\link[rstan]{stanfit}} model object.
#'
#' @return Invisbly returns a list of sampling formulas.
#'
#' @export
extract_model = function(object) {
    model_code = get_model_code(object)

    parsed_model = parse_model(model_code)
    parsed_vars = get_variables(parsed_model)
    parsed_samp = get_sampling_stmts(parsed_model)

    samp_vars = map_chr(parsed_samp, ~ as.character(f_lhs(.)))
    type = map_chr(samp_vars, function(var) {
        if (stringr::str_ends(parsed_vars[var], "data")) "data" else "parameter"
    })
    print_order = order(type, samp_vars, decreasing=c(T, F))

    cat(paste0("Sampling statements for model ", object@model_name, ":\n"))
    purrr::walk(print_order, ~ cat(sprintf("  %-9s   %s\n", type[.], as.character(parsed_samp[.]))))
    invisible(parsed_samp)
}

# Check that the model object is correct, and extract its Stan code
get_model_code = function(object) {
    if (is(object, "stanfit")) {
        return(object@stanmodel@model_code)
    } else {
        stop("`object` must be a stanfit.")
    }
}

# Change prior specification to list of lists with prior values
prep_specs = function(specs) {
    if (is.null(specs)) {
        return(list(list()))
    } else {
        if (is(specs, "list")) {
            purrr::walk(specs, function(el) {
                if (!is(el, "list"))
                    stop("`priors` must be a formula or a list of formulas")
            })
            return(specs)
        } else if (is(specs, "data.frame")) {
            return(as.list(specs) %>% purrr::transpose())
        } else {
            stop("`specs` must be NULL, a data frame, or a list of lists")
        }
    }
}