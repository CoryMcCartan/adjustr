# Pass this function a density function.
# It will return a curried version which first takes the points at
# which to evaluate the density, and then takes the distribution parameters
make_dens = function(f) {
    function(x) {
        function(...) {
            f(x, ..., log=TRUE)
        }
    }
}

# Given a sampling formula and a `vars_data` object containing (1) model draws,
# (2) model input data (if applicable), and (3) prior specification data (if
# applicable), compute the log probability of the prior distribution for each
# MCMC draw of the parameter of interest
calc_lp = function(samp, vars_data) {
    # plug in RHS sampling distribution
    distr = tryCatch(call_fn(f_rhs(samp), distr_env),
                error=function(e)
                    stop("Distribution ", as.character(samp)[3], " not supported."))
    # plug in LHS, then RHS values. eval_tidy will throw error if no matches found
    distr = distr(eval_tidy(f_lhs(samp), vars_data))
    params = map(call_args(f_rhs(samp)), eval_tidy, vars_data)

    apply(exec(distr, !!!params), 1:2, sum)
}

# Given a list of sampling formulas, compile the necessary data from
# `object` MCMC draws and `data` model data
# NOTE: fills data first with MCMC data, so these values will override any
# passed in data
get_base_data = function(object, samps, parsed_vars, data, extra_names=NULL) {
    iter = object@sim$iter - object@sim$warmup
    chains = object@sim$chains
    reshape_data = function(x) {
        x = as.array(x)
        new_dim = c(iter, chains)
        new_x = array(rep(0, prod(new_dim)), dim=new_dim)
        apply(new_x, 1:2, function(y) x) %>%
            aperm(c(length(dim(x)) + 1:2, 1:length(dim(x))))
    }

    map(samps, function(samp) {
        vars = get_stmt_vars(samp)
        # vars stored in MCMC draws
        vars_inmodel = intersect(vars, names(parsed_vars))
        vars_indraws = vars_inmodel[!stringr::str_ends(parsed_vars[vars_inmodel], "data")]
        vars_indata = vars_inmodel[stringr::str_ends(parsed_vars[vars_inmodel], "data")]
        # check data vars provided
        found = vars_indata %in% names(data)
        if (!all(found)) stop(paste(vars_indata[!found], collapse=", "), " not found")
        # combine draws and data
        base_data = append(
            map(vars_indraws, ~ rstan::extract(object, ., permuted=FALSE)) %>%
                set_names(vars_indraws),
            map(vars_indata, ~ reshape_data(data[[.]])) %>%
                set_names(vars_indata),
        )
        # check all data found
        found = vars %in% c(names(base_data), extra_names)
        if (!all(found)) stop(paste(vars[!found], collapse=", "), " not found")
        base_data
    })
}

# Given a `samps` list of samp formula, compute the total log probability
# for each draw of the parameters of interest
calc_original_lp = function(object, samps, parsed_vars, data) {
    # figure out what data we need and calculate and sum lp
    base_data = get_base_data(object, samps, parsed_vars, data)
    purrr::reduce(map2(samps, base_data, calc_lp), `+`)
}

# Given a `samps` list of samp formula, compute the total log probability for
# each draw of the parameters of interest, for each samp specification given
# in `specs`
calc_specs_lp = function(object, samps, parsed_vars, data, specs) {
    base_data = get_base_data(object, samps, parsed_vars, data, names(specs[[1]]))
    map(specs, function(spec) {
        purrr::reduce(map2(samps, base_data, ~ calc_lp(.x, append(.y, spec))), `+`)
    })
}


# Mapping of Stan distribution names to R functions
 distrs = list(
    bernoulli = function(x, p, ...) dbinom(x, 1, p, ...),
    bernoulli_logit = function(x, p, ...) dbinom(x, 1, plogis(p), ...),
    binomial = dbinom,
    binomial_logit = function(x, n, p, ...) dbinom(x, n, plogis(p), ...),
    hypergeometric = dhyper,
    poisson = dpois,
    poisson_log = function(x, alpha, ...) dpois(x, exp(alpha), ...),
    normal = dnorm,
    std_normal = dnorm,
    student_t = function(x, df, loc, scale, ...) dt((x - loc)/scale, df, ...),
    cauchy = dcauchy,
    laplace = function(x, mu, sigma, ...) dexp(abs(x - mu), 1/sigma, ...),
    logistic = dlogis,
    lognormal = dlnorm,
    chi_square = dchisq,
    exponential = dexp,
    gamma = function(x, alpha, beta, ...) dgamma(x, shape=alpha, rate=beta, ...),
    weibull = dweibull,
    beta = dbeta,
    beta_proportion = function(x, mu, k, ...) dbeta(x, mu*k, (1-mu)*k, ...),
    uniform = dunif
)
# Grab even more distributions from `extraDistr` if available (called from .onLoad())
distrs_onload = function() {
    if (requireNamespace("extraDistr", quietly=TRUE)) {
        distrs <<- append(distrs, list(
            beta_binomial = extraDistr::dbbinom,
            categorical = extraDistr::dcat,
            gumbel = extraDistr::dgumbel,
            inv_chi_square = extraDistr::dinvchisq,
            scaled_inv_chi_square = extraDistr::dinvchisq,
            inv_gamma = extraDistr::dinvgamma,
            frechet = function(x, lambda, sigma, ...)
                extraDistr::dfrechet(x, lambda, sigma=sigma, ...),
            rayleigh = extraDistr::drayleigh,
            pareto = extraDistr::dpareto
        ))
    } else {
        message("`extraDistr` package not found. Install to access more ",
                "distributions, like inverse chi-square and beta-binomial.")
    }
    # Turn mapping into an environment suitable for metaprogramming,
    # and turn each density into its curried form (see `make_dens` above)
    distr_env <<- new_environment(purrr::map(distrs, make_dens))
}