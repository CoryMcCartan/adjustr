# Pass this function a density function.
# It will return a curried version which first takes the points at
# which to evaluate the density, and then takes the distribution parameters
make_dens = function(f) {
    function(x) {
        function(...) {
            f(x, ..., log=T)
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
    params = call_args(f_rhs(samp)) %>% map(eval_tidy, vars_data)

    exec(distr, !!!params) %>%
        apply(1:2, sum) # handle multivariate cases
}

# Given a list of sampling formulas, compile the necessary data from
# `object` MCMC draws and `data` model data
# NOTE: fills data first with MCMC data, so these values will override any
# passed in data
get_base_data = function(object, samps, parsed_vars, data, extra_names=NULL) {
    map(samps, function(samp) {
        vars = c(as.character(f_lhs(samp)), map_chr(call_args(f_rhs(samp)), as.character))
        # vars stored in MCMC draws
        vars_inmodel = intersect(vars, names(parsed_vars))
        vars_indraws = vars_inmodel[!stringr::str_ends(parsed_vars[vars_inmodel], "data")]
        base_data = map(vars_indraws, ~ rstan::extract(object, ., permuted=F)) %>%
            set_names(vars_indraws) %>%
            append(data) # add rest of the `data` (+ `specs`)
        # check all data found
        found = vars %in% c(names(base_data), extra_names) |
            !is.na(suppressWarnings(as.numeric(vars)))
        if (!all(found)) stop(paste(vars[!found], collapse=", "), " not found")
        base_data
    })
}

# Given a `samps` list of samp formula, compute the total log probability
# for each draw of the parameters of interest
calc_original_lp = function(object, samps, parsed_vars, data) {
    # figure out what data we need and calculate and sum lp
    base_data = get_base_data(object, samps, parsed_vars, data)
    map2(samps, base_data, calc_lp) %>%
        purrr::reduce(`+`)
}

# Given a `samps` list of samp formula, compute the total log probability for
# each draw of the parameters of interest, for each samp specification given
# in `specs`
calc_specs_lp = function(object, samps, parsed_vars, data, specs) {
    base_data = get_base_data(object, samps, parsed_vars, data, names(specs[[1]]))
    map(specs, function(spec) {
        # append this spec to base_data and calculate LP
        map2(samps, base_data, ~ calc_lp(.x, append(.y, spec))) %>%
            purrr::reduce(`+`)
    })
}


# Mapping of Stan distribution names to R functions
# Grab even more distributions from `extraDistr` if avaialable
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
if (requireNamespace("extraDistr", quietly=T)) {
    distrs = append(distrs, list(
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
distr_env = new_environment(purrr::map(distrs, make_dens), parent=empty_env())
