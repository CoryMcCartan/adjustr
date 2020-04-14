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

# Given a `prior` formula and a `vars_data` object containing (1) model draws,
# (2) model input data (if applicable), and (3) prior specification data (if
# applicable), compute the log probability of the prior distribution for each
# MCMC draw of the parameter of interest
calc_lp = function(prior, vars_data) {
    distr = tryCatch(call_fn(f_rhs(prior), distr_env),
                error = function(e) stop(paste("Distribution", as.character(prior)[3], "not supported.")))
    distr = distr(eval_tidy(f_lhs(prior), vars_data)) # plug in LHS values
    params = call_args(f_rhs(prior)) %>% map(eval_tidy, vars_data) # plug in RHS values
    exec(distr, !!!params) %>%
        apply(1:2, sum) # handle multivariate cases
}

# Given a `priors` list of prior formula, compile the necessary data from
# `object` MCMC draws and `data` model data
get_base_data = function(object, priors, parsed_vars, data) {
    map(priors, function(prior) {
        vars = c(as.character(f_lhs(prior)),
                 map_chr(call_args(f_rhs(prior)), as.character))
        # vars stored in MCMC draws
        vars_inmodel = intersect(vars, names(parsed_vars))
        vars_indraws = vars_inmodel[!stringr::str_ends(parsed_vars[vars_inmodel], "data")]
        vars_data = map(vars_indraws, ~ extract(object, ., permuted=F)) %>%
            set_names(vars_indraws) %>%
            append(data) # add rest of the data (`specs` + `data`)
    })
}

# Given a `priors` list of prior formula, compute the total log probability
# for each draw of the parameters of interest
calc_original_lp = function(object, priors, parsed_vars, data) {
    # figure out what data we need and calculate and sum lp
    map2(priors, get_base_data(object, priors, parsed_vars, data), calc_lp) %>%
        reduce(`+`)
}

# Given a `priors` list of prior formula, compute the total log probability for
# each draw of the parameters of interest, for each prior specification given
# in `specs`
calc_specs_lp = function(object, priors, parsed_vars, data, specs) {
    base_data = get_base_data(object, priors, parsed_vars, data)
    map(specs, function(spec) {
        # append this spec to base_data and calculate LP
        map2(priors, base_data, ~ calc_lp(.x, append(.y, spec))) %>%
            reduce(`+`)
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
}
# Turn mapping into an environment suitable for metaprogramming,
# and turn each density into its curried form (see `make_dens` above)
distr_env = new_environment(purrr::map(distrs, make_dens), parent=empty_env())
