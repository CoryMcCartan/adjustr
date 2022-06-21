context("Log probability calculations")

test_that("`make_dens` curries correctly", {
    f = make_dens(dexp)
    expect_error(f(1, 2), "unused argument") # right number of arguments
    test_val = runif(1, 0, 10)
    test_x = runif(1, 0.5, 2)
    expect_equal(f(test_x)(test_val), dexp(test_x, test_val, log=TRUE))
})

test_that("Constant parameter log probabilities are calculated correctly", {
    draws = matrix(-3:3, ncol=1)
    lp = calc_lp(y ~ std_normal(), list(y=draws))
    expect_equal(lp, matrix(dnorm(-3:3, log=TRUE), ncol=1))
})

test_that("Model parameter log probabilities are calculated correctly", {
    draws = matrix(-3:3, ncol=1)
    mus = -4:2
    sigmas = c(2, 1, 1, 1, 2, 1, 2)
    lp = calc_lp(y ~ normal(mu, sigma), list(y=draws, mu=mus, sigma=sigmas))
    expect_equal(lp, matrix(dnorm(-3:3, mus, sigmas, log=TRUE), ncol=1))
})

test_that("Data is assembled correctly", {
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    bd = get_base_data(eightschools_m, list(eta ~ student_t(df, 0, tau)),
                       parsed$vars, list(df=1:2), "df")

    expect_length(bd, 1)
    expect_named(bd[[1]], c("eta", "tau"), ignore.order=TRUE)
    expect_equal(dim(bd[[1]]$eta), c(10, 2, 8))
    expect_equal(dim(bd[[1]]$tau), c(10, 2, 1))
})

test_that("MCMC draws are preferred over provided data", {
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    bd = get_base_data(eightschools_m, list(eta ~ student_t(2, 0, tau)),
                       parsed$vars, list(tau=3))

    expect_length(bd, 1)
    expect_named(bd[[1]], c("eta", "tau"), ignore.order=TRUE)
    expect_equal(dim(bd[[1]]$eta), c(10, 2, 8))
    expect_equal(dim(bd[[1]]$tau), c(10, 2, 1))
})

test_that("Parameter-less specification data is correctly assembled", {
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    bd = get_base_data(eightschools_m, list(y ~ std_normal()), parsed$vars,
                       list(y=c(28,  8, -3,  7, -1,  1, 18, 12), J=8))

    expect_length(bd, 1)
    expect_named(bd[[1]], "y")
    expect_equal(dim(bd[[1]]$y), c(10, 2, 8))
})


test_that("Error thrown for missing data", {
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)

    expect_error(get_base_data(eightschools_m, list(eta ~ normal(gamma, sigma)),
                               parsed$vars, list()), "sigma not found")
    expect_error(get_base_data(eightschools_m, list(eta ~ normal(gamma, 2)),
                               parsed$vars, list()), "gamma not found")
})

test_that("Model log probability is correctly calculated", {
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    form = eta ~ normal(0, 1)
    draws = rstan::extract(eightschools_m, "eta", permuted=FALSE)
    exp_lp =  2*apply(dnorm(draws, 0, 1, log=TRUE), 1:2, sum)
    lp = calc_original_lp(eightschools_m, list(form, form), parsed$vars, list())
    expect_equal(exp_lp, lp)
})

test_that("Alternate specifications log probabilities are correctly calculated", {
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    form = eta ~ normal(0, s)
    draws = rstan::extract(eightschools_m, "eta", permuted=FALSE)
    exp_lp =  2*apply(dnorm(draws, 0, 1, log=TRUE), 1:2, sum)
    lp = calc_specs_lp(eightschools_m, list(form, form), parsed$vars, list(), list(list(s=1)))
    expect_equal(exp_lp, lp[[1]])
})


