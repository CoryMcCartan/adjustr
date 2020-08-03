context("Weight generation and model helpers")

setup({
    pkg_env$model_d = list(J = 8,
                           y = c(28,  8, -3,  7, -1,  1, 18, 12),
                           sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
})

test_that("Identical specification gives warning", {
    spec = make_spec(y ~ normal(theta, sigma))
    expect_warning(adjust_weights(spec, eightschools_m, pkg_env$model_d), "equal to old")
})

test_that("High Pareto k values lead to discarded weights", {
    spec = make_spec(y ~ normal(theta, 1.1*sigma))
    obj = adjust_weights(spec, eightschools_m, pkg_env$model_d)
    expect_true(is.na(obj$.weights[[1]]))
})

test_that("Weights calculated correctly (normal/inflated)", {
    theta_draws = rstan::extract(eightschools_m, "theta", permuted=F)
    y = pkg_env$model_d$y
    sigma = pkg_env$model_d$sigma

    ref_lp = apply(theta_draws, 1:2, function(theta) sum(dnorm(y, theta, sigma, log=T)))
    new_lp = apply(theta_draws, 1:2, function(theta) sum(dnorm(y, theta, 1.1*sigma, log=T)))
    lratio = new_lp - ref_lp
    dim(lratio) = c(dim(lratio), 1)
    r_eff = loo::relative_eff(as.array(exp(-lratio)))
    psis_wgt = suppressWarnings(loo::psis(lratio, r_eff=r_eff))
    pareto_k = loo::pareto_k_values(psis_wgt)
    weights = as.numeric(loo::weights.importance_sampling(psis_wgt, log=F))

    spec = make_spec(y ~ normal(theta, 1.1*sigma))
    obj = adjust_weights(spec, eightschools_m, pkg_env$model_d, keep_bad=T, incl_orig=F)

    expect_s3_class(obj, "adjustr_weighted")
    expect_s3_class(obj, "tbl_df")
    expect_true(is.adjustr_weighted(obj))
    expect_true("draws" %in% names(attributes(obj)))
    expect_true("data" %in% names(attributes(obj)))
    expect_equal(weights, obj$.weights[[1]])
    expect_equal(pareto_k, obj$.pareto_k)
})

test_that("Weights calculated correctly (normal/student_t)", {
    theta_draws = rstan::extract(eightschools_m, "theta", permuted=F)
    y = pkg_env$model_d$y
    sigma = pkg_env$model_d$sigma

    ref_lp = apply(theta_draws, 1:2, function(theta) sum(dnorm(y, theta, sigma, log=T)))
    new_lp = apply(theta_draws, 1:2, function(theta) sum(dt((y-theta)/sigma, 6, log=T)))
    lratio = new_lp - ref_lp
    dim(lratio) = c(dim(lratio), 1)
    r_eff = loo::relative_eff(exp(-lratio))
    psis_wgt = suppressWarnings(loo::psis(lratio, r_eff=r_eff))
    pareto_k = loo::pareto_k_values(psis_wgt)
    weights = as.numeric(loo::weights.importance_sampling(psis_wgt, log=F))

    spec = make_spec(y ~ student_t(df, theta, sigma), df=5:6)
    obj = adjust_weights(spec, eightschools_m, pkg_env$model_d, keep_bad=T, incl_orig=F)

    expect_equal(weights, obj$.weights[[2]])
    expect_equal(pareto_k, obj$.pareto_k[2])
})

test_that("Weights calculated correctly (no data normal/student_t)", {
    eta_draws = rstan::extract(eightschools_m, "eta", permuted=F)

    ref_lp = apply(eta_draws, 1:2, function(eta) sum(dnorm(eta, log=T)))
    new_lp = apply(eta_draws, 1:2, function(eta) sum(dt(eta, 4, log=T)))
    lratio = new_lp - ref_lp
    dim(lratio) = c(dim(lratio), 1)
    r_eff = loo::relative_eff(exp(-lratio))
    psis_wgt = suppressWarnings(loo::psis(lratio, r_eff=r_eff))
    pareto_k = loo::pareto_k_values(psis_wgt)
    weights = as.numeric(loo::weights.importance_sampling(psis_wgt, log=F))

    spec = make_spec(eta ~ student_t(4, 0, 1))
    obj = adjust_weights(spec, eightschools_m, keep_bad=T, incl_orig=F)

    expect_equal(weights, obj$.weights[[1]])
    expect_equal(pareto_k, obj$.pareto_k)
})


test_that("Weights extracted correctly", {
    spec = make_spec(y ~ student_t(df, theta, sigma), df=5)
    obj = adjust_weights(spec, eightschools_m, pkg_env$model_d, keep_bad=T, incl_orig=F)
    pulled = pull(obj)

    expect_is(pulled, "numeric")
    expect_length(pulled, 20)

    spec2 = make_spec(y ~ student_t(df, theta, sigma), df=5:6)
    obj = adjust_weights(spec2, eightschools_m, pkg_env$model_d, keep_bad=T, incl_orig=F)
    pulled = pull(obj)

    expect_is(pulled, "list")
    expect_length(pulled, 2)
    expect_equal(purrr::map_int(pulled, length), c(20, 20))
})

test_that("Sampling statements printed correctly", {
    expect_output(extract_samp_stmts(eightschools_m),
"Sampling statements for model 2c8d1d8a30137533422c438f23b83428:
  parameter   tau ~ uniform(-1e+100, 1e+100)
  parameter   mu ~ uniform(-1e+100, 1e+100)
  parameter   eta ~ std_normal()
  data        y ~ normal(theta, sigma)", fixed=T)
})

test_that("Fit objects extracted correctly", {
    obj = list(stanfit="stanreg", fit="brmsfit")

    class(obj) = "stanreg"
    expect_equal(get_fit_obj(obj), "stanreg")

    class(obj) = "brmsfit"
    expect_equal(get_fit_obj(obj), "brmsfit")

    class(obj) = "list"
    expect_error(get_fit_obj(obj), "must be of class")
})
