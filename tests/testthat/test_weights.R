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
    theta_draws = rstan::extract(eightschools_m, "theta", permuted=FALSE)
    y = pkg_env$model_d$y
    sigma = pkg_env$model_d$sigma

    ref_lp = apply(theta_draws, 1:2, function(theta) sum(dnorm(y, theta, sigma, log=TRUE)))
    new_lp = apply(theta_draws, 1:2, function(theta) sum(dnorm(y, theta, 1.1*sigma, log=TRUE)))
    lratio = new_lp - ref_lp
    dim(lratio) = c(dim(lratio), 1)
    r_eff = loo::relative_eff(as.array(exp(-lratio)))
    psis_wgt = suppressWarnings(loo::psis(lratio, r_eff=r_eff))
    pareto_k = loo::pareto_k_values(psis_wgt)
    weights = as.numeric(loo::weights.importance_sampling(psis_wgt, log=FALSE))

    spec = make_spec(y ~ normal(theta, 1.1*sigma))
    obj = adjust_weights(spec, eightschools_m, pkg_env$model_d, keep_bad=TRUE, incl_orig=FALSE)

    expect_s3_class(obj, "adjustr_weighted")
    expect_s3_class(obj, "tbl_df")
    expect_true(is.adjustr_weighted(obj))
    expect_true("draws" %in% names(attributes(obj)))
    expect_true("data" %in% names(attributes(obj)))
    expect_equal(weights, obj$.weights[[1]])
    expect_equal(pareto_k, obj$.pareto_k)
})

test_that("Weights calculated correctly (normal/student_t)", {
    theta_draws = rstan::extract(eightschools_m, "theta", permuted=FALSE)
    y = pkg_env$model_d$y
    sigma = pkg_env$model_d$sigma

    ref_lp = apply(theta_draws, 1:2, function(theta) sum(dnorm(y, theta, sigma, log=TRUE)))
    new_lp = apply(theta_draws, 1:2, function(theta) sum(dt((y-theta)/sigma, 6, log=TRUE)))
    lratio = new_lp - ref_lp
    dim(lratio) = c(dim(lratio), 1)
    r_eff = loo::relative_eff(exp(-lratio))
    psis_wgt = suppressWarnings(loo::psis(lratio, r_eff=r_eff))
    pareto_k = loo::pareto_k_values(psis_wgt)
    weights = as.numeric(loo::weights.importance_sampling(psis_wgt, log=FALSE))

    spec = make_spec(y ~ student_t(df, theta, sigma), df=5:6)
    obj = adjust_weights(spec, eightschools_m, pkg_env$model_d, keep_bad=TRUE, incl_orig=FALSE)

    expect_equal(weights, obj$.weights[[2]])
    expect_equal(pareto_k, obj$.pareto_k[2])
})

test_that("Weights calculated correctly (no data normal/student_t)", {
    eta_draws = rstan::extract(eightschools_m, "eta", permuted=FALSE)

    ref_lp = apply(eta_draws, 1:2, function(eta) sum(dnorm(eta, log=TRUE)))
    new_lp = apply(eta_draws, 1:2, function(eta) sum(dt(eta, 4, log=TRUE)))
    lratio = new_lp - ref_lp
    dim(lratio) = c(dim(lratio), 1)
    r_eff = loo::relative_eff(exp(-lratio))
    psis_wgt = suppressWarnings(loo::psis(lratio, r_eff=r_eff))
    pareto_k = loo::pareto_k_values(psis_wgt)
    weights = as.numeric(loo::weights.importance_sampling(psis_wgt, log=FALSE))

    spec = make_spec(eta ~ student_t(4, 0, 1))
    obj = adjust_weights(spec, eightschools_m, keep_bad=TRUE, incl_orig=FALSE)

    expect_equal(weights, obj$.weights[[1]])
    expect_equal(pareto_k, obj$.pareto_k)
})


test_that("Weights extracted correctly", {
    spec = make_spec(y ~ student_t(df, theta, sigma), df=5)
    obj = adjust_weights(spec, eightschools_m, pkg_env$model_d, keep_bad=TRUE, incl_orig=FALSE)
    pulled = pull(obj)

    expect_is(pulled, "numeric")
    expect_length(pulled, 20)

    spec2 = make_spec(y ~ student_t(df, theta, sigma), df=5:6)
    obj = adjust_weights(spec2, eightschools_m, pkg_env$model_d, keep_bad=TRUE, incl_orig=FALSE)
    pulled = pull(obj)

    expect_is(pulled, "list")
    expect_length(pulled, 2)
    expect_equal(vapply(pulled, length, 0L), c(20L, 20L))
})

test_that("Sampling statements printed correctly", {
    expect_output(extract_samp_stmts(eightschools_m),
"Sampling statements for model 2c8d1d8a30137533422c438f23b83428:
  parameter   eta ~ std_normal()
  parameter   mu ~ uniform(-1e+100, 1e+100)
  parameter   tau ~ uniform(-1e+100, 1e+100)
  data        y ~ normal(theta, sigma)", fixed=TRUE)
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

test_that("adjust_weights includes original model row by default", {
    spec = make_spec(eta ~ student_t(4, 0, 1))
    obj = adjust_weights(spec, eightschools_m, keep_bad = TRUE, incl_orig = TRUE)
    expect_true(any(obj$.samp == "<original model>"))
    expect_equal(obj$.pareto_k[nrow(obj)], -Inf)
    # Original model weights should all be 1
    orig_wgts = obj$.weights[[nrow(obj)]]
    expect_true(all(orig_wgts == 1))
})

test_that("adjust_weights without original model row", {
    spec = make_spec(eta ~ student_t(4, 0, 1))
    obj = adjust_weights(spec, eightschools_m, keep_bad = TRUE, incl_orig = FALSE)
    expect_false(any(obj$.samp == "<original model>", na.rm = TRUE))
    expect_equal(nrow(obj), 1)
})

test_that("adjust_weights with multiple specifications", {
    spec = make_spec(eta ~ student_t(df, 0, 1), df = c(3, 5, 7))
    obj = adjust_weights(spec, eightschools_m, keep_bad = TRUE, incl_orig = FALSE)
    expect_equal(nrow(obj), 3)
    expect_length(obj$.weights, 3)
    expect_length(obj$.pareto_k, 3)
    # Each weight vector should have the right length
    expect_true(all(vapply(obj$.weights, length, 0L) == attr(obj, "iter")))
})

test_that("adjust_weights preserves spec parameters in output", {
    spec = make_spec(eta ~ student_t(df, 0, 1), df = c(3, 5, 7))
    obj = adjust_weights(spec, eightschools_m, keep_bad = TRUE, incl_orig = FALSE)
    expect_true("df" %in% names(obj))
    expect_equal(obj$df, c(3, 5, 7))
})

test_that("adjust_weights errors when data needed but not provided", {
    spec = make_spec(y ~ normal(theta, sigma))
    expect_error(adjust_weights(spec, eightschools_m))
})
