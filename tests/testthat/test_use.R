context("Summarizing and using computed weights")

test_that("Resampling indices can be made from a vector", {
    expect_equal(get_resampling_idxs(c(0, 0, 1)), c(3, 3, 3))
    expect_equal(get_resampling_idxs(c(0, 1, 0), frac=1/3), 2)
    expect_true(is.na(get_resampling_idxs(NA)))
})

test_that("Resampling indices can be made from a list", {
    expect_equal(get_resampling_idxs(list(c(0, 0, 1), c(0, 1, 0))),
                                     list(c(3, 3, 3), c(2, 2, 2)))
    expect_equal(get_resampling_idxs(list(NA, c(0, 1, 0))),
                                     list(NA_integer_, c(2, 2, 2)))
})

test_that("Resampling indices can be made from an adjustr_weighted object", {
    spec = make_spec(eta ~ student_t(7, 0, 1))
    obj = adjust_weights(spec, eightschools_m, keep_bad=TRUE)
    obj = get_resampling_idxs(obj)
    expect_is(obj$.idxs, "list")
    expect_length(obj$.idxs[[1]], 20)
})

test_that("Resampling indices fail with negative `frac`", {
    expect_error(get_resampling_idxs(c(0, 0, 1), frac=-0.5), "must be nonnegative")
})

test_that("Weighted array functions compute correctly", {
    y = as.array(1:5)
    dim(y) = c(dim(y), 1)
    wgt = c(1, 1, 2, 5, 1)
    wtd_mean = weighted.mean(y, wgt)

    expect_equal(wtd_array_mean(y, wgt), wtd_mean)
    expect_equal(wtd_array_var(y, wgt), weighted.mean((y - wtd_mean)^2, wgt))
    expect_equal(wtd_array_sd(y, wgt), sqrt(weighted.mean((y - wtd_mean)^2, wgt)))
    expect_equal(wtd_array_quantile(y, rep(1, 5), 0.2), 1)
    expect_equal(wtd_array_median(y, rep(1, 5)), 2.5)
})

test_that("Empty call to `summarize` should change nothing", {
    obj = tibble(.weights=list(c(1,1,1), c(1,1,4)))
    class(obj) = c("adjustr_weighted", class(obj))
    expect_identical(summarize(obj), obj)
})

test_that("Non-summary call to `summarize` should throw error", {
    obj = tibble(.weights=list(c(1,1,1), c(1,1,4)))
    attr(obj, "draws") = list(theta=matrix(c(3,5,7,1,1,1), ncol=2))
    attr(obj, "iter") = 3
    class(obj) = c("adjustr_weighted", class(obj))

    expect_error(summarize(obj, theta), "must summarize posterior draws")
})

test_that("Basic summaries are computed correctly", {
    obj = tibble(.weights=list(c(1,1,1), c(1,1,4)))
    attr(obj, "draws") = list(theta=matrix(c(3,5,7,1,1,1), ncol=2))
    attr(obj, "iter") = 3
    class(obj) = c("adjustr_weighted", class(obj))

    sum1 = summarize(obj, mean(theta[1]))
    expect_is(sum1, "adjustr_weighted")
    expect_equal(sum1$`mean(theta[1])`, 5:6)

    sum2 = summarize(obj, th = mean(theta))
    expect_is(sum2, "adjustr_weighted")
    expect_equal(sum2$th, list(c(5, 1), c(6, 1)))

    sum3 = summarize(obj, W = wasserstein(theta[1]))
    expect_is(sum3, "adjustr_weighted")
    expect_equal(sum3$W[1], 0)

    expect_error(summarise.adjustr_weighted(as_tibble(obj)), "is not TRUE")
})

test_that("`summarize` uses data correctly", {
    obj = tibble(.weights=list(c(1,1,1), c(1,1,4)))
    attr(obj, "draws") = list(eta=matrix(c(3,5,7), nrow=3, ncol=1))
    attr(obj, "iter") = 3
    class(obj) = c("adjustr_weighted", class(obj))
    model_d = list(theta=4)

    expect_equal(summarize(obj, th=mean(theta+eta), .model_data=model_d)$th, c(9,10))

    attr(obj, "data") = model_d
    expect_equal(summarize(obj, th=mean(theta+eta))$th, c(9,10))
})


test_that("Resampling-based summaries are computed correctly", {
    obj = tibble(.weights=list(c(1,0,0), c(0,0,1)))
    attr(obj, "draws") = list(theta=matrix(c(3,5,7), nrow=3, ncol=1))
    attr(obj, "iter") = 3
    class(obj) = c("adjustr_weighted", class(obj))

    sum1 = summarize(obj, th=mean(theta), .resampling=TRUE)
    expect_equal(sum1$th, c(3,7))

    sum2 = summarize(obj, th=quantile(theta, 0.05), .resampling=TRUE)
    expect_equal(sum2$th, c(3,7))
})


test_that("Plotting function handles arguments correctly", {
    obj = tibble(.weights=list(c(1,0,0), c(0,0,1), c(1,1,1)),
                 .samp=c("y ~ normal(0, 1)", "y ~ normal(0, 2)", "<original model>"))
    attr(obj, "draws") = list(theta=matrix(c(3,5,7), nrow=3, ncol=1))
    attr(obj, "iter") = 3
    class(obj) = c("adjustr_weighted", class(obj))

    expect_is(spec_plot(obj, 1, theta), "ggplot")
    expect_is(spec_plot(obj, 1, theta, only_mean=TRUE), "ggplot")

    expect_error(spec_plot(obj, 1, theta, outer_level=0.4), "should be less than")
})
