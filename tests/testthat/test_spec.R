context("Specification creation")

test_that("Specifications can be created out of any number of formulas with no data", {
    spec1 = make_spec(y ~ normal(mu, sigma))
    expect_equal(spec1$params, list(list()))
    expect_length(spec1$samp, 1)
    expect_is(spec1$samp[[1]], "formula")

    spec2 = make_spec(y ~ normal(mu, sigma), sigma ~ gamma(alpha, beta))
    expect_equal(spec2$params, list(list()))
    expect_length(spec2$samp, 2)
    expect_is(spec2$samp[[1]], "formula")
    expect_is(spec2$samp[[2]], "formula")
})

test_that("Empty specifications generate warnigns", {
    expect_warning(make_spec(), "No sampling statements provided.")
    expect_warning(make_spec(df=1:5), "No sampling statements provided.")
})

test_that("Specifications can be created with named vectors", {
    spec = make_spec(y ~ std_normal(), df=1:5)
    expect_equal(spec$params, purrr::transpose(list(df=1:5)))
    expect_error(make_spec(y ~ std_normal(), 1:5), "named vector")
})

test_that("Specifications can be created with data frames", {
    dat = tibble(df=1:5)
    spec = make_spec(y ~ std_normal(), dat)
    expect_equal(spec$params, purrr::transpose(as.list(dat)))
})

test_that("Specifications can be created with lists", {
    dat = tibble(df=1:5)
    ldat = purrr::transpose(as.list(dat))
    expect_equal(make_spec(y ~ std_normal(), ldat)$params, ldat)
    expect_equal(make_spec(y ~ std_normal(), as.list(dat))$params, ldat)

    expect_error(make_spec(y ~ std_normal(), list(list(a=1:3), 7)),
                 "List-of-list arguments must be coercible to data frames")
    expect_error(make_spec(y ~ std_normal(), list(list(a=1:3), list(c=2:3))),
                 "NAs found. Check input parameters and format.")
    expect_error(make_spec(y ~ std_normal(), list(a=1:3, b=2:3)),
                 "List-of-vector arguments must be coercible to data frames")
    expect_error(make_spec(y ~ std_normal(), list(y~x)), "must be lists of lists")
})

test_that("Specifications generics work correctly", {
    spec = make_spec(y ~ std_normal(), df=1:5)

    expect_output(print(spec), "Specification parameters:\n df\n  1")
    expect_true(is.adjustr_spec(spec))
    expect_equal(length(spec), 5)

    spec2 = make_spec(y ~ std_normal(), df=1:5, mu=0:4)
    expect_equal(length(filter(spec2, df < 4)), 3)
    expect_equal(slice(arrange(spec2, desc(df)), 1)$params[[1]]$df, 5)
    expect_equal(names(rename(spec2, b=df)$params[[1]]), c("b", "mu"))
    expect_equal(length(select(spec2, df)$params[[1]]), 1)
    expect_equal(filter(make_spec(y ~ normal())), make_spec(y ~ normal()))
    expect_is(as.data.frame(spec2), "data.frame")
    expect_is(as.data.frame(make_spec(y ~ normal(), x ~ normal())), "data.frame")
})

test_that("Specifications with multiple parameters have correct structure", {
    spec = make_spec(eta ~ student_t(df, 0, scale), df = 1:3, scale = c(1, 2, 3))
    expect_length(spec, 3)
    expect_equal(spec$params[[1]]$df, 1)
    expect_equal(spec$params[[1]]$scale, 1)
    expect_equal(spec$params[[3]]$df, 3)
    expect_equal(spec$params[[3]]$scale, 3)
})

test_that("Specifications with crossing parameters via data frame", {
    params = tidyr::crossing(df = 1:3, scale = c(1, 2))
    spec = make_spec(eta ~ student_t(df, 0, scale), params)
    expect_length(spec, 6)
    expect_length(spec$samp, 1)
})

test_that("Specifications with multiple sampling statements", {
    spec = make_spec(eta ~ student_t(df, 0, 1),
                     y ~ normal(theta, infl * sigma),
                     df = 1:3, infl = c(1, 1.5, 2))
    expect_length(spec$samp, 2)
    expect_length(spec, 3)
})

test_that("as.data.frame for single vs multiple sampling statements", {
    spec1 = make_spec(y ~ normal(0, s), s = 1:3)
    df1 = as.data.frame(spec1)
    expect_true(".samp" %in% names(df1))
    expect_false(".samp_1" %in% names(df1))

    spec2 = make_spec(y ~ normal(0, s), x ~ normal(0, 1), s = 1:3)
    df2 = as.data.frame(spec2)
    expect_true(".samp_1" %in% names(df2))
    expect_true(".samp_2" %in% names(df2))
    expect_false(".samp" %in% names(df2))
})

test_that("Specifications with no parameters convert to data frame", {
    spec = make_spec(y ~ normal(0, 1))
    df = as.data.frame(spec)
    expect_is(df, "data.frame")
    expect_equal(nrow(df), 1)
})

test_that("slice preserves adjustr_spec class", {
    spec = make_spec(eta ~ student_t(df, 0, 1), df = 1:10)
    sliced = slice(spec, 3:5)
    expect_true(is.adjustr_spec(sliced))
    expect_length(sliced, 3)
    expect_equal(sliced$params[[1]]$df, 3)
})

test_that("NAs in parameters throw error", {
    expect_error(make_spec(y ~ normal(0, s), s = c(1, NA, 3)),
                 "NAs found")
})