context("Specification creation")

test_that("Specifications can be created out of any number of formulas with no data", {
    spec1 = make_spec(y ~ normal(mu, sigma))
    expect_equal(spec1$params, list())
    expect_length(spec1$samp, 1)
    expect_is(spec1$samp[[1]], "formula")

    spec2 = make_spec(y ~ normal(mu, sigma), sigma ~ gamma(alpha, beta))
    expect_equal(spec1$params, list())
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
    expect_is(as.data.frame(spec2), "data.frame")
    expect_is(as.data.frame(make_spec(y ~ normal(), x ~ normal())), "data.frame")
})