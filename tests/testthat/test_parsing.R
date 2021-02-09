context("Stan model parsing")

test_that("Empty model handled correctly", {
    parsed = parse_model("")
    expect_equal(length(parsed$vars), 0)
    expect_equal(length(parsed$samps), 0)
})

test_that("Correct parsed variables", {
    correct_vars = c(J= "data", y= "data", sigma="data", mu="parameters",
                     tau="parameters", eta="parameters",
                     theta="transformed parameters")
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    expect_equal(parsed$vars, correct_vars)
})

test_that("Correct parsed sampling statements", {
    correct_samp = list(eta ~ std_normal(), y ~ normal(theta, sigma),
                        mu ~ uniform(-1e+100, 1e+100), tau ~ uniform(-1e+100, 1e+100))
    code = eightschools_m@stanmodel@model_code
    parsed = parse_model(code)
    expect_equal(parsed$samp, correct_samp)
})

test_that("Provided sampling statements can be matched to model", {
    model_samp = list(eta ~ std_normal(), y ~ normal(theta, sigma))
    prov_samp = list(eta ~ exponential(5))
    matched = match_sampling_stmts(prov_samp, model_samp)
    expect_false(identical(matched, prov_samp))
    expect_length(matched, 1)
    expect_equal(matched[[1]], model_samp[[1]])
})

test_that("Extra sampling statements not in model throw an error", {
    model_samp = list(eta ~ std_normal(), y ~ normal(theta, sigma))
    prov_samp = list(eta ~ exponential(5), x ~ normal(theta, sigma))
    expect_error(match_sampling_stmts(prov_samp, model_samp),
                 "No matching sampling statement found for x ~ normal\\(theta, sigma\\)")
})

test_that("Variables are correctly extracted from sampling statements", {
    expect_equal(get_stmt_vars(y ~ normal(theta %*% eta/2, 2%%4 + sigma)),
                 c("y", "theta", "eta", "sigma"))
    expect_error(get_stmt_vars(y ~ .), "y ~ \\. does not")
})
