context("Stan model parsing")


test_that("Model parses and returns parsing table", {
    code = test_env$eightschools_m@stanmodel@model_code
    parsed_model = parse_model(code)
    expect_equal(names(parsed_model), c("name", "value", "pos", "depth", "i"))
})

test_that("Correct parsed variables", {
    correct_vars = c(J= "data", y= "data", sigma="data", mu="parameters",
                     tau="parameters", eta="parameters",
                     theta="transformed parameters")
    code = test_env$eightschools_m@stanmodel@model_code
    parsed_model = parse_model(code)
    expect_equal(get_variables(parsed_model), correct_vars)
})

test_that("Correct parsed sampling statements", {
    code = test_env$eightschools_m@stanmodel@model_code
    parsed_model = parse_model(code)
    expect_equal(get_sampling_stmts(parsed_model), correct_samp)
})