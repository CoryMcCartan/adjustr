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
    obj = adjust_weights(spec, eightschools_m, keep_bad=T)
    obj = get_resampling_idxs(obj)
    expect_is(obj$.idxs, "list")
    expect_length(obj$.idxs[[1]], 20)
})

test_that("Resampling indices fail with negative `frac`", {
    expect_error(get_resampling_idxs(c(0, 0, 1), frac=-0.5), "must be nonnegative")
})