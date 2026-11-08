test_that("refitted coefficients accept the reported MKL rounding difference", {
  observed <- 2.04085310337
  reference <- 2.04085310316
  expect_failure(expect_equal(observed, reference, tolerance = 1e-10))
  expect_equal(observed, reference, tolerance = model_fit_tolerance)
  expect_equal(1e-9, 0, tolerance = model_fit_tolerance)
})

test_that("refit tolerance still rejects larger differences and invalid results", {
  reference <- 2.04085310316
  expect_failure(expect_equal(reference + 1e-5, reference,
    tolerance = model_fit_tolerance))
  expect_failure(expect_equal(c(0.2, 0.3, 0.5), c(0.2 + 1e-5, 0.3, 0.5 - 1e-5),
    tolerance = model_fit_tolerance))
  expect_failure(expect_equal(1e-5, 0, tolerance = model_fit_tolerance))
  for (invalid in c(NA_real_, NaN, Inf, -Inf)) {
    expect_failure(expect_equal(invalid, reference, tolerance = model_fit_tolerance))
  }
  expect_failure(expect_equal(list(coefficient = reference, status = "failed"),
    list(coefficient = reference, status = "convergence"),
    tolerance = model_fit_tolerance))
})
