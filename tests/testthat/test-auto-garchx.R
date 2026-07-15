simulate_garch11_data <- function(
    n = 300L,
    omega = 0.08,
    alpha = 0.08,
    beta = 0.85,
    seed = 123
) {
  set.seed(seed)
  
  innovations <- stats::rnorm(n)
  variance <- numeric(n)
  y <- numeric(n)
  
  variance[1L] <- omega / (1 - alpha - beta)
  y[1L] <- sqrt(variance[1L]) * innovations[1L]
  
  for (i in 2:n) {
    variance[i] <-
      omega +
      alpha * y[i - 1L]^2 +
      beta * variance[i - 1L]
    
    y[i] <- sqrt(variance[i]) * innovations[i]
  }
  
  y
}


test_that("auto_garchx returns a valid BIC-selected model", {
  y <- simulate_garch11_data()
  
  out <- auto_garchx(
    y = y,
    max.p = 1L,
    max.q = 1L,
    search_o = 0L,
    final_refit = FALSE,
    vcov.type = "ordinary",
    verbose = FALSE
  )
  
  # Check the complete return structure
  expect_named(
    out,
    c(
      "bic",
      "aic",
      "fit",
      "q",
      "p",
      "o",
      "refit_used"
    ),
    ignore.order = FALSE
  )
  
  # Check that a usable GARCH-X model was returned
  expect_s3_class(out$fit, "garchx")
  
  # Check that the selected orders are valid
  expect_true(out$q %in% 0:1)
  expect_true(out$p %in% 0:1)
  expect_identical(out$o, 0L)
  
  # The trivial GARCH(0, 0, 0) model should never be selected
  expect_false(
    out$q == 0L &&
      out$p == 0L &&
      out$o == 0L
  )
  
  # Check information criteria
  expect_type(out$bic, "double")
  expect_type(out$aic, "double")
  expect_length(out$bic, 1L)
  expect_length(out$aic, 1L)
  expect_true(is.finite(out$bic))
  expect_true(is.finite(out$aic))
  
  # No refit should be used when final_refit is FALSE
  expect_identical(out$refit_used, FALSE)
})


test_that("auto_garchx reports information criteria consistently", {
  y <- simulate_garch11_data(seed = 456)
  
  out <- auto_garchx(
    y = y,
    max.p = 1L,
    max.q = 1L,
    search_o = 0L,
    final_refit = FALSE,
    vcov.type = "ordinary",
    verbose = FALSE
  )
  
  log_likelihood <- as.numeric(
    stats::logLik(out$fit)
  )
  
  parameter_count <- length(
    stats::coef(out$fit)
  )
  
  expected_bic <-
    -2 * log_likelihood +
    parameter_count * log(length(y))
  
  expected_aic <-
    -2 * log_likelihood +
    2 * parameter_count
  
  expect_equal(
    out$bic,
    expected_bic,
    tolerance = 1e-8
  )
  
  expect_equal(
    out$aic,
    expected_aic,
    tolerance = 1e-8
  )
})


test_that("auto_garchx returns a model that can generate forecasts", {
  y <- simulate_garch11_data(seed = 789)
  
  out <- auto_garchx(
    y = y,
    max.p = 1L,
    max.q = 1L,
    search_o = 0L,
    final_refit = FALSE,
    vcov.type = "ordinary",
    verbose = FALSE
  )
  
  forecast <- stats::predict(
    out$fit,
    n.ahead = 3L
  )
  
  forecast_values <- if (is.list(forecast)) {
    as.numeric(forecast$pred)
  } else {
    as.numeric(forecast)
  }
  
  expect_length(forecast_values, 3L)
  expect_true(all(is.finite(forecast_values)))
  expect_true(all(forecast_values > 0))
})


test_that("auto_garchx supports an exogenous regressor", {
  y <- simulate_garch11_data(seed = 321)
  n <- length(y)
  
  shock_indicator <- as.numeric(
    seq_len(n) > 220L
  )
  
  # Introduce a moderate variance increase associated with the regressor
  y_with_xreg <- y * sqrt(
    1 + 0.5 * shock_indicator
  )
  
  xreg <- matrix(
    shock_indicator,
    ncol = 1L,
    dimnames = list(
      NULL,
      "shock_indicator"
    )
  )
  
  out <- auto_garchx(
    y = y_with_xreg,
    xreg = xreg,
    max.p = 1L,
    max.q = 1L,
    search_o = 0L,
    final_refit = FALSE,
    vcov.type = "ordinary",
    verbose = FALSE
  )
  
  expect_s3_class(out$fit, "garchx")
  expect_true(is.finite(out$bic))
  expect_true(is.finite(out$aic))
  expect_true(out$q %in% 0:1)
  expect_true(out$p %in% 0:1)
})


test_that("auto_garchx handles the optional refit branch", {
  y <- simulate_garch11_data(seed = 654)
  
  out <- auto_garchx(
    y = y,
    max.p = 1L,
    max.q = 1L,
    search_o = 0L,
    final_refit = TRUE,
    clamp_factor = c(0.1, 10),
    vcov.type = "ordinary",
    verbose = FALSE
  )
  
  expect_s3_class(out$fit, "garchx")
  expect_type(out$refit_used, "logical")
  expect_length(out$refit_used, 1L)
  expect_true(is.finite(out$bic))
  expect_true(is.finite(out$aic))
  
  # The refit may or may not improve BIC.
  # Both TRUE and FALSE are valid outcomes.
  expect_true(out$refit_used %in% c(TRUE, FALSE))
})


test_that("auto_garchx rejects invalid response inputs", {
  expect_error(
    auto_garchx(
      y = "not numeric"
    ),
    regexp = "'y' must be a numeric vector",
    fixed = TRUE
  )
  
  expect_error(
    auto_garchx(
      y = matrix(
        stats::rnorm(100),
        ncol = 1L
      )
    ),
    regexp = "'y' must be a numeric vector",
    fixed = TRUE
  )
  
  y_with_na <- stats::rnorm(100)
  y_with_na[10L] <- NA_real_
  
  expect_error(
    auto_garchx(
      y = y_with_na
    ),
    regexp = "'y' contains NA/Inf/NaN",
    fixed = TRUE
  )
  
  y_with_inf <- stats::rnorm(100)
  y_with_inf[10L] <- Inf
  
  expect_error(
    auto_garchx(
      y = y_with_inf
    ),
    regexp = "'y' contains NA/Inf/NaN",
    fixed = TRUE
  )
})


test_that("auto_garchx validates exogenous regressors", {
  y <- simulate_garch11_data(
    n = 200L
  )
  
  expect_error(
    auto_garchx(
      y = y,
      xreg = matrix(
        stats::rnorm(199),
        ncol = 1L
      )
    ),
    regexp = "nrow\\(xreg\\) must equal length\\(y\\)"
  )
  
  xreg_with_na <- matrix(
    stats::rnorm(200),
    ncol = 1L
  )
  
  xreg_with_na[25L, 1L] <- NA_real_
  
  expect_error(
    auto_garchx(
      y = y,
      xreg = xreg_with_na
    ),
    regexp = "'xreg' contains NA/Inf/NaN",
    fixed = TRUE
  )
})


test_that("auto_garchx rejects invalid order bounds", {
  y <- simulate_garch11_data(
    n = 200L
  )
  
  expect_error(
    auto_garchx(
      y = y,
      max.p = -1L,
      max.q = 1L
    ),
    regexp = "must be nonnegative integers"
  )
  
  expect_error(
    auto_garchx(
      y = y,
      max.p = 1L,
      max.q = -1L
    ),
    regexp = "must be nonnegative integers"
  )
  
  expect_error(
    auto_garchx(
      y = y,
      max.p = c(1L, 2L),
      max.q = 1L
    ),
    regexp = "must be nonnegative integers"
  )
})


test_that("auto_garchx validates clamp factors and covariance type", {
  y <- simulate_garch11_data(
    n = 200L
  )
  
  expect_error(
    auto_garchx(
      y = y,
      clamp_factor = 1
    ),
    regexp = "positive numeric vector of length 2"
  )
  
  expect_error(
    auto_garchx(
      y = y,
      clamp_factor = c(-0.1, 10)
    ),
    regexp = "positive numeric vector of length 2"
  )
  
  expect_error(
    auto_garchx(
      y = y,
      clamp_factor = c(NA_real_, 10)
    ),
    regexp = "positive numeric vector of length 2"
  )
  
  expect_error(
    auto_garchx(
      y = y,
      vcov.type = "invalid"
    ),
    regexp = "should be one of"
  )
})


test_that("auto_garchx reports failure when the search grid has no model", {
  y <- simulate_garch11_data(
    n = 200L
  )
  
  # The only candidate is GARCH(0, 0, 0), which the function skips.
  expect_error(
    auto_garchx(
      y = y,
      max.p = 0L,
      max.q = 0L,
      search_o = 0L,
      final_refit = FALSE,
      verbose = FALSE
    ),
    regexp = "no converged model in grid search"
  )
})


test_that("auto_garch is an alias of auto_garchx", {
  expect_identical(
    auto_garch,
    auto_garchx
  )
})