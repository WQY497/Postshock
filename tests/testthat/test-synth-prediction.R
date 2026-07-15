make_synthprediction_toy_data <- function() {
  set.seed(123)
  
  n <- 90L
  shock_time <- 70L
  shock_length <- 5L
  
  # Generate stationary target and donor series
  target <- as.numeric(
    stats::arima.sim(
      model = list(ar = 0.4),
      n = n,
      sd = 0.3
    )
  ) + 10
  
  donor1 <- as.numeric(
    stats::arima.sim(
      model = list(ar = 0.3),
      n = n,
      sd = 0.3
    )
  ) + 10
  
  donor2 <- as.numeric(
    stats::arima.sim(
      model = list(ar = 0.2),
      n = n,
      sd = 0.3
    )
  ) + 10
  
  # Add known post-shock level shifts to the donor series
  post_rows <- (shock_time + 1L):(shock_time + shock_length)
  
  donor1[post_rows] <- donor1[post_rows] + 1.5
  donor2[post_rows] <- donor2[post_rows] + 3.0
  
  # Construct pre-shock covariates for donor matching
  x_target <- cbind(
    x1 = seq_len(n) / n,
    x2 = sin(seq_len(n) / 10)
  )
  
  # Donor 1 is designed to be similar to the target
  x_donor1 <- x_target +
    matrix(
      rnorm(2L * n, sd = 0.02),
      nrow = n,
      ncol = 2L
    )
  
  colnames(x_donor1) <- c("x1", "x2")
  
  # Donor 2 is designed to be less similar to the target
  x_donor2 <- cbind(
    x1 = seq_len(n) / n + 2,
    x2 = cos(seq_len(n) / 7)
  )
  
  list(
    Y = list(
      target = target,
      donor1 = donor1,
      donor2 = donor2
    ),
    X = list(
      target = x_target,
      donor1 = x_donor1,
      donor2 = x_donor2
    ),
    shock_time = rep(shock_time, 3L),
    shock_length = rep(shock_length, 3L)
  )
}


test_that("SynthPrediction returns valid k-step forecasts", {
  dat <- make_synthprediction_toy_data()
  
  out <- SynthPrediction(
    Y_series_list = dat$Y,
    covariates_series_list = dat$X,
    shock_time_vec = dat$shock_time,
    shock_length_vec = dat$shock_length,
    k = 3L,
    use_dbw = FALSE,
    arima_order = c(1L, 0L, 0L),
    seasonal = FALSE,
    plots = FALSE
  )
  
  # Check the top-level output structure
  expect_named(
    out,
    c("linear_combinations", "predictions", "meta"),
    ignore.order = FALSE
  )
  
  # Check the forecast output structure
  expect_named(
    out$predictions,
    c("unadjusted", "adjusted", "arithmetic_mean"),
    ignore.order = FALSE
  )
  
  # Check that all forecast vectors have length k
  expect_length(out$predictions$unadjusted, 3L)
  expect_length(out$predictions$adjusted, 3L)
  expect_length(out$predictions$arithmetic_mean, 3L)
  
  # Check that all forecast values are finite
  expect_true(all(is.finite(out$predictions$unadjusted)))
  expect_true(all(is.finite(out$predictions$adjusted)))
  expect_true(all(is.finite(out$predictions$arithmetic_mean)))
  
  # Check donor metadata
  expect_identical(out$meta$n_donors, 2L)
  expect_length(out$meta$omega_vec, 2L)
})


test_that("SynthPrediction uses equal weights when DBW is disabled", {
  dat <- make_synthprediction_toy_data()
  
  out <- SynthPrediction(
    Y_series_list = dat$Y,
    covariates_series_list = dat$X,
    shock_time_vec = dat$shock_time,
    shock_length_vec = dat$shock_length,
    k = 2L,
    use_dbw = FALSE,
    arima_order = c(1L, 0L, 0L),
    seasonal = FALSE,
    plots = FALSE
  )
  
  # Two donors should receive equal weights
  expect_equal(
    out$linear_combinations,
    c(0.5, 0.5),
    tolerance = 1e-10
  )
  
  # Donor weights should sum to one
  expect_equal(
    sum(out$linear_combinations),
    1,
    tolerance = 1e-10
  )
  
  # The adjusted forecast should equal the baseline plus
  # the combined donor shock effect
  expect_equal(
    out$predictions$adjusted,
    out$predictions$unadjusted +
      out$meta$combined_omega,
    tolerance = 1e-8
  )
  
  # Under equal donor weights, the weighted adjustment
  # should equal the arithmetic-mean adjustment
  expect_equal(
    out$predictions$adjusted,
    out$predictions$arithmetic_mean,
    tolerance = 1e-8
  )
})


test_that("SynthPrediction returns feasible DBW weights", {
  dat <- make_synthprediction_toy_data()
  
  out <- SynthPrediction(
    Y_series_list = dat$Y,
    covariates_series_list = dat$X,
    shock_time_vec = dat$shock_time,
    shock_length_vec = dat$shock_length,
    k = 2L,
    use_dbw = TRUE,
    dbw_indices = 1:2,
    dbw_center = FALSE,
    dbw_scale = FALSE,
    penalty_normchoice = "l2",
    penalty_lambda = 1e-2,
    arima_order = c(1L, 0L, 0L),
    seasonal = FALSE,
    plots = FALSE
  )
  
  # Check the donor weight vector
  expect_type(out$linear_combinations, "double")
  expect_length(out$linear_combinations, 2L)
  expect_true(all(is.finite(out$linear_combinations)))
  
  # DBW weights should satisfy the simplex constraint
  expect_equal(
    sum(out$linear_combinations),
    1,
    tolerance = 1e-6
  )
  
  # DBW weights should satisfy the lower and upper bounds
  expect_true(
    all(out$linear_combinations >= -1e-7)
  )
  
  expect_true(
    all(out$linear_combinations <= 1 + 1e-7)
  )
  
  # The DBW optimization should report successful convergence
  expect_identical(
    out$meta$dbw_status,
    "convergence"
  )
  
  # The adjusted forecast should equal the baseline plus
  # the DBW-weighted donor shock effect
  expect_equal(
    out$predictions$adjusted,
    out$predictions$unadjusted +
      out$meta$combined_omega,
    tolerance = 1e-8
  )
})

test_that("SynthPrediction maps Date shock times to the correct rows", {
  dat <- make_synthprediction_toy_data()
  
  idx <- seq.Date(
    from = as.Date("2023-01-01"),
    by = "day",
    length.out = 90L
  )
  
  # Give every series and covariate matrix the same Date index
  Y_xts <- lapply(
    dat$Y,
    function(x) xts::xts(x, order.by = idx)
  )
  
  X_xts <- lapply(
    dat$X,
    function(x) xts::xts(x, order.by = idx)
  )
  
  # Reference result using the integer row number
  out_integer <- SynthPrediction(
    Y_series_list = Y_xts,
    covariates_series_list = X_xts,
    shock_time_vec = rep(70L, 3L),
    shock_length_vec = dat$shock_length,
    k = 2L,
    use_dbw = FALSE,
    arima_order = c(1L, 0L, 0L),
    seasonal = FALSE,
    plots = FALSE
  )
  
  # Same shock, now specified using the Date label at row 70
  out_date <- SynthPrediction(
    Y_series_list = Y_xts,
    covariates_series_list = X_xts,
    shock_time_vec = rep(idx[70L], 3L),
    shock_length_vec = dat$shock_length,
    k = 2L,
    use_dbw = FALSE,
    arima_order = c(1L, 0L, 0L),
    seasonal = FALSE,
    plots = FALSE
  )
  
  # Date row 70 and integer row 70 must produce the same result
  expect_equal(
    out_date$linear_combinations,
    out_integer$linear_combinations,
    tolerance = 1e-10
  )
  
  expect_equal(
    out_date$predictions$unadjusted,
    out_integer$predictions$unadjusted,
    tolerance = 1e-8
  )
  
  expect_equal(
    out_date$predictions$adjusted,
    out_integer$predictions$adjusted,
    tolerance = 1e-8
  )
  
  expect_equal(
    out_date$meta$omega_vec,
    out_integer$meta$omega_vec,
    tolerance = 1e-8
  )
})