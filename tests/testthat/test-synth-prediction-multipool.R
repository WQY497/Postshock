make_multipool_toy_data <- function() {
  set.seed(321)
  
  n <- 100L
  shock_time <- 75L
  shock_length <- 6L
  
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
  
  # Add distinct post-shock level shifts to the donor series
  post_rows <- (shock_time + 1L):(shock_time + shock_length)
  
  donor1[post_rows] <- donor1[post_rows] + 1.0
  donor2[post_rows] <- donor2[post_rows] + 2.0
  
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
    x1 = seq_len(n) / n + 1.5,
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


run_multipool_toy <- function(
    dat,
    pools,
    pool_agg = "sum",
    pool_weights = NULL,
    k = 2L,
    base_use_dbw = FALSE,
    pool_use_dbw = TRUE
) {
  SynthPredictionMultiPool(
    Y_series_list = dat$Y,
    covariates_series_list = dat$X,
    shock_time_vec = dat$shock_time,
    shock_length_vec = dat$shock_length,
    k = k,
    pools = pools,
    pool_agg = pool_agg,
    pool_weights = pool_weights,
    base_use_dbw = base_use_dbw,
    pool_use_dbw = pool_use_dbw,
    arima_order = c(1L, 0L, 0L),
    seasonal = FALSE,
    plots = FALSE
  )
}


test_that("SynthPredictionMultiPool returns a valid sum-aggregated forecast", {
  dat <- make_multipool_toy_data()
  
  pools <- list(
    channel_one = "donor1",
    channel_two = "donor2"
  )
  
  out <- expect_no_warning(
    run_multipool_toy(
      dat = dat,
      pools = pools,
      pool_agg = "sum"
    )
  )
  
  # Check the top-level output structure
  expect_named(
    out,
    c("linear_combinations", "predictions", "meta"),
    ignore.order = FALSE
  )
  
  # Check that the MultiPool forecast was added to the prediction object
  expect_named(
    out$predictions,
    c(
      "unadjusted",
      "adjusted",
      "arithmetic_mean",
      "adjusted_multipool"
    ),
    ignore.order = FALSE
  )
  
  # Check the MultiPool metadata structure
  expect_named(
    out$meta$multi_pool,
    c(
      "pools",
      "target_name",
      "pool_effects",
      "pool_agg",
      "pool_weights",
      "aggregated_effect",
      "base_use_dbw",
      "pool_use_dbw"
    ),
    ignore.order = FALSE
  )
  
  # Check that one effect was estimated for each pool
  expect_length(
    out$meta$multi_pool$pool_effects,
    2L
  )
  
  expect_identical(
    names(out$meta$multi_pool$pool_effects),
    names(pools)
  )
  
  expect_true(
    all(is.finite(out$meta$multi_pool$pool_effects))
  )
  
  # Under sum aggregation, the total effect should equal
  # the sum of the individual pool effects
  expect_equal(
    out$meta$multi_pool$aggregated_effect,
    sum(out$meta$multi_pool$pool_effects),
    tolerance = 1e-8
  )
  
  # The MultiPool forecast should equal the raw forecast
  # plus the aggregated pool effect
  expect_equal(
    out$predictions$adjusted_multipool,
    out$predictions$unadjusted +
      out$meta$multi_pool$aggregated_effect,
    tolerance = 1e-8
  )
  
  expect_length(
    out$predictions$adjusted_multipool,
    2L
  )
  
  expect_true(
    all(is.finite(out$predictions$adjusted_multipool))
  )
  
  expect_identical(
    out$meta$multi_pool$pool_agg,
    "sum"
  )
})


test_that("SynthPredictionMultiPool applies user-specified pool weights", {
  dat <- make_multipool_toy_data()
  
  pools <- list(
    channel_one = "donor1",
    channel_two = "donor2"
  )
  
  pool_weights <- c(
    channel_one = 0.25,
    channel_two = 0.75
  )
  
  out <- expect_no_warning(
    run_multipool_toy(
      dat = dat,
      pools = pools,
      pool_agg = "weighted",
      pool_weights = pool_weights
    )
  )
  
  pool_effects <- out$meta$multi_pool$pool_effects
  
  expected_effect <- sum(
    pool_weights[names(pool_effects)] * pool_effects
  )
  
  # Check the weighted aggregation formula
  expect_equal(
    out$meta$multi_pool$aggregated_effect,
    expected_effect,
    tolerance = 1e-8
  )
  
  # Check the resulting MultiPool forecast
  expect_equal(
    out$predictions$adjusted_multipool,
    out$predictions$unadjusted + expected_effect,
    tolerance = 1e-8
  )
  
  expect_identical(
    out$meta$multi_pool$pool_agg,
    "weighted"
  )
  
  expect_equal(
    out$meta$multi_pool$pool_weights,
    pool_weights
  )
})


test_that("SynthPredictionMultiPool supports a single-donor pool without warnings", {
  dat <- make_multipool_toy_data()
  
  pools <- list(
    single_channel = "donor1"
  )
  
  out <- expect_no_warning(
    run_multipool_toy(
      dat = dat,
      pools = pools,
      pool_agg = "sum",
      pool_use_dbw = TRUE
    )
  )
  
  pool_effects <- out$meta$multi_pool$pool_effects
  
  # A single pool should produce exactly one pool effect
  expect_length(pool_effects, 1L)
  expect_true(is.finite(pool_effects[[1L]]))
  
  # With one pool and sum aggregation, the aggregated effect
  # should equal that pool's estimated effect
  expect_equal(
    out$meta$multi_pool$aggregated_effect,
    pool_effects[[1L]],
    tolerance = 1e-8
  )
  
  expect_equal(
    out$predictions$adjusted_multipool,
    out$predictions$unadjusted + pool_effects[[1L]],
    tolerance = 1e-8
  )
})


test_that("SynthPredictionMultiPool rejects unnamed series lists", {
  dat <- make_multipool_toy_data()
  
  dat$Y <- unname(dat$Y)
  
  pools <- list(
    channel_one = "donor1"
  )
  
  expect_error(
    run_multipool_toy(
      dat = dat,
      pools = pools
    ),
    regexp = "must be a.*named.*list"
  )
})


test_that("SynthPredictionMultiPool rejects unknown donor names", {
  dat <- make_multipool_toy_data()
  
  pools <- list(
    invalid_channel = "unknown_donor"
  )
  
  expect_error(
    run_multipool_toy(
      dat = dat,
      pools = pools,
      pool_use_dbw = FALSE
    ),
    regexp = "unknown Y series"
  )
})


test_that("SynthPredictionMultiPool rejects pools containing only the target", {
  dat <- make_multipool_toy_data()
  
  pools <- list(
    invalid_channel = "target"
  )
  
  expect_error(
    run_multipool_toy(
      dat = dat,
      pools = pools,
      pool_use_dbw = FALSE
    ),
    regexp = "at least one donor"
  )
})


test_that("SynthPredictionMultiPool requires weights for weighted aggregation", {
  dat <- make_multipool_toy_data()
  
  pools <- list(
    channel_one = "donor1",
    channel_two = "donor2"
  )
  
  expect_error(
    run_multipool_toy(
      dat = dat,
      pools = pools,
      pool_agg = "weighted",
      pool_weights = NULL,
      pool_use_dbw = FALSE
    ),
    regexp = "pool_weights is required"
  )
})


test_that("SynthPredictionMultiPool requires weights for every pool", {
  dat <- make_multipool_toy_data()
  
  pools <- list(
    channel_one = "donor1",
    channel_two = "donor2"
  )
  
  incomplete_weights <- c(
    channel_one = 1
  )
  
  expect_error(
    run_multipool_toy(
      dat = dat,
      pools = pools,
      pool_agg = "weighted",
      pool_weights = incomplete_weights,
      pool_use_dbw = FALSE
    ),
    regexp = "must include all pool names"
  )
})