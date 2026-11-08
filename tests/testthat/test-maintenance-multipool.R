test_that("MultiPool forwards shared non-default settings into every internal call", {
  arguments <- maintenance_input()
  arguments$use_dbw <- NULL
  arguments$pools <- list(pair = c("donorA", "donorB"), single = "donorC")
  settings <- list(dbw_scale = FALSE, dbw_center = FALSE, dbw_indices = 2L,
    princ_comp_input = 1L, covariate_indices = 1L, geometric_sets = list("reserved"),
    days_before_shocktime_vec = 3L, arima_order = c(2L, 0L, 1L),
    seasonal = TRUE, seasonal_order = c(1L, 0L, 0L), seasonal_period = 7L,
    user_ic_choice = "bic", stepwise = FALSE, approximation = FALSE,
    penalty_lambda = 0.4, penalty_normchoice = "l2",
    control_shocks = list(list(list(time = 9L, length = 1L, shape = "point"))))
  arguments[names(settings)] <- settings
  arguments$base_use_dbw <- FALSE
  arguments$pool_use_dbw <- TRUE
  outer <- SynthPrediction
  calls <- list()
  testthat::local_mocked_bindings(SynthPrediction = function(...) {
    values <- list(...)
    calls[[length(calls) + 1L]] <<- values
    list(predictions = list(unadjusted = c(10, 11)), meta = list(combined_omega = 2))
  }, .package = "postshock")
  direct <- do.call(SynthPredictionMultiPool, arguments)
  direct_calls <- calls
  calls <- list()
  wrapped <- do.call(outer, arguments)
  expect_equal(wrapped, direct)
  expect_equal(lapply(calls, function(call) call[sort(names(call))]),
               lapply(direct_calls, function(call) call[sort(names(call))]))
  expect_length(calls, 3L)
  for (call in calls) {
    expect_equal(call[names(settings)], settings)
    expect_null(call$pools)
  }
  expect_identical(vapply(calls, function(call) call$use_dbw, logical(1)), c(FALSE, TRUE, TRUE))
  expect_equal(wrapped$predictions$adjusted_multipool, c(14, 15))
})

test_that("MultiPool fixed ARIMA and L2 settings agree for real direct and wrapper fits", {
  arguments <- maintenance_input()
  arguments$use_dbw <- NULL
  arguments$pools <- list(pair = c("donorA", "donorB"), single = "donorC")
  arguments$dbw_scale <- FALSE
  arguments$dbw_center <- FALSE
  arguments$dbw_indices <- 1:2
  arguments$days_before_shocktime_vec <- 2L
  arguments$penalty_lambda <- 0.35
  arguments$penalty_normchoice <- "l2"
  direct <- suppressWarnings(do.call(SynthPredictionMultiPool, arguments))
  wrapped <- suppressWarnings(do.call(SynthPrediction, arguments))
  expect_equal(wrapped, direct, tolerance = model_fit_tolerance)
  expect_equal(wrapped$meta$arima_order[[2L]][1:3], c(1, 0, 0))
})

test_that("MultiPool rejects conflicting duplicated internal settings explicitly", {
  arguments <- maintenance_input()
  arguments$pools <- list(all = c("donorA", "donorB", "donorC"))
  expect_error(do.call(SynthPredictionMultiPool, arguments), "use_dbw.*base_use_dbw.*pool_use_dbw")
  arguments$use_dbw <- NULL
  expect_error(do.call(SynthPredictionMultiPool, c(arguments, list(dbw_scale = TRUE, dbw_scale = FALSE))),
               "duplicate.*dbw_scale")
})
