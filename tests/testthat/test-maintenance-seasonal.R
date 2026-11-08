test_that("automatic mean models receive explicit frequencies without altering rows", {
  for (period in c(7L, 12L)) {
    arguments <- maintenance_input()
    arguments$arima_order <- NULL
    arguments$seasonal <- TRUE
    arguments$seasonal_period <- period
    arguments$covariate_indices <- 1L
    arguments$k <- 3L
    calls <- list()
    testthat::local_mocked_bindings(auto.arima = function(y, xreg = NULL, seasonal, ...) {
      calls[[length(calls) + 1L]] <<- list(values = as.numeric(y), frequency = stats::frequency(y),
        seasonal = seasonal, xreg = xreg)
      forecast::Arima(as.numeric(y), order = c(1, 0, 0), seasonal = FALSE, xreg = xreg)
    }, .package = "forecast")
    fit <- suppressWarnings(do.call(SynthPrediction, arguments))
    expect_equal(vapply(calls, function(call) call$frequency, numeric(1)), rep(period, 4))
    expect_true(all(vapply(calls, function(call) call$seasonal, logical(1))))
    for (donor in 1:3) {
      expect_identical(calls[[donor]]$values, as.numeric(arguments$Y_series_list[[donor + 1L]][1:73]))
      expect_equal(as.numeric(calls[[donor]]$xreg[, 1L]), as.numeric(arguments$covariates_series_list[[donor + 1L]][1:73, 1L]))
    }
    expect_identical(calls[[4L]]$values, as.numeric(arguments$Y_series_list[[1L]][1:70]))
    expect_equal(dim(calls[[4L]]$xreg), c(70, 1))
    expect_length(fit$predictions$adjusted, 3L)
  }
})

test_that("detected frequencies use each training window and safely disable failed detection", {
  for (detected in list(7L, 1L, NA_real_, structure("failed", class = "try-error"))) {
    arguments <- maintenance_input()
    arguments$arima_order <- NULL
    arguments$seasonal <- TRUE
    windows <- list()
    flags <- logical()
    frequencies <- numeric()
    testthat::local_mocked_bindings(findfrequency = function(x) {
      windows[[length(windows) + 1L]] <<- x
      if (inherits(detected, "try-error")) stop("detection failed")
      detected
    }, auto.arima = function(y, xreg = NULL, seasonal, ...) {
      flags <<- c(flags, seasonal)
      frequencies <<- c(frequencies, stats::frequency(y))
      forecast::Arima(as.numeric(y), order = c(1, 0, 0), seasonal = FALSE, xreg = xreg)
    }, .package = "forecast")
    suppressWarnings(do.call(SynthPrediction, arguments))
    active <- identical(detected, 7L)
    expect_equal(flags, rep(active, 4L))
    expect_equal(frequencies, rep(if (active) 7 else 1, 4L))
    expect_equal(lengths(windows), c(73L, 73L, 73L, 70L))
    expect_identical(windows[[4L]], as.numeric(arguments$Y_series_list[[1L]][1:70]))
  }
})

test_that("seasonal periods are strict integers but ignored with seasonal FALSE", {
  arguments <- maintenance_input()
  arguments$seasonal <- TRUE
  for (period in list(0, -1, 1.5, NA_real_, NaN, Inf, numeric(), c(7, 12), "7", TRUE, 7 + 0i)) {
    arguments["seasonal_period"] <- list(period)
    expect_error(do.call(SynthPrediction, arguments), "seasonal_period.*single finite positive integer")
  }
  arguments <- maintenance_input()
  original <- suppressWarnings(do.call(SynthPrediction, arguments))
  arguments$seasonal_period <- "unused"
  expect_equal(suppressWarnings(do.call(SynthPrediction, arguments)), original)
  arguments$seasonal <- TRUE
  arguments$seasonal_period <- 1L
  expect_equal(suppressWarnings(do.call(SynthPrediction, arguments)), original)
})

test_that("real seasonal detection and fixed seasonal orders remain supported", {
  set.seed(8712)
  values <- sin(2 * pi * seq_len(144) / 12) + stats::rnorm(144, sd = 0.04)
  detected <- forecast::findfrequency(values)
  expect_gt(detected, 1)
  arguments <- maintenance_input()
  arguments$seasonal <- TRUE
  arguments$seasonal_period <- 12L
  arguments$seasonal_order <- c(1L, 0L, 0L)
  fixed <- suppressWarnings(do.call(SynthPrediction, arguments))
  expect_equal(fixed$meta$arima_order[[2L]][3:5], c(1, 0, 12))
  arguments$arima_order <- NULL
  arguments$seasonal_order <- NULL
  arguments$stepwise <- TRUE
  arguments$approximation <- FALSE
  arguments$use_dbw <- NULL
  arguments$pools <- list(all = c("donorA", "donorB", "donorC"))
  direct <- suppressWarnings(do.call(SynthPredictionMultiPool, arguments))
  wrapped <- suppressWarnings(do.call(SynthPrediction, arguments))
  expect_equal(wrapped, direct, tolerance = model_fit_tolerance)
})

test_that("real per-window frequency detection reaches the automatic engine", {
  arguments <- maintenance_input()
  arguments$arima_order <- NULL
  arguments$seasonal <- TRUE
  set.seed(741)
  dates <- as.Date("2020-01-01") + 0:143
  values <- sin(2 * pi * seq_len(144) / 12) + stats::rnorm(144, sd = 0.04)
  arguments$Y_series_list <- lapply(seq_len(4L), function(series)
    xts::xts(values + series, dates))
  arguments$shock_time_vec <- rep(120L, 4L)
  frequencies <- numeric()
  testthat::local_mocked_bindings(auto.arima = function(y, xreg = NULL, seasonal, ...) {
    frequencies <<- c(frequencies, stats::frequency(y))
    expect_true(seasonal)
    forecast::Arima(as.numeric(y), order = c(1, 0, 0), seasonal = FALSE, xreg = xreg)
  }, .package = "forecast")
  first <- suppressWarnings(do.call(SynthPrediction, arguments))
  expected <- c(rep(forecast::findfrequency(values[1:123]), 3), forecast::findfrequency(values[1:120]))
  expect_equal(frequencies, expected)
  expect_true(all(expected > 1))
  arguments$Y_series_list[[1L]][121:144] <- 1e9
  second <- suppressWarnings(do.call(SynthPrediction, arguments))
  expect_equal(second$predictions, first$predictions, tolerance = model_fit_tolerance)
  expect_equal(tail(frequencies, 4L), expected)
})

test_that("donor fallback explicitly remains nonseasonal", {
  arguments <- maintenance_input()
  arguments$arima_order <- NULL
  arguments$seasonal <- TRUE
  arguments$seasonal_period <- 12L
  original_fixed <- forecast::Arima
  count <- 0L
  fallback_frequencies <- numeric()
  testthat::local_mocked_bindings(auto.arima = function(y, ...) {
    count <<- count + 1L
    if (count <= 3L) stop("deliberate failure")
    original_fixed(as.numeric(y), order = c(1, 0, 0), seasonal = FALSE)
  }, Arima = function(y, ...) {
    fallback_frequencies <<- c(fallback_frequencies, stats::frequency(y))
    original_fixed(y, ...)
  }, .package = "forecast")
  warnings <- character()
  fit <- withCallingHandlers(do.call(SynthPrediction, arguments), warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  })
  expect_length(warnings, 3L)
  expect_true(all(grepl("non-seasonal Arima", warnings)))
  expect_equal(head(fallback_frequencies, 3L), rep(1, 3L))
  expect_length(fit$predictions$adjusted, 2L)
})
