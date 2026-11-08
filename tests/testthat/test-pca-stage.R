pca_stage_input <- function() {
  set.seed(20261004)
  count <- 70L
  dates <- as.Date("2020-01-01") + seq_len(count) - 1L
  responses <- lapply(seq_len(4L), function(series) {
    values <- 5 + series / 3 + as.numeric(stats::arima.sim(list(ar = 0.6), count, sd = 0.4))
    if (series > 1L) values[61:63] <- values[61:63] + series
    xts::xts(values, dates)
  })
  features <- rbind(c(1, 2, -1), c(-2, 1, 3), c(2, 3, 1), c(0.4, -2, 0.8))
  covariates <- lapply(seq_len(4L), function(series) {
    values <- matrix(stats::rnorm(count * 3L), count, 3L)
    values[60L, ] <- features[series, ]
    colnames(values) <- c("market", "rate", "index")
    xts::xts(values, dates)
  })
  names(responses) <- names(covariates) <- c("target", "donor1", "donor2", "donor3")
  list(Y_series_list = responses, covariates_series_list = covariates,
    shock_time_vec = rep(60L, 4L), shock_length_vec = rep(3L, 4L), k = 2L,
    dbw_indices = 1:3, arima_order = c(1L, 0L, 0L), seasonal = FALSE)
}

pca_stage_dbw_arguments <- function(arguments, components = NULL, lookback = 1L,
                                    center = TRUE, scale = TRUE,
                                    penalty_lambda = 0, penalty_normchoice = "l1") {
  list(X = lapply(arguments$covariates_series_list, as.matrix),
    dbw_indices = arguments$dbw_indices, shock_time_vec = arguments$shock_time_vec,
    center = center, scale = scale, sum_to_1 = TRUE,
    bounded_below_by = 0, bounded_above_by = 1, normchoice = "l2",
    princ_comp_count = components, X_lookback_indices = lookback,
    penalty_lambda = penalty_lambda, penalty_normchoice = penalty_normchoice)
}

pca_stage_mock <- function(arguments, output = NULL, error_text = NULL) {
  testthat::local_mocked_bindings(dbw = function(...) {
    if (!is.null(error_text)) stop(error_text)
    output
  }, .package = "postshock")
  do.call(SynthPrediction, arguments)
}

test_that("L2-penalized PCA weights match independent DBW and constructed PC scores", {
  arguments <- pca_stage_input()
  arguments$penalty_lambda <- 1e-4
  arguments$penalty_normchoice <- "l2"
  plain <- do.call(SynthPrediction, arguments)
  for (settings in list(c(TRUE, TRUE), c(FALSE, FALSE), c(TRUE, FALSE), c(FALSE, TRUE))) {
    arguments$princ_comp_input <- 1L
    arguments$dbw_center <- settings[1L]
    arguments$dbw_scale <- settings[2L]
    reduced <- do.call(SynthPrediction, arguments)
    reference <- do.call(dbw, pca_stage_dbw_arguments(arguments, 1L,
      center = settings[1L], scale = settings[2L],
      penalty_lambda = arguments$penalty_lambda,
      penalty_normchoice = arguments$penalty_normchoice))
    feature_matrix <- rbind(c(1, 2, -1), c(-2, 1, 3), c(2, 3, 1), c(0.4, -2, 0.8))
    processed <- base::scale(feature_matrix, center = settings[1L], scale = settings[2L])
    scores <- processed %*% svd(processed)$v[, 1L, drop = FALSE]
    score_fit <- dbw(X = lapply(seq_len(4L), function(series) scores[series, , drop = FALSE]),
      dbw_indices = 1L, shock_time_vec = rep(1L, 4L), center = FALSE, scale = FALSE,
      penalty_lambda = arguments$penalty_lambda,
      penalty_normchoice = arguments$penalty_normchoice)
    expect_identical(reduced$meta$dbw_status, "convergence")
    expect_identical(reference$convergence, "convergence")
    expect_identical(score_fit$convergence, "convergence")
    expect_equal(reduced$linear_combinations, reference$opt_params, tolerance = model_fit_tolerance)
    expect_equal(reduced$linear_combinations, score_fit$opt_params, tolerance = model_fit_tolerance)
    expect_equal(reduced$meta$omega_vec, plain$meta$omega_vec, tolerance = model_fit_tolerance)
    expect_equal(reduced$predictions$unadjusted, plain$predictions$unadjusted, tolerance = model_fit_tolerance)
    expect_equal(reduced$meta$combined_omega,
      sum(reference$opt_params * reduced$meta$omega_vec), tolerance = model_fit_tolerance)
    expect_equal(reduced$predictions$adjusted,
      reduced$predictions$unadjusted + reduced$meta$combined_omega, tolerance = 1e-10)
    expect_equal(sum(reduced$linear_combinations), 1, tolerance = 1e-6)
    expect_true(all(reduced$linear_combinations >= -1e-7 &
      reduced$linear_combinations <= 1 + 1e-7))
  }
  arguments$dbw_center <- arguments$dbw_scale <- TRUE
  reduced <- do.call(SynthPrediction, arguments)
  expect_gt(max(abs(reduced$linear_combinations - plain$linear_combinations)), 1e-3)
  expect_identical(names(reduced), names(plain))
})

test_that("unpenalized PCA retains a unique smooth matching example", {
  arguments <- pca_stage_input()
  features <- rbind(c(0, -1), c(-1, 0), c(1, 0), c(0, 1))
  for (series in seq_along(arguments$covariates_series_list)) {
    arguments$covariates_series_list[[series]][60L, 1:2] <- features[series, ]
  }
  arguments$dbw_indices <- 1:2
  arguments$princ_comp_input <- 2L
  arguments$dbw_center <- TRUE
  arguments$dbw_scale <- FALSE
  result <- do.call(SynthPrediction, arguments)
  reference <- do.call(dbw, pca_stage_dbw_arguments(arguments, 2L, scale = FALSE))
  expect_identical(formals(SynthPrediction)$penalty_lambda, 0)
  expect_identical(formals(dbw)$penalty_lambda, 0)
  expect_identical(result$meta$dbw_status, "convergence")
  expect_identical(reference$convergence, "convergence")
  expect_equal(result$linear_combinations, reference$opt_params, tolerance = model_fit_tolerance)
  expect_equal(result$linear_combinations, c(0.5, 0.5, 0), tolerance = 1e-6)
  expect_equal(result$predictions$adjusted, result$predictions$unadjusted +
    sum(result$linear_combinations * result$meta$omega_vec), tolerance = 1e-10)
})

test_that("the original one-component unpenalized example has multiple optima", {
  arguments <- pca_stage_input()
  features <- do.call(rbind, lapply(arguments$covariates_series_list,
    function(series) as.numeric(series[60L, ])))
  processed <- base::scale(features, center = TRUE, scale = FALSE)
  scores <- as.numeric(processed %*% svd(processed)$v[, 1L, drop = FALSE])
  first_share <- (scores[1L] - scores[3L]) / (scores[2L] - scores[3L])
  last_share <- (scores[1L] - scores[3L]) / (scores[4L] - scores[3L])
  solutions <- rbind(c(first_share, 1 - first_share, 0),
    c(0, 1 - last_share, last_share))
  expect_true(all(is.finite(solutions)))
  expect_true(all(solutions >= 0 & solutions <= 1))
  expect_equal(rowSums(solutions), c(1, 1), tolerance = 1e-12)
  expect_equal(as.numeric(solutions %*% scores[-1L]), rep(scores[1L], 2L),
    tolerance = 1e-12)
  expect_gt(max(abs(solutions[1L, ] - solutions[2L, ])), 0.01)
})

test_that("both PCA parameters reject invalid scalar counts before indexing", {
  arguments <- pca_stage_input()
  invalid_values <- list(0, -1, 1.5, NA_real_, NaN, Inf, -Inf, numeric(),
    c(1, 2), "one", TRUE, list(1), 1 + 1i)
  for (invalid in invalid_values) {
    dbw_arguments <- pca_stage_dbw_arguments(arguments)
    dbw_arguments["princ_comp_count"] <- list(invalid)
    arguments["princ_comp_input"] <- list(invalid)
    expect_error(suppressWarnings(do.call(dbw, dbw_arguments)),
      "princ_comp_count.*single finite positive integer")
    expect_error(suppressWarnings(do.call(SynthPrediction, arguments)),
      "princ_comp_input.*single finite positive integer")
  }
  dbw_arguments <- pca_stage_dbw_arguments(arguments, 1e30)
  arguments$princ_comp_input <- 1e30
  expect_error(do.call(dbw, dbw_arguments), "available.*directions")
  expect_error(do.call(SynthPrediction, arguments), "available.*directions")
})

test_that("lookback PCA uses available SVD directions rather than feature count or rank", {
  arguments <- pca_stage_input()
  arguments$days_before_shocktime_vec <- 3L
  arguments$princ_comp_input <- 5L
  dbw_arguments <- pca_stage_dbw_arguments(arguments, 5L, lookback = 1:3)
  expect_error(suppressWarnings(do.call(dbw, dbw_arguments)),
    "requested 5.*available 4.*directions")
  expect_error(suppressWarnings(do.call(SynthPrediction, arguments)),
    "requested 5.*available 4.*directions")
  dbw_arguments$princ_comp_count <- arguments$princ_comp_input <- 4L
  reference <- suppressWarnings(do.call(dbw, dbw_arguments))
  result <- suppressWarnings(do.call(SynthPrediction, arguments))
  expect_equal(result$linear_combinations, reference$opt_params, tolerance = model_fit_tolerance)
  expect_length(result$linear_combinations, 3L)
  feature_matrix <- do.call(rbind, lapply(arguments$covariates_series_list,
    function(series) as.numeric(t(as.matrix(series[c(60L, 59L, 58L), ])))))
  expect_identical(dim(feature_matrix), c(4L, 9L))
  expect_lte(qr(base::scale(feature_matrix))$rank, 3L)
  expect_equal(result$predictions$adjusted, result$predictions$unadjusted +
    sum(result$linear_combinations * result$meta$omega_vec), tolerance = 1e-10)
})

test_that("constant removal and the existing single-feature restriction are respected", {
  arguments <- pca_stage_input()
  for (series in seq_along(arguments$covariates_series_list)) {
    arguments$covariates_series_list[[series]][, 3L] <- 1
  }
  dbw_arguments <- pca_stage_dbw_arguments(arguments, 3L)
  arguments$princ_comp_input <- 3L
  expect_error(suppressWarnings(do.call(dbw, dbw_arguments)),
    "requested 3.*available 2.*directions")
  expect_error(suppressWarnings(do.call(SynthPrediction, arguments)),
    "requested 3.*available 2.*directions")
  arguments$princ_comp_input <- 1L
  dbw_arguments$princ_comp_count <- 1L
  expect_equal(suppressWarnings(do.call(SynthPrediction, arguments))$linear_combinations,
    suppressWarnings(do.call(dbw, dbw_arguments))$opt_params, tolerance = model_fit_tolerance)
  arguments$dbw_indices <- 1L
  dbw_arguments$dbw_indices <- 1L
  expect_error(do.call(dbw, dbw_arguments), "PCA is not meaningful with 1 column")
  expect_error(do.call(SynthPrediction, arguments), "PCA is not meaningful with 1 column")
})

test_that("collinearity does not introduce a new rank-based PCA rejection", {
  arguments <- pca_stage_input()
  for (series in seq_along(arguments$covariates_series_list)) {
    values <- as.numeric(arguments$covariates_series_list[[series]][, 1L])
    arguments$covariates_series_list[[series]][, 2L] <- 2 * values
    arguments$covariates_series_list[[series]][, 3L] <- -values
  }
  result <- suppressWarnings(do.call(dbw, pca_stage_dbw_arguments(arguments, 3L)))
  expect_named(result, c("opt_params", "convergence", "loss"))
  expect_length(result$opt_params, 3L)
  expect_true(all(is.finite(result$opt_params)))
})

test_that("explicit PCA stops on preprocessing and DBW errors without equal fallback", {
  arguments <- pca_stage_input()
  arguments$princ_comp_input <- 1L
  expect_error(pca_stage_mock(arguments, error_text = "controlled SVD or solver error"),
    "explicit PCA.*controlled SVD or solver error")
  for (series in seq_along(arguments$covariates_series_list)) {
    arguments$covariates_series_list[[series]][, ] <- 1
  }
  arguments$dbw_indices <- NULL
  expect_error(do.call(SynthPrediction, arguments), "explicit PCA.*no informative")
  arguments$dbw_indices <- 1:3
  expect_error(suppressWarnings(do.call(SynthPrediction, arguments)),
    "explicit PCA.*All columns were constant")
  arguments <- pca_stage_input()
  arguments$princ_comp_input <- 1L
  arguments$covariates_series_list[[2L]][60L, 2L] <- NA_real_
  expect_error(do.call(SynthPrediction, arguments), "explicit PCA.*NA/Inf")
})

test_that("explicit PCA validates DBW string status and usable simplex weights", {
  arguments <- pca_stage_input()
  arguments$princ_comp_input <- 1L
  for (status in list("failed_convergence", "unknown", NULL, 0L, NA_character_,
                      c("convergence", "convergence"))) {
    expect_error(pca_stage_mock(arguments,
      list(opt_params = c(0.2, 0.3, 0.5), convergence = status, loss = 0)),
      "explicit PCA.*convergence")
  }
  expect_error(pca_stage_mock(arguments, 42), "explicit PCA.*convergence")
  invalid_weights <- list(NULL, numeric(), c(0.5, 0.5), rep(0.25, 4),
    c(NA_real_, 0.5, 0.5), c(Inf, 0, 0), c("0", "0", "1"),
    c(0 + 0i, 0 + 0i, 1 + 0i), matrix(c(0.2, 0.3, 0.5), 1),
    c(-0.01, 0.51, 0.5), c(1.01, -0.005, -0.005), c(0.2, 0.3, 0.6))
  for (weights in invalid_weights) {
    expect_error(pca_stage_mock(arguments,
      list(opt_params = weights, convergence = "convergence", loss = 0)),
      "explicit PCA.*weights")
  }
  for (weights in list(c(-5e-8, 0.4, 0.60000005), c(0.2, 0.3, 0.5000005))) {
    result <- pca_stage_mock(arguments,
      list(opt_params = weights, convergence = "convergence", loss = 0))
    expect_identical(result$linear_combinations, weights)
    expect_equal(result$predictions$adjusted, result$predictions$unadjusted +
      sum(weights * result$meta$omega_vec), tolerance = 1e-10)
  }
  for (weights in list(c(-2e-7, 0.4, 0.6000002), c(0.2, 0.3, 0.500002))) {
    expect_error(pca_stage_mock(arguments,
      list(opt_params = weights, convergence = "convergence", loss = 0)),
      "explicit PCA.*weights")
  }
})

test_that("NULL retains warning fallback and legacy returned-status semantics", {
  arguments <- pca_stage_input()
  expect_warning(result <- pca_stage_mock(arguments, error_text = "legacy solver error"),
    "\\[DBW\\] failed: legacy solver error")
  expect_equal(result$linear_combinations, rep(1 / 3, 3), tolerance = 1e-10)
  returned <- list(opt_params = c(0.2, 0.3, 0.5),
    convergence = "failed_convergence", loss = 1)
  result <- pca_stage_mock(arguments, returned)
  expect_identical(result$linear_combinations, returned$opt_params)
  expect_identical(result$meta$dbw_status, "failed_convergence")
  for (series in seq_along(arguments$covariates_series_list)) {
    arguments$covariates_series_list[[series]][, ] <- 1
  }
  arguments$dbw_indices <- NULL
  expect_warning(result <- do.call(SynthPrediction, arguments), "no informative")
  expect_equal(result$linear_combinations, rep(1 / 3, 3), tolerance = 1e-10)
})

test_that("unused PCA is ignored and the single-donor path remains legal", {
  arguments <- pca_stage_input()
  arguments$use_dbw <- FALSE
  for (request in list(99, "unused", numeric())) {
    arguments["princ_comp_input"] <- list(request)
    expect_equal(do.call(SynthPrediction, arguments)$linear_combinations,
      rep(1 / 3, 3), tolerance = 1e-10)
  }
  arguments <- pca_stage_input()
  arguments$Y_series_list <- arguments$Y_series_list[1:2]
  arguments$covariates_series_list <- arguments$covariates_series_list[1:2]
  arguments$shock_time_vec <- arguments$shock_time_vec[1:2]
  arguments$shock_length_vec <- arguments$shock_length_vec[1:2]
  arguments$princ_comp_input <- 1L
  result <- do.call(SynthPrediction, arguments)
  expect_identical(result$linear_combinations, 1)
  expect_identical(result$meta$dbw_status, "convergence")
  expect_equal(result$predictions$adjusted,
    result$predictions$unadjusted + result$meta$omega_vec, tolerance = 1e-10)
})

test_that("MultiPool PCA forwarding honors each existing DBW switch", {
  arguments <- pca_stage_input()[c("Y_series_list", "covariates_series_list",
    "shock_time_vec", "shock_length_vec", "k")]
  arguments$pools <- list(one = c("donor1", "donor2"), two = c("donor2", "donor3"))
  original_dbw <- getFromNamespace("dbw", "postshock")
  for (switches in list(c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE), c(FALSE, FALSE))) {
    arguments$base_use_dbw <- switches[1L]
    arguments$pool_use_dbw <- switches[2L]
    arguments$princ_comp_input <- if (any(switches)) 1L else "unused"
    direct <- suppressWarnings(do.call(SynthPredictionMultiPool, arguments))
    recorded <- list()
    routed <- local({
      testthat::local_mocked_bindings(dbw = function(...) {
        supplied <- list(...)
        recorded[[length(recorded) + 1L]] <<- supplied$princ_comp_count
        do.call(original_dbw, supplied)
      }, .package = "postshock")
      suppressWarnings(do.call(SynthPrediction, arguments))
    })
    expect_equal(routed, direct, tolerance = model_fit_tolerance)
    expect_length(recorded, as.integer(switches[1L]) + 2L * as.integer(switches[2L]))
    if (any(switches)) expect_true(all(vapply(recorded, identical, logical(1), 1L)))
    expect_equal(routed$predictions$adjusted_multipool,
      routed$predictions$unadjusted + sum(routed$meta$multi_pool$pool_effects),
      tolerance = 1e-10)
  }
  for (switches in list(c(FALSE, TRUE), c(TRUE, FALSE))) {
    arguments$base_use_dbw <- switches[1L]
    arguments$pool_use_dbw <- switches[2L]
    arguments$princ_comp_input <- 5L
    expect_error(do.call(SynthPrediction, arguments), "available.*directions")
    expect_error(do.call(SynthPredictionMultiPool, arguments), "available.*directions")
  }
})

