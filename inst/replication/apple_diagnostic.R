check_apple_adjustment <- function(meta) {
  reference_effects <- c(0.425272469331, 0.218634577538, 1.990999985436)
  reference_weights <- c(0.583683226784, 0.0000000562775710497, 0.416316716939)
  reference_adjustment <- 1.077110996827
  effect_tolerance <- 1e-9
  weight_tolerance <- 1e-7
  aggregation_tolerance <- 1e-10
  effects <- meta$omega_vec
  weights <- meta$weights
  combined <- meta$combined_omega
  stopifnot(is.numeric(effects), !is.complex(effects), length(effects) == 3L,
            is.numeric(weights), !is.complex(weights), length(weights) == 3L,
            is.numeric(combined), !is.complex(combined), length(combined) == 1L,
            all(is.finite(c(effects, weights, combined))),
            all(abs(effects - reference_effects) < effect_tolerance),
            all(abs(weights - reference_weights) < weight_tolerance))
  recomputed <- sum(weights * effects)
  aggregation_error <- abs(combined - recomputed)
  reference_rounding <- abs(sum(reference_weights * reference_effects) - reference_adjustment)
  propagated_bound <- sum(abs(reference_effects)) * weight_tolerance +
    sum(abs(weights)) * effect_tolerance + reference_rounding + aggregation_tolerance
  reference_error <- abs(combined - reference_adjustment)
  stopifnot(aggregation_error < aggregation_tolerance,
            reference_error < propagated_bound)
  invisible(list(reference_error = reference_error, propagated_bound = propagated_bound,
                 aggregation_error = aggregation_error,
                 aggregation_tolerance = aggregation_tolerance))
}

apple_diagnostics <- function(verbose = FALSE) {
  suppressPackageStartupMessages({
    library(postshock)
    library(xts)
    library(zoo)
  })
  source(system.file("replication", "check_runtime.R", package = "postshock", mustWork = TRUE), local = TRUE)
  check_postshock_runtime(verbose = verbose)
  data("aapl_s1_input", package = "postshock")
  data("aapl_onestep_input", package = "postshock")
  data("aapl_control_input", package = "postshock")
  apple <- aapl_onestep_input
  market_names <- c("GSPC", "IRX")
  panels <- list(apple$Z_tgt_base,
                cbind(apple$Y_series_list[[1L]], apple$X_for_dbw[[1L]]),
                cbind(aapl_control_input$Y[[1L]], aapl_control_input$X[[1L]]))
  for (panel in panels) stopifnot(identical(colnames(panel), c("AAPL", market_names)))
  overlap_checks <- list()
  for (first in seq_len(length(panels) - 1L)) {
    for (second in seq.int(first + 1L, length(panels))) {
      overlap <- intersect(as.character(index(panels[[first]])), as.character(index(panels[[second]])))
      difference <- as.matrix(panels[[first]][match(overlap, as.character(index(panels[[first]]))), ]) -
        as.matrix(panels[[second]][match(overlap, as.character(index(panels[[second]]))), ])
      stopifnot(max(abs(difference)) == 0)
      overlap_checks[[length(overlap_checks) + 1L]] <- data.frame(
        first = first, second = second, rows = length(overlap), max_difference = max(abs(difference)))
    }
  }
  target_panel <- do.call(rbind, panels)
  target_panel <- target_panel[!duplicated(index(target_panel)), ]
  target_panel <- target_panel[order(index(target_panel)), ]
  stopifnot(!anyNA(target_panel), !anyDuplicated(index(target_panel)))
  take_window <- function(origin, count) {
    available <- target_panel[as.Date(index(target_panel)) <= as.Date(origin), ]
    stopifnot(NROW(available) >= count)
    tail(available, count)
  }
  window_s1 <- take_window("2020-11-10", 25L)
  window_s2 <- take_window("2020-11-11", 25L)
  window_s3 <- take_window("2020-11-09", 24L)
  stopifnot(identical(as.matrix(window_s1), as.matrix(apple$Z_tgt_base)),
            sum(as.Date(index(target_panel)) <= as.Date("2020-11-09")) == 24L)
  donor_response <- apple$Y_series_list[-1L]
  donor_market <- apple$X_for_dbw[-1L]
  effect_rows <- vapply(seq_along(donor_response), function(donor_number) {
    match(as.Date(apple$donor_effect_days[donor_number]), as.Date(index(donor_response[[donor_number]])))
  }, integer(1))
  stopifnot(identical(effect_rows, rep(13L, 3L)))
  boundaries <- effect_rows - 1L
  capture_run <- function(expression) {
    conditions <- character()
    value <- withCallingHandlers(expression, warning = function(condition) {
      conditions <<- c(conditions, conditionMessage(condition))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = unique(conditions))
  }
  point_control <- function(dates) {
    position <- match(as.Date(apple$control_shock_day), as.Date(dates))
    stopifnot(!is.na(position))
    list(time = position, length = 1L, shape = "point")
  }
  matching_diagnostic <- function(features, weights) {
    standardized <- scale(features, center = TRUE, scale = TRUE)
    donor_features <- standardized[-1L, , drop = FALSE]
    affine <- sweep(donor_features[-NROW(donor_features), , drop = FALSE],
                    2L, donor_features[NROW(donor_features), ], "-")
    c(rows = NROW(features), columns = NCOL(features), centered_rank = qr(standardized)$rank,
      donor_rank = qr(donor_features)$rank, affine_rank = qr(t(affine))$rank,
      weight_sum = sum(weights), weight_min = min(weights), weight_max = max(weights),
      matching_loss = sqrt(sum((standardized[1L, ] - as.numeric(weights %*% donor_features))^2)))
  }
  run_standard <- function(window, horizon, with_control = TRUE) {
    responses <- c(list(window[, "AAPL", drop = FALSE]), donor_response)
    markets <- c(list(window[, market_names, drop = FALSE]), donor_market)
    cutoff <- c(NROW(window), boundaries)
    controls <- NULL
    if (with_control) {
      controls <- vector("list", length(responses))
      controls[[1L]] <- list(point_control(index(window)))
    }
    arguments <- list(
      Y_series_list = responses, covariates_series_list = markets,
      shock_time_vec = cutoff, shock_length_vec = rep(1L, 4L), k = horizon,
      use_dbw = TRUE, dbw_indices = 1:2, covariate_indices = NULL,
      arima_order = NULL, seasonal = FALSE, user_ic_choice = "aicc",
      stepwise = TRUE, approximation = TRUE, control_shocks = controls,
      dbw_center = TRUE, dbw_scale = TRUE, days_before_shocktime_vec = 10L,
      penalty_lambda = 0, penalty_normchoice = "l1")
    if (identical(window, window_s1) && with_control) {
      fields <- c("Y_series_list", "covariates_series_list", "shock_time_vec",
                  "shock_length_vec", "control_shocks")
      frozen_arguments <- c(aapl_s1_input[fields], aapl_s1_input$settings)
      frozen_arguments$k <- horizon
      stopifnot(isTRUE(all.equal(arguments, frozen_arguments[names(arguments)], tolerance = 0)))
      arguments <- frozen_arguments
    }
    run <- capture_run(do.call(SynthPrediction, arguments))
    fit <- run$value
    control_matrix <- build_control_xreg(NROW(window), if (with_control) controls[[1L]] else NULL)
    model_run <- capture_run(forecast::auto.arima(
      responses[[1L]], xreg = control_matrix, ic = "aicc",
      seasonal = FALSE, stepwise = TRUE, approximation = TRUE))
    target_model <- model_run$value
    future_control <- if (is.null(control_matrix)) NULL else
      matrix(0, horizon, NCOL(control_matrix), dimnames = list(NULL, colnames(control_matrix)))
    direct <- tryCatch(
      if (is.null(future_control)) stats::predict(target_model, n.ahead = horizon)
      else stats::predict(target_model, n.ahead = horizon, newxreg = future_control),
      error = identity)
    fallback <- inherits(direct, "error")
    if (fallback) {
      verified <- if (is.null(future_control)) forecast::forecast(target_model, h = horizon) else
        forecast::forecast(target_model, h = horizon, xreg = future_control)
      verified_raw <- as.numeric(verified$mean)
    } else verified_raw <- as.numeric(direct$pred)
    stopifnot(max(abs(verified_raw - fit$predictions$unadjusted)) < 1e-9,
              max(abs(fit$predictions$adjusted - fit$predictions$unadjusted -
                        fit$meta$combined_omega)) < 1e-10,
              max(abs(fit$predictions$arithmetic_mean - fit$predictions$unadjusted -
                        mean(fit$meta$omega_vec))) < 1e-10,
              all(is.finite(unlist(fit$predictions))),
              abs(sum(fit$meta$weights) - 1) < 1e-6,
              all(fit$meta$weights >= -1e-7), all(fit$meta$weights <= 1 + 1e-7))
    features <- do.call(rbind, lapply(seq_along(markets), function(series_number)
      as.numeric(t(as.matrix(markets[[series_number]][cutoff[series_number] - 0:9, , drop = FALSE])))))
    list(fit = fit, arguments = arguments, window = window, target_model = target_model,
         control_matrix = control_matrix, future_control = future_control,
         controls = controls, warnings = unique(c(run$warnings, model_run$warnings)),
         fallback = fallback, fallback_reason = if (fallback) conditionMessage(direct) else "",
         matching = matching_diagnostic(features, fit$meta$weights))
  }
  standard <- list(S1 = run_standard(window_s1, 2L),
                   S1_no_control = run_standard(window_s1, 2L, FALSE),
                   S2 = run_standard(window_s2, 1L),
                   S3 = run_standard(window_s3, 2L))
  s1_first_step_check <- run_standard(window_s1, 1L)
  for (field in names(standard$S1$fit$predictions))
    stopifnot(abs(s1_first_step_check$fit$predictions[[field]][1L] -
                    standard$S1$fit$predictions[[field]][1L]) < 1e-10)
  for (result in standard) stopifnot(identical(result$fit$meta$omega_vec, standard$S1$fit$meta$omega_vec))
  stopifnot(identical(standard$S1$fit$meta$weights, standard$S1_no_control$fit$meta$weights))
  original_matching <- apple$X_for_dbw[[1L]][as.Date(index(apple$X_for_dbw[[1L]])) <= as.Date("2020-11-10"), ]
  stopifnot(identical(as.matrix(tail(original_matching, 10L)),
                      as.matrix(tail(window_s1[, market_names], 10L))))
  regression_data <- function(window) {
    values <- as.matrix(window)
    data.frame(response = values[-1L, "AAPL"], lag_AAPL = head(values[, "AAPL"], -1L),
               lag_GSPC = head(values[, "GSPC"], -1L), lag_IRX = head(values[, "IRX"], -1L))
  }
  fit_local <- function(window) {
    training <- regression_data(window)
    training$control <- as.numeric(build_control_xreg(
      NROW(training), point_control(index(window)[-1L]))[, 1L])
    lm(response ~ lag_AAPL + lag_GSPC + lag_IRX + control, data = training)
  }
  model_m <- fit_local(window_s1)
  future_m <- data.frame(lag_AAPL = as.numeric(tail(window_s1[, "AAPL"], 1L)),
                         lag_GSPC = as.numeric(tail(window_s1[, "GSPC"], 1L)),
                         lag_IRX = as.numeric(tail(window_s1[, "IRX"], 1L)), control = 0)
  raw_m <- as.numeric(predict(model_m, future_m))
  theta_m <- standard$S1$fit$meta$combined_omega
  equal_m <- mean(standard$S1$fit$meta$omega_vec)
  actual_on <- function(dates) as.numeric(target_panel[match(as.Date(dates), as.Date(index(target_panel))), "AAPL"])
  audit_values <- c(119.489998, 119.977803, 1.077111, 0.878302, 121.054914, 120.856105)
  computed_values <- c(actual_on("2020-11-11"), raw_m, theta_m, equal_m, raw_m + theta_m, raw_m + equal_m)
  stopifnot(all(abs(computed_values - audit_values) < 5e-7))
  donor_models_s4 <- lapply(seq_along(donor_response), function(donor_number) {
    window <- cbind(donor_response[[donor_number]], donor_market[[donor_number]])
    training <- regression_data(window)
    training$effect <- as.numeric(as.Date(index(window)[-1L]) == apple$donor_effect_days[donor_number])
    stopifnot(sum(training$effect) == 1L)
    lm(response ~ lag_AAPL + lag_GSPC + lag_IRX + effect, data = training)
  })
  effects_s4 <- vapply(donor_models_s4, function(model) unname(coef(model)["effect"]), numeric(1))
  fit_market_ar <- function(values, order) {
    count <- length(values)
    training <- data.frame(response = values[seq.int(order + 1L, count)])
    for (lag_number in seq_len(order))
      training[[paste0("lag", lag_number)]] <- values[seq.int(order + 1L, count) - lag_number]
    model <- lm(response ~ ., data = training)
    future <- as.data.frame(as.list(setNames(rev(tail(values, order)), paste0("lag", seq_len(order)))))
    list(model = model, prediction = as.numeric(predict(model, future)))
  }
  market_models_s4 <- list(GSPC = fit_market_ar(as.numeric(window_s3[, "GSPC"]), 3L),
                           IRX = fit_market_ar(as.numeric(window_s3[, "IRX"]), 1L))
  predicted_market <- vapply(market_models_s4, function(result) result$prediction, numeric(1))
  model_s4 <- fit_local(window_s3)
  first_input <- data.frame(lag_AAPL = as.numeric(tail(window_s3[, "AAPL"], 1L)),
                            lag_GSPC = as.numeric(tail(window_s3[, "GSPC"], 1L)),
                            lag_IRX = as.numeric(tail(window_s3[, "IRX"], 1L)), control = 0)
  raw_s4_first <- as.numeric(predict(model_s4, first_input))
  second_input <- data.frame(lag_AAPL = raw_s4_first, lag_GSPC = predicted_market["GSPC"],
                             lag_IRX = predicted_market["IRX"], control = 0)
  raw_s4 <- c(raw_s4_first, as.numeric(predict(model_s4, second_input)))
  matching_s4 <- c(list(matrix(predicted_market, nrow = 1L, dimnames = list(NULL, market_names))),
    lapply(seq_along(donor_market), function(donor_number)
      as.matrix(donor_market[[donor_number]][boundaries[donor_number], , drop = FALSE])))
  weights_run_s4 <- capture_run(dbw(
    X = matching_s4, dbw_indices = 1:2, shock_time_vec = rep(1L, 4L),
    center = TRUE, scale = TRUE, sum_to_1 = TRUE,
    bounded_below_by = 0, bounded_above_by = 1, princ_comp_count = NULL,
    normchoice = "l2", penalty_lambda = 0, penalty_normchoice = "l1",
    X_lookback_indices = NULL, verbose = FALSE))
  weights_s4 <- weights_run_s4$value$opt_params
  theta_s4 <- sum(weights_s4 * effects_s4)
  equal_s4 <- mean(effects_s4)
  stopifnot(all(is.finite(c(raw_s4, effects_s4, weights_s4, theta_s4))),
            abs(sum(weights_s4) - 1) < 1e-6, all(weights_s4 >= -1e-7), all(weights_s4 <= 1 + 1e-7))
  result_rows <- function(specification, origin, dates, raw, dbw_forecast, equal_forecast) {
    actual <- actual_on(dates)
    theta <- dbw_forecast - raw
    raw_error <- raw - actual
    dbw_error <- dbw_forecast - actual
    equal_error <- equal_forecast - actual
    stopifnot(max(abs(dbw_error - (raw_error + theta))) < 1e-10)
    data.frame(specification = specification, origin = origin, evaluation = dates,
               horizon = seq_along(raw), actual = actual, raw = raw, theta_DBW = theta,
               DBW = dbw_forecast, equal = equal_forecast,
               raw_error = raw_error, DBW_error = dbw_error, equal_error = equal_error,
               raw_AE = abs(raw_error), DBW_AE = abs(dbw_error), equal_AE = abs(equal_error),
               AE_change_DBW = abs(dbw_error) - abs(raw_error),
               AE_change_equal = abs(equal_error) - abs(raw_error))
  }
  rows <- list(M = result_rows("M", "2020-11-10", "2020-11-11", raw_m, raw_m + theta_m, raw_m + equal_m))
  evaluations <- list(S1 = c("2020-11-11", "2020-11-12"),
                      S1_no_control = c("2020-11-11", "2020-11-12"),
                      S2 = "2020-11-12", S3 = c("2020-11-10", "2020-11-11"))
  for (specification in names(standard)) {
    result <- standard[[specification]]
    predictions <- result$fit$predictions
    rows[[specification]] <- result_rows(specification, as.character(tail(index(result$window), 1L)),
      evaluations[[specification]], predictions$unadjusted, predictions$adjusted, predictions$arithmetic_mean)
  }
  rows$S4 <- result_rows("S4", "2020-11-09", c("2020-11-10", "2020-11-11"),
                         raw_s4, raw_s4 + c(0, theta_s4), raw_s4 + c(0, equal_s4))
  results <- do.call(rbind, rows)
  rownames(results) <- NULL

  main <- s1_first_step_check
  main_values <- c(actual_on("2020-11-11"),
                   main$fit$predictions$unadjusted,
                   main$fit$predictions$adjusted,
                   main$fit$predictions$arithmetic_mean)
  main_table <- data.frame(
    Forecast = c("Realized close", "Unadjusted baseline", "DBW-adjusted", "Equal-weight adjusted"),
    Value = main_values,
    "Absolute error" = c(NA_real_, abs(main_values[-1L] - main_values[1L])),
    check.names = FALSE)
  verified <- c(main_values, main_table[-1L, "Absolute error"])
  expected <- c(119.489998, 116.049774, 117.126885, 116.928076,
                3.440224, 2.363113, 2.561922)
  adjustment_verification <- check_apple_adjustment(main$fit$meta)
  stopifnot(all(abs(verified - expected) < 5e-7),
            as.character(apple$control_shock_day) == "2020-10-09",
            identical(which(main$control_matrix[, 1L] == 1), 3L),
            all(main$future_control == 0),
            identical(as.character(main$fit$meta$dbw_status), "convergence"))
  donor_dates <- do.call(rbind, lapply(seq_along(donor_response), function(donor_number) {
    dates <- as.Date(index(donor_response[[donor_number]]))
    data.frame(donor = donor_number, full_start = dates[1L], full_end = tail(dates, 1L),
               boundary = dates[boundaries[donor_number]], effect = dates[effect_rows[donor_number]],
               matching_start = dates[boundaries[donor_number] - 9L],
               inferred_Tstar = dates[effect_rows[donor_number] - 2L])
  }))
  training_m <- regression_data(window_s1)
  baseline_models <- list(
    "OLS AR(1)" = lm(response ~ lag_AAPL, data = training_m),
    "OLS with markets" = lm(response ~ lag_AAPL + lag_GSPC + lag_IRX, data = training_m))
  baseline_values <- c("S1: package" = main_values[2L],
    vapply(baseline_models, function(model) as.numeric(predict(model, future_m)), numeric(1)),
    "M: markets and control" = raw_m)
  baseline_comparison <- data.frame(
    Baseline = names(baseline_values), Raw = unname(baseline_values),
    DBW = unname(baseline_values) + theta_m,
    Raw_AE = abs(unname(baseline_values) - main_values[1L]),
    DBW_AE = abs(unname(baseline_values) + theta_m - main_values[1L]))
  weight_sets <- list("Current S1" = main$fit$meta$weights,
                      "Equal" = rep(1 / length(effects_s4), length(effects_s4)),
                      "Published rounded" = c(0.712, 0, 0.288))
  decomposition <- do.call(rbind, lapply(names(weight_sets), function(label) {
    weights <- weight_sets[[label]]
    adjustments <- c(sum(weights * main$fit$meta$omega_vec), sum(weights * effects_s4))
    data.frame(Weights = label, Donor_model = c("D1: package", "D2: full-window OLS"),
               Adjustment = adjustments, M_adjusted = raw_m + adjustments,
               M_AE = abs(raw_m + adjustments - main_values[1L]))
  }))
  expected_diagnostics <- c(
    119.977803, 121.054914, 120.856105,
    116.103617, 117.180728, 116.981919,
    116.030901, 117.108012, 116.909203,
    118.614426, 119.652038, 119.492728,
    116.297291, 117.390794, 117.175593,
    116.281986, 117.375489, 117.160289,
    121.515079, 118.878242, 119.631361, 120.150985)
  current_diagnostics <- c(
    unlist(rows$M[c("raw", "DBW", "equal")], use.names = FALSE),
    unlist(rows$S1[2L, c("raw", "DBW", "equal")], use.names = FALSE),
    unlist(rows$S1_no_control[1L, c("raw", "DBW", "equal")], use.names = FALSE),
    unlist(rows$S2[c("raw", "DBW", "equal")], use.names = FALSE),
    unlist(rows$S3[1L, c("raw", "DBW", "equal")], use.names = FALSE),
    unlist(rows$S3[2L, c("raw", "DBW", "equal")], use.names = FALSE),
    raw_s4_first, unlist(rows$S4[2L, c("raw", "DBW", "equal")], use.names = FALSE))
  stopifnot(all(abs(current_diagnostics - expected_diagnostics) < 5e-7))

  if (verbose) {
    cat("\nFROZEN-DATA OVERLAP AND NUMERICAL CHECKS: PASS\n\nOVERLAP CHECKS\n")
    cat("S1 k=1 versus first step of k=2: PASS\n")
    print(do.call(rbind, overlap_checks))
    cat("\nNUMERICAL RESULTS\n")
    print(results)
    cat("\nADJUSTMENT CHECK: PROPAGATED COMPONENT ERROR BOUND\n")
    print(adjustment_verification)
    cat("\nSTANDARD MODEL DIAGNOSTICS\n")
    for (specification in names(standard)) {
      result <- standard[[specification]]
      cat("\n", specification, "\n")
      print(list(window = range(index(result$window)), raw_rows = NROW(result$window),
                 fitting_nobs = stats::nobs(result$target_model),
                 target_order = forecast::arimaorder(result$target_model),
                 target_coefficients = coef(result$target_model),
                 optimizer_code = result$target_model$code,
                 control_rows = if (is.null(result$control_matrix)) integer() else which(result$control_matrix[, 1L] == 1),
                 future_control = result$future_control,
                 target_matching_dates = as.character(index(tail(result$window, 10L))),
                 donor_orders = result$fit$meta$arima_order[-1L],
                 donor_effects = result$fit$meta$omega_vec, weights = result$fit$meta$weights,
                 DBW_status = result$fit$meta$dbw_status, theta = result$fit$meta$combined_omega,
                 equal_theta = mean(result$fit$meta$omega_vec),
                 matching_diagnostic = result$matching, fallback = result$fallback,
                 fallback_reason = result$fallback_reason, warnings = result$warnings))
    }
    cat("\nDONOR DATES\n")
    for (donor_number in seq_along(donor_response)) {
      dates <- as.Date(index(donor_response[[donor_number]]))
      print(data.frame(donor = donor_number, full_start = dates[1L], full_end = tail(dates, 1L),
        package_fit_start = dates[1L], package_fit_end = dates[effect_rows[donor_number]],
        package_fit_rows = effect_rows[donor_number],
        boundary = dates[boundaries[donor_number]], effect = dates[effect_rows[donor_number]],
        package_matching_start = dates[boundaries[donor_number] - 9L],
        package_matching_end = dates[boundaries[donor_number]],
        inferred_paper_Tstar = dates[effect_rows[donor_number] - 2L]))
    }
    cat("\nM AND S4 DIAGNOSTICS\n")
    print(list(M = list(coefficients = coef(model_m), nobs = nobs(model_m), rank = model_m$rank),
      S4 = list(target_coefficients = coef(model_s4), target_nobs = nobs(model_s4), target_rank = model_s4$rank,
        market_coefficients = lapply(market_models_s4, function(result) coef(result$model)),
        market_nobs = vapply(market_models_s4, function(result) nobs(result$model), numeric(1)),
        market_ranks = vapply(market_models_s4, function(result) result$model$rank, numeric(1)),
        forecast_market_date = "2020-11-10", predicted_market = predicted_market,
        first_input = first_input, second_input = second_input,
        donor_nobs = vapply(donor_models_s4, nobs, numeric(1)),
        donor_ranks = vapply(donor_models_s4, function(model) model$rank, numeric(1)),
        donor_coefficients = lapply(donor_models_s4, coef), donor_effects = effects_s4,
        weights = weights_s4, theta = theta_s4, equal_theta = equal_s4,
        DBW_status = weights_run_s4$value$convergence, DBW_loss = weights_run_s4$value$loss,
        matching_features = do.call(rbind, matching_s4),
        matching_diagnostic = matching_diagnostic(do.call(rbind, matching_s4), weights_s4),
        warnings = weights_run_s4$warnings)))
    cat("\nDIRECT COMPARISONS\n")
    print(c(S1_raw_minus_M = rows$S1$raw[1L] - rows$M$raw,
            S1_DBW_AE_minus_M_raw_AE = rows$S1$DBW_AE[1L] - rows$M$raw_AE,
            control_change_h1 = rows$S1$raw[1L] - rows$S1_no_control$raw[1L],
            control_change_h2 = rows$S1$raw[2L] - rows$S1_no_control$raw[2L]))
    cat("\nSESSION INFORMATION\n")
    print(sessionInfo())
    cat("\nR HOME\n"); print(R.home())
    cat("\nLIBRARY PATHS\n"); print(.libPaths())
    cat("\nALL ASSERTIONS PASSED\n")

  }

  list(main = main, main_table = main_table, results = results, rows = rows,
       standard = standard, target_panel = target_panel, apple = apple,
       donor_dates = donor_dates, model_m = model_m, model_s4 = model_s4,
       donor_models_s4 = donor_models_s4, effects_s4 = effects_s4,
       weights_s4 = weights_s4, theta_s4 = theta_s4, equal_s4 = equal_s4,
       predicted_market = predicted_market, market_models_s4 = market_models_s4,
       weights_run_s4 = weights_run_s4, baseline_comparison = baseline_comparison,
       decomposition = decomposition, overlap_checks = overlap_checks,
       verification = data.frame(Audit = expected, Executable = verified,
                                 Status = rep("PASS", length(verified))),
       adjustment_verification = adjustment_verification,
       session = sessionInfo())
}

if (sys.nframe() == 0L) {
  options(digits = 12, width = 120)
  apple_results <- apple_diagnostics(verbose = TRUE)
  print(apple_results$verification)
}
