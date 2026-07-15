# These tests isolate SynthVolForecast() from numerical optimizers.
# auto_garchx(), dbw(), and garchx() are tested separately; this file verifies
# input validation, design construction, argument forwarding, aggregation, and
# the public return contract of SynthVolForecast().

# Test helpers -------------------------------------------------------------

new_svf_mock_fit <- function(
    pred = 0.5,
    omega = NULL,
    omega_se = 0.01,
    fitted_values = numeric(),
    missing_coefficient = FALSE
) {
  structure(
    list(
      pred = as.numeric(pred),
      omega = omega,
      omega_se = omega_se,
      fitted_values = as.numeric(fitted_values),
      missing_coefficient = isTRUE(missing_coefficient)
    ),
    class = "svf_mock_fit"
  )
}

local_svf_mocks <- function(
    omegas = numeric(),
    weights = numeric(),
    target_pred = 0.5,
    target_order = c(1L, 1L, 0L),
    donor_orders = NULL,
    omega_se = rep(0.01, length(omegas)),
    dbw_convergence = 0L,
    missing_coefficient = FALSE,
    capture = NULL,
    frame = parent.frame()
) {
  if (is.null(donor_orders)) {
    donor_orders <- rep(list(c(1L, 1L, 0L)), length(omegas))
  }
  
  if (length(omega_se) == 1L && length(omegas) > 1L) {
    omega_se <- rep(omega_se, length(omegas))
  }
  
  if (!is.null(capture)) {
    capture$auto_calls <- list()
    capture$garchx_calls <- list()
    capture$dbw_calls <- list()
    capture$predict_n_ahead <- integer()
  }
  
  donor_index <- 0L
  
  next_mock_fit <- function(y, xreg) {
    if (is.null(xreg)) {
      return(new_svf_mock_fit(
        pred = target_pred,
        fitted_values = rep(mean(target_pred), length(y))
      ))
    }
    
    donor_index <<- donor_index + 1L
    
    if (donor_index > length(omegas)) {
      stop("The test requested more donor fits than the mock provides.")
    }
    
    new_svf_mock_fit(
      pred = 0.25,
      omega = omegas[[donor_index]],
      omega_se = omega_se[[donor_index]],
      fitted_values = rep(0.25, length(y)),
      missing_coefficient = missing_coefficient
    )
  }
  
  auto_mock <- function(y, xreg = NULL, max.p = 3, max.q = 3,
                        vcov.type = "robust", ...) {
    if (!is.null(capture)) {
      capture$auto_calls[[length(capture$auto_calls) + 1L]] <- list(
        y = y,
        xreg = xreg,
        max.p = max.p,
        max.q = max.q,
        vcov.type = vcov.type,
        dots = list(...)
      )
    }
    
    fit <- next_mock_fit(y, xreg)
    
    if (is.null(xreg)) {
      order <- target_order
    } else {
      order <- donor_orders[[donor_index]]
    }
    
    list(
      q = order[[1L]],
      p = order[[2L]],
      o = order[[3L]],
      fit = fit
    )
  }
  
  garchx_mock <- function(y, order = c(1L, 1L), xreg = NULL,
                          vcov.type = "robust", backcast.values = NULL, ...) {
    if (!is.null(capture)) {
      capture$garchx_calls[[length(capture$garchx_calls) + 1L]] <- list(
        y = y,
        order = order,
        xreg = xreg,
        vcov.type = vcov.type,
        backcast.values = backcast.values,
        dots = list(...)
      )
    }
    
    next_mock_fit(y, xreg)
  }
  
  dbw_mock <- function(...) {
    args <- list(...)
    
    if (!is.null(capture)) {
      capture$dbw_calls[[length(capture$dbw_calls) + 1L]] <- args
    }
    
    list(
      opt_params = weights,
      convergence = dbw_convergence
    )
  }
  
  predict_mock <- function(object, n.ahead = 1L, ...) {
    if (!is.null(capture)) {
      capture$predict_n_ahead <- c(
        capture$predict_n_ahead,
        as.integer(n.ahead)
      )
    }
    
    list(pred = rep_len(object$pred, as.integer(n.ahead)))
  }
  
  fitted_mock <- function(object, ...) {
    object$fitted_values
  }
  
  coeftest_mock <- function(x, ...) {
    coefficient_name <- if (isTRUE(x$missing_coefficient)) {
      "unrelated_coefficient"
    } else {
      "post_shock_indicator"
    }
    
    matrix(
      c(x$omega, x$omega_se),
      nrow = 1L,
      dimnames = list(
        coefficient_name,
        c("Estimate", "Std. Error")
      )
    )
  }
  
  testthat::local_mocked_bindings(
    auto_garchx = auto_mock,
    dbw = dbw_mock,
    garchx = garchx_mock,
    predict = predict_mock,
    fitted = fitted_mock,
    coeftest = coeftest_mock,
    .env = frame
  )
  
  invisible(capture)
}

make_svf_args <- function(n_donors = 2L, n = 8L) {
  n_donors <- as.integer(n_donors)
  n <- as.integer(n)
  n_total <- n_donors + 1L
  time_index <- seq_len(n)
  
  target_y <- 0.05 * time_index
  donor_y <- lapply(seq_len(n_donors), function(i) {
    0.10 * time_index + i / 100
  })
  
  target_x <- data.frame(
    x1 = time_index,
    x2 = time_index^2
  )
  
  donor_x <- lapply(seq_len(n_donors), function(i) {
    data.frame(
      x1 = time_index + i,
      x2 = (time_index + i)^2
    )
  })
  
  list(
    Y_series_list = c(list(target_y), donor_y),
    covariates_series_list = donor_x,
    target_covariates = target_x,
    shock_time_vec = rep(4L, n_total),
    shock_length_vec = rep(2L, n_total),
    k = 2L,
    plots = FALSE
  )
}

call_svf <- function(args, ...) {
  overrides <- list(...)
  
  if (length(overrides)) {
    args[names(overrides)] <- overrides
  }
  
  do.call(SynthVolForecast, args)
}


# Core output contract -----------------------------------------------------

test_that("SynthVolForecast combines donor effects using DBW weights", {
  args <- make_svf_args(n_donors = 2L)
  
  local_svf_mocks(
    omegas = c(0.20, -0.10),
    weights = c(0.75, 0.25),
    target_pred = c(0.40, 0.50),
    target_order = c(2L, 1L, 0L),
    donor_orders = list(c(1L, 1L, 0L), c(2L, 1L, 1L)),
    omega_se = c(0.03, 0.04),
    dbw_convergence = 0L
  )
  
  out <- call_svf(args)
  
  expect_identical(
    names(out),
    c("linear_combinations", "predictions", "meta")
  )
  expect_equal(out$linear_combinations, c(0.75, 0.25))
  expect_identical(
    names(out$predictions),
    c("unadjusted", "adjusted", "arithmetic_mean")
  )
  
  expect_equal(out$predictions$unadjusted, c(0.40, 0.50))
  expect_equal(out$predictions$adjusted, c(0.525, 0.625))
  expect_equal(out$predictions$arithmetic_mean, c(0.45, 0.55))
  
  expect_equal(out$meta$n_donors, 2L)
  expect_equal(out$meta$weights, c(0.75, 0.25))
  expect_equal(out$meta$omega_vec, c(0.20, -0.10))
  expect_equal(out$meta$omega_se, c(0.03, 0.04))
  expect_equal(out$meta$combined_omega, 0.125)
  expect_equal(out$meta$dbw_status, 0L)
  expect_equal(out$meta$donor_orders[[1L]], c(1L, 1L, 0L))
  expect_equal(out$meta$donor_orders[[2L]], c(2L, 1L, 1L))
  expect_equal(out$meta$target_order, c(2L, 1L, 0L))
  
  expect_false("donor_fits" %in% names(out))
  expect_false("target_fit" %in% names(out))
})


test_that("SynthVolForecast returns fitted objects only when requested", {
  args <- make_svf_args(n_donors = 2L)
  
  local_svf_mocks(
    omegas = c(0.10, 0.20),
    weights = c(0.40, 0.60),
    target_pred = c(0.30, 0.35)
  )
  
  out <- call_svf(args, return_fits = TRUE)
  
  expect_length(out$donor_fits, 2L)
  expect_s3_class(out$donor_fits[[1L]], "svf_mock_fit")
  expect_s3_class(out$donor_fits[[2L]], "svf_mock_fit")
  expect_s3_class(out$target_fit, "svf_mock_fit")
})


test_that("SynthVolForecast no-donor path returns the unadjusted forecast", {
  args <- make_svf_args(n_donors = 0L)
  capture <- new.env(parent = emptyenv())
  
  local_svf_mocks(
    target_pred = c(-1, 0, 0.30),
    target_order = c(1L, 2L, 0L),
    capture = capture
  )
  
  out <- call_svf(args, k = 3L)
  
  expected <- c(1e-8, 1e-8, 0.30)
  
  expect_equal(out$linear_combinations, numeric(0))
  expect_equal(out$predictions$unadjusted, expected)
  expect_equal(out$predictions$adjusted, expected)
  expect_equal(out$predictions$arithmetic_mean, expected)
  
  expect_equal(out$meta$n_donors, 0L)
  expect_equal(out$meta$weights, numeric(0))
  expect_equal(out$meta$omega_vec, numeric(0))
  expect_equal(out$meta$omega_se, numeric(0))
  expect_equal(out$meta$combined_omega, 0)
  expect_equal(out$meta$target_order, c(1L, 2L, 0L))
  expect_s3_class(out$target_fit, "svf_mock_fit")
  
  expect_length(capture$auto_calls, 1L)
  expect_length(capture$dbw_calls, 0L)
  expect_equal(capture$predict_n_ahead, 3L)
})


test_that("variance forecasts are bounded below by the numerical floor", {
  args <- make_svf_args(n_donors = 1L)
  
  local_svf_mocks(
    omegas = -2,
    weights = 1,
    target_pred = c(-1, 1e-12)
  )
  
  out <- call_svf(args)
  
  expect_equal(out$predictions$unadjusted, rep(1e-8, 2L))
  expect_equal(out$predictions$adjusted, rep(1e-8, 2L))
  expect_equal(out$predictions$arithmetic_mean, rep(1e-8, 2L))
  expect_true(all(out$predictions$unadjusted > 0))
  expect_true(all(out$predictions$adjusted > 0))
  expect_true(all(out$predictions$arithmetic_mean > 0))
})


# Model fitting and argument forwarding -----------------------------------

test_that("automatic GARCH selection receives the correct data and controls", {
  args <- make_svf_args(n_donors = 2L)
  capture <- new.env(parent = emptyenv())
  
  local_svf_mocks(
    omegas = c(0.10, 0.20),
    weights = c(0.50, 0.50),
    target_pred = c(0.30, 0.40),
    capture = capture
  )
  
  out <- call_svf(
    args,
    max.p = 2L,
    max.q = 1L,
    vcov.type = "ROBUST"
  )
  
  expect_length(capture$auto_calls, 3L)
  expect_length(capture$garchx_calls, 0L)
  
  expect_equal(capture$auto_calls[[1L]]$y, args$Y_series_list[[2L]])
  expect_equal(capture$auto_calls[[2L]]$y, args$Y_series_list[[3L]])
  expect_equal(
    capture$auto_calls[[3L]]$y,
    args$Y_series_list[[1L]][seq_len(args$shock_time_vec[[1L]])]
  )
  expect_null(capture$auto_calls[[3L]]$xreg)
  
  for (one_call in capture$auto_calls) {
    expect_equal(one_call$max.p, 2L)
    expect_equal(one_call$max.q, 1L)
    expect_equal(one_call$vcov.type, "robust")
  }
  
  expect_equal(out$meta$n_donors, 2L)
})


test_that("fixed GARCH order bypasses auto_garchx and forwards fit arguments", {
  args <- make_svf_args(n_donors = 2L)
  capture <- new.env(parent = emptyenv())
  fixed_order <- c(2L, 1L, 1L)
  
  local_svf_mocks(
    omegas = c(0.10, 0.20),
    weights = c(0.25, 0.75),
    target_pred = c(0.30, 0.40),
    capture = capture
  )
  
  out <- call_svf(
    args,
    garch_order = fixed_order,
    backcast.initial = 0.75,
    vcov.type = "ordinary"
  )
  
  expect_length(capture$auto_calls, 0L)
  expect_length(capture$garchx_calls, 3L)
  
  for (one_call in capture$garchx_calls) {
    expect_equal(one_call$order, fixed_order)
    expect_equal(one_call$vcov.type, "ordinary")
    expect_equal(one_call$backcast.values, 0.75)
  }
  
  expect_equal(capture$garchx_calls[[1L]]$y, args$Y_series_list[[2L]])
  expect_equal(capture$garchx_calls[[2L]]$y, args$Y_series_list[[3L]])
  expect_equal(
    capture$garchx_calls[[3L]]$y,
    args$Y_series_list[[1L]][seq_len(args$shock_time_vec[[1L]])]
  )
  expect_null(capture$garchx_calls[[3L]]$xreg)
  
  expect_equal(out$meta$target_order, fixed_order)
  expect_equal(out$meta$donor_orders, rep(list(fixed_order), 2L))
})


test_that("DBW receives user controls and fixed simplex constraints", {
  args <- make_svf_args(n_donors = 2L)
  capture <- new.env(parent = emptyenv())
  
  local_svf_mocks(
    omegas = c(0.10, 0.20),
    weights = c(0.30, 0.70),
    target_pred = c(0.30, 0.40),
    capture = capture
  )
  
  call_svf(
    args,
    dbw_scale = FALSE,
    dbw_center = FALSE,
    dbw_indices = 2L,
    princ_comp_input = 1L,
    penalty_lambda = 0.25,
    penalty_normchoice = "l2"
  )
  
  expect_length(capture$dbw_calls, 1L)
  dbw_call <- capture$dbw_calls[[1L]]
  
  expect_equal(dbw_call$dbw_indices, 2L)
  expect_equal(dbw_call$shock_time_vec, args$shock_time_vec)
  expect_false(dbw_call$scale)
  expect_false(dbw_call$center)
  expect_true(dbw_call$sum_to_1)
  expect_equal(dbw_call$bounded_below_by, 0)
  expect_equal(dbw_call$bounded_above_by, 1)
  expect_equal(dbw_call$princ_comp_count, 1L)
  expect_equal(dbw_call$normchoice, "l2")
  expect_equal(dbw_call$penalty_normchoice, "l2")
  expect_equal(dbw_call$penalty_lambda, 0.25)
})


test_that("numeric and matrix covariates are coerced before DBW", {
  args <- make_svf_args(n_donors = 1L)
  capture <- new.env(parent = emptyenv())
  
  args$target_covariates <- seq_len(8L)
  args$covariates_series_list <- list(
    matrix(seq_len(8L), ncol = 1L, dimnames = list(NULL, "z"))
  )
  
  local_svf_mocks(
    omegas = 0.10,
    weights = 1,
    target_pred = c(0.30, 0.40),
    capture = capture
  )
  
  call_svf(args, dbw_indices = NULL)
  
  dbw_call <- capture$dbw_calls[[1L]]
  
  expect_equal(dbw_call$dbw_indices, 1L)
  expect_length(dbw_call$X, 2L)
  expect_s3_class(dbw_call$X[[1L]], "data.frame")
  expect_s3_class(dbw_call$X[[2L]], "data.frame")
  expect_identical(names(dbw_call$X[[1L]]), "V1")
  expect_identical(names(dbw_call$X[[2L]]), "z")
})


# Donor design construction ------------------------------------------------

test_that("selected donor covariates and shock indicator are constructed correctly", {
  args <- make_svf_args(n_donors = 1L)
  capture <- new.env(parent = emptyenv())
  
  donor_y <- 101:108
  donor_x <- data.frame(
    ignored = 11:18,
    selected = c(21, 22, NA, 24, 25, 26, 27, 28)
  )
  
  args$Y_series_list[[2L]] <- donor_y
  args$covariates_series_list[[1L]] <- donor_x
  args$shock_time_vec <- c(4L, 4L)
  args$shock_length_vec <- c(2L, 2L)
  
  local_svf_mocks(
    omegas = 0.15,
    weights = 1,
    target_pred = c(0.30, 0.40),
    capture = capture
  )
  
  call_svf(args, covariate_indices = 2L)
  
  donor_call <- capture$auto_calls[[1L]]
  
  expect_equal(donor_call$y, donor_y[-3L])
  expect_equal(
    colnames(donor_call$xreg),
    c("selected", "post_shock_indicator")
  )
  expect_equal(unname(donor_call$xreg[, "selected"]), donor_x$selected[-3L])
  expect_equal(
    unname(donor_call$xreg[, "post_shock_indicator"]),
    c(0, 0, 0, 1, 1, 0, 0)
  )
})


test_that("NULL covariate_indices uses only the post-shock indicator", {
  args <- make_svf_args(n_donors = 1L)
  capture <- new.env(parent = emptyenv())
  
  args$covariates_series_list[[1L]] <- data.frame(
    x1 = rep(NA_real_, 8L),
    x2 = rep(NA_real_, 8L)
  )
  
  local_svf_mocks(
    omegas = 0.15,
    weights = 1,
    target_pred = c(0.30, 0.40),
    capture = capture
  )
  
  call_svf(args, covariate_indices = NULL)
  
  donor_call <- capture$auto_calls[[1L]]
  
  expect_equal(donor_call$y, args$Y_series_list[[2L]])
  expect_equal(dim(donor_call$xreg), c(8L, 1L))
  expect_identical(colnames(donor_call$xreg), "post_shock_indicator")
  expect_equal(
    as.numeric(donor_call$xreg[, 1L]),
    c(0, 0, 0, 0, 1, 1, 0, 0)
  )
})


test_that("forecast horizon is coerced to an integer and bounded below by one", {
  args <- make_svf_args(n_donors = 0L)
  capture <- new.env(parent = emptyenv())
  
  local_svf_mocks(
    target_pred = c(0.20, 0.30, 0.40),
    capture = capture
  )
  
  out_zero <- call_svf(args, k = 0)
  out_fractional <- call_svf(args, k = 2.9)
  
  expect_length(out_zero$predictions$unadjusted, 1L)
  expect_length(out_fractional$predictions$unadjusted, 2L)
  expect_equal(capture$predict_n_ahead, c(1L, 2L))
})


# Validation and failure modes --------------------------------------------

test_that("top-level input lengths are validated", {
  empty_args <- make_svf_args(n_donors = 0L)
  empty_args$Y_series_list <- list()
  empty_args$shock_time_vec <- integer()
  empty_args$shock_length_vec <- integer()
  
  expect_error(
    call_svf(empty_args),
    "must contain at least the target series"
  )
  
  args <- make_svf_args(n_donors = 1L)
  
  expect_error(
    call_svf(args, shock_time_vec = 4L),
    "shock_time_vec"
  )
  
  expect_error(
    call_svf(args, shock_length_vec = 2L),
    "shock_length_vec"
  )
  
  expect_error(
    call_svf(args, covariates_series_list = list()),
    "covariates_series_list"
  )
})


test_that("invalid target and donor shock indices are rejected", {
  no_donor_args <- make_svf_args(n_donors = 0L)
  
  expect_error(
    call_svf(no_donor_args, shock_time_vec = 0L),
    "invalid target shock index"
  )
  
  expect_error(
    call_svf(no_donor_args, shock_time_vec = 9L),
    "invalid target shock index"
  )
  
  donor_args <- make_svf_args(n_donors = 1L)
  
  expect_error(
    call_svf(donor_args, shock_time_vec = c(4L, 0L)),
    "invalid shock index for donor 1"
  )
  
  expect_error(
    call_svf(donor_args, shock_time_vec = c(4L, 9L)),
    "invalid shock index for donor 1"
  )
})


test_that("unsupported covariate types are rejected", {
  args <- make_svf_args(n_donors = 1L)
  
  expect_error(
    call_svf(args, target_covariates = "not numeric"),
    "covariates must be numeric, matrix, or data.frame"
  )
  
  expect_error(
    call_svf(args, covariates_series_list = list("not numeric")),
    "covariates must be numeric, matrix, or data.frame"
  )
})


test_that("scalar tuning arguments are validated", {
  args <- make_svf_args(n_donors = 0L)
  
  bad_pca_values <- list(0, -1, Inf, c(1, 2), "one")
  
  for (bad_value in bad_pca_values) {
    expect_error(
      call_svf(args, princ_comp_input = bad_value),
      "princ_comp_input"
    )
  }
  
  bad_backcast_values <- list(0, -1, Inf, NA_real_, c(1, 2), "one")
  
  for (bad_value in bad_backcast_values) {
    expect_error(call_svf(args, backcast.initial = bad_value))
  }
  
  expect_error(
    call_svf(args, vcov.type = "unsupported"),
    "should be one of"
  )
})


test_that("DBW weight length must match the number of donor effects", {
  args <- make_svf_args(n_donors = 2L)
  
  local_svf_mocks(
    omegas = c(0.10, 0.20),
    weights = 1,
    target_pred = c(0.30, 0.40)
  )
  
  expect_error(
    call_svf(args),
    "DBW weight length does not match number of donor effects"
  )
})


test_that("a missing post-shock coefficient produces an informative error", {
  args <- make_svf_args(n_donors = 1L)
  
  local_svf_mocks(
    omegas = 0.10,
    weights = 1,
    target_pred = c(0.30, 0.40),
    missing_coefficient = TRUE
  )
  
  expect_error(
    call_svf(args),
    "post_shock_indicator coefficient not found for donor 1"
  )
})