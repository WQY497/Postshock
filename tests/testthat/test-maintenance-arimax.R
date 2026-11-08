test_that("donor ARIMAX rejects nonfinite selected inputs before either fitting engine", {
  testthat::local_mocked_bindings(Arima = function(...) stop("ENGINE ENTERED"),
    auto.arima = function(...) stop("ENGINE ENTERED"), .package = "forecast")
  for (fixed in c(TRUE, FALSE)) {
    for (field in c("response", "first")) {
      for (position in c(1L, 37L, 71L, 73L)) {
        for (bad in list(NA_real_, NaN, Inf, -Inf)) {
          arguments <- maintenance_input()
          arguments$covariate_indices <- 1L
          if (!fixed) arguments$arima_order <- NULL
          if (field == "response") arguments$Y_series_list[[2L]][position] <- bad else
            arguments$covariates_series_list[[2L]][position, 1L] <- bad
          error <- tryCatch(do.call(SynthPrediction, arguments), error = identity)
          expect_s3_class(error, "error")
          expect_match(conditionMessage(error), paste0("donorA.*", field, ".*row ", position, " "))
          expect_match(conditionMessage(error), as.character(as.Date("2020-01-01") + position - 1L), fixed = TRUE)
          expect_match(conditionMessage(error), if (is.nan(bad)) "NaN" else if (is.na(bad)) "NA" else if (bad > 0) "\\+Inf" else "-Inf")
        }
      }
    }
  }
})

test_that("donor ARIMAX validates indexed alignment, window lengths and selected columns", {
  arguments <- maintenance_input()
  arguments$covariate_indices <- 1L
  original <- arguments
  arguments$covariates_series_list[[2L]] <- xts::xts(
    as.matrix(arguments$covariates_series_list[[2L]]),
    order.by = zoo::index(arguments$covariates_series_list[[2L]]) + 1L)
  expect_error(do.call(SynthPrediction, arguments), "donorA.*index mismatch.*row 1")
  arguments <- original
  arguments$covariates_series_list[[2L]] <- arguments$covariates_series_list[[2L]][1:70, ]
  expect_error(do.call(SynthPrediction, arguments), "donorA.*covariates.*73.*70")
  arguments <- original
  arguments$Y_series_list[[2L]] <- arguments$Y_series_list[[2L]][1:72]
  expect_error(do.call(SynthPrediction, arguments), "donorA.*response.*73.*72")
  arguments <- original
  arguments$covariate_indices <- 99L
  expect_error(do.call(SynthPrediction, arguments), "donorA.*covariate.*99")
  arguments <- original
  names(arguments$Y_series_list) <- NULL
  arguments$Y_series_list[[2L]][2L] <- NA_real_
  expect_error(do.call(SynthPrediction, arguments), "donor 1.*response.*row 2")
})

test_that("complete donor ARIMAX preserves raw rows, coefficients and unused data", {
  arguments <- maintenance_input()
  arguments$covariate_indices <- 1L
  reference <- suppressWarnings(do.call(SynthPrediction, arguments))
  for (donor in 2:4) {
    indicator <- c(rep(0, 70L), rep(1, 3L))
    design <- cbind(as.matrix(arguments$covariates_series_list[[donor]][1:73, 1L, drop = FALSE]), indicator)
    independent <- forecast::Arima(arguments$Y_series_list[[donor]][1:73],
      order = c(1, 0, 0), seasonal = FALSE, xreg = design)
    expect_equal(reference$meta$omega_vec[donor - 1L], unname(tail(stats::coef(independent), 1L)),
      tolerance = model_fit_tolerance)
  }
  arguments$covariates_series_list[[2L]][, 2L] <- NA_real_
  arguments$covariates_series_list[[2L]][74:90, 1L] <- Inf
  arguments$Y_series_list[[2L]][74:90] <- NA_real_
  unused <- suppressWarnings(do.call(SynthPrediction, arguments))
  expect_equal(unused, reference, tolerance = model_fit_tolerance)
  arguments$Y_series_list[[2L]] <- as.numeric(arguments$Y_series_list[[2L]])
  arguments$covariates_series_list[[2L]] <- as.matrix(arguments$covariates_series_list[[2L]])
  expect_equal(suppressWarnings(do.call(SynthPrediction, arguments)), reference,
    tolerance = model_fit_tolerance)
})

test_that("donor ARIMAX errors survive wrapper and direct MultiPool dispatch", {
  arguments <- maintenance_input()
  arguments$covariate_indices <- 1L
  arguments$use_dbw <- NULL
  arguments$pools <- list(all = c("donorA", "donorB", "donorC"))
  arguments$Y_series_list[[2L]][72L] <- NA_real_
  expect_error(do.call(SynthPrediction, arguments), "donorA.*response.*row 72")
  expect_error(do.call(SynthPredictionMultiPool, arguments), "donorA.*response.*row 72")
})

