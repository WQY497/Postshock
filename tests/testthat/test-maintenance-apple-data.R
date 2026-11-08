test_that("Apple S1 data reconstructs exactly from unchanged frozen sources", {
  input <- new.env()
  utils::data("aapl_s1_input", "aapl_onestep_input", "aapl_control_input", package = "postshock", envir = input)
  builder <- new.env()
  sys.source(system.file("replication", "apple_data_builder.R", package = "postshock", mustWork = TRUE), builder)
  rebuilt <- builder$build_aapl_s1_input(input$aapl_onestep_input, input$aapl_control_input)
  expect_identical(rebuilt, input$aapl_s1_input)
  expect_equal(rebuilt$shock_time_vec, c(25, 12, 12, 12))
  expect_equal(rebuilt$dates$fitted_rows, c(25, 13, 13, 13))
  expect_equal(as.character(rebuilt$dates$effect), c("2020-11-11", "2017-06-07", "2018-10-31", "2019-06-04"))
  expect_true(all(as.Date(zoo::index(rebuilt$Y_series_list[[1L]])) < rebuilt$evaluation$date))
  expect_true(all(as.Date(zoo::index(rebuilt$covariates_series_list[[1L]])) < rebuilt$evaluation$date))
  expect_equal(which(build_control_xreg(25L, rebuilt$control_shocks[[1L]])[, 1L] == 1), 3L)
  expect_false(any(c("predictions", "forecast", "fit") %in% names(rebuilt)))
})

test_that("installed shared Apple analysis preserves article checks and independent evaluation", {
  input <- new.env()
  utils::data("aapl_s1_input", package = "postshock", envir = input)
  frozen <- input$aapl_s1_input
  fields <- c("Y_series_list", "covariates_series_list", "shock_time_vec", "shock_length_vec", "control_shocks")
  arguments <- c(frozen[fields], frozen$settings)
  direct <- suppressWarnings(do.call(SynthPrediction, arguments))
  shared <- new.env()
  sys.source(system.file("replication", "apple_diagnostic.R", package = "postshock", mustWork = TRUE), shared)
  diagnostic <- shared$apple_diagnostics()
  expect_equal(direct$predictions, diagnostic$main$fit$predictions, tolerance = model_fit_tolerance)
  expect_equal(direct$meta$weights, diagnostic$main$fit$meta$weights, tolerance = model_fit_tolerance)
  expect_equal(frozen$evaluation$realized, diagnostic$main_table$Value[1L], tolerance = 0)
  frozen$evaluation$realized <- -99999
  expect_identical(c(frozen[fields], frozen$settings), arguments)
  expect_true(all(diagnostic$main$future_control == 0))
})

