maintenance_garch_fit <- function(quality = 1, status = 0, coefficients = c(intercept = 0.2, arch1 = 0.1, garch1 = 0.5), tag = "fit") {
  structure(list(convergence = status, par = coefficients, y.coredata = rep(0, 30),
    y.n = 30, maxpqrpluss1 = 1, fitted = zoo::zoo(rep(exp(-quality), 30)),
    objective.fun = 1, tag = tag, vcov = matrix(NA_real_, 3, 3)), class = "garchx")
}

test_that("GARCH grid ignores every unqualified lower-BIC candidate", {
  invalid <- list(
    maintenance_garch_fit(20, 1), maintenance_garch_fit(20, NULL),
    maintenance_garch_fit(20, NA_real_), maintenance_garch_fit(20, NaN),
    maintenance_garch_fit(20, Inf), maintenance_garch_fit(20, -Inf),
    maintenance_garch_fit(20, "0"), maintenance_garch_fit(20, FALSE),
    maintenance_garch_fit(20, c(0, 0)), maintenance_garch_fit(20, 0 + 0i),
    maintenance_garch_fit(20, 0, numeric()), maintenance_garch_fit(20, 0, c(intercept = Inf)),
    maintenance_garch_fit(20, 0, c(intercept = NaN)), maintenance_garch_fit(20, 0, "bad"),
    maintenance_garch_fit(Inf, 0), maintenance_garch_fit(NA_real_, 0),
    structure("fit error", class = "try-error"))
  for (bad in invalid) {
    count <- 0L
    testthat::local_mocked_bindings(garchx = function(...) {
      count <<- count + 1L
      if (count == 1L) bad else maintenance_garch_fit(tag = paste0("grid", count))
    }, .package = "garchx")
    fit <- auto_garchx(rep(1, 30), max.p = 1L, max.q = 0L, final_refit = FALSE, verbose = FALSE)
    expect_identical(fit$fit$tag, "grid2")
    expect_equal(c(fit$q, fit$p, fit$o), c(0, 1, 0))
  }
})

test_that("GARCH grid preserves exact BIC tie order and reports all rejection categories", {
  count <- 0L
  testthat::local_mocked_bindings(garchx = function(...) {
    count <<- count + 1L
    maintenance_garch_fit(tag = paste0("grid", count))
  }, .package = "garchx")
  fit <- auto_garchx(rep(1, 30), max.p = 1L, max.q = 0L, final_refit = FALSE, verbose = FALSE)
  expect_identical(fit$fit$tag, "grid1")
  expect_equal(c(fit$q, fit$p, fit$o), c(0, 0, 1))
  for (case in c("fit error", "convergence", "coefficients", "BIC")) {
    testthat::local_mocked_bindings(garchx = function(...) {
      switch(case, "fit error" = stop("deliberate engine failure"),
        convergence = maintenance_garch_fit(status = NULL),
        coefficients = maintenance_garch_fit(coefficients = c(intercept = Inf)),
        BIC = maintenance_garch_fit(quality = NA_real_))
    }, .package = "garchx")
    expect_error(auto_garchx(rep(1, 30), max.p = 1L, max.q = 0L, verbose = FALSE),
                 paste0("no converged model.*", case))
  }
})

test_that("backcast refit uses exactly the grid eligibility rules and strict improvement", {
  grid <- maintenance_garch_fit(tag = "grid")
  refits <- list(
    improve = maintenance_garch_fit(2, tag = "refit"),
    equal = maintenance_garch_fit(1, tag = "refit"),
    worse = maintenance_garch_fit(0, tag = "refit"),
    unconverged = maintenance_garch_fit(20, 1), missing = maintenance_garch_fit(20, NULL),
    bad_parameters = maintenance_garch_fit(20, 0, c(intercept = Inf)),
    bad_bic = maintenance_garch_fit(NA_real_), error = structure("failed", class = "try-error"))
  for (case in names(refits)) {
    backcasts <- list()
    testthat::local_mocked_bindings(garchx = function(..., backcast.values = NULL) {
      if (is.null(backcast.values)) return(grid)
      backcasts[[1L]] <<- backcast.values
      refits[[case]]
    }, .package = "garchx")
    fit <- auto_garchx(seq_len(30) / 30, max.p = 1L, max.q = 0L, search_o = 0L, verbose = FALSE)
    expect_identical(fit$refit_used, case == "improve")
    expect_identical(fit$fit, if (case == "improve") refits[[case]] else grid)
    expect_equal(backcasts[[1L]], 0.2 / (1 - 0.1 - 0.5), tolerance = 1e-12)
  }
})

test_that("backcast skip conditions preserve a qualified grid model", {
  for (coefficients in list(c(arch1 = 0.1), c(intercept = 0.2, arch1 = 2))) {
    grid <- maintenance_garch_fit(coefficients = coefficients)
    testthat::local_mocked_bindings(garchx = function(..., backcast.values = NULL) {
      if (!is.null(backcast.values)) stop("UNEXPECTED REFIT")
      grid
    }, .package = "garchx")
    fit <- auto_garchx(seq_len(30), max.p = 1L, max.q = 0L, search_o = 0L, verbose = FALSE)
    expect_identical(fit$fit, grid)
    expect_false(fit$refit_used)
  }
})

test_that("real garchx exposes top-level numeric convergence and engine order mapping", {
  set.seed(51026)
  fit <- suppressWarnings(auto_garchx(stats::rnorm(180), max.p = 1, max.q = 1,
    search_o = 0L, final_refit = TRUE, verbose = FALSE))
  expect_true(is.numeric(fit$fit$convergence))
  expect_equal(fit$fit$convergence, 0)
  expect_equal(fit$fit$order, c(fit$q, fit$p, fit$o))
  expect_equal(fit$fit$garchOrder, fit$q)
  expect_equal(fit$fit$archOrder, fit$p)
})
