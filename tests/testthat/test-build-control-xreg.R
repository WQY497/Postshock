test_that("build_control_xreg returns NULL when no control shocks are supplied", {
  expect_null(
    build_control_xreg(
      last = 10L,
      control_spec = NULL
    )
  )
  
  expect_null(
    build_control_xreg(
      last = 10L,
      control_spec = list()
    )
  )
})


test_that("build_control_xreg creates a point-shock indicator", {
  out <- build_control_xreg(
    last = 8L,
    control_spec = list(
      time = 3L,
      length = 1L,
      shape = "point"
    )
  )
  
  expected <- matrix(
    c(0, 0, 1, 0, 0, 0, 0, 0),
    ncol = 1L
  )
  
  colnames(expected) <- "ctrl_shock_1"
  
  expect_equal(out, expected)
  expect_identical(storage.mode(out), "double")
})


test_that("build_control_xreg creates a finite-window shock indicator", {
  out <- build_control_xreg(
    last = 8L,
    control_spec = list(
      time = 3L,
      length = 3L,
      shape = "window"
    )
  )
  
  expected <- matrix(
    c(0, 0, 1, 1, 1, 0, 0, 0),
    ncol = 1L
  )
  
  colnames(expected) <- "ctrl_shock_1"
  
  expect_equal(out, expected)
})


test_that("build_control_xreg creates a step-shock indicator", {
  out <- build_control_xreg(
    last = 8L,
    control_spec = list(
      time = 5L,
      length = 1L,
      shape = "step"
    )
  )
  
  expected <- matrix(
    c(0, 0, 0, 0, 1, 1, 1, 1),
    ncol = 1L
  )
  
  colnames(expected) <- "ctrl_shock_1"
  
  expect_equal(out, expected)
})


test_that("build_control_xreg creates one column for each control shock", {
  control_spec <- list(
    list(
      time = 2L,
      length = 1L,
      shape = "point"
    ),
    list(
      time = 4L,
      length = 2L,
      shape = "window"
    ),
    list(
      time = 6L,
      length = 1L,
      shape = "step"
    )
  )
  
  out <- build_control_xreg(
    last = 8L,
    control_spec = control_spec
  )
  
  expected <- cbind(
    c(0, 1, 0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 1, 1, 0, 0, 0),
    c(0, 0, 0, 0, 0, 1, 1, 1)
  )
  
  colnames(expected) <- c(
    "ctrl_shock_1",
    "ctrl_shock_2",
    "ctrl_shock_3"
  )
  
  storage.mode(expected) <- "double"
  
  expect_equal(out, expected)
  expect_equal(dim(out), c(8L, 3L))
  expect_identical(
    colnames(out),
    c("ctrl_shock_1", "ctrl_shock_2", "ctrl_shock_3")
  )
})


test_that("build_control_xreg accepts a data-frame specification", {
  control_spec <- data.frame(
    time = c(2L, 5L),
    length = c(1L, 2L),
    shape = c("point", "window")
  )
  
  out <- build_control_xreg(
    last = 7L,
    control_spec = control_spec
  )
  
  expected <- cbind(
    c(0, 1, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 1, 0)
  )
  
  colnames(expected) <- c(
    "ctrl_shock_1",
    "ctrl_shock_2"
  )
  
  storage.mode(expected) <- "double"
  
  expect_equal(out, expected)
})


test_that("build_control_xreg supplies data-frame defaults", {
  control_spec <- data.frame(
    time = c(2L, 5L)
  )
  
  out <- build_control_xreg(
    last = 7L,
    control_spec = control_spec
  )
  
  expected <- cbind(
    c(0, 1, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 0, 0)
  )
  
  colnames(expected) <- c(
    "ctrl_shock_1",
    "ctrl_shock_2"
  )
  
  storage.mode(expected) <- "double"
  
  expect_equal(out, expected)
})


test_that("build_control_xreg rejects invalid training lengths", {
  expect_error(
    build_control_xreg(
      last = 0L,
      control_spec = NULL
    ),
    regexp = "single positive integer"
  )
  
  expect_error(
    build_control_xreg(
      last = NA_integer_,
      control_spec = NULL
    ),
    regexp = "single positive integer"
  )
  
  expect_error(
    build_control_xreg(
      last = c(5L, 10L),
      control_spec = NULL
    ),
    regexp = "single positive integer"
  )
})


test_that("build_control_xreg rejects invalid control specifications", {
  expect_error(
    build_control_xreg(
      last = 10L,
      control_spec = "invalid"
    ),
    regexp = "NULL, a list, or a data.frame"
  )
  
  expect_error(
    build_control_xreg(
      last = 10L,
      control_spec = data.frame(
        length = 2L,
        shape = "window"
      )
    ),
    regexp = "must contain a `time` column"
  )
  
  expect_error(
    build_control_xreg(
      last = 10L,
      control_spec = list(
        length = 2L,
        shape = "window"
      )
    ),
    regexp = "valid `time`"
  )
})


test_that("build_control_xreg rejects invalid shock lengths and shapes", {
  expect_error(
    build_control_xreg(
      last = 10L,
      control_spec = list(
        time = 3L,
        length = 0L,
        shape = "window"
      )
    ),
    regexp = "`length` must be >= 1"
  )
  
  expect_error(
    build_control_xreg(
      last = 10L,
      control_spec = list(
        time = 3L,
        length = 2L,
        shape = "triangle"
      )
    ),
    regexp = "Unknown shock shape"
  )
})