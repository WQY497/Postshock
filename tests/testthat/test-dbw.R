make_dbw_toy_data <- function() {
  set.seed(1)
  
  T <- 50L
  
  target <- data.frame(
    x1 = rnorm(T),
    x2 = rnorm(T, 0.2)
  )
  
  donor1 <- data.frame(
    x1 = rnorm(T, 0.1),
    x2 = rnorm(T, -0.1)
  )
  
  donor2 <- data.frame(
    x1 = rnorm(T, -0.2),
    x2 = rnorm(T, 0.4)
  )
  
  list(
    X = list(target, donor1, donor2),
    shock_time = c(40L, 40L, 40L)
  )
}


test_that("dbw returns feasible simplex weights on toy data", {
  dat <- make_dbw_toy_data()
  
  res <- dbw(
    X = dat$X,
    dbw_indices = 1:2,
    shock_time_vec = dat$shock_time,
    center = TRUE,
    scale = TRUE,
    sum_to_1 = TRUE,
    bounded_below_by = 0,
    bounded_above_by = 1,
    normchoice = "l2",
    penalty_normchoice = "l2",
    penalty_lambda = 1e-2
  )
  
  # 检查整体返回结构
  expect_named(
    res,
    c("opt_params", "convergence", "loss"),
    ignore.order = FALSE
  )
  
  # 两个 donors 应该返回两个权重
  expect_type(res$opt_params, "double")
  expect_length(res$opt_params, 2L)
  expect_true(all(is.finite(res$opt_params)))
  
  # simplex 约束：权重和为 1
  expect_equal(
    sum(res$opt_params),
    1,
    tolerance = 1e-6
  )
  
  # box constraints：权重应当位于 [0, 1]
  # 留一点数值优化误差空间
  expect_true(all(res$opt_params >= -1e-7))
  expect_true(all(res$opt_params <= 1 + 1e-7))
  
  # 优化必须真正收敛
  expect_identical(res$convergence, "convergence")
  
  # loss 必须是一个有限的非负标量
  expect_type(res$loss, "double")
  expect_length(res$loss, 1L)
  expect_true(is.finite(res$loss))
  expect_true(res$loss >= 0)
})


test_that("dbw rejects infeasible simplex bounds", {
  dat <- make_dbw_toy_data()
  
  # 两个 donor，每个权重至少为 0.6，
  # 那么权重和至少为 1.2，不可能同时满足 sum(W) = 1
  expect_error(
    dbw(
      X = dat$X,
      dbw_indices = 1:2,
      shock_time_vec = dat$shock_time,
      sum_to_1 = TRUE,
      bounded_below_by = 0.6,
      bounded_above_by = 1
    ),
    regexp = "Infeasible constraints"
  )
})


test_that("dbw rejects input without donors", {
  dat <- make_dbw_toy_data()
  
  expect_error(
    dbw(
      X = dat$X[1],
      dbw_indices = 1:2,
      shock_time_vec = 40L
    ),
    regexp = "no donors|empty",
    ignore.case = TRUE
  )
})


test_that("dbw rejects invalid norm choices", {
  dat <- make_dbw_toy_data()
  
  expect_error(
    dbw(
      X = dat$X,
      dbw_indices = 1:2,
      shock_time_vec = dat$shock_time,
      normchoice = "invalid"
    ),
    regexp = "normchoice must be"
  )
  
  expect_error(
    dbw(
      X = dat$X,
      dbw_indices = 1:2,
      shock_time_vec = dat$shock_time,
      penalty_normchoice = "invalid"
    ),
    regexp = "penalty_normchoice must be"
  )
})



test_that("dbw returns weight one for a single donor under simplex constraints", {
  set.seed(10)
  
  target <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40)
  )
  
  donor <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40)
  )
  
  expect_no_warning({
    res <- dbw(
      X = list(target, donor),
      dbw_indices = 1:2,
      shock_time_vec = c(30L, 30L),
      center = TRUE,
      scale = TRUE,
      sum_to_1 = TRUE,
      bounded_below_by = 0,
      bounded_above_by = 1
    )
  })
  
  expect_identical(res$opt_params, 1)
  expect_identical(res$convergence, "convergence")
  expect_true(is.finite(res$loss))
})