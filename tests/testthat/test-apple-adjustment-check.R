apple_adjustment_checker <- function() {
  shared <- new.env(parent = globalenv())
  sys.source(system.file("replication", "apple_diagnostic.R", package = "postshock",
                        mustWork = TRUE), shared)
  shared$check_apple_adjustment
}

apple_adjustment_reference <- function() {
  effects <- c(0.425272469331, 0.218634577538, 1.990999985436)
  weights <- c(0.583683226784, 0.0000000562775710497, 0.416316716939)
  list(omega_vec = effects, weights = weights, combined_omega = sum(weights * effects))
}

test_that("Apple aggregate reference bounds follow unchanged component tolerances", {
  checker <- apple_adjustment_checker()
  reference <- apple_adjustment_reference()
  result <- checker(reference)
  expect_lt(result$aggregation_error, 1e-10)
  expect_equal(result$propagated_bound,
    sum(abs(reference$omega_vec)) * 1e-7 + sum(abs(reference$weights)) * 1e-9 +
      abs(reference$combined_omega - 1.077110996827) + 1e-10, tolerance = 1e-15)
  shifted <- reference
  shifted$weights <- shifted$weights + c(-5e-8, 0, 5e-8)
  shifted$combined_omega <- sum(shifted$weights * shifted$omega_vec)
  expect_gt(abs(shifted$combined_omega - 1.077110996827), 1e-9)
  expect_lt(checker(shifted)$reference_error, result$propagated_bound)
})

test_that("Apple diagnostic rejects aggregation errors despite acceptable references", {
  checker <- apple_adjustment_checker()
  reference <- apple_adjustment_reference()
  corrupted <- reference
  corrupted$combined_omega <- corrupted$combined_omega + 5e-10
  expect_lt(abs(corrupted$combined_omega - 1.077110996827), 1e-9)
  expect_error(checker(corrupted), "aggregation_error")
  corrupted$combined_omega <- corrupted$combined_omega + 1e-8
  expect_error(checker(corrupted), "aggregation_error")
})

test_that("recorded Accelerate output reproduces the old failure and passes the derived checks", {
  observed <- list(
    omega_vec = c(0.42527246933059615, 0.21863457753676729, 1.9909999854370486),
    weights = c(0.58368322017212315, 5.6277568609149234e-08, 0.41631672355030935),
    combined_omega = 1.0771110071795209)
  expect_gt(abs(observed$combined_omega - 1.077110996827), 1e-9)
  verification <- apple_adjustment_checker()(observed)
  expect_lt(verification$aggregation_error, 1e-10)
  expect_lt(verification$reference_error, verification$propagated_bound)
})

test_that("Apple component reference checks and finite-value checks remain enforced", {
  checker <- apple_adjustment_checker()
  reference <- apple_adjustment_reference()
  wrong_weights <- reference
  wrong_weights$weights <- wrong_weights$weights + c(-2e-7, 0, 2e-7)
  wrong_weights$combined_omega <- sum(wrong_weights$weights * wrong_weights$omega_vec)
  expect_error(checker(wrong_weights), "weight_tolerance")
  wrong_effects <- reference
  wrong_effects$omega_vec[1L] <- wrong_effects$omega_vec[1L] + 2e-9
  wrong_effects$combined_omega <- sum(wrong_effects$weights * wrong_effects$omega_vec)
  expect_error(checker(wrong_effects), "effect_tolerance")
  for (field in names(reference)) {
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
      invalid <- reference
      invalid[[field]][1L] <- value
      expect_error(checker(invalid))
    }
    invalid <- reference
    invalid[[field]] <- numeric()
    expect_error(checker(invalid))
  }
})
