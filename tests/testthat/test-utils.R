context("Utility Functions")

test_that("brewer_ramp returns correct output", {

  expect_true(length(brewer_ramp(100, plt = "Spectral")) == 100)
  expect_false(length(brewer_ramp(100, plt = "Spectral")) == 50)
  expect_error(brewer_ramp(100, plt = "SamBob"))

})
