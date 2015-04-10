context("Utility Functions")

test_that("brewer_ramp returns correct output", {

  expect_true(length(brewer_ramp(100, plt = "Spectral")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Spectral", delete.lights = FALSE)) == 100)
  expect_false(length(brewer_ramp(100, plt = "Spectral")) == 50)
  expect_false(length(brewer_ramp(100, plt = "Spectral", delete.lights = FALSE)) == 50)

  expect_true(length(brewer_ramp(100, plt = "Accent")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Accent", delete.lights = FALSE)) == 100)

  expect_true(length(brewer_ramp(100, plt = "Oranges")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Oranges", delete.lights = FALSE)) == 100)

  expect_true(length(brewer_ramp(100, plt = "Set1")) == 100)

  expect_error(brewer_ramp(100, plt = "Jimmy"))
  expect_error(brewer_ramp(100, plt = "Jimmy", delete.lights = FALSE))
  expect_error(brewer_ramp(-1, plt = "Spectral"))

})


test_that("deleteAttr returns correct output", {

  l <- list(a = 1:5, b = 6:10)
  expect_is(deleteAttr(l, 5), "list")
  expect_true(length(unique(sapply(deleteAttr(l, 2:3), length))) == 1)

  l2 <- list(a = 1:3, b = 5:20)
  expect_error(deleteAttr(l2, 2:4))
  expect_error(deleteAttr(as.data.frame(l), 1))

  expect_equal(l, deleteAttr(l, NULL))

})
