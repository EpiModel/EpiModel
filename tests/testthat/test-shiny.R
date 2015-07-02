context("Shiny app")

test_that("epiweb works", {
  expect_error(epiweb(class = "foo"))
})
