context("geom_bands ggplot helper")

test_that("geom_bands returns a ggplot2 layer", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    time = rep(1:10, 3),
    sim = rep(1:3, each = 10),
    i.num = stats::runif(30, 0, 100)
  )

  layer <- geom_bands(mapping = ggplot2::aes(time, i.num),
                      lower = 0.1, upper = 0.9, fill = "firebrick")
  expect_s3_class(layer, "Layer")
  expect_s3_class(layer, "ggproto")
})

test_that("geom_bands integrates into a ggplot object", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    time = rep(1:10, 3),
    sim = rep(1:3, each = 10),
    i.num = stats::runif(30, 0, 100)
  )

  p <- ggplot2::ggplot() +
    geom_bands(data = df,
               mapping = ggplot2::aes(time, i.num),
               lower = 0.25, upper = 0.75, fill = "steelblue")
  expect_s3_class(p, "ggplot")

  # Building the plot exercises the underlying stat_summary call.
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})
