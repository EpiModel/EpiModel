context("plot.epi.data.frame and as.epi.data.frame")

build_small_netsim <- function() {
  nw <- network_initialize(n = 50)
  est <- netest(nw,
                formation = ~edges,
                target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3, rec.rate = 0.05)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS",
                         nsims = 2,
                         nsteps = 8,
                         verbose = FALSE)
  netsim(est, param, init, control)
}

test_that("as.epi.data.frame attaches the class on valid input", {
  df <- data.frame(time = rep(1:5, 2),
                   sim = rep(1:2, each = 5),
                   value = stats::runif(10))
  out <- as.epi.data.frame(df)
  expect_s3_class(out, "epi.data.frame")
  expect_s3_class(out, "data.frame")
})

test_that("as.epi.data.frame validates input shape", {
  expect_error(as.epi.data.frame(list(a = 1)),
               "must be a `data.frame`")
  expect_error(as.epi.data.frame(data.frame(x = 1)),
               "must contain a `time` and a `sim` column")
  # Uneven rows per sim.
  df_bad <- data.frame(time = c(1, 2, 1), sim = c(1, 1, 2))
  expect_error(as.epi.data.frame(df_bad),
               "same number of time step")
})

test_that("plot.epi.data.frame renders without error from a netsim roundtrip", {
  skip_on_cran()
  mod <- build_small_netsim()
  df <- as.data.frame(mod)
  expect_s3_class(df, "epi.data.frame")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(df, y = "i.num"))
  expect_silent(plot(df, y = "i.num", mean.line = FALSE,
                     qnts = 0.9, sim.lines = TRUE))
})

test_that("plot.epi.data.frame works on a manually built data.frame", {
  skip_on_cran()
  df <- data.frame(time = rep(1:6, 3),
                   sim = rep(1:3, each = 6),
                   i.num = c(10, 11, 12, 13, 14, 15,
                             10, 12, 14, 16, 18, 20,
                             10, 13, 16, 19, 22, 25))
  df <- as.epi.data.frame(df)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(df, y = "i.num"))
})
