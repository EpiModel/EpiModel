context("DCM extended models")

# SI Models ---------------------------------------------------------------

test_that("SI, 1G, CL: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$control$nruns, 1)
  expect_equal(x$epi$i.num[150, 1], 388.1552, tol = 0.0001)
  expect_equal(max(x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num + x$epi$i.num), 501)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 1G, CL: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1),
                     act.rate = 0.25)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$control$nruns, 9)
  expect_equal(x$epi$i.num[150, 1], 38.3715, tol = 0.0001)
  expect_equal(x$epi$i.num[150, 4], 500.9153, tol = 0.0001)
  expect_equal(x$epi$i.num[150, 9], 501)
  expect_equal(max(x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num + x$epi$i.num), 501)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 1G, CL: varying act.rate", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.1,
                     act.rate = seq(0.25, 2, 0.25))
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$control$nruns, 8)
  expect_equal(x$epi$i.num[150, 1], 38.3715, tol = 0.0001)
  expect_equal(x$epi$i.num[150, 2], 388.1552, tol = 0.0001)
  expect_equal(x$epi$i.num[150, 5], 500.998, tol = 0.0001)
  expect_equal(max(x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num + x$epi$i.num), 501)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 1G, CL: varying inf.prob and act.rate", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.8, 0.1),
                     act.rate = seq(0.25, 2, 0.25))
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$control$nruns, 8)
  expect_equal(x$epi$i.num[150, 1], 38.3715, tol = 0.0001)
  expect_equal(x$epi$i.num[150, 2], 500.9153, tol = 0.0001)
  expect_equal(x$epi$i.num[10, 5], 178.344, tol = 0.0001)
  expect_equal(max(x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num + x$epi$i.num), 501)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 2G, CL: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1")
  init <- init.dcm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[200, 1], 278.8645, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[200, 1], 212.65, tol = 0.0001)
  expect_equal(max(x$epi$i.num), 501, tol = 0.0001)
  expect_equal(max(x$epi$s.num + x$epi$i.num), 501)
  expect_equal(max(x$epi$i.num + x$epi$i.num.g2 +
                     x$epi$s.num + x$epi$s.num.g2), 1001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 2G, CL: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1),
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1")
  init <- init.dcm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[200, 2], 278.8645, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[200, 2], 212.65, tol = 0.0001)
  expect_equal(x$epi$i.num[200, 4], 494.4093, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[200, 4], 403.4559, tol = 0.0001)
  expect_equal(max(x$epi$i.num), 501)
  expect_equal(max(x$epi$s.num + x$epi$i.num), 501)
  expect_equal(max(x$epi$i.num + x$epi$i.num.g2 +
                     x$epi$s.num + x$epi$s.num.g2), 1001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 1G, OP: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2, act.rate = 0.25,
                     a.rate = 1 / 100, ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[150, 1], 177.8674, tol = 0.0001)
  expect_equal(x$epi$s.num[200, 1], 154.4402, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 1G, OP: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[50, 1], 1.9681, tol = 0.0001)
  expect_equal(x$epi$i.num[50, 5], 21.7698, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 2G, OP: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     a.rate = 1 / 100, a.rate.g2 = NA, ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100, di.rate = 1 / 90,
                     di.rate.g2 = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[250, 1], 134.9654, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[250, 1], 98.9927, tol = 0.0001)
  expect_equal(max(x$epi$i.num), 327, tol = 0.001)
  expect_equal(max(x$epi$i.num.g2), 277, tol = 0.001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SI, 2G, OP: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     a.rate = 1 / 100, a.rate.g2 = NA, ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100, di.rate = 1 / 90,
                     di.rate.g2 = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SI", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[250, 1], 15.08, tol = 0.0001)
  expect_equal(x$epi$i.num[250, 5], 320.8842, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[250, 1], 15.0784, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[250, 5], 222.7064, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})


# SIR Models --------------------------------------------------------------

test_that("SIR, 1G, CL: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2, act.rate = 0.25,
                     rec.rate = 1 / 50)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[250, 1], 107.788, tol = 0.0001)
  expect_equal(max(x$epi$s.num + x$epi$i.num + x$epi$r.num), 501)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIR, 1G, CL: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[100, 1], 1.6207, tol = 0.0001)
  expect_equal(x$epi$i.num[100, 5], 122.3146, tol = 0.0001)
  expect_equal(max(x$epi$s.num + x$epi$i.num + x$epi$r.num), 501)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIR, 2G, CL: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     rec.rate = 1 / 100, rec.rate.g2 = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 0, r.num.g2 = 0)
  control <- control.dcm(type = "SIR", nsteps = 1000, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[250, 1], 137.0072, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[250, 1], 102.8639, tol = 0.0001)
  expect_equal(max(x$epi$s.num + x$epi$i.num + x$epi$r.num), 501)
  expect_equal(max(x$epi$s.num + x$epi$i.num + x$epi$r.num +
                     x$epi$s.num.g2 + x$epi$i.num.g2 + x$epi$r.num.g2), 1001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIR, 2G, CL: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     rec.rate = 1 / 100, rec.rate.g2 = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 0, r.num.g2 = 0)
  control <- control.dcm(type = "SIR", nsteps = 1000, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[250, 1], 18.9723, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[250, 5], 179.5086, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIR, 1G, OP: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 2, rec.rate = 1 / 50, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90, dr.rate = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[50, 1], 272.9788, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIR, 1G, OP: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 1, rec.rate = 1 / 50, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90, dr.rate = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[50, 1], 26.55524, tol = 0.0001)
  expect_equal(x$epi$i.num[50, 5], 297.3653, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIR, 2G, OP: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2, inf.prob.g2 = 0.1,
                     act.rate = 1, balance = "g1",
                     rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100, a.rate.g2 = NA,
                     ds.rate = 1 / 100, ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90, di.rate.g2 = 1 / 90,
                     dr.rate = 1 / 100, dr.rate.g2 = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 1, r.num.g2 = 0)
  control <- control.dcm(type = "SIR", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[50, 1], 154.4498, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[50, 1], 114.8753, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIR, 2G, OP: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05), inf.prob.g2 = 0.1,
                     act.rate = 1, balance = "g1",
                     rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100, a.rate.g2 = NA,
                     ds.rate = 1 / 100, ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90, di.rate.g2 = 1 / 90,
                     dr.rate = 1 / 100, dr.rate.g2 = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 1, r.num.g2 = 0)
  control <- control.dcm(type = "SIR", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[50, 1], 26.55524, tol = 0.0001)
  expect_equal(x$epi$i.num[50, 5], 298.6245, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[50, 5], 213.1567, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})


# SIS Models --------------------------------------------------------------

test_that("SIS, 1G, CL: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[250, 1], 256.7584, tol = 0.0001)
  expect_equal(max(x$epi$s.num + x$epi$i.num), 501)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIS, 1G, CL: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1),
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[50, 1], 1.2741, tol = 0.0001)
  expect_equal(x$epi$i.num[50, 5], 122.0875, tol = 0.0001)
  expect_equal(x$epi$i.num[50, 9], 447.6173, tol = 0.0001)
  expect_equal(max(x$epi$s.num + x$epi$i.num), 501)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIS, 2G, CL: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     rec.rate = 1 / 100, rec.rate.g2 = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SIS", nsteps = 1000, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[250, 1], 164.0793, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[250, 1], 121.0725, tol = 0.0001)
  expect_equal(max(x$epi$s.num + x$epi$i.num +
                     x$epi$s.num.g2 + x$epi$i.num.g2), 1001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIS, 2G, CL: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1),
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     rec.rate = 1 / 100, rec.rate.g2 = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SIS", nsteps = 1000, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[100, 1], 2.2167, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[100, 1], 2.1833, tol = 0.0001)
  expect_equal(x$epi$i.num[100, 9], 199.7184, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[100, 9], 74.7278, tol = 0.0001)
  expect_equal(max(x$epi$s.num + x$epi$i.num +
                     x$epi$s.num.g2 + x$epi$i.num.g2), 1001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIS, 1G, OP: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.5, rec.rate = 1 / 50,
                     a.rate = 1 / 100, ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[100, 1], 249.2884, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIS, 1G, OP: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1),
                     act.rate = 0.5, rec.rate = 1 / 50,
                     a.rate = 1 / 100, ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[100, 1], 6.3054, tol = 0.0001)
  expect_equal(x$epi$i.num[100, 9], 428.3781, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIS, 2G, OP: 1 run", {
  skip_on_cran()
  param <- param.dcm(inf.prob = 0.2,
                     inf.prob.g2 = 0.1, act.rate = 0.5,
                     balance = "g1", rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100, a.rate.g2 = NA, ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100, di.rate = 1 / 90,
                     di.rate.g2 = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 1)
  control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[150, 1], 182.0338, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[150, 1], 136.3023, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})

test_that("SIS, 2G, OP: varying inf.prob", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1), inf.prob.g2 = 0.1,
                     act.rate = 0.5, balance = "g1",
                     rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100, a.rate.g2 = NA,
                     ds.rate = 1 / 100, ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90, di.rate.g2 = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 1)
  control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_equal(x$epi$i.num[100, 1], 6.3054, tol = 0.0001)
  expect_equal(x$epi$i.num[100, 9], 424.9165, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[100, 1], 6.3054, tol = 0.0001)
  expect_equal(x$epi$i.num.g2[100, 9], 279.5831, tol = 0.0001)
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"), run = 2)
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})


# Other cases -------------------------------------------------------------

test_that("SIS, 2G, OP: balance = g2", {
  skip_on_cran()
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1), inf.prob.g2 = 0.1,
                     act.rate.g2 = 0.5, balance = "g2",
                     rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100, a.rate.g2 = NA,
                     ds.rate = 1 / 100, ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90, di.rate.g2 = 1 / 90)
  init <- init.dcm(s.num = 5000, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 1)
  control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
  plot(x)
  plot(x, y = "i.num")
  plot(x, y = c("i.num", "s.num"))
  expect_output(summary(x, at = 2), regexp = "EpiModel Summary")
})
