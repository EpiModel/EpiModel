context("DCM standard models")

# SI Models ---------------------------------------------------------------

test_that("SI, 1G, CL: 1 run", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SI, 1G, CL: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1),
                     act.rate = 0.25)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})


test_that("SI, 1G, CL: varying inf.prob and act.rate", {
  param <- param.dcm(inf.prob = seq(0.1, 0.8, 0.1),
                     act.rate = seq(0.25, 2, 0.25))
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})


test_that("SI, 1G, OP: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, b.rate = 1/100,
                     ds.rate = 1/100, di.rate = 1/90)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})


test_that("SI, 2G, OP: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     b.rate = 1/100, b.rate.g2 = NA, ds.rate = 1/100,
                     ds.rate.g2 = 1/100, di.rate = 1/90, di.rate.g2 = 1/90)
  init <- init.dcm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SI", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})


# SIR Models --------------------------------------------------------------

test_that("SIR, 1G, CL: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, rec.rate = 1/50,)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SIR, 1G, OP: 1 run", {
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 2, rec.rate = 1/50, b.rate = 1/100,
                     ds.rate = 1/100, di.rate = 1/90, dr.rate = 1/100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SIR, 2G, OP: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05), inf.prob.g2 = 0.1,
                     act.rate = 1, balance = "g1",
                     rec.rate = 1/50, rec.rate.g2 = 1/50,
                     b.rate = 1/100, b.rate.g2 = NA,
                     ds.rate = 1/100, ds.rate.g2 = 1/100,
                     di.rate = 1/90, di.rate.g2 = 1/90,
                     dr.rate = 1/100, dr.rate.g2 = 1/100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 1, r.num.g2 = 0)
  control <- control.dcm(type = "SIR", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})



# SIS Models --------------------------------------------------------------

test_that("SIS, 1G, CL: 1 run", {
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.25, rec.rate = 1/50)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SIS", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SIS, 2G, CL: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1),
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     rec.rate = 1/100, rec.rate.g2 = 1/100)
  init <- init.dcm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SIS", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SIS, 1G, OP: 1 run", {
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.5, rec.rate = 1/50,
                     b.rate = 1/100, ds.rate = 1/100, di.rate = 1/90)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SIS", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})


