context("ICM standard models")

# SI Models ---------------------------------------------------------------

test_that("SI, 1G, CL: 1 sim", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SI", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SI, 1G, CL: 2 sims", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SI, 1G, OP: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SI", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SI, 2G, OP: 2 sims", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     a.rate = 1 / 100, a.rate.g2 = NA, ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100, di.rate = 1 / 90,
                     di.rate.g2 = 1 / 90)
  init <- init.icm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
  control <- control.icm(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})


# SIR Models --------------------------------------------------------------

test_that("SIR, 1G, CL: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.icm(type = "SIR", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SIR, 1G, OP: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 2, rec.rate = 1 / 50, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90, dr.rate = 1 / 100)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.icm(type = "SIR", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SIR, 2G, OP: 2 sims", {
  param <- param.icm(inf.prob = 0.2, inf.prob.g2 = 0.1,
                     act.rate = 1, balance = "g1",
                     rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100, a.rate.g2 = NA,
                     ds.rate = 1 / 100, ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90, di.rate.g2 = 1 / 90,
                     dr.rate = 1 / 100, dr.rate.g2 = 1 / 100)
  init <- init.icm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 1, r.num.g2 = 0)
  control <- control.icm(type = "SIR", nsteps = 5, nsims = 2, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})



# SIS Models --------------------------------------------------------------

test_that("SIS, 1G, CL: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SIS", nsteps = 5, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SIS, 2G, CL: 2 sims", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     rec.rate = 1 / 100, rec.rate.g2 = 1 / 100)
  init <- init.icm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 0)
  control <- control.icm(type = "SIS", nsteps = 5, nsims = 2, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})

test_that("SIS, 1G, OP: 1 sim", {
  param <- param.icm(inf.prob = 0.2,
                     act.rate = 0.5, rec.rate = 1 / 50,
                     a.rate = 1 / 100, ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SIS", nsteps = 10, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})



# Other -------------------------------------------------------------------

test_that("Act rate balance works when specified to g2", {
  param <- param.icm(inf.prob = 0.2,
                     inf.prob.g2 = 0.1,
                     act.rate.g2 = 0.5,
                     balance = "g2",
                     rec.rate = 1 / 50,
                     rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100,
                     a.rate.g2 = NA,
                     ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90,
                     di.rate.g2 = 1 / 90)
  init <- init.icm(s.num = 50,
                   i.num = 10,
                   s.num.g2 = 100,
                   i.num.g2 = 1)
  control <- control.icm(type = "SIS",
                         nsteps = 10,
                         nsims = 1,
                         verbose = FALSE)
  x <- icm(param, init, control)
  expect_is(x, "icm")
})
