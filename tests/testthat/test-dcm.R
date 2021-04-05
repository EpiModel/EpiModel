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
                     act.rate = 0.25, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})


test_that("SI, 2G, OP: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     a.rate = 1 / 100, a.rate.g2 = NA, ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100, di.rate = 1 / 90,
                     di.rate.g2 = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SI", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})


# SIR Models --------------------------------------------------------------

test_that("SIR, 1G, CL: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05),
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SIR, 1G, OP: 1 run", {
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 2, rec.rate = 1 / 50, a.rate = 1 / 100,
                     ds.rate = 1 / 100, di.rate = 1 / 90, dr.rate = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SIR, 2G, OP: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.5, 0.05), inf.prob.g2 = 0.1,
                     act.rate = 1, balance = "g1",
                     rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     a.rate = 1 / 100, a.rate.g2 = NA,
                     ds.rate = 1 / 100, ds.rate.g2 = 1 / 100,
                     di.rate = 1 / 90, di.rate.g2 = 1 / 90,
                     dr.rate = 1 / 100, dr.rate.g2 = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1, r.num = 0,
                   s.num.g2 = 500, i.num.g2 = 1, r.num.g2 = 0)
  control <- control.dcm(type = "SIR", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})



# SIS Models --------------------------------------------------------------

test_that("SIS, 1G, CL: 1 run", {
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.25, rec.rate = 1 / 50)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SIS", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SIS, 2G, CL: varying inf.prob", {
  param <- param.dcm(inf.prob = seq(0.1, 0.9, 0.1),
                     act.rate = 0.25, inf.prob.g2 = 0.1, balance = "g1",
                     rec.rate = 1 / 100, rec.rate.g2 = 1 / 100)
  init <- init.dcm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 0)
  control <- control.dcm(type = "SIS", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})

test_that("SIS, 1G, OP: 1 run", {
  param <- param.dcm(inf.prob = 0.2,
                     act.rate = 0.5, rec.rate = 1 / 50,
                     a.rate = 1 / 100, ds.rate = 1 / 100, di.rate = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SIS", nsteps = 10, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_is(x, "dcm")
})



# Check flows -------------------------------------------------------------

test_that("si.flow correct for closed SI model, RK4 method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 2)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(df$i.num[2], df$i.num[1] + df$si.flow[1])
  expect_equal(df$s.num[2], df$s.num[1] - df$si.flow[1])
  expect_equal(df$si.flow[1], 96.58919, tol = 0.0001)
})


test_that("si.flow correct for closed SI model, RK4 method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4,
                     a.rate = 0.02, ds.rate = 0.01, di.rate = 0.01)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 2)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(df$num[2],
               df$num[1] + df$a.flow[1] - df$ds.flow[1] - df$di.flow[1])
  expect_equal(df$i.num[2], df$i.num[1] + df$si.flow[1] - df$di.flow[1])
  expect_equal(df$si.flow[1], 96.06876, tol = 0.0001)
})

test_that("si.flow correct for closed SI model, Euler method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 2, odemethod = "euler")
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(df$i.num[2], df$i.num[1] + df$si.flow[1])
  expect_equal(df$s.num[2], df$s.num[1] - df$si.flow[1])
  expect_equal(df$si.flow[1], 67.76348, tol = 0.0001)
})


test_that("si.flow correct for closed SI model, Euler method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4,
                     a.rate = 0.02, ds.rate = 0.01, di.rate = 0.01)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 2, odemethod = "euler")
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(df$num[2],
               df$num[1] + df$a.flow[1] - df$ds.flow[1] - df$di.flow[1])
  expect_equal(df$i.num[2], df$i.num[1] + df$si.flow[1] - df$di.flow[1])
  expect_equal(df$si.flow[1], 67.76348, tol = 0.0001)
})


# Check dt fractional -----------------------------------------------------

test_that("fractional dt returns, RK4 method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 10, dt = 0.5)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(tail(df$i.num, 1), df$i.num[1] + sum(df$si.flow, na.rm = TRUE))
  expect_equal(df$time, seq(1, 10, 0.5))
})

test_that("fractional dt returns, Euler method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 10, dt = 0.5,
                         odemethod = "euler")
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(tail(df$i.num, 1), df$i.num[1] + sum(df$si.flow, na.rm = TRUE))
  expect_equal(df$time, seq(1, 10, 0.5))
})


# DCM interventions -------------------------------------------------------

test_that("DCM interventions, SI model", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4,
                     inter.eff = 1)
  expect_equal(param$inter.start, 1)

  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4,
                     inter.eff = 1, inter.start = 5)
  expect_equal(param$inter.start, 5)

  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 10)
  mod <- dcm(param, init, control)

  df <- as.data.frame(mod)
  expect_equal(sum(df$si.flow[5:10], na.rm = TRUE), 0)
  expect_true(length(unique(df$i.num[5:10])) == 1)
  expect_true(length(unique(df$s.num[5:10])) == 1)
})

test_that("DCM interventions, SIS model", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4, rec.rate = 0.1,
                     inter.eff = 1, inter.start = 5)
  expect_equal(param$inter.start, 5)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SIS", nsteps = 10)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(sum(df$si.flow[5:10], na.rm = TRUE), 0)
  expect_equal(df$s.num[10], 28221.24, tol = 0.0001)
  expect_equal(df$i.num[10], 528.7624, tol = 0.0001)
})

test_that("DCM interventions, SIR model", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4, rec.rate = 0.1,
                     inter.eff = 1, inter.start = 5)
  expect_equal(param$inter.start, 5)
  init <- init.dcm(s.num = 28650, i.num = 100, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 10)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(sum(df$si.flow[5:10], na.rm = TRUE), 0)
  expect_true(length(unique(df$s.num[5:10])) == 1)
  expect_equal(df$i.num[10], 526.6549, tol = 0.0001)
})



# param, init, control ------------------------------------------------------

test_that("control checks", {
  control <- control.dcm(type = "SI", nsteps = 10, foo = "boo")
  expect_true(control$foo == "boo")
  expect_error(control.dcm(type = "SI"), "Specify nsteps")
  expect_error(control.dcm(type = "SEIR", nsteps = 10), "Specify type as")
})



# dcm extension check -----------------------------------------------------

test_that("test extension model, NA output", {
  eMod2 <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {

      # Total population size
      num <- s.num + i.num + r.num

      # Define lambdas
      lambda.direct <- beta.direct * i.num
      lambda.envirn <- beta.envirn * e.num

      # Differential equations
      dS <- -lambda.direct*s.num - lambda.envirn*s.num + mu*num - mu*s.num
      dI <- lambda.direct*s.num + lambda.envirn*s.num - mu*i.num - gamma*i.num
      dR <- gamma*i.num - mu*r.num
      dE <- xi*i.num - delta*e.num

      # Outputs
      list(c(dS, dI, dR, dE,
             si.direct.flow = lambda.direct*s.num,
             si.envirn.flow = lambda.envirn*s.num))
    })
  }

  # Parameterize the model with the following parameters this time
  param <- param.dcm(beta.direct = 0.0001,
                     beta.envirn = 0.00001,
                     mu = 0.01,
                     gamma = 0.2,
                     xi = 2,
                     delta = 0.5)

  # Use similar initial conditions and control settings as the first model, but
  # make sure to update the new.mod function name!
  init <- init.dcm(s.num = 10000, i.num = 1, r.num = 0, e.num = 0,
                   si.direct.flow = 0, si.envirn.flow = 0)
  control <- control.dcm(nsteps = 10, dt = 0.1, new.mod = eMod2)

  # Run the simulation and examine the output as a data frame
  sim <- dcm(param, init, control)
  expect_s3_class(sim, "dcm")

  df <- as.data.frame(sim)
  expect_true(all(!is.na(df)))
})
