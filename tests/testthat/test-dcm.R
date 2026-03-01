context("DCM Standard Models")

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

test_that("si.flow correct for closed SI model, lsoda method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 2)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(df$i.num[2], df$i.num[1] + df$si.flow[1])
  expect_equal(df$s.num[2], df$s.num[1] - df$si.flow[1])
  expect_equal(df$si.flow[1], 96.72141, tol = 0.0001)
})


test_that("si.flow correct for open SI model, lsoda method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4,
                     a.rate = 0.02, ds.rate = 0.01, di.rate = 0.01)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 2)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(df$num[2],
               df$num[1] + df$a.flow[1] - df$ds.flow[1] - df$di.flow[1])
  expect_equal(df$i.num[2], df$i.num[1] + df$si.flow[1] - df$di.flow[1])
  expect_equal(df$si.flow[1], 96.19311, tol = 0.0001)
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
  expect_equal(df$time, seq(1, 10, 0.5))
})

test_that("fractional dt returns, Euler method", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SI", nsteps = 10, dt = 0.5,
                         odemethod = "euler")
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(df$time, seq(1, 10, 0.5))
})

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
  expect_equal(sum(df$si.flow[5:10], na.rm = TRUE), 0, tol = 0.01)
  expect_equal(max(df$i.num[5:10]) - min(df$i.num[5:10]), 0, tol = 0.01)
  expect_equal(max(df$s.num[5:10]) - min(df$s.num[5:10]), 0, tol = 0.01)
})

test_that("DCM interventions, SIS model", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4, rec.rate = 0.1,
                     inter.eff = 1, inter.start = 5)
  expect_equal(param$inter.start, 5)
  init <- init.dcm(s.num = 28650, i.num = 100)
  control <- control.dcm(type = "SIS", nsteps = 10)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(sum(df$si.flow[5:10], na.rm = TRUE), 0, tol = 0.01)
  expect_equal(df$s.num[10], 28155.07472, tol = 0.0001)
  expect_equal(df$i.num[10], 594.92528, tol = 0.0001)
})

test_that("DCM interventions, SIR model", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 3.4, rec.rate = 0.1,
                     inter.eff = 1, inter.start = 5)
  expect_equal(param$inter.start, 5)
  init <- init.dcm(s.num = 28650, i.num = 100, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 10)
  mod <- dcm(param, init, control)
  df <- as.data.frame(mod)
  expect_equal(sum(df$si.flow[5:10], na.rm = TRUE), 0, tol = 0.01)
  expect_equal(max(df$s.num[5:10]) - min(df$s.num[5:10]), 0, tol = 0.01)
  expect_equal(df$i.num[10], 592.14131, tol = 0.5)
})



# param, init, control ------------------------------------------------------

test_that("control checks", {
  control <- control.dcm(type = "SI", nsteps = 10, foo = "boo")
  expect_true(control$foo == "boo")
  expect_error(control.dcm(type = "SI"), "Specify nsteps")
  expect_error(control.dcm(type = "SEIR", nsteps = 10), "Specify type as")
})



# dcm extension check -----------------------------------------------------


context("DCM Delayed Differential Equation Models")

test_that("Delayed differntial equation models function", {

  lagsi <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {

      num <- s.num + i.num

      if (t < 2) {
        i.num.foi <- 0
      } else {
        i.num.foi <- lagvalue(t - 1, 2)
      }

      lambda <- inf.prob * act.rate * i.num.foi / num

      dS <- -lambda * s.num
      dI <- lambda * s.num

      list(c(dS, dI),
           num = num,
           si.flow = lambda * s.num,
           i.num.foi = i.num.foi)
    })
  }

  param <- param.dcm(inf.prob = 0.5, act.rate = 0.25)
  init <- init.dcm(s.num = 100, i.num = 10)
  control <- control.dcm(nsteps = 25, new.mod = lagsi, dede = TRUE)
  mod <- dcm(param, init, control)

  expect_is(mod, "dcm")
  expect_true(mod$control$dede)

})


context("New DCM Models")


test_that("New DCMs Example 1", {

  ## SI function
  intSI <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {

      ## Dynamic Calculations
      # Population size
      num <- s.num + i.num

      # Intervention start time: if time > start,
      #   then multiply lambda by relative hazard
      if (t < start.time) {
        lambda <- inf.prob * act.rate * i.num / num
      } else {
        lambda <- (inf.prob * act.rate * i.num / num) * rel.haz
      }

      ## Flows
      si.flow <- lambda * s.num

      ## Differential Equations
      dS <- -lambda * s.num
      dI <- lambda * s.num

      ## Output
      list(c(dS, dI, si.flow),
           num = num)
    })
  }


  init <- init.dcm(s.num = 999, i.num = 1, si.flow = 0)
  control <- control.dcm(nsteps = 250, dt = 1, new.mod = intSI, verbose = FALSE)
  param1 <- param.dcm(inf.prob = 0.5, act.rate = 0.1,
                      start.time = seq(100, 200, 100), rel.haz = 1)
  param05 <- param.dcm(inf.prob = 0.5, act.rate = 0.1,
                       start.time = seq(100, 200, 100), rel.haz = 0.5)
  mod1 <- dcm(param1, init, control)
  mod2 <- dcm(param05, init, control)

  expect_is(mod1, "dcm")
  expect_is(mod2, "dcm")
})


test_that("New DCMs Example 2", {

  ## Q Mod Function
  Qmod <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {

      ## Dynamic Calculations ##

      # Popsize and prevalence
      h.num <- sh.num + ih.num
      l.num <- sl.num + il.num
      num <- h.num + l.num
      prev <- (ih.num + il.num) / num

      # Contact rates for high specified as a function of
      #   mean and low rates
      c.high <- (c.mean * num - c.low * l.num) / h.num

      # Mixing matrix calculations based on variable Q statistic
      g.hh <- ((c.high * h.num) + (Q * c.low * l.num)) /
        ((c.high * h.num) + (c.low * l.num))
      g.lh <- 1 - g.hh
      g.hl <- (1 - g.hh) * ((c.high * h.num) / (c.low * l.num))
      g.ll <- 1 - g.hl

      # Probability that contact is infected based on mixing probabilities
      p.high <- (g.hh * ih.num / h.num) + (g.lh * il.num / l.num)
      p.low <- (g.ll * il.num / l.num) + (g.hl * ih.num / h.num)

      # Force of infection for high and low groups
      lambda.high <- rho * c.high * p.high
      lambda.low <- rho * c.low * p.low


      ## Differential Equations ##
      dS.high <- -lambda.high * sh.num + nu * ih.num
      dI.high <- lambda.high * sh.num - nu * ih.num

      dS.low <- -lambda.low * sl.num + nu * il.num
      dI.low <- lambda.low * sl.num - nu * il.num


      ## Output ##
      list(c(dS.high, dI.high,
             dS.low, dI.low),
           num = num,
           prev = prev)
    })
  }


  param <- param.dcm(c.mean = 2,
                     c.low = 1.4,
                     rho = 0.75,
                     nu = 6,
                     Q = c(-0.45, 0.5, 1))
  init <- init.dcm(sh.num = 2e7 * 0.02,
                   ih.num = 1,
                   sl.num = 2e7 * 0.98,
                   il.num = 1)
  control <- control.dcm(nsteps = 6, dt = 0.02, new.mod = Qmod, verbose = FALSE)

  mod <- dcm(param, init, control)
  expect_is(mod, "dcm")

})


test_that("DCM inital conditions ordering correct", {

  SEIR <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {
      num <- s.num + e.num + i.num + r.num
      lambda <- inf.prob * act.rate * i.num / num
      se.flow <-  lambda * s.num
      ei.flow <-  sx.rate * e.num
      ir.flow <-  rec.rate * i.num
      dS <- -lambda * s.num
      dE <- lambda * s.num - sx.rate * e.num
      dI <- sx.rate * e.num - rec.rate * i.num
      dR <- rec.rate * i.num
      list(c(dS, dE, dI, dR, se.flow, ei.flow, ir.flow),
           num = num)
    })
  }

  init <- init.dcm(s.num = 980, e.num = 10, i.num = 10, r.num = 0,
                   se.flow = 0, ei.flow = 0, ir.flow = 0)
  expect_identical(names(init), c("s.num", "e.num", "i.num", "r.num",
                                  "se.flow", "ei.flow", "ir.flow"))

  param <- param.dcm(inf.prob = 0.2, act.rate = 0.5, sx.rate = 0.1,
                     rec.rate = 0.05)
  control <- control.dcm(nsteps = 10, dt = 1, new.mod = SEIR)

  mod <- dcm(param, init, control)
  expect_is(mod, "dcm")
  expect_identical(names(as.data.frame(mod))[3], "e.num")
})


test_that("Non-sensitivity parameter vector", {

  ## SI function
  intSI <- function(t, t0, parms) {
    with(as.list(c(t0, parms)), {

      ## Dynamic Calculations
      # Population size
      num <- s.num + i.num

      if (t < start.time) {
        lambda <- inf.prob[1] * act.rate * i.num / num
      } else {
        lambda <- inf.prob[2] * act.rate * i.num / num
      }

      ## Flows
      si.flow <- lambda * s.num

      ## Differential Equations
      dS <- -si.flow
      dI <- si.flow

      ## Output
      list(c(dS, dI, si.flow),
           num = num)
    })
  }

  param <- param.dcm(inf.prob = c(0.5, 0.05), act.rate = 0.1, start.time = 100)
  init <- init.dcm(s.num = 999, i.num = 1, si.flow = 0)
  control <- control.dcm(nsteps = 250, new.mod = intSI, sens.param = FALSE)
  mod <- dcm(param, init, control)
  expect_is(mod, "dcm")
  expect_true(any(!is.na(as.data.frame(mod))))
})
