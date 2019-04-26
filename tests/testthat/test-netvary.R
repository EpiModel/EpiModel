context("Time Varying Network Parameters")

test_that("time varying parameters for one-mode", {
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration  = 20)

  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  probs <- c(0.25, 0.02)
  durs <- c(10, 1)
  inf.probs <- rep(probs, durs)

  acts <- c(0.5, 2)
  act.rates <- rep(acts, durs)

  rec.rates <- rep(0:1, c(20, 1))

  # Parameters, initial conditions, and controls for model
  param <- param.net(inf.prob = inf.probs, act.rate = act.rates,
                     rec.rate = rec.rates)
  init <- init.net(i.num = 50)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1,
                         save.transmat = TRUE, verbose = FALSE)
  mod <- netsim(est1, param, init, control)
  expect_is(mod, "netsim")
})


test_that("time varying parameters for two-group models", {
  nw <- network.initialize(n = 100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "group", rep(c(1,2), each = 50))
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration  = 20)

  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  probs <- c(0.25, 0.01)
  durs <- c(10, 1)
  inf.probs <- rep(probs, durs)
  inf.probs.g2 <- inf.probs*2

  acts <- c(1, 2)
  act.rates <- rep(acts, durs)

  rec.rates <- rep(0:1, c(20, 1))
  rec.rates.g2 <- rep(0:1, c(10, 1))

  param <- param.net(inf.prob = inf.probs,
                     inf.prob.g2 = inf.probs.g2,
                     act.rate = act.rates,
                     rec.rate = rec.rates,
                     rec.rate.g2 = rec.rates.g2)
  init <- init.net(i.num = 10, i.num.g2 = 10,
                   r.num = 0, r.num.g2 = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1,
                         save.transmat = TRUE, verbose = FALSE)

  mod <- netsim(est1, param, init, control)
  expect_is(mod, "netsim")
})
