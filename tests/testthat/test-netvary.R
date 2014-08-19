## Testing time-varying trans.rate and act.rate in network models

context("Time Varying Network Parameters")

test_that("time varying trans.rate and act.rate for one-mode", {
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~ edges
  target.stats <- 50
  dissolution <- ~ offset(edges)
  duration <- 20
  coef.diss <- dissolution_coefs(dissolution, duration)

  est1 <- netest(nw,
                 formation,
                 dissolution,
                 target.stats,
                 coef.diss,
                 verbose = FALSE)


  probs <- c(0.25, 0.02)
  durs <- c(10, 1)
  trans.rates <- rep(probs, durs)

  acts <- c(0.5, 2)
  act.rates <- rep(acts, durs)

  # Parameters, initial conditions, and controls for model
  param <- param.net(trans.rate = trans.rates, act.rate = act.rates)
  init <- init.net(i.num = 50)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         save.transmat = TRUE, verbose = FALSE)
  mod <- netsim(est1, param, init, control)
  expect_is(mod, "netsim")
})


test_that("time varying trans.rate and act.rate for bipartite", {
  nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
  formation <- ~ edges
  target.stats <- 50
  dissolution <- ~ offset(edges)
  duration <- 20
  coef.diss <- dissolution_coefs(dissolution, duration)

  est1 <- netest(nw,
                 formation,
                 dissolution,
                 target.stats,
                 coef.diss,
                 verbose = FALSE)


  probs <- c(0.25, 0.01)
  durs <- c(10, 1)
  trans.rates <- rep(probs, durs)
  trans.rates.m2 <- trans.rates*2

  acts <- c(1, 2)
  act.rates <- rep(acts, durs)

  param <- param.net(trans.rate = trans.rates,
                     trans.rate.m2 = trans.rates.m2,
                     act.rate = act.rates)
  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         save.transmat = TRUE, verbose = FALSE)

  mod <- netsim(est1, param, init, control)
  expect_is(mod, "netsim")
})
