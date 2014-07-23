context("net")


# netest ------------------------------------------------------------------

test_that("netest works for edges only model", {
  nw <- network.initialize(n = 50, directed = FALSE)
  est <- netest(
    nw,
    formation = ~ edges,
    dissolution = ~offset(edges),
    target.stats = 25,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE)
    expect_is(est, "netest")
})

test_that("netest works for edges + nodematch model", {
  nw <- network.initialize(n = 50, directed = FALSE)
  nw <- set.vertex.attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(
    nw,
    formation = ~ edges + nodematch("race"),
    dissolution = ~offset(edges),
    target.stats = c(25, 10),
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE)
  expect_is(est, "netest")
})



# netsim ------------------------------------------------------------------

nw <- network.initialize(n = 50, directed = FALSE)
nw <- set.vertex.attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(
  nw,
  formation = ~ edges + nodematch("race"),
  dissolution = ~offset(edges),
  target.stats = c(25, 10),
  coef.diss = dissolution_coefs(~offset(edges), 10, 0),
  verbose = FALSE)

# Edges + nodematch, one-mode, closed

test_that("netsim for edges only, SI, one-mode, closed, 1 sim", {
  param <- param.net(trans.rate = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
})

test_that("netsim for edges only, SI, one-mode, closed, 2 sim", {
  param <- param.net(trans.rate = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 2, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
})

test_that("netsim for edges only, SIS, one-mode, closed, 1 sim", {
  param <- param.net(trans.rate = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsims = 1, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
})

test_that("netsim for edges only, SIS, one-mode, closed, 2 sim", {
  param <- param.net(trans.rate = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsims = 2, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
})

test_that("netsim for edges only, SIR, one-mode, closed, 1 sim", {
  param <- param.net(trans.rate = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsims = 1, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
})

test_that("netsim for edges only, SIR, one-mode, closed, 2 sim", {
  param <- param.net(trans.rate = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsims = 2, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
})



# Cases -------------------------------------------------------------------


test_that("Updating attributes in open populations", {
  nw <- network.initialize(n = 50, bipartite = 25, directed = FALSE)
  nw <- set.vertex.attribute(nw, attrname = "risk",
                             value = rbinom(50, 1, 0.5))

  formation <- ~edges + nodefactor("risk")
  target.stats <- c(25, 36)
  dissolution <- ~offset(edges)
  coef.diss <- dissolution_coefs(dissolution, 38, d.rate = 0.002)

  est1 <- netest(nw,
                 formation,
                 dissolution,
                 target.stats,
                 coef.diss,
                 verbose = FALSE)

  probs <- c(0.2055, 0.0088, 0.0614, 0)
  durs <- c(3, 100, 9, 10)
  trans.rates <- rep(probs, durs)
  trans.ratesf <- trans.rates*2
  param <- param.net(trans.rate = trans.ratesf, act.rate = 1,
                     trans.rate.m2 = trans.rates,
                     b.rate = 0.05, b.rate.m2 = NA,
                     ds.rate = 0.002, ds.rate.m2 = 0.002,
                     di.rate = 0.008, di.rate.m2 = 0.008)

  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         epi.by = "risk", verbose = FALSE)

  sim1 <- netsim(est1, param, init, control)
  expect_is(sim1, "netsim")
})
