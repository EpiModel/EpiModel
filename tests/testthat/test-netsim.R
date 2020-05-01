context("Standard network models")

nw <- network_initialize(n = 50, directed = FALSE)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(nw, formation = ~edges + nodematch("race"), target.stats = c(25, 10),
              coef.diss = dissolution_coefs(~offset(edges), 10, 0),
              verbose = FALSE)

# Edges + nodematch, one-mode, closed

test_that("netsim for edges only, SI, one-mode, closed, 1 sim", {
  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "network")
  test_net(mod)
})

test_that("netsim for edges only, SI, one-mode, closed, 2 sim", {
  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 2, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "network")
  test_net(mod)
})

test_that("netsim for edges only, SIS, one-mode, closed, 1 sim", {
  param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsims = 1, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "network")
  test_net(mod)
})

test_that("netsim for edges only, SIS, one-mode, closed, 2 sim", {
  param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsims = 2, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
})

test_that("netsim for edges only, SIR, one-mode, closed, 1 sim", {
  param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsims = 1, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "network")
  test_net(mod)
})

test_that("netsim for edges only, SIR, one-mode, closed, 2 sim", {
  param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsims = 2, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "network")
  test_net(mod)
})
