context("test netsim parallel")

test_that("netsim par, 1 sim, 1 core", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 24,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, ncores = 1,
                         verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  expect_true(mod$control$nsims == 1)
  expect_true(mod$control$ncores == 1)

})


test_that("netsim par, 2 sims, 2 cores", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 24,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 2, nsteps = 5, ncores = 2,
                         verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  expect_true(mod$control$nsims == 2)
  expect_true(mod$control$ncores == 2)

})

test_that("netsim par, 1 sim, 2 cores", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 24,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, ncores = 2,
                         verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  expect_true(mod$control$nsims == 1)
  expect_true(mod$control$ncores == 1)

})

test_that("netsim par, 2 sims, 1 cores", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 24,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 2, nsteps = 5, ncores = 1,
                         verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  expect_true(mod$control$nsims == 2)
  expect_true(mod$control$ncores == 1)

})
