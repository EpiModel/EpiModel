context("Netsim Parallel")

test_that("1 sim on 1 core", {
  skip_on_cran()
  nw <- network.initialize(n = 50, directed = FALSE)
  formation <- ~ edges
  target.stats <- 25
  dissolution <- ~ offset(edges)
  duration <- 50
  coef.diss <- dissolution_coefs(dissolution, duration)
  est <- netest(nw, formation, dissolution,
                target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.25)
  init <- init.net(i.num = 50)
  control <- control.net(type = "SI", nsteps = 5, verbose = FALSE,
                         par.type = "single", nsims = 1, ncores = 1)
  sims <- netsim_parallel(est, param, init, control)
  expect_is(sims, "netsim")
  expect_output(sims, "simulations: 1")
})

test_that("2 sims on 1 core", {
  skip_on_cran()
  nw <- network.initialize(n = 50, directed = FALSE)
  formation <- ~ edges
  target.stats <- 25
  dissolution <- ~ offset(edges)
  duration <- 50
  coef.diss <- dissolution_coefs(dissolution, duration)
  est <- netest(nw, formation, dissolution,
                target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.25)
  init <- init.net(i.num = 50)
  control <- control.net(type = "SI", nsteps = 5, verbose = FALSE,
                         par.type = "single", nsims = 2, ncores = 1)
  sims <- netsim_parallel(est, param, init, control)
  expect_is(sims, "netsim")
  expect_output(sims, "simulations: 2")
})

test_that("1 sim on (not really) 2 cores",{
  #skip_on_cran()
  nw <- network.initialize(n = 50, directed = FALSE)
  formation <- ~ edges
  target.stats <- 25
  dissolution <- ~ offset(edges)
  duration <- 50
  coef.diss <- dissolution_coefs(dissolution, duration)
  est <- netest(nw, formation, dissolution,
                target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.25)
  init <- init.net(i.num = 50)
  control <- control.net(type = "SI", nsteps = 5, verbose = FALSE,
                         par.type = "single", nsims = 1, ncores = 2)
  sims <- netsim_parallel(est, param, init, control)
  expect_is(sims, "netsim")
  expect_output(sims, "simulations: 1")
})

test_that("2 sims on 2 cores", {
  skip_on_cran()
  nw <- network.initialize(n = 50, directed = FALSE)
  formation <- ~ edges
  target.stats <- 25
  dissolution <- ~ offset(edges)
  duration <- 50
  coef.diss <- dissolution_coefs(dissolution, duration)
  est <- netest(nw, formation, dissolution,
                target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.25)
  init <- init.net(i.num = 50)
  control <- control.net(type = "SI", nsteps = 5, verbose = FALSE,
                         par.type = "single", nsims = 2, ncores = 2)
  sims <- netsim_parallel(est, param, init, control)
  expect_is(sims, "netsim")
  expect_output(sims, "simulations: 2")
})




