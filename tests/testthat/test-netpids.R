context("Use PIDs for Network Simulations")

test_that("PIDs for one-mode models", {

  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50,
                                 d.rate = 0.01)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.1, rec.rate = 0.02,
                     b.rate = 0.01, ds.rate = 0.01, di.rate = 0.01)
  init <- init.net(i.num = 10)

  # Use pids
  control <- control.net(type = "SIS", nsteps = 10, verbose = FALSE,
                         save.network = TRUE, nsims = 1)
  sim <- netsim(est, param, init, control)
  expect_true(sim$network$sim1$gal$vertex.pid == "vertex.pid")

  # Do not use pids
  control <- control.net(type = "SIS", nsteps = 10, verbose = FALSE,
                         save.network = TRUE, nsims = 1, use.pids = FALSE)
  sim <- netsim(est, param, init, control)
  expect_true(sim$network$sim1$gal$vertex.pid == "tergm_pid")

})


test_that("PIDs for bipartite models", {

  nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50,
                                 d.rate = 0.01)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.1, inf.prob.m2 = 0.1,
                     b.rate = 0.01, ds.rate = 0.01, di.rate = 0.01,
                     b.rate.m2 = 0.01, ds.rate.m2 = 0.01, di.rate.m2 = 0.01)
  init <- init.net(i.num = 10, i.num.m2 = 10)

  # Use pids
  control <- control.net(type = "SI", nsteps = 10, verbose = FALSE,
                         save.network = TRUE, nsims = 1)
  simb <- netsim(est, param, init, control)
  expect_true(simb$network$sim1$gal$vertex.pid == "vertex.names")

  # Do not use pids
  control <- control.net(type = "SI", nsteps = 10, verbose = FALSE,
                         save.network = TRUE, nsims = 1, use.pids = FALSE)
  simb <- netsim(est, param, init, control)
  expect_true(simb$network$sim1$gal$vertex.pid == "tergm_pid")

})
