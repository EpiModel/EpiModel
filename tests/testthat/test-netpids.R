context("Use PIDs for Network Simulations")

test_that("PIDs for one-mode models", {

  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50,
                                 d.rate = 0.01)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.1, rec.rate = 0.02,
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01)
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


test_that("PIDs for two-group models", {

  nw <- network.initialize(n = 100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "group", rep(1:2, each = 50))
  formation <- ~edges + nodematch("group")
  target.stats <- c(50, 0)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50,
                                 d.rate = 0.01)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.1,
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01,
                     a.rate.g2 = 0.01, ds.rate.g2 = 0.01, di.rate.g2 = 0.01)
  init <- init.net(i.num = 10, i.num.g2 = 10)

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
