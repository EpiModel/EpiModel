context("Standard network models")

nw <- network_initialize(n = 50)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(nw, formation = ~edges + nodematch("race"),
              target.stats = c(25, 10),
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

test_that("netsim implicit save.network option", {
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # Epidemic model
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  mod1 <- netsim(est1, param, init, control)
  expect_s3_class(get_network(mod1), "networkDynamic")

  control <- control.net(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE,
                         save.network = FALSE)
  mod2 <- netsim(est1, param, init, control)
  expect_error(get_network(mod2))
})

test_that("netsim for edges only, SIR, one-mode, closed, 2 sim, set.control.stergm", {
  param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(i.num = 10, r.num = 0)
  expect_warning(control <- control.net(type = "SIR", nsims = 2, nsteps = 5, verbose = FALSE,
                                        set.control.stergm = control.simulate.network()),
                 "set.control.stergm is deprecated")
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "network")
  test_net(mod)
})

test_that("netsim duration 1", {
  estd1 <- netest(nw, formation = ~edges + nodematch("race"),
                  target.stats = c(25, 10),
                  coef.diss = dissolution_coefs(~offset(edges), 1, 0),
                  verbose = FALSE)
  param <- param.net(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.05)
  init <- init.net(r.num = 0, status.vector = rep("s", 50))
  control <- control.net(type = "SIR", nsims = 1, nsteps = 5, 
                         resimulate.network = TRUE, verbose = FALSE,
                         nwupdate.FUN = NULL)
  set.seed(0)
  mod <- netsim(estd1, param, init, control)
  expect_is(mod, "netsim")
  plot(mod)
  plot(mod, type = "formation")
  plot(mod, type = "network")
  test_net(mod)
  
  # compare to manually produced networkDynamic
  set.seed(0)
  sim <- simulate(estd1$formation,
                  coef = coef(estd1$fit),
                  basis = estd1$fit$newnetwork,
                  control = control.simulate.formula(MCMC.burnin = 2e5),
                  dynamic = FALSE)
  for(i in 1:5) {
    suppressWarnings(sim <- simulate(estd1$formation,
                                     basis = sim,
                                     time.slices = 1,
                                     time.start = i,
                                     time.offset = 0,
                                     coef = c(estd1$coef.form),
                                     dynamic = TRUE))
  }
  expect_identical(sim$mel, mod$network$sim1[[1]]$mel)  
})
