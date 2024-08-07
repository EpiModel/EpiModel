context("Network attributes with arrivals")

test_that("Updating attributes in open populations", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, attrname = "group", rep(1:2, each = 25))

  formation <- ~edges + nodefactor("group")
  target.stats <- c(25, 36)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 38,
                                 d.rate = 0.002)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  probs <- c(0.2055, 0.0088, 0.0614, 0)
  durs <- c(3, 100, 9, 10)
  inf.probs <- rep(probs, durs)
  inf.probsf <- inf.probs * 2
  param <- param.net(inf.prob = inf.probs, act.rate = 1,
                     inf.prob.g2 = inf.probs,
                     a.rate = 0.05, a.rate.g2 = NA,
                     ds.rate = 0.05, ds.rate.g2 = 0.05,
                     di.rate = 0.05, di.rate.g2 = 0.05)

  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SI", nsteps = 20, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)

  sim1 <- netsim(est1, param, init, control)
  expect_is(sim1, "netsim")
})

test_that("SIR model with epi.by parameter", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, attrname = "race", rep(0:1, each = 25))
  formation <- ~edges + nodefactor("race")
  target.stats <- c(25, 25)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 50)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.1, act.rate = 1, rec.rate = 0.005)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1,
                         epi.by = "race", verbose = FALSE, verbose.int = 0)
  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")
  expect_true(all(c("s.num.race0", "s.num.race1", "i.num.race0", "i.num.race1",
                    "r.num.race0", "r.num.race1") %in% names(sim$epi)))
})


test_that("Serosorting model in open population", {
  skip_on_cran()
  n <- 100
  nw <- network_initialize(n = n)

  prev <- 0.2
  infIds <- sample(1:n, n * prev)
  nw <- set_vertex_attribute(nw, "status", "s")
  nw <- set_vertex_attribute(nw, "status", "i", infIds)
  nw <- set_vertex_attribute(nw, "race", rbinom(n, 1, 0.5))

  formation <- ~edges + nodefactor("status", levels = -1) +
                nodematch("status") + nodematch("race")
  target.stats <- c(36, 55, 25, 18)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 5, d.rate = 0.01)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.8, a.rate = 0.05,
                     ds.rate = 0.01, di.rate = 0.01)
  init <- init.net()
  control <- control.net(type = "SI", nsteps = 20, nsims = 1,
                         nwstats.formula = ~edges +
                                            meandeg +
                                            nodefactor("status",
                                                       levels = NULL) +
                                            nodematch("status"),
                         tergmLite = FALSE,
                         resimulate.network = TRUE,
                         save.run = TRUE,
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")

  nD <- get_network(sim)
  tea1 <- get.vertex.attribute.active(nD, "testatus", at = 1)
  expect_true(sum(!is.na(tea1)) == n)

  tea20 <- get.vertex.attribute.active(nD, "testatus", at = 20)
  expect_true(sum(is.na(tea20)) == 0)

  fstat.nw <- get_vertex_attribute(nD, "status")
  fstat.attr <- sim$run[[1]]$attr$status

  expect_identical(tea20, fstat.nw)
  expect_identical(fstat.nw, fstat.attr)

})

test_that("Serosorting model in closed population", {
  skip_on_cran()
  n <- 100
  nw <- network_initialize(n = n)

  prev <- 0.2
  infIds <- sample(1:n, n * prev)
  nw <- set_vertex_attribute(nw, "status", "s")
  nw <- set_vertex_attribute(nw, "status", "i", infIds)
  nw <- set_vertex_attribute(nw, "race", rbinom(n, 1, 0.5))

  formation <- ~edges + nodefactor("status", levels = -1) +
    nodematch("status") + nodematch("race")
  target.stats <- c(36, 55, 25, 18)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 5)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.8)
  init <- init.net()
  control <- control.net(type = "SI", nsteps = 20, nsims = 1,
                         nwstats.formula = ~edges +
                           meandeg +
                           nodefactor("status", levels = NULL) +
                           nodematch("status"),
                         tergmLite = FALSE,
                         resimulate.network = TRUE,
                         save.run = TRUE,
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")

  nD <- get_network(sim)
  nD
  tea1 <- get.vertex.attribute.active(nD, "testatus", at = 1)

  expect_true(sum(!is.na(tea1)) == n)

  tea20 <- get.vertex.attribute.active(nD, "testatus", at = 20)
  expect_true(sum(is.na(tea20)) == 0)

  fstat.nw <- get_vertex_attribute(nD, "status")
  fstat.attr <- sim$run[[1]]$attr$status

  expect_identical(tea20, fstat.nw)
  expect_identical(fstat.nw, fstat.attr)

})


test_that("Serosorting model in open population, with tergmLite", {
  skip_on_cran()
  n <- 100
  nw <- network_initialize(n = n)

  prev <- 0.2
  infIds <- sample(1:n, n * prev)
  nw <- set_vertex_attribute(nw, "status", "s")
  nw <- set_vertex_attribute(nw, "status", "i", infIds)
  nw <- set_vertex_attribute(nw, "race", rbinom(n, 1, 0.5))

  formation <- ~edges + nodefactor("status", levels = -1) +
    nodematch("status") + nodematch("race")
  target.stats <- c(36, 55, 25, 18)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 5, d.rate = 0.01)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.8, a.rate = 0.05,
                     ds.rate = 0.01, di.rate = 0.01)
  init <- init.net()
  control <- control.net(type = "SI", nsteps = 20, nsims = 1,
                         nwstats.formula = ~edges +
                           meandeg +
                           nodefactor("status", levels = NULL) +
                           nodematch("status"),
                         tergmLite = TRUE,
                         resimulate.network = TRUE,
                         save.run = TRUE,
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")

})

test_that("Save attributes to output", {
  skip_on_cran()
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 25))
  formation <- ~edges + nodematch("group")
  target.stats <- c(25, 0)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 38,
                                 d.rate = 0.01)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.2, act.rate = 1,
                     inf.prob.g2 = 0.2,
                     a.rate = 0.01, a.rate.g2 = NA,
                     ds.rate = 0.01, ds.rate.g2 = 0.01,
                     di.rate = 0.01, di.rate.g2 = 0.01)

  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2,
                         save.run = TRUE, resimulate.network = TRUE,
                         verbose = FALSE)

  sim1 <- netsim(est1, param, init, control)
  expect_is(sim1, "netsim")
  expect_is(sim1$run[[1]]$attr, "list")
  expect_true(all(c("entrTime", "exitTime") %in% names(sim1$run[[1]]$attr)))
})


test_that("Check TE Status Variable Against Epi Stats", {

  skip_on_cran()

  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 38)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # SIR
  param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.01)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsims = 1, nsteps = 100, verbose = FALSE)
  sim <- netsim(est, param, init, control)
  times <- sample(1:100, 10)
  for (at in times) {
    df <- as.data.frame(sim)[at, ]
    nwd <- get_network(sim, collapse = TRUE, at = at)
    attr <- get_vertex_attribute(nwd, "testatus")
    expect_true(sum(attr == "s") == df$s.num)
    expect_true(sum(attr == "i") == df$i.num)
    expect_true(sum(attr == "r") == df$r.num)
  }

  # SIS
  param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.01)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsims = 1, nsteps = 100, verbose = FALSE)
  sim <- netsim(est, param, init, control)
  times <- sample(1:100, 10)
  for (at in times) {
    df <- as.data.frame(sim)[at, ]
    nwd <- get_network(sim, collapse = TRUE, at = at)
    attr <- get_vertex_attribute(nwd, "testatus")
    expect_true(sum(attr == "s") == df$s.num)
    expect_true(sum(attr == "i") == df$i.num)
  }

})
