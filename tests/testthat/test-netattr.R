context("Network attributes with arrivals")

test_that("Updating attributes in open populations", {
  nw <- network.initialize(n = 50, directed = FALSE)
  nw <- set.vertex.attribute(nw, attrname = "group",
                             value = rbinom(50, 1, 0.5)+1) #FLAG

  formation <- ~edges + nodefactor("group")
  target.stats <- c(25, 36)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 38, d.rate = 0.002)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  probs <- c(0.2055, 0.0088, 0.0614, 0)
  durs <- c(3, 100, 9, 10)
  inf.probs <- rep(probs, durs)
  inf.probsf <- inf.probs*2
  param <- param.net(inf.prob = inf.probs, act.rate = 1,
                     inf.prob.g2 = inf.probs,
                     a.rate = 0.05, a.rate.g2 = NA, 
                     ds.rate = 0.002, ds.rate.g2 = 0.002,
                     di.rate = 0.008, di.rate.g2 = 0.008)

  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         epi.by = "group", verbose = FALSE)

  sim1 <- netsim(est1, param, init, control)
  expect_is(sim1, "netsim")
})

test_that("SIR model with epi.by parameter", {
  skip_on_cran()
  nw <- network.initialize(n = 50, directed = FALSE)
  nw <- set.vertex.attribute(nw, attrname = "race", rep(0:1, each = 25))
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
  n <- 100
  nw <- network.initialize(n, directed = FALSE)

  prev <- 0.2
  infIds <- sample(1:n, n*prev)
  nw <- set.vertex.attribute(nw, "status", "s")
  nw <- set.vertex.attribute(nw, "status", "i", infIds)
  nw <- set.vertex.attribute(nw, "race", rbinom(n, 1, 0.5))

  formation <- ~edges + nodefactor("status", base = 1) +
                nodematch("status") + nodematch("race")
  target.stats <- c(36, 55, 25, 18)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 50, d.rate = 0.01)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.03, a.rate = 0.01,
                     ds.rate = 0.01, di.rate = 0.01)
  init <- init.net()
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         nwstats.formula = ~edges +
                                            meandeg +
                                            nodefactor("status", base = 0) +
                                            nodematch("status"),
                         save.network = FALSE,
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")
})


test_that("Save attributes to output", {
  skip_on_cran()
  nw <- network.initialize(n = 50, directed = FALSE)
  nw <- set.vertex.attribute(nw, "group", rep(c(1,2), each = 25))
  formation <- ~edges
  target.stats <- 25
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 38, d.rate = 0.01)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.2, act.rate = 1,
                     inf.prob.g2 = 0.2,
                     a.rate = 0.01, a.rate.g2 = NA,
                     ds.rate = 0.01, ds.rate.g2 = 0.01,
                     di.rate = 0.01, di.rate.g2 = 0.01)

  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2,
                         save.other = "attr", verbose = FALSE)

  sim1 <- netsim(est1, param, init, control)
  expect_is(sim1, "netsim")
  expect_is(sim1$attr[[1]], "list")
  expect_true(all(c("entrTime", "exitTime") %in% names(sim1$attr[[1]])))
})
