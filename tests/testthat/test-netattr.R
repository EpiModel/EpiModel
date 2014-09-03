context("Network attributes with births")

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
  inf.probs <- rep(probs, durs)
  inf.probsf <- inf.probs*2
  param <- param.net(inf.prob = inf.probsf, act.rate = 1,
                     inf.prob.m2 = inf.probs,
                     b.rate = 0.05, b.rate.m2 = NA,
                     ds.rate = 0.002, ds.rate.m2 = 0.002,
                     di.rate = 0.008, di.rate.m2 = 0.008)

  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         epi.by = "risk", verbose = FALSE)

  sim1 <- netsim(est1, param, init, control)
  expect_is(sim1, "netsim")
})


test_that("Serosorting model in open population", {
  n <- 50
  nw <- network.initialize(n, directed = FALSE)

  prev <- 0.2
  infIds <- sample(1:n, n*prev)
  nw <- set.vertex.attribute(nw, "status", "s")
  nw <- set.vertex.attribute(nw, "status", "i", infIds)
  nw <- set.vertex.attribute(nw, "race", rbinom(n, 1, 0.5))

  formation <- ~ edges + nodefactor("status", base = 1) +
                 nodematch("status") + nodematch("race")
  target.stats <- c(18, 3, 15, 10)

  dissolution <- ~offset(edges)
  coef.diss <- dissolution_coefs(dissolution, 50, d.rate = 0.01)

  est <- netest(nw,
                formation,
                dissolution,
                target.stats,
                coef.diss,
                verbose = FALSE)

  param <- param.net(inf.prob = 0.03, b.rate = 0.01,
                     ds.rate = 0.01, di.rate = 0.01)
  init <- init.net()
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         nwstats.formula = ~ edges +
                                             meandeg +
                                             nodefactor("status", base = 0) +
                                             nodematch("status"),
                         save.network = FALSE,
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")
})
