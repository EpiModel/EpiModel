context("Network models with tergmLite (All SOC)")

# tergmLite = TRUE -----

test_that("tergmLite: 1G, Closed", {
  skip_on_cran()

  # SI
  num <- 1000
  nw <- network_initialize(n = num)
  formation <- ~edges
  target.stats <- 400
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
  est <- netest(nw, formation, target.stats, coef.diss)

  param <- param.net(inf.prob = 0.1, act.rate = 5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 20, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = 1)
  capture_output(
    summary(sim, at = 20)
  )

  # SIS
  param <- param.net(inf.prob = 0.5, act.rate = 1, rec.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = 1)
  capture_output(
    summary(sim, at = 10)
  )

  # SIR
  param <- param.net(inf.prob = 0.5, act.rate = 1, rec.rate = 0.02)
  init <- init.net(i.num = 10, r.num = 5)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = 1)
  capture_output(
    summary(sim, at = 10)
  )
})

test_that("tergmLite: 2G, Closed", {
  skip_on_cran()

  # SI
  num1 <- num2 <- 500
  nw <- network_initialize(n = num1 + num2)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = num1))
  formation <- ~edges + nodematch("group")
  target.stats <- c(400, 0)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss)

  # Parameters
  param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2)
  init <- init.net(i.num = 20, i.num.g2 = 20)
  control <- control.net(type = "SI", nsteps = 20, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim)
  capture_output(
    summary(sim, at = 20)
  )

  # SIS
  param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2,
                     rec.rate = 0.02, rec.rate.g2 = 0.02)
  init <- init.net(i.num = 20, i.num.g2 = 20)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim)
  capture_output(
    summary(sim, at = 10)
  )

  # SIR
  param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2,
                     rec.rate = 0.02, rec.rate.g2 = 0.02)
  init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = FALSE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, mean.col = c(1, 2, 3, 1, 2, 3))
  capture_output(
    summary(sim, at = 10)
  )

})


test_that("tergmLite: 1G, Open", {
  skip_on_cran()

  # SI
  num <- 1000
  nw <- network_initialize(n = num)
  formation <- ~edges
  target.stats <- 400
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10,
                                 d.rate = 0.005)
  est <- netest(nw, formation, target.stats, coef.diss)

  # Parameters
  param <- param.net(inf.prob = 0.4, act.rate = 1,
                     a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE,
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, y = "si.flow", ylim = c(0, 20))
  capture_output(
    summary(sim, at = 10)
  )

  # SIS
  param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                     a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, y = c("si.flow", "is.flow"), legend = TRUE)
  capture_output(
    summary(sim, at = 10)
  )

  # SIR
  param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                     a.rate = 0.005, di.rate = 0.005, ds.rate = 0.005,
                     dr.rate = 0.005)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, qnts = 1)
  plot(sim, y = "num", ylim = c(800, 1200))
})


test_that("tergmLite: 2G, Open", {
  skip_on_cran()

  # SI
  num1 <- num2 <- 500
  nw <- network_initialize(n = num1 + num2)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = num1))
  formation <- ~edges + nodematch("group")
  target.stats <- c(400, 0)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25,
                                 d.rate = 0.005)
  est <- netest(nw, formation, target.stats, coef.diss)

  # Parameters
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                     act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                     di.rate = 0.005, ds.rate = 0.005,
                     di.rate.g2 = 0.005, ds.rate.g2 = 0.005)
  init <- init.net(i.num = 50, i.num.g2 = 50)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, qnts = 1, ylim = c(0, 500))

  # SIS
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                     act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                     di.rate = 0.005, ds.rate = 0.005,
                     di.rate.g2 = 0.005, ds.rate.g2 = 0.005,
                     rec.rate = 0.02, rec.rate.g2 = 0.02)
  init <- init.net(i.num = 50, i.num.g2 = 50)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, qnts = 1, ylim = c(0, 500))

  # SIR
  param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                     act.rate = 5, a.rate = 0.005, a.rate.g2 = 0.005,
                     di.rate = 0.005, ds.rate = 0.005,
                     di.rate.g2 = 0.005, ds.rate.g2 = 0.005,
                     dr.rate = 0.005, dr.rate.g2 = 0.005,
                     rec.rate = 0.005, rec.rate.g2 = 0.005)
  init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
  plot(sim, qnts = FALSE, sim.lines = TRUE)
  plot(sim, qnts = 1, ylim = c(0, 500))
})


test_that("Models with duration = 1", {
  skip_on_cran()

  num <- 1000
  nw <- network_initialize(n = num)
  formation <- ~edges
  target.stats <- 400
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
  est <- netest(nw, formation, target.stats, coef.diss)

  param <- param.net(inf.prob = 0.1, act.rate = 5)
  init <- init.net(i.num = 10)

  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         tergmLite = FALSE, resimulate.network = TRUE, verbose = FALSE)
  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")

  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         tergmLite = TRUE, resimulate.network = TRUE, verbose = FALSE)
  sim <- netsim(est, param, init, control)
  expect_s3_class(sim, "netsim")
})

test_that("edgecov works with tergmLite", {
  nw <- network_initialize(n = 20)
  nw %v% "race" <- rep(c("a", "b"), length.out = 20)

  # create matrix for edgecov term
  ec <- matrix(runif(400), 20, 20)
  ec <- (ec + t(ec))/3
  nw %n% "ec" <- ec

  fit <- netest(nw,
                formation = ~edges + nodefactor("race") + edgecov("ec") ,
                target.stats = c(30, 15, 10),
                coef.diss = dissolution_coefs(dissolution = ~offset(edges), duration = 1),
                edapprox = TRUE)

  coefsign <- sign(fit$coef.form.crude[3])

  resim_fun <- function(dat, at) {
    n <- attr(dat$el[[1]], "n")
    m <- matrix(if(at == 3) +Inf else -Inf, n, n)*coefsign
    attr(dat$el[[1]], "ec") <- m
    resim_nets(dat, at)
  }

  param <- param.net(inf.prob = 0.3, act.rate = 1)
  init <- init.net(i.num = 5, r.num = 0)
  control <- control.net(type = NULL, nsteps = 4, nsims = 1, ncores = 1,
                         resim_nets.FUN = resim_fun, tergmLite = TRUE,
                         resimulate.network = TRUE, nwstats.formula = ~edges,
                         set.control.ergm = control.simulate.formula(MCMC.burnin = 1e5),
                         verbose = FALSE)
  sim <- netsim(fit, param, init, control)
  expect_s3_class(sim, "netsim")

  expect_equal(unname(sim$stats$nwstats$sim1[[1]][2,1]), 0)
  expect_equal(unname(sim$stats$nwstats$sim1[[1]][3,1]), network.dyadcount(nw))
  expect_equal(unname(sim$stats$nwstats$sim1[[1]][4,1]), 0)
})
