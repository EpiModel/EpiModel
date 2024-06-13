context("Network model restart")

test_that("network models can be restarted", {
  skip_on_cran()

  nw <- network_initialize(n = 100)
  est.vit <- netest(nw, formation = ~edges, target.stats = 25,
                    coef.diss = dissolution_coefs(~offset(edges), 10, 0.02),
                    verbose = FALSE)

  param <- param.net(inf.prob = 0.5, act.rate = 2, a.rate = 0.02,
                     ds.rate = 0.02, di.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE,
                         save.run = TRUE,
                         save.other = c())
  x <- netsim(est.vit, param, init, control)

  control <- control.net(type = "SI", nsteps = 10, start = 6,
                         nsims = 1, verbose = FALSE)
  x2 <- netsim(x, param, init, control)

  expect_is(x, "netsim")
  expect_is(x2, "netsim")
  expect_true(x$control$nsteps == 5)
  expect_true(x2$control$nsteps == 10)

  plot(x)
  plot(x2)

})


test_that("restart error flags", {
  skip_on_cran()

  nw <- network_initialize(n = 100)
  est.vit <- netest(nw, formation = ~edges, target.stats = 25,
                    coef.diss = dissolution_coefs(~offset(edges), 10, 0.02),
                    verbose = FALSE)

  param <- param.net(inf.prob = 0.5, act.rate = 2, a.rate = 0.02,
                     ds.rate = 0.02, di.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5,
                         nsims = 1, resimulate.network = TRUE,
                         verbose = FALSE,
                         save.run = TRUE)
  x <- netsim(est.vit, param, init, control)

  control <- control.net(type = "SI", nsteps = 5, start = 10,
                         nsims = 1, verbose = FALSE)
  expect_error(netsim(x, param, init, control),
               "control setting nsteps must be >")

  control <- control.net(type = "SI", nsteps = 10, start = 7,
                         nsims = 1, verbose = FALSE)
  expect_error(netsim(x, param, init, control),
               "control setting start must be 1")

  control <- control.net(type = "SI", nsteps = 10, start = 6,
                         nsims = 1, verbose = FALSE)
  x$run <- NULL
  expect_error(netsim(x, param, init, control), "x must contain `run` to restart simulation, see `save.run` control setting")

})

test_that("reinitialization works with open population, nwterms, and epi.by", {
  nw <- network_initialize(n = 50)
  nw %v% "race" <- rep(0:1, length.out = 50)
  est <- netest(nw, formation = ~edges + nodematch("race"),
                target.stats = c(25, 15),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0.05),
                verbose = FALSE)

  param <- param.net(inf.prob = 0.5, act.rate = 2, a.rate = 0.05,
                     ds.rate = 0.05, di.rate = 0.05)

  init <- init.net(i.num = 10)

  for (tergmLite in c(FALSE, TRUE)) {
    control <- control.net(type = "SI", nsteps = 5,
                           nsims = 2, resimulate.network = TRUE,
                           verbose = FALSE, tergmLite = tergmLite,
                           epi.by = "race",
                           save.run = TRUE,
                           save.other = c(if (tergmLite) c("el", "net_attr")))

    x <- netsim(est, param, init, control)
    expect_is(x, "netsim")
    control$start <- 6
    control$nsteps <- 11
    y <- netsim(x, param, init, control)
    expect_is(y, "netsim")
  }
})
