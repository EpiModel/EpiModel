context("Netsim Parallel")

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
  control <- control.net(type = "SI", nsims = 2, ncores = 2, nsteps = 5,
                         verbose = FALSE, future.use.plan = TRUE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  expect_true(mod$control$nsims == 2)
  expect_true(mod$control$ncores == 2)
  expect_true(inherits(future::plan(), "sequential"))
})

test_that("netsim par, 2 sims, 2 cores - future.use.plan", {
  skip_on_cran()
  skip_on_os("windows")
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 24,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 2, nsteps = 5,
                         verbose = FALSE, future.use.plan = TRUE)
  with(future::plan("multicore", workers = 2), local = TRUE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")
  expect_true(mod$control$nsims == 2)

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

test_that("user code runs with its global objects on multisession workers", {
  skip_on_cran()

  # Functions and formulas from a script look up their objects in the global
  # environment, which multisession workers do not share (#1104)
  genv <- globalenv()
  assign("epimodel_test_lv", 2, envir = genv)
  assign("epimodel_test_rate", function() 0.5, envir = genv)
  on.exit(rm(list = c("epimodel_test_lv", "epimodel_test_rate"), envir = genv))

  my_mod <- function(dat, at) {
    EpiModel::set_epi(dat, "rate", epimodel_test_rate())
  }
  environment(my_mod) <- genv
  formation <- ~edges + nodefactor("grp", levels = epimodel_test_lv)
  environment(formation) <- genv

  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "grp", rep(1:2, 25))
  coef.diss <- dissolution_coefs(~offset(edges), 10, 0)
  est <- netest(nw, formation, target.stats = c(24, 12), coef.diss,
                verbose = FALSE)
  est_nn <- netest(nw, formation, target.stats = c(24, 12), coef.diss,
                   nested.edapprox = FALSE, verbose = FALSE)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 2, nsteps = 3, ncores = 2,
                         infection.FUN = infection.net, my.FUN = my_mod,
                         resimulate.network = TRUE, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_equal(unname(unlist(mod$epi$rate[3, ])), c(0.5, 0.5))

  mod <- netsim(est_nn, param, init, control)
  expect_s3_class(mod, "netsim")

  dx <- netdx(est_nn, nsims = 2, nsteps = 3, ncores = 2, verbose = FALSE)
  expect_s3_class(dx, "netdx")
})
