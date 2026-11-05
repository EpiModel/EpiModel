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

test_that("worker_globals() finds what user code needs from this session", {
  # Functions and formulas from a script look up their objects in the global
  # environment and on the search path, which multisession workers do not
  # share (#1104)
  genv <- globalenv()
  skip_if(exists("edges", envir = genv, inherits = FALSE))
  assign("epimodel_test_lv", 2, envir = genv)
  assign("epimodel_test_rate", function() 0.5, envir = genv)
  assign("edges", 100, envir = genv)
  on.exit(rm(list = c("epimodel_test_lv", "epimodel_test_rate", "edges"),
             envir = genv))

  # EpiModel and stats are attached on the workers, testthat is not, and
  # qualified calls and locals are fine
  my_mod <- function(dat, at) {
    expect_true(TRUE)
    testthat::expect_false(FALSE)
    lv <- rbinom(1, 1, 0.5)
    set_epi(dat, "rate", epimodel_test_rate() + lv)
  }
  environment(my_mod) <- genv

  # a closure carries its own objects
  carry_env <- new.env(parent = genv)
  carry_env$epimodel_test_lv <- 3
  carried_mod <- function(dat, at) epimodel_test_lv
  environment(carried_mod) <- carry_env

  # term names are not R objects, term arguments are; trim_env() carries them
  formation <- ~edges + nodefactor("grp", levels = epimodel_test_lv) +
    offset(edges)
  environment(formation) <- genv
  carried <- trim_env(formation, keep = "epimodel_test_lv")
  nwstats <- ~edges + offset(nodematch("grp"))
  environment(nwstats) <- genv

  found <- worker_globals(list(
    x = list(formation = formation, carried = carried),
    nwstats.formula = nwstats,
    control = list(my.FUN = my_mod, carried.FUN = carried_mod)
  ))
  expect_named(found, c("x$formation", "control$my.FUN"))
  expect_identical(found[["x$formation"]], "epimodel_test_lv")
  expect_setequal(found[["control$my.FUN"]],
                  c("epimodel_test_rate()", "expect_true() from testthat"))

  # through the trimmed formula environment of a fitted model
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "grp", rep(1:2, 25))
  formation <- ~edges + nodefactor("grp", levels = epimodel_test_lv)
  environment(formation) <- genv
  coef.diss <- dissolution_coefs(~offset(edges), 10, 0)
  for (nested in c(TRUE, FALSE)) {
    est <- netest(nw, formation, target.stats = c(24, 12), coef.diss,
                  nested.edapprox = nested, verbose = FALSE)
    expect_identical(worker_globals(list(x = est)),
                     list("x$formation" = "epimodel_test_lv"))
  }
})

test_that("netsim() and netdx() warn about what multisession workers lack", {
  skip_on_cran()

  genv <- globalenv()
  assign("epimodel_test_lv", 2, envir = genv)
  assign("epimodel_test_rate", function() 0.5, envir = genv)
  on.exit(rm(list = c("epimodel_test_lv", "epimodel_test_rate"), envir = genv))

  # a global helper in a branch that never runs, so the simulation completes
  my_mod <- function(dat, at) {
    if (at < 0) dat <- set_epi(dat, "rate", epimodel_test_rate())
    dat
  }
  environment(my_mod) <- genv

  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "grp", rep(1:2, 25))
  coef.diss <- dissolution_coefs(~offset(edges), 10, 0)
  est <- netest(nw, ~edges, target.stats = 24, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)

  control <- control.net(type = NULL, nsims = 2, nsteps = 3, ncores = 2,
                         infection.FUN = infection.net, my.FUN = my_mod,
                         verbose = FALSE)
  expect_warning(mod <- netsim(est, param, init, control),
                 "control\\$my.FUN: epimodel_test_rate\\(\\)")
  expect_s3_class(mod, "netsim")

  control <- control.net(type = NULL, nsims = 2, nsteps = 3, ncores = 1,
                         infection.FUN = infection.net, my.FUN = my_mod,
                         verbose = FALSE)
  expect_no_warning(netsim(est, param, init, control),
                    message = "separate R sessions")

  # nothing is exported, so a formula that needs a global object fails on the
  # workers, and one that carries it with trim_env() runs
  formation <- ~edges + nodefactor("grp", levels = epimodel_test_lv)
  environment(formation) <- genv
  est <- netest(nw, formation, target.stats = c(24, 12), coef.diss,
                verbose = FALSE)
  expect_error(
    expect_warning(netdx(est, nsims = 2, nsteps = 3, ncores = 2, verbose = FALSE),
                   "x\\$formation: epimodel_test_lv"),
    "epimodel_test_lv"
  )

  formation <- trim_env(formation, keep = "epimodel_test_lv")
  est <- netest(nw, formation, target.stats = c(24, 12), coef.diss,
                verbose = FALSE)
  expect_no_warning(
    dx <- netdx(est, nsims = 2, nsteps = 3, ncores = 2, verbose = FALSE),
    message = "separate R sessions"
  )
  expect_s3_class(dx, "netdx")
})
