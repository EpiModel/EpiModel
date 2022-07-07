context("Stochastic model merging")


# merge.icm ---------------------------------------------------------------

test_that("merge for ICM", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 3, verbose = FALSE)
  x <- icm(param, init, control)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  y <- icm(param, init, control)
  z <- merge(x, y)
  expect_is(z, "icm")
  expect_true(z$control$nsims == 4)
  expect_true(dim(z$epi$i.num)[2] == 4)

  param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  y <- icm(param, init, control)
  z <- merge(x, y)
  expect_is(z, "icm")
})

test_that("merge 1 sim each", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  y <- icm(param, init, control)
  z <- merge(x, y)
  expect_is(z, "icm")
})

test_that("merge errors", {
  param <- param.icm(inf.prob = 0.4, act.rate = 0.8)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 3, verbose = FALSE)
  x <- icm(param, init, control)
  param <- param.icm(inf.prob = 0.2, act.rate = 0.8)
  control <- control.icm(type = "SI", nsteps = 10,
                         nsims = 1, verbose = FALSE)
  y <- icm(param, init, control)
  expect_error(merge(x, y), "x and y have different parameters")
})

# merge.netsim ------------------------------------------------------------

test_that("merge for netsim", {
  nw <- network_initialize(n = 100)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 1)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsteps = 20, nsims = 2,
                         save.nwstats = TRUE,
                         nwstats.formula = ~edges + degree(0),
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)
  z <- merge(x, y)
  expect_is(z, "netsim")
  expect_true(z$control$nsims == 4)
  expect_true(dim(z$epi$i.num)[2] == 4)
})

test_that("merge for netsim", {
  nw <- network_initialize(n = 100)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 1)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsteps = 20, nsims = 2,
                         save.nwstats = TRUE,
                         nwstats.formula = ~edges + degree(0),
                         verbose = FALSE, save.other = "attr")
  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)
  z <- merge(x, y, keep.other = TRUE)
  expect_is(z, "netsim")
  expect_true(length(z$attr) == 4)
  expect_true(length(z$attr[[1]]) == 6)
  z <- merge(x, y, keep.other = FALSE)
  expect_true(any(names(z) == "attr") == FALSE)
})

test_that("merge works for open sims saving nw stats", {
  nw <- network_initialize(n = 100)
  est <- netest(nw, formation = ~edges, target.stats = 20,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0.01),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.9, a.rate = 0.01, ds.rate = 0.01,
                     di.rate = 0.01)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsteps = 5, save.nwstats = TRUE,
                         nwstats.formula =
                           ~edges + meandeg + degree(0) + concurrent,
                         resimulate.network = TRUE, verbose = FALSE)

  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)
  z <- merge(x, y)

  nws <- get_nwstats(z)
  expect_true(nrow(nws) == 10)
  expect_true(length(unique(nws$sim)) == 2)

})

test_that("merge.netsim works as expected for transmat", {
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # Epidemic model
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 5, nsims = 2, verbose = FALSE)
  mod <- netsim(est, param, init, control)

  expect_equal(length(mod$stats$transmat), 2)

  mod2 <- merge(mod, mod)
  expect_equal(length(mod2$stats$transmat), 4)
  
  mod3 <- merge(mod, mod, keep.transmat = FALSE)
  expect_true(is.null(mod3$stats$transmat))
  
  mod4 <- merge(mod2, mod3)
  expect_true(is.null(mod4$stats$transmat))
})
