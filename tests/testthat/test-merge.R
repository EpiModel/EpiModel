context("Stochastic model merging")

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
})

test_that("merge works for open sims saving nw stats", {
  nw <- network.initialize(n = 50, directed = FALSE)
  est <- netest(nw, formation = ~edges, target.stats = 10,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  param <- param.net(inf.prob = 0.9, b.rate = 0.01,
                     ds.rate = 0.01,
                     di.rate = 0.01)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsteps = 5, save.stats = TRUE,
                         nwstats.formula = ~edges + meandeg + degree(0) + concurrent,
                         verbose = FALSE)

  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)
  z <- merge(x, y)

  expect_equal(length(z$stats), 2)
  expect_true(all(sapply(z$stats$nwstats, dim)[1,] == 5) &
              all(sapply(z$stats$nwstats, dim)[2,] == 4))
})

