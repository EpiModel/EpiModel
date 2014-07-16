context("merge")

test_that("merge: works for open sims saving nw stats", {
  nw <- network.initialize(n = 50, directed = FALSE)
  est <- netest(nw,
                formation = ~ edges,
                dissolution = ~offset(edges),
                target.stats = 10,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  param <- param.net(trans.rate = 0.9, b.rate = 0.01,
                     ds.rate = 0.01,
                     di.rate = 0.01)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI",
                         nsteps = 5,
                         save.stats = TRUE,
                         stats.formula = ~edges + meandeg + degree(0) + concurrent,
                         verbose = FALSE)

  x <- netsim(est, param, init, control)
  y <- netsim(est, param, init, control)
  z <- merge(x, y)

  expect_equal(length(z$stats), 2)
  expect_true(all(sapply(z$stats$nwstats, dim)[1,] == 5) &
              all(sapply(z$stats$nwstats, dim)[2,] == 1))
})

