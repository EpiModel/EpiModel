context("Get functions for network models")


# Model simulation --------------------------------------------------------

nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
formation <- ~ edges
target.stats <- 50
dissolution <- ~ offset(edges)
duration <- 20
coef.diss <- dissolution_coefs(dissolution, duration)
est <- netest(nw,
               formation,
               dissolution,
               target.stats,
               coef.diss,
               verbose = FALSE)
param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
init <- init.net(i.num = 10, i.num.m2 = 10)
control <- control.net(type = "SI", nsteps = 10, nsims = 3,
                       verbose = FALSE)
mod <- netsim(est, param, init, control)


# get network -------------------------------------------------------------

test_that("get_network extracts to nD object", {
  a <- get_network(mod, sim = 2)
  expect_is(a, "networkDynamic")
})

test_that("get_network extracts and collapses to network object", {
  a <- get_network(mod, collapse = TRUE, at = 3)
  expect_is(a, "network")
})

test_that("get_network yields warning when missing network", {
  mod$network <- NULL
  expect_error(get_network(mod, sim = 1))
})

test_that("get_network yields warning for incorrect sim", {
  expect_error(get_network(mod, sim = 5))
})

test_that("get_network error flags", {
  expect_error(get_network(list(a = 1), 1), "must be of class netsim")
  expect_error(get_network(mod, 4), "Specify sim between 1 and 3")
  expect_error(get_network(mod, 1, collapse = TRUE), "Specify collapse time")
  expect_error(get_network(mod, 1, collapse = TRUE, at = 200), "Specify collapse time")
  expect_error(get_network(mod, 1, 2), "Specify network")
})


# get transmat ------------------------------------------------------------

test_that("get_transmat extracts data frame", {
  a <- get_transmat(mod, sim = 1)
  expect_is(a, "data.frame")
})

test_that("get_transmat yields warning for incorrect sim", {

})

test_that("get_transmat error flags", {
  expect_error(get_transmat(mod, sim = 5))
  expect_error(get_transmat(list(a = 1)))
  mod$stats$transmat <- NULL
  expect_error(get_transmat(mod, 1))
})


# get nwstats -------------------------------------------------------------

test_that("get_nwstats extracts data frame", {
  a <- get_nwstats(mod, sim = 1)
  b <- a[[1]]
  expect_is(a, "list")
  expect_is(b, "data.frame")
})

test_that("get_nwstats yields warning for incorrect sim", {
  expect_error(get_nwstats(mod, sim = 5))
})
