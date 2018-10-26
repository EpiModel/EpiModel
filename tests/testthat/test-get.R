context("Get functions for network models")


# Model simulation --------------------------------------------------------

nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
formation <- ~edges
target.stats <- 50
dissolution <- ~offset(edges)
duration <- 20
coef.diss <- dissolution_coefs(dissolution, duration)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
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

test_that("get_transmat error flags", {
  expect_error(get_transmat(mod, sim = 5), "Specify sim between 1 and 3")
  expect_error(get_transmat(list(a = 1)), "x must be of class netsim")
  mod$stats$transmat <- NULL
  expect_error(get_transmat(mod, 1), "transmat not saved")
})


# get nwstats -------------------------------------------------------------

test_that("get_nwstats extracts data frame", {
  a <- get_nwstats(mod, sim = 1)
  expect_is(a, "data.frame")
  expect_equal(get_nwstats(mod, sim = 1:3), get_nwstats(mod))
})

test_that("get_nwstats error flags", {
  expect_error(get_nwstats(list(a = 1)), "x must be of class netsim")
  expect_error(get_nwstats(mod, sim = 5))
  expect_error(get_nwstats(mod, sim = 1, network = 2), "Specify network between 1")
  mod$stats$nwstats <- NULL
  expect_error(get_nwstats(mod), "Network statistics not saved")
})


# get sims ----------------------------------------------------------------

test_that("get_sims extracts simulations", {
  expect_is(get_sims(mod, sims = 1), "netsim")
  expect_is(get_sims(mod, sims = 2:3), "netsim")
  expect_is(get_sims(mod, sims = "mean", var = "i.num"), "netsim")
  expect_is(get_sims(mod, sims = 1:3), "netsim")
  expect_is(get_sims(mod, sims = 1, var = c("i.num", "s.num")), "netsim")
  expect_is(get_sims(mod, sims = c(1, 3), var = c("i.num", "s.num")), "netsim")
})


test_that("get_sims error flags", {
  expect_error(get_sims(list(a = 1)), "x must be of class netsim")
  expect_error(get_sims(mod), "Specify sims as a vector")
})

