context("Get functions for network models")


# Model simulation --------------------------------------------------------

nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
dissolution <- ~offset(edges)
duration <- 20
coef.diss <- dissolution_coefs(dissolution, duration)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
init <- init.net(i.num = 10, i.num.g2 = 10)
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
  expect_error(get_network(mod, 4), "Specify a single sim between 1 and 3")
  expect_error(get_network(mod, 1, collapse = TRUE), "Specify collapse time")
  expect_error(get_network(mod, 1, collapse = TRUE, at = 200),
               "Specify collapse time")
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
  expect_error(get_nwstats(mod, sim = 1, network = 2),
               "Specify network between 1")
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
  expect_equal(length(get_sims(mod, sims = 1:3)$diss.stats), 3)
  expect_equal(names(get_sims(mod, sims = 1:3)$diss.stats), c("sim1", "sim2", "sim3"))
  expect_equal(length(get_sims(mod, sims = c(1, 3))$diss.stats), 2)
  expect_equal(names(get_sims(mod, sims = c(1, 3))$diss.stats), c("sim1", "sim2"))
  expect_equal(length(get_sims(mod, sims = 2)$diss.stats), 1)
  expect_equal(names(get_sims(mod, sims = 2)$diss.stats), c("sim1"))
  expect_equal(length(merge(get_sims(mod, sims = 2:3),
                            get_sims(mod, sims = c(1,3)))$diss.stats), 4)
  plot(get_sims(mod, sims = 2), type = "duration")
  plot(get_sims(mod, sims = c(1, 3)), type = "dissolution")
})


test_that("get_sims error flags", {
  expect_error(get_sims(list(a = 1)), "x must be of class netsim")
  expect_error(get_sims(mod), "Specify sims as a vector")

# get parameter set ------------------------------------------------------------
  nw <- network_initialize(n = 50)

  est <- netest(
    nw, formation = ~edges,
    target.stats = c(25),
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )

  init <- init.net(i.num = 10)

  my.randoms <- list(
    act.rate = param_random(c(0.25, 0.5, 0.75)),
    dummy.param = function() rbeta(1, 1, 2),
    dummy.strat.param = function() c(
      rnorm(1, 0, 10),
      rnorm(1, 10, 1)
    )
  )

  param <- param.net(
    inf.prob = 0.3,
    dummy = c(0, 1, 2),
    random.params = my.randoms
  )

  control <- control.net(type = "SI", nsims = 3, nsteps = 5, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  d.set <- get_param_set(mod)

  set.colnames <- c(
    "sim",
    "inf.prob",
    "dummy_1",
    "dummy_2",
    "dummy_3",
    "act.rate",
    "vital",
    "dummy.param",
    "dummy.strat.param_1",
    "dummy.strat.param_2",
    "groups"
  )

  expect_is(d.set, "data.frame")
  expect_true(setequal(names(d.set), set.colnames))
  expect_error(get_param_set(control), "`sims` must be of class netsim")
  expect_equal(dim(get_param_set(mod)), c(3, length(set.colnames)))
})

dxs <- netdx(est, dynamic = FALSE, nsims = 5,
             nwstats.formula = ~edges + nodemix("group", levels2 = TRUE), verbose = FALSE)
dxd <- netdx(est, dynamic = TRUE, nsims = 5, nsteps = 3, verbose = FALSE)

control1 <- control.net(type = "SI", nsteps = 2, nsims = 3,
                        verbose = FALSE, tergmLite = FALSE,
                        resimulate.network = TRUE,
                        nwstats.formula = ~edges + triangle)
mod1 <- netsim(est, param, init, control1)

control2 <- control.net(type = "SI", nsteps = 3, nsims = 4,
                        verbose = FALSE, tergmLite = FALSE,
                        resimulate.network = FALSE)
mod2 <- netsim(est, param, init, control2)

control3 <- control.net(type = "SI", nsteps = 4, nsims = 2,
                        verbose = FALSE, tergmLite = TRUE,
                        resimulate.network = TRUE,
                        nwstats.formula = ~edges +
                                           nodematch("group", diff = TRUE))
mod3 <- netsim(est, param, init, control3)

test_that("get_nwstats with mode = list behaves as expected", {
  expect_equal(unique(lapply(get_nwstats(dxs, mode = "list"), class)),
               list(c("matrix", "array")))
  expect_equal(unique(lapply(get_nwstats(dxd, mode = "list"), class)),
               list(c("matrix", "array")))
  expect_equal(unique(lapply(get_nwstats(mod1, mode = "list"), class)),
               list(c("matrix", "array")))
  expect_equal(unique(lapply(get_nwstats(mod2, mode = "list"), class)),
               list(c("matrix", "array")))
  expect_equal(unique(lapply(get_nwstats(mod3, mode = "list"), class)),
               list(c("matrix", "array")))

  expect_equal(unique(lapply(get_nwstats(dxs, mode = "list"), dim)),
               list(c(5, 4)))
  expect_equal(unique(lapply(get_nwstats(dxd, mode = "list"), dim)),
               list(c(3, 1)))
  expect_equal(unique(lapply(get_nwstats(mod1, mode = "list"), dim)),
               list(c(2, 2)))
  expect_equal(unique(lapply(get_nwstats(mod2, mode = "list"), dim)),
               list(c(3, 1)))
  expect_equal(unique(lapply(get_nwstats(mod3, mode = "list"), dim)),
               list(c(4, 3)))

  expect_equal(length(get_nwstats(dxs, mode = "list")), 1)
  expect_equal(length(get_nwstats(dxd, mode = "list")), 5)
  expect_equal(length(get_nwstats(mod1, mode = "list")), 3)
  expect_equal(length(get_nwstats(mod2, mode = "list")), 4)
  expect_equal(length(get_nwstats(mod3, mode = "list")), 2)

  expect_equal(length(get_nwstats(dxs, sim = c(1), mode = "list")), 1)
  expect_equal(length(get_nwstats(dxd, sim = c(5,3,1), mode = "list")), 3)
  expect_equal(length(get_nwstats(mod1, sim = c(2,3), mode = "list")), 2)
  expect_equal(length(get_nwstats(mod2, sim = c(1,3,2), mode = "list")), 3)
  expect_equal(length(get_nwstats(mod3, sim = c(1,2), mode = "list")), 2)
})

test_that("get_network_attributes functions as intended", {
  nw <- network.initialize(10, directed = FALSE, bipartite = 4)
  nw %n% "newattr" <- "string"
  expect_equal(get_network_attributes(nw), list(bipartite = 4,
                                                directed = FALSE,
                                                hyper = FALSE,
                                                loops = FALSE,
                                                multiple = FALSE,
                                                n = 10,
                                                newattr = "string"))
})
