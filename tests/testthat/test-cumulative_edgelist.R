context("cumulative edgelist")

nw <- network_initialize(n = 100)
est <- netest(
  nw,
  formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~ offset(edges), 10, 0),
  verbose = FALSE
)

param <- param.net(
  inf.prob = 0.3,
  act.rate = 0.1
)

init <- init.net(i.num = 10)

test_that("netsim, SI, Cumulative Edgelist", {
  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 100,
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    truncate.el.cuml = NULL,
    cumulative.edgelist = TRUE,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)
  d <- get_cumulative_edgelists_df(mod[[1]])

  expect_is(d, "data.frame")
  expect_equal(colnames(d), c("head", "tail", "start", "stop", "network"))

  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 100,
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    truncate.el.cuml = 40,
    cumulative.edgelist = TRUE,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)
  d <- get_cumulative_edgelists_df(mod[[1]])
  expect_gte(min(d$stop, na.rm = TRUE), 40)

  d <- get_partners(mod[[1]], 100, 1:4)
  expect_equal(colnames(d), c("index", "partner", "start", "stop", "network"))
})
