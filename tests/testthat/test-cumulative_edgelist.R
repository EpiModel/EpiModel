context("cumulative edgelist")

nw <- network_initialize(n = 50)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(
  nw,
  formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~ offset(edges), 10, 0),
  verbose = FALSE
)

init <- init.net(i.num = 10)

test_that("netsim, SI, Cumulative Edgelist", {
  control <- control.net(
    type = NULL, # must be NULL as we use a custom module
    nsims = 1,
    nsteps = 50,
    verbose = FALSE,
    infection.FUN = infection.net,
    cumulative_edgelist.FUN = cumulative_edgelist.net,
    truncate.el_cuml = NULL,
    raw.output = TRUE
  )

  param <- param.net(
    inf.prob = 0.3,
    act.rate = 0.1
  )

  mod <- netsim(est, param, init, control)
  d <- get_cumulative_edgelists_df(mod[[1]])

  expect_is(d, "data.frame")

  control <- control.net(
    type = NULL, # must be NULL as we use a custom module
    nsims = 1,
    nsteps = 50,
    verbose = FALSE,
    tergmLite = TRUE,
    infection.FUN = infection.net,
    cumulative_edgelist.FUN = cumulative_edgelist.net,
    truncate.el_cuml = 10,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)
  d <- get_cumulative_edgelists_df(mod[[1]])
  expect_lte(min(d$stop, na.rm = TRUE), 40)
})
