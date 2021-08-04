context("Records: attr_history and raw objects")

test_that("Time varying elements", {
  test_logger <- function(dat, at) {
    nodes <- get_posit_ids(dat)

    some_nodes <- nodes[runif(length(nodes)) < 0.1]
    dat <- record_attr_history(
      dat, at,
      "attr_norm",
      some_nodes,
      rnorm(length(some_nodes))
    )

    some_nodes <- nodes[runif(length(nodes)) < 0.1]
    dat <- record_attr_history(
      dat, at,
      "attr_unif",
      some_nodes,
      runif(length(some_nodes))
    )

    some_nodes <- nodes[runif(length(nodes)) < 0.1]
    dat <- record_attr_history(
      dat, at,
      "attr_fix",
      some_nodes,
      at
    )

    return(dat)
  }

  param <- param.net(
    inf.prob = 0.1,
    act.rate = 0.1
  )

  # Enable the module in `control`
  control <- control.net(
    type = NULL, # must be NULL as we use a custom module
    nsims = 1,
    nsteps = 20,
    verbose = FALSE,
    infection.FUN = infection.net,
    logger.FUN = test_logger
  )

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

  mod <- netsim(est, param, init, control)
  attr_history <- get_attr_history(mod)
  expect_is(attr_history, "list")
  expect_is(attr_history[[1]], "data.frame")
  expect_equal(names(attr_history), c("attr_norm", "attr_unif", "attr_fix"))
})
