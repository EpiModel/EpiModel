context("Network model with param updater")

test_that("netsim with param updater", {
  # Create the parame.updater.list
  param.updater.list <- list(
    # this is one updater
    list(
      at = 100,
      verbose = TRUE,
      param = list(
        inf.prob = 0.3,
        act.rate = 0.3
      )
    ),
    # this is another updater
    list(
      at = 125,
      verbose = TRUE,
      param = list(
        # inf.prob = function(x) plogis(qlogis(x) - log(10)),
        # act.rate = function(x) plogis(qlogis(x) - log(10))
        inf.prob = 0.01
      )
    )
  )

  # Do not forget to add it to `param`
  param <- param.net(
    inf.prob = 0.1,
    act.rate = 0.1,
    param.updater.list = param.updater.list
  )

  # Enable the module in `control`
  control <- control.net(
    type = NULL, # must be NULL as we use a custom module
    nsims = 1,
    nsteps = 200,
    verbose = FALSE,
    updater.FUN = updater.net,
    infection.FUN = infection.net
  )

  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(
    nw,
    formation = ~edges,
    target.stats = 25,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )

  init <- init.net(i.num = 10)

  expect_message(mod <- netsim(est, param, init, control))
})
