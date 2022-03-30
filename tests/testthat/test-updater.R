context("Network model with param updater")

test_that("netsim with param updater", {
  # Create the list.param.updaters
  list.param.updaters <- list(
    # this is one updater
    list(
      at = 10,
      verbose = TRUE,
      param = list(
        inf.prob = 0.3,
        act.rate = 0.3
      )
    ),
    # this is another updater
    list(
      at = 20,
      verbose = TRUE,
      param = list(
        # inf.prob = function(x) plogis(qlogis(x) - log(10)),
        # act.rate = function(x) plogis(qlogis(x) - log(10))
        inf.prob = 0.01
      )
    )
  )

  # Create the list.control.updaters
  list.control.updaters <- list(
    # this is one updater
    list(
      at = 15,
      verbose = TRUE,
      control = list(
        verbose = TRUE
      )
    ),
    # this is another updater
    list(
      at = 25,
      verbose = TRUE,
      control = list(
        verbose = FALSE
      )
    ),
    list(
      at = 30,
      verbose = TRUE,
      control = list(
        resimulate.network = FALSE
      )
    )
  )

  # Do not forget to add it to `param`
  param <- param.net(
    inf.prob = 0.1,
    act.rate = 0.1,
    .param.updater.list = list.param.updaters
  )

  # Enable the module in `control`
  control <- control.net(
    type = NULL, # must be NULL as we use a custom module
    nsims = 1,
    nsteps = 50,
    verbose = FALSE,
    infection.FUN = infection.net,
    .control.updater.list = list.control.updaters,
    resimulate.network = TRUE
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

  # `resimulate.network` is turned of at step 30. We check that the number of
  # observations in the "networkDynamic" object is < than 30 and not 50 (the
  # number of timestep in the simulation)
  n_obs <- length(
    get.network.attribute(mod$network[[1]][[1]], 'net.obs.period')$observations
  )
  expect_lt(n_obs, 30)
})
