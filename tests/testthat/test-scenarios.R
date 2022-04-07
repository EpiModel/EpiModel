context("Network model with scenarios")

test_that("SIS with scenarios", {
  set.seed(10)

  nw <- network_initialize(n = 200)
  est <- netest(nw,
    formation = ~edges, target.stats = 60,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )

  param <- param.net(inf.prob = 0.9, rec.rate = 0.01, act.rate = 2)
  control <- control.net(type = "SIS", nsims = 1, nsteps = 50, verbose = FALSE)
  init <- init.net(i.num = 10)

  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf.prob, ~rec.rate,
    "base", 0, 0.9, 0.01,
    "multiple_changes", 0, 0.1, 0.04,
    "multiple_changes", 20, 0.9, 0.01,
    "multiple_changes", 40, 0.1, 0.1
  )

  scenarios.list <- create_scenario_list(scenarios.df)
  expect_length(scenarios.list, 2)

  sc.param <- use_scenario(param, scenarios.list[[1]])
  expect_silent(netsim(est, sc.param, init, control))
  sc.param <- use_scenario(param, scenarios.list[[2]])
  expect_message(netsim(est, sc.param, init, control))

  # .at not a integer
  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf.prob, ~rec.rate,
    "multiple_changes", "text", 0.1, 0.1
  )
  expect_error(scenarios.list <- create_scenario_list(scenarios.df))

  # inf_prob with an underscore
  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf_prob, ~rec.rate,
    "multiple_changes", 0, 0.1, 0.1
  )
  expect_error(scenarios.list <- create_scenario_list(scenarios.df))

  # rec.rate2 not in param
  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf.prob, ~rec.rate2,
    "multiple_changes", 0, 0.1, 0.1
  )
  scenarios.list <- create_scenario_list(scenarios.df)
  expect_error(sc.param <- use_scenario(param, scenarios.list[[1]]))
})


