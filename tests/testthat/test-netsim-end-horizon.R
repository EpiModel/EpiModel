test_that("netsim with checkpoint", {
  nw <- network_initialize(n = 50)
  est <- netest(
    nw,
    formation = ~edges,
    target.stats = 24,
    coef.diss = dissolution_coefs(~ offset(edges), 10, 0),
    verbose = FALSE
  )

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)

  end_horizon <- list(at = 5, modules = c("resim_nets.FUN", "recovery.FUN"))

  control <- control.net(
    type = "SI", nsims = 1, nsteps = 10, ncores = 1, verbose = FALSE,
    end.horizon = end_horizon
  )

  netsim(est, param, init, control)

  end_horizon_err_at <- list(at = 0, modules = "resim_nets.FUN")
  control$end.horizon <- end_horizon_err_at
  expect_error(netsim(est, param, init, control))

  end_horizon_err_at <- list(at = "a", modules = "resim_nets.FUN")
  control$end.horizon <- end_horizon_err_at
  expect_error(netsim(est, param, init, control))

  end_horizon_err_mod <- list(at = "a", modules = "not_there.FUN")
  control$end.horizon <- end_horizon_err_mod
  expect_error(netsim(est, param, init, control))

  # test for actual removal of 2 modules at once
  #   - resimulation module (nwstats will stop being produced)
  #   - infection module (prevalence will stay constant; see notes)
  end_horizon <- list(at = 5, modules = c("resim_nets.FUN", "infection.FUN"))
  control <- control.net(
    type = "SI", nsims = 1, nsteps = 20, ncores = 1, resimulate.network = TRUE,
    verbose = FALSE, end.horizon = end_horizon
  )
  sim <- netsim(est, param, init, control)
  expect_true(nrow(get_nwstats(sim)) == 4)

  # SI module without arrival or departure
  # prevalence stays constant if the infection module is disabled
  expect_length(unique(sim$epi$i.num[5:20, 1]), 1)
})
