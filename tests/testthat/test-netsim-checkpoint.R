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
  control <- control.net(
    type = "SI", nsims = 3, nsteps = 10, ncores = 2, verbose = FALSE,
    .checkpoint.dir = "cp", .checkpoint.steps = 3, .checkpoint.keep = TRUE
  )

  mod <- netsim(est, param, init, control)

  # check that the "checkpoint dir" contains 3 files (.checkpoint.keep = TRUE)
  expect_length(list.files("cp"), 3)

  control <- control.net(
    type = "SI", nsims = 3, nsteps = 10, ncores = 2, verbose = FALSE,
    .checkpoint.dir = "cp", .checkpoint.steps = 3, .checkpoint.keep = FALSE
  )

  mod <- netsim(est, param, init, control)
  # check that the "checkpoint dir" has been removed (.checkpoint.keep = FALSE)
  expect_false(dir.exists("cp"))

  control <- control.net(
    type = "SI", nsims = 1, nsteps = 5, ncores = 1, verbose = FALSE,
    .checkpoint.dir = "cp", .checkpoint.steps = -3
  )
  mod <- netsim(est, param, init, control)

  # .checkpoint.steps cannot be negative, thus the ".checkpointed" flag is
  # set to FALSE
  expect_false(mod$control[[".checkpointed"]])
})
