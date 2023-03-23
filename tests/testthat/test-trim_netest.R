context("trim_netest functionality")

test_that("trim_netest keep argument behaves as expected", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "sex", rbinom(50, 1, 0.5))

  attrname <- "sex"

  est <- netest(nw, formation = ~edges + nodematch(attrname),
                target.stats = c(25, 10),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)

  dxs <- netdx(est, nsims = 10, dynamic = FALSE)
  dxd <- netdx(est, nsims = 2, nsteps = 10, dynamic = TRUE)

  param <- param.net(inf.prob = 0.3, act.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
  sim <- netsim(est, param, init, control)

  expect_is(dxs, "netdx")
  expect_is(dxd, "netdx")
  expect_is(sim, "netsim")

  trim_est <- trim_netest(est)

  expect_error(netdx(trim_est, nsims = 10, dynamic = FALSE),
               "object 'attrname' not found")
  expect_error(netdx(trim_est, nsims = 2, nsteps = 10, dynamic = TRUE),
               "object 'attrname' not found")
  expect_error(netsim(trim_est, param, init, control),
               "object 'attrname' not found")

  trim_est_keep <- trim_netest(est, keep = "attrname")

  trim_dxs_keep <- netdx(trim_est_keep, nsims = 10, dynamic = FALSE)
  trim_dxd_keep <- netdx(trim_est_keep, nsims = 2, nsteps = 10, dynamic = TRUE)
  trim_sim_keep <- netsim(trim_est_keep, param, init, control)

  expect_is(trim_dxs_keep, "netdx")
  expect_is(trim_dxd_keep, "netdx")
  expect_is(trim_sim_keep, "netsim")
})
