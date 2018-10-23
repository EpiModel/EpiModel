context("Network extended models")

################################################################################

test_that("edges models", {
  skip_on_cran()

  nw <- network.initialize(n = 100, directed = FALSE)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  expect_is(est, "netest")

  ## "SI, 1M, CL: 1 sim"
  param <- param.net(inf.prob = 0.5)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsims = 1, nsteps = 25, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) >= 1)
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SI, 1M, CL: 1 sim, inf.prob = 0"
  param <- param.net(inf.prob = 0)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsims = 1, nsteps = 25, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) == 1)
  expect_true(max(x$epi$si.flow) == 0)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SI, 1M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsims = 2, nsteps = 25, verbose = FALSE,
                         nwstats.formula = ~edges + meandeg + concurrent)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) >= 1)
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "formation", plots.joined = FALSE)
  plot(x, type = "formation", stats = "edges")
  plot(x, type = "formation", stats = c("edges", "meandeg"))
  plot(x, type = "formation", sim.lines = TRUE, qnts.smooth = FALSE)
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SI, 1M, CL: 2 sim, use TEAs"
  param <- param.net(inf.prob = 1)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", nsims = 2, nsteps = 25,
                         verbose = FALSE, tea.status = TRUE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) >= 1)
  expect_true(max(x$epi$i.num) <= 100)
  expect_true(x$control$tea.status, TRUE)
  expect_true(sum(get.vertex.attribute.active(x$network[[1]],
                                              prefix = "testatus", at = 1) == "i") >= 0)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SI, 1M, CL: 2 sim, set status.vector"
  param <- param.net(inf.prob = 0.5)
  init <- init.net(status.vector = c(rep("i", 10),
                                     rep("s", 90)))
  control <- control.net(type = "SI", nsims = 2, nsteps = 25,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  # this tests that starting infected is same across sims
  expect_true(all(x$epi$i.num[1,] == 10))
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  expect_error(plot(x, type = "network", col.status = TRUE))
  test_net(x)
  rm(x)

  ## "SIR, 1M, CL: 1 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.02)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsims = 1, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) >= 1)
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SIR, 1M, CL: 1 sim, inf.prob=0"
  param <- param.net(inf.prob = 0, rec.rate = 0.02)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsims = 1,
                         nsteps = 25, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) == 10)
  expect_true(max(x$epi$si.flow) == 0)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SIR, 1M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.1)
  init <- init.net(i.num = 1, r.num = 0)
  control <- control.net(type = "SIR", nsims = 2, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) >= 1)
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SIR, 1M, CL: 2 sim, set status.vector"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.1)
  init <- init.net(status.vector = rep(c("s", "i"), each = 50))
  control <- control.net(type = "SIR", nsims = 2, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  # this tests that starting infected is same across sims
  expect_true(all(x$epi$i.num[1,] == 50))
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SIS, 1M, CL: 1 sim"
  param <- param.net(inf.prob = 0.9, rec.rate = 0.01)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SIS", nsims = 1, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) >= 1)
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SIS, 1M, CL: 1 sim, inf.prob=0"
  param <- param.net(inf.prob = 0, rec.rate = 0.01)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SIS", nsims = 1, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) == 1)
  expect_true(max(x$epi$si.flow) == 0)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SIS, 1M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.01)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SIS", nsims = 2, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) >= 1)
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIS, 1M, CL: 2 sim, test TEAs"
  param <- param.net(inf.prob = 1, rec.rate = 0.01)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SIS", nsims = 2, nsteps = 25,
                         tea.status = TRUE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) <= 100)
  expect_true(x$control$tea.status, TRUE)
  expect_true(sum(get.vertex.attribute.active(x$network[[1]],
                                              prefix = "testatus", at = 1) == "i") >= 0)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SIS, 1M, CL: 2 sim, set status.vector"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.01)
  init <- init.net(status.vector = c(rep("i", 10), rep("s", 90)))
  control <- control.net(type = "SIS", nsims = 2, nsteps = 25,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  # this tests that starting infected is same across sims
  expect_true(all(x$epi$i.num[1,] == 10))
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  expect_error(plot(x, type = "network", sims = "mean", col.status = TRUE))
  test_net(x)
})


################################################################################

test_that("High departure rate models", {

  ## "netsim: 1M, ds.rate = 0.5"
  nw <- network.initialize(n = 25, directed = FALSE)
  est <- netest(nw, formation = ~edges, target.stats = 12,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0.01),
                edapprox = TRUE, verbose = FALSE)
  param <- param.net(inf.prob = 0.5, act.rate = 2,
                     b.rate = 0.01, ds.rate = 0.5,
                     di.rate = 0.25)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 25,
                         nsims = 1, tea.status = FALSE,
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_equal(unique(sapply(x$epi, nrow)), 25)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "netsim: 1M, di.rate = 0.5"
  param <- param.net(inf.prob = 0.1, act.rate = 2, b.rate = 0.01,
                     ds.rate = 0.01, di.rate = 0.5)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 25,
                         nsims = 1, tea.status = FALSE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_equal(unique(sapply(x$epi, nrow)), 25)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)

})



################################################################################

test_that("edges bipartite models", {
  skip_on_cran()

  nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
  est5 <- netest(nw, formation = ~edges, target.stats = 25,
                 coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                 edapprox = TRUE, verbose = FALSE)
  expect_is(est5, "netest")

  ## "SI, 2M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.1)
  init <- init.net(i.num = 10, i.num.m2 = 0)
  control <- control.net(type = "SI", nsims = 2, nsteps = 25,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est5, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$modes, 2)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", shp.bip = "triangle")
  expect_error(plot(x, type = "network", shp.bip = TRUE))
  test_net(x)
  rm(x)

  ## "SIR, 2M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.1,
                     rec.rate = 0.1, rec.rate.m2 = 0.1)
  init <- init.net(i.num = 10, i.num.m2 = 10,
                   r.num = 0, r.num.m2 = 0)
  control <- control.net(type = "SIR", nsims = 2, nsteps = 25,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est5, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$modes, 2)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIR, 2M, CL: rec.rate = 0"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.25,
                     rec.rate = 0, rec.rate.m2 = 0, act.rate = 2)
  init <- init.net(i.num = 10, i.num.m2 = 10,
                   r.num = 0, r.num.m2 = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 2,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est5, param, init, control)
  expect_equal(max(x$epi$ir.flow), 0)
  expect_equal(max(x$epi$ir.flow.m2), 0)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIS, 2M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.25,
                     rec.rate = 0.01, rec.rate.m2 = 0.01)
  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SIS", nsims = 2, nsteps = 25,
                         tea.status = FALSE, verbose = FALSE)
  x <- netsim(est5, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$modes, 2)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIS, 2M, CL: rec.rate = 0"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.25,
                     rec.rate = 0, rec.rate.m2 = 0, act.rate = 2)
  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 2,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est5, param, init, control)
  expect_equal(max(x$epi$is.flow), 0)
  expect_equal(max(x$epi$is.flow.m2), 0)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

})

################################################################################

test_that("Open population 1 mode models", {
  skip_on_cran()

  nw <- network.initialize(n = 100, directed = FALSE)
  est.vit <- netest(nw, formation = ~edges, target.stats = 25,
                    coef.diss = dissolution_coefs(~offset(edges), 10, 0.02),
                    verbose = FALSE)

  ## "SI, 1M, OP: 1 sim"
  param <- param.net(inf.prob = 0.5, act.rate = 2, b.rate = 0.02,
                     ds.rate = 0.02, di.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SI, 1M, OP, deterministic arrivals and departures: 1 sim"
  param <- param.net(inf.prob = 0.5, act.rate = 2, a.rate = 0.02,
                     ds.rate = 0.02, di.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 25, nsims = 1,
                         verbose = FALSE, tea.status = FALSE,
                         b.rand = FALSE, d.rand = FALSE)
  x <- netsim(est.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  test_net(x)
  rm(x)

  ## "SI, 1M, OP: 2 sim"
  param <- param.net(inf.prob = 0.5, act.rate = 2, b.rate = 0.02,
                     ds.rate = 0.02, di.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 10,
                         nsims = 2,
                         verbose = FALSE,
                         tea.status = FALSE)
  x <- netsim(est.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIR, 1M OP: 1 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.1, act.rate = 2,
                     b.rate = 0.02, ds.rate = 0.02, di.rate = 0.02,
                     dr.rate = 0.02)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIR, 1M OP, deterministic arrivals, recoveries and departures: 1 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.1, act.rate = 2,
                     a.rate = 0.02, ds.rate = 0.02, di.rate = 0.02,
                     dr.rate = 0.02)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsteps = 25, nsims = 1,
                         verbose = FALSE, tea.status = FALSE,
                         b.rand = FALSE, d.rand = FALSE, rec.rand = FALSE)
  x <- netsim(est.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  test_net(x)
  rm(x)

  ## "SIR, 1M, OP: 2 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.1, act.rate = 2,
                     a.rate = 0.02, ds.rate = 0.02, di.rate = 0.02,
                     dr.rate = 0.02)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 2,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIS, 1M, OP: 1 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.01, act.rate = 2,
                     a.rate = 0.02, ds.rate = 0.02, di.rate = 0.02)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIS, 1M, OP: 2 sim"
  control <- control.net(type = "SIS", nsteps = 10, nsims = 2,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

})

################################################################################

test_that("Open-population bipartite models", {
  skip_on_cran()

  nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
  est5.vit <- netest(nw, formation = ~edges, target.stats = 25,
                     coef.diss = dissolution_coefs(~offset(edges), 10, 0.02),
                     edapprox = TRUE, verbose = FALSE)

  ## "SI, 2M, OP: 1 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.1, act.rate = 2,
                     a.rate = 0.02, ds.rate = 0.02, di.rate = 0.02,
                     a.rate.m2 = 0.02, ds.rate.m2 = 0.02, di.rate.m2 = 0.02)
  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est5.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SI, 2M, OP, deterministic arrivals and departures: 1 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.1, act.rate = 2,
                     a.rate = 0.02, ds.rate = 0.02, di.rate = 0.02,
                     a.rate.m2 = 0.02, ds.rate.m2 = 0.02, di.rate.m2 = 0.02)
  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SI", nsteps = 25, nsims = 1,
                         verbose = FALSE, tea.status = FALSE,
                         b.rand = FALSE, d.rand = FALSE)
  x <- netsim(est5.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  test_net(x)
  rm(x)

  ## "SI, 2M, OP: 2 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.1, act.rate = 2,
                     a.rate = 0.02, ds.rate = 0.02, di.rate = 0.02,
                     a.rate.m2 = 0.02, ds.rate.m2 = 0.02, di.rate.m2 = 0.02)
  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est5.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIR, 2M, OP: 1 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.1, rec.rate = 0.1,
                     rec.rate.m2 = 0.1, act.rate = 2, a.rate = 0.02,
                     a.rate.m2 = NA, ds.rate = 0.02, ds.rate.m2 = 0.02,
                     di.rate = 0.02, di.rate.m2 = 0.02, dr.rate = 0.02,
                     dr.rate.m2 = 0.02)
  init <- init.net(i.num = 10, i.num.m2 = 0,
                   r.num = 0, r.num.m2 = 10)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est5.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)


  ## "SIR, 2M, OP, deterministic arrivals, departures, and recoveries: 1 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.m2 = 0.1, rec.rate = 0.1,
                     rec.rate.m2 = 0.1, act.rate = 2, a.rate = 0.02,
                     a.rate.m2 = NA, ds.rate = 0.02, ds.rate.m2 = 0.02,
                     di.rate = 0.02, di.rate.m2 = 0.02, dr.rate = 0.02,
                     dr.rate.m2 = 0.02)
  init <- init.net(i.num = 10, i.num.m2 = 0,
                   r.num = 0, r.num.m2 = 10)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1,
                         verbose = FALSE, tea.status = FALSE,
                         b.rand = FALSE, d.rand = FALSE, rec.rand = FALSE)
  x <- netsim(est5.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  test_net(x)
  rm(x)

  ## "SIR, 2M, OP: 3 sim"
  control <- control.net(type = "SIR", nsteps = 10, nsims = 2,
                         verbose = FALSE, tea.status = FALSE)
  x <- netsim(est5.vit, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$vital, TRUE)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

})

################################################################################

test_that("Extinction open-population models", {
  skip_on_cran()

  nw <- network.initialize(n = 25, bipartite = 10, directed = FALSE)
  est <- netest(nw, formation = ~edges, target.stats = 15,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0.02),
                edapprox = TRUE, verbose = FALSE)

  ## "netsim: 2M, ds.rate = 0.5"
  param <- param.net(inf.prob = 0.1, inf.prob.m2 = 0.1, act.rate = 2,
                     a.rate = 0.02, ds.rate = 0.5, di.rate = 0.5,
                     a.rate.m2 = 0.02, ds.rate.m2 = 0.02, di.rate.m2 = 0.02)
  init <- init.net(i.num = 1, i.num.m2 = 0)
  control <- control.net(type = "SI", nsteps = 30, nsims = 1,
                         tea.status = FALSE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "netsim: 2M, ds.rate.m2 = 0.5"
  param <- param.net(inf.prob = 0.1, inf.prob.m2 = 0.1, act.rate = 2,
                     a.rate = 0.02, ds.rate = 0.02, di.rate = 0.02,
                     a.rate.m2 = 0.02, ds.rate.m2 = 0.5, di.rate.m2 = 0.5)
  init <- init.net(i.num = 1, i.num.m2 = 0)
  control <- control.net(type = "SI", nsteps = 30, nsims = 1,
                         tea.status = FALSE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

})


################################################################################

test_that("Extended post-simulation diagnosntic tests", {
  skip_on_cran()

  nw <- network.initialize(100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "risk", rep(1:5, each = 20))

  est <- netest(nw,
                formation = ~edges + nodefactor("risk"),
                target.stats = c(50, 20, 20, 20, 20),
                coef.diss = dissolution_coefs(~offset(edges), 25))

  dx <- netdx(est, nsims = 2, nsteps = 10)
  plot(dx)

  param <- param.net(inf.prob = 0.5)
  init <- init.net(i.num = 10)

  control <- control.net(type = "SI", nsteps = 10, nsims = 2)
  sim <- netsim(est, param, init, control)
  plot(sim, type = "formation")

  control <- control.net(type = "SI", nsteps = 10, nsims = 1)
  sim <- netsim(est, param, init, control)
  plot(sim, type = "formation")

  est <- netest(nw,
                formation = ~edges + concurrent,
                target.stats = c(50, 20),
                coef.diss = dissolution_coefs(~offset(edges), 25))

  dx <- netdx(est, nsims = 5, nsteps = 100)
  plot(dx)

  param <- param.net(inf.prob = 0.5)
  init <- init.net(i.num = 10)

  control <- control.net(type = "SI", nsteps = 10, nsims = 2)
  sim <- netsim(est, param, init, control)
  plot(sim, type = "formation")

  control <- control.net(type = "SI", nsteps = 10, nsims = 1)
  sim <- netsim(est, param, init, control)
  plot(sim, type = "formation")

})


################################################################################

test_that("status.vector and infTime.vector", {

  n <- 100
  nw <- network.initialize(n = n, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  # Epidemic model
  param <- param.net(inf.prob = 0.3, rec.rate = 0.01)

  status <- sample(c("s", "i"), size = n, replace = TRUE, prob = c(0.8, 0.2))
  infTime <- rep(NA, n)
  infTime[which(status == "i")] <- -rgeom(sum(status == "i"), prob = 0.01) + 2
  init <- init.net(status.vector = status, infTime.vector = infTime)

  control <- control.net(type = "SIS", nsteps = 100, nsims = 5, verbose.int = 0)
  mod1 <- netsim(est1, param, init, control)
  expect_is(mod1, "netsim")

  control <- control.net(type = "SIR", nsteps = 100, nsims = 5, verbose.int = 0)
  mod2 <- netsim(est1, param, init, control)
  expect_is(mod2, "netsim")

})


