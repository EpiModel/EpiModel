context("Network extended models")

test_that("edges models", {
  skip_on_cran()

  nw <- network_initialize(n = 100)
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
  expect_true(max(x$epi$si.flow, na.rm = TRUE) == 0)
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
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) >= 1)
  expect_true(max(x$epi$i.num) <= 100)
  expect_true(sum(get.vertex.attribute.active(x$network[[1]][[1]],
                                              prefix = "testatus",
                                              at = 1) == "i") >= 0)
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
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  # this tests that starting infected is same across sims
  expect_true(all(x$epi$i.num[1, ] == 10))
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network", col.status = TRUE)
  test_net(x)
  rm(x)

  ## "SIR, 1M, CL: 1 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.01)
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
  param <- param.net(inf.prob = 0, rec.rate = 0.01)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsims = 1,
                         nsteps = 25, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) == 10)
  expect_true(max(x$epi$si.flow, na.rm = TRUE) == 0)
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
  expect_true(all(x$epi$i.num[1, ] == 50))
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
  expect_true(max(x$epi$si.flow, na.rm = TRUE) == 0)
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
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_true(max(x$epi$i.num) <= 100)
  expect_true(sum(get.vertex.attribute.active(x$network[[1]][[1]],
                                              prefix = "testatus",
                                              at = 1) == "i") >= 0)
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
                         verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  # this tests that starting infected is same across sims
  expect_true(all(x$epi$i.num[1, ] == 10))
  expect_true(max(x$epi$i.num) <= 100)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", sims = "mean", col.status = TRUE)
  test_net(x)
})


################################################################################

test_that("High departure rate models", {
  skip_on_cran()
  ## "netsim: 1M, ds.rate = 0.5"
  nw <- network_initialize(n = 25)
  est <- netest(nw, formation = ~edges, target.stats = 12,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0.01),
                edapprox = TRUE, verbose = FALSE)
  param <- param.net(inf.prob = 0.5, act.rate = 2,
                     a.rate = 0.01, ds.rate = 0.25,
                     di.rate = 0.1)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 25, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_equal(unique(sapply(x$epi, nrow)), 25)
  expect_output(summary(x, at = 25), "EpiModel Summary")

  capture_output(
    get_nwstats(x)
  )

  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "netsim: 1M, di.rate = 0.5"
  param <- param.net(inf.prob = 0.1, act.rate = 2, a.rate = 0.01,
                     ds.rate = 0.01, di.rate = 0.25)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 25, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_equal(unique(sapply(x$epi, nrow)), 25)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)

})

################################################################################

test_that("erroneous two-group models", {

  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "group", rep(0:1, each = 50))
  est <- netest(nw, formation = ~edges, target.stats = 25,
                 coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                 edapprox = TRUE, verbose = FALSE)

  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.1)
  init <- init.net(i.num = 10, i.num.g2 = 0)
  control <- control.net(type = "SI", nsims = 2, nsteps = 25, verbose = FALSE)
  expect_error(netsim(est, param, init, control))
})


################################################################################

test_that("edges two-group models", {
  skip_on_cran()

  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
  est5 <- netest(nw, formation = ~edges, target.stats = 25,
                 coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                 edapprox = TRUE, verbose = FALSE)
  expect_is(est5, "netest")

  ## "SI, 2M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.1)
  init <- init.net(i.num = 10, i.num.g2 = 0)
  control <- control.net(type = "SI", nsims = 2, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est5, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$groups, 2)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  plot(x, type = "network")
  plot(x, type = "network", shp.g2 = "triangle")
  #FLAG: Adding "two-group" attribute to network type
  expect_error(plot(x, type = "network", shp.g2 = TRUE))
  test_net(x)
  rm(x)

  ## "SIR, 2M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.1,
                     rec.rate = 0.1, rec.rate.g2 = 0.1)
  init <- init.net(i.num = 10, i.num.g2 = 10,
                   r.num = 0, r.num.g2 = 0)
  control <- control.net(type = "SIR", nsims = 2, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est5, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$groups, 2)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIR, 2M, CL: rec.rate = 0"
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.25,
                     rec.rate = 0, rec.rate.g2 = 0, act.rate = 2)
  init <- init.net(i.num = 10, i.num.g2 = 10,
                   r.num = 0, r.num.g2 = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 2,
                         verbose = FALSE)
  x <- netsim(est5, param, init, control)
  expect_equal(max(x$epi$ir.flow, na.rm = TRUE), 0)
  expect_equal(max(x$epi$ir.flow.g2, na.rm = TRUE), 0)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIS, 2M, CL: 2 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.25,
                     rec.rate = 0.01, rec.rate.g2 = 0.01)
  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SIS", nsims = 2, nsteps = 25,
                         verbose = FALSE)
  x <- netsim(est5, param, init, control)
  expect_is(x, "netsim")
  expect_is(as.data.frame(x), "data.frame")
  expect_equal(x$param$groups, 2)
  expect_output(summary(x, at = 25), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "SIS, 2M, CL: rec.rate = 0"
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.25,
                     rec.rate = 0, rec.rate.g2 = 0, act.rate = 2)
  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 2,
                         verbose = FALSE)
  x <- netsim(est5, param, init, control)
  expect_equal(max(x$epi$is.flow, na.rm = TRUE), 0)
  expect_equal(max(x$epi$is.flow.g2, na.rm = TRUE), 0)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

})

################################################################################

test_that("Open population 1 group models", {
  skip_on_cran()

  nw <- network_initialize(n = 100)
  est.vit <- netest(nw, formation = ~edges, target.stats = 25,
                    coef.diss = dissolution_coefs(~offset(edges), 10, 0.01),
                    verbose = FALSE)

  ## "SI, 1M, OP: 1 sim"
  param <- param.net(inf.prob = 0.5, act.rate = 2, a.rate = 0.01,
                     ds.rate = 0.01, di.rate = 0.01)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
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

  ## "SI, 1M, OP: 2 sim"
  param <- param.net(inf.prob = 0.5, act.rate = 2, a.rate = 0.01,
                     ds.rate = 0.01, di.rate = 0.01)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2,
                         resimulate.network = TRUE, verbose = FALSE)
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
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01,
                     dr.rate = 0.01)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
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

  ## "SIR, 1M, OP: 2 sim"
  param <- param.net(inf.prob = 0.5, rec.rate = 0.1, act.rate = 2,
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01,
                     dr.rate = 0.01)
  init <- init.net(i.num = 10, r.num = 0)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 2,
                         resimulate.network = TRUE, verbose = FALSE)
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
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
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
                         resimulate.network = TRUE, verbose = FALSE)
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

test_that("Open-population two-group models", {
  skip_on_cran()

  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
  est5.vit <- netest(nw, formation = ~edges + nodematch("group"),
                     target.stats = c(25, 0),
                     coef.diss = dissolution_coefs(~offset(edges), 10, 0.01),
                     edapprox = TRUE, verbose = FALSE)

  ## "SI, 2M, OP: 1 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.1, act.rate = 2,
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01,
                     a.rate.g2 = 0.01, ds.rate.g2 = 0.01, di.rate.g2 = 0.01)
  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
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


  ## "SI, 2M, OP: 2 sim"
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.1, act.rate = 2,
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01,
                     a.rate.g2 = 0.01, ds.rate.g2 = 0.01, di.rate.g2 = 0.01)
  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2,
                         resimulate.network = TRUE, verbose = FALSE)
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
  param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.1, rec.rate = 0.1,
                     rec.rate.g2 = 0.1, act.rate = 2, a.rate = 0.01,
                     a.rate.g2 = NA, ds.rate = 0.01, ds.rate.g2 = 0.01,
                     di.rate = 0.01, di.rate.g2 = 0.01, dr.rate = 0.01,
                     dr.rate.g2 = 0.01)
  init <- init.net(i.num = 10, i.num.g2 = 0,
                   r.num = 0, r.num.g2 = 10)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
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

  ## "SIR, 2M, OP: 3 sim"
  control <- control.net(type = "SIR", nsteps = 10, nsims = 2,
                         resimulate.network = TRUE, verbose = FALSE)
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

  nw <- network_initialize(n = 25)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, c(15, 10)))
  est <- netest(nw, formation = ~edges + nodematch("group"),
                target.stats = c(15, 0),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0.01),
                edapprox = TRUE, verbose = FALSE)

  ## "netsim: 2M, ds.rate = 0.5"
  param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.1, act.rate = 2,
                     a.rate = 0.01, ds.rate = 0.5, di.rate = 0.5,
                     a.rate.g2 = 0.01, ds.rate.g2 = 0.01, di.rate.g2 = 0.01)
  init <- init.net(i.num = 5, i.num.g2 = 0)
  control <- control.net(type = "SI", nsteps = 30, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

  ## "netsim: 2M, ds.rate.g2 = 0.5"
  param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.1, act.rate = 2,
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01,
                     a.rate.g2 = 0.01, ds.rate.g2 = 0.5, di.rate.g2 = 0.5)
  init <- init.net(i.num = 5, i.num.g2 = 0)
  control <- control.net(type = "SI", nsteps = 30, nsims = 1,
                         resimulate.network = TRUE, verbose = FALSE)
  x <- netsim(est, param, init, control)
  expect_output(summary(x, at = 10), "EpiModel Summary")
  plot(x)
  plot(x, y = "si.flow", mean.smooth = TRUE)
  plot(x, type = "formation")
  test_net(x)
  rm(x)

})


################################################################################

test_that("status.vector and infTime.vector", {

  n <- 100
  nw <- network_initialize(n = n)
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

  control <- control.net(type = "SIS", nsteps = 100, nsims = 5, verbose = FALSE)
  mod1 <- netsim(est1, param, init, control)
  expect_is(mod1, "netsim")

  control <- control.net(type = "SIR", nsteps = 100, nsims = 5, verbose = FALSE)
  mod2 <- netsim(est1, param, init, control)
  expect_is(mod2, "netsim")

})
