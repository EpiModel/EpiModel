context("Print functions")

test_that("print.dcm", {
  param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.dcm(s.num = 500, i.num = 1)
  control <- control.dcm(type = "SI", nsteps = 500)
  mod <- dcm(param, init, control)
  expect_output(print(mod), "EpiModel Simulation")
  expect_output(print(mod), "Model class: dcm")
  expect_output(print(mod), "Model type: SI")
  expect_output(print(mod), "inf.prob = 0.2")
  expect_output(print(mod), "act.rate = 0.25")
  expect_output(print(mod), "Variables: s.num i.num si.flow num")
})

test_that("print.icm", {
  param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SI", nsteps = 10, nsims = 2, verbose = FALSE)
  mod <- icm(param, init, control)
  expect_output(print(mod), "EpiModel Simulation")
  expect_output(print(mod), "Model class: icm")
  expect_output(print(mod), "Model type: SI")
  expect_output(print(mod), "inf.prob = 0.2")
  expect_output(print(mod), "act.rate = 0.25")
  expect_output(print(mod), "Variables: s.num i.num num si.flow")
})

test_that("print.netest", {
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  expect_output(print(est), "EpiModel Network Estimation")
  expect_output(print(est), "Model class: netest")
  expect_output(print(est), "Estimation Method: ERGM with Edges Approximation")
  expect_output(print(est), "Formation: ~edges")
  expect_output(print(est), "Target Statistics: 50")
  expect_output(print(est), "Target Statistics: 10")
})

test_that("print.netsim", {
  nw <- network_initialize(n = 100)
  nw <- set_vertex_attribute(nw, "group", rep(c(1, 2), each = 50))
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15)
  init <- init.net(i.num = 10, i.num.g2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_output(print(mod), "EpiModel Simulation")
  expect_output(print(mod), "Model class: netsim")
  expect_output(print(mod), "Model type: SI")
  expect_output(print(mod), "No. NW groups: 2")
  expect_output(print(mod),
                "Variables: s.num i.num num s.num.g2 i.num.g2 num.g2")
})

test_that("print.disscoefs", {
  o <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
  expect_output(print(o), "Dissolution Coefficients")
  expect_output(print(o), "Crude Coefficient: 3.178054")
  expect_output(print(o), "Mortality/Exit Rate: 0")

  o <- dissolution_coefs(dissolution = ~offset(edges), duration = 25,
                         d.rate = 0.001)
  expect_is(o, "disscoef")
  expect_output(print(o), "Mortality/Exit Rate: 0.001")

  o <- dissolution_coefs(dissolution =
                           ~offset(edges) + offset(nodematch("race")),
                         duration = c(20, 10))
  expect_output(print(o), "2.944439 -0.7472144")
  expect_output(print(o), "Mortality/Exit Rate: 0")
  expect_output(print(o), "Target Statistics: 20 10")

  o <- dissolution_coefs(dissolution =
                           ~offset(edges) + offset(nodematch("race")),
                         duration = c(20, 10), d.rate = 0.001)
  expect_output(print(o), "Crude Coefficient: 2.944439 -0.7472144")
  expect_output(print(o), "Adjusted Coefficient: 2.98524 -0.7678231")
})

test_that("print.param", {
  p <- param.dcm(inf.prob = 0.1, rec.rate = 0.1)
  expect_output(print(p), "DCM Parameters")

  p <- param.icm(inf.prob = 0.1, rec.rate = 0.1)
  expect_output(print(p), "ICM Parameters")

  p <- param.net(inf.prob = 0.1, rec.rate = 0.1)
  expect_output(print(p), "Network Model Parameters")
  expect_output(print(p), "act.rate = 1")
})

test_that("print.init", {
  i <- init.dcm(s.num = 10, i.num = 10)
  expect_output(print(i), "DCM Initial Conditions")

  i <- init.icm(s.num = 10, i.num = 10)
  expect_output(print(i), "ICM Initial Conditions")

  i <- init.net(s.num = 10, i.num = 10)
  expect_output(print(i), "Network Model Initial Conditions")
})

test_that("print.control", {
  co <- control.dcm(type = "SI", nsteps = 10)
  expect_output(print(co), "DCM Control Settings")
  expect_output(print(co), "odemethod = rk4")

  co <- control.icm(type = "SI", nsteps = 10)
  expect_output(print(co), "ICM Control Settings")
  expect_output(print(co), "Base Modules: initialize.FUN")

  co <- control.net(type = "SI", nsteps = 10)
  expect_output(print(co), "Network Model Control Settings")
  #FLAG 4/23
  #expect_output(print(co), "Base Modules: initialize.FUN")
})
