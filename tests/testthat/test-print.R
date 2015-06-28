context("print functions")

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
  expect_output(print(mod), "Compartments: s.num i.num num")
  expect_output(print(mod), "Flows: si.flow")
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
  expect_output(print(mod), "Compartments: s.num i.num num")
  expect_output(print(mod), "Flows: si.flow")
})

test_that("print.netest", {
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  expect_output(print(est), "EpiModel Network Estimation")
  expect_output(print(est), "Model class: netest")
  expect_output(print(est), "Estimation Method: ERGM with Edges Approximation")
  expect_output(print(est), "Formation: ~edges")
  expect_output(print(est), "Formation Targets: 50")
  expect_output(print(est), "Edge Duration Target: 10")
})

test_that("print.netsim", {
  nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
  init <- init.net(i.num = 10, i.num.m2 = 10)
  control <- control.net(type = "SI", nsteps = 10, nsims = 1, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_output(print(mod), "EpiModel Simulation")
  expect_output(print(mod), "Model class: netsim")
  expect_output(print(mod), "Model type: SI")
  expect_output(print(mod), "No. NW modes: 2")
  expect_output(print(mod), "Compartments: s.num i.num num s.num.m2 i.num.m2 num.m2")
  expect_output(print(mod), "Flows: si.flow si.flow.m2")
})

test_that("print.disscoefs", {
  o <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
  expect_output(print(o), "Dissolution Coefficients")
  expect_output(print(o), "Crude Coefficient: 3.178054")
  expect_output(print(o), "Death rate: 0")

  o <- dissolution_coefs(dissolution = ~offset(edges), duration = 25,
                         d.rate = 0.001)
  expect_is(o, "disscoef")
  expect_output(print(o), "Death rate: 0.001")

  o <- dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
                         duration = c(20, 10))
  expect_output(print(o), "2.944439 -0.7472144")
  expect_output(print(o), "Death rate: 0")
  expect_output(print(o), "Edge Duration: 20 10")

  o <- dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("race")),
                         duration = c(20, 10), d.rate = 0.001)
  expect_output(print(o), "Crude Coefficient: 2.944439 -0.7472144")
  expect_output(print(o), "Adjusted Coefficient: 2.98524 -0.7678231")
})
