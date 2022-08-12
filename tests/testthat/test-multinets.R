context("multiple network models")

test_that("netsim runs with multiple networks", {
  nw <- network_initialize(n = 50)
  nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(nw, formation = ~edges + nodematch("race"),
                target.stats = c(25, 10),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3, act.rate = 0.1)
  init <- init.net(i.num = 1)
  control <- control.net(type = "SI", 
                         nsims = 5, 
                         nsteps = 15, 
                         nwstats.formula = ~edges+nodematch("race")+degree(0:10),
                         verbose = TRUE)

  set.seed(0)
  mod <- netsim(list(est, est), param, init, control)  

  print(mod)

  plot(mod, type = "formation")
  plot(mod, type = "dissolution")
  plot(mod, type = "duration")

  plot(mod, type = "formation", network = 2)
  plot(mod, type = "dissolution", network = 2)
  plot(mod, type = "duration", network = 2)
})
