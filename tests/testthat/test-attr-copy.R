context("Attribute copying between network and dat$attr")

################################################################################

test_that("Copying attributes from network to attribute list",{

  num1 <- num2 <- 500
  nw <- network_initialize(num1 + num2)
  nw <- set_vertex_attribute(nw, "group", rep(1:2, each = num1))
  nw <- set_vertex_attribute(nw, "race", sample(c("B","W"), num1 + num2, replace = TRUE))
  nw <- set_vertex_attribute(nw, "region", sample(1:4, num1 + num2, replace = TRUE))
  formation <- ~edges + nodematch("group")
  target.stats <- c(400, 0)
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
  est <- netest(nw, formation, target.stats, coef.diss)

  param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                     act.rate = 5)
  init <- init.net(i.num = 50, i.num.g2 = 50)
  control <- control.net(type = "SI", nsteps = 10, nsims = 2, tergmLite = FALSE,
                         raw.output = TRUE, verbose = FALSE)

  sim <- netsim(est, param, init, control)

  # Character attribute
  dat.attr <- prop.table(table(sim[[1]]$attr$race))
  nw.attr <- prop.table(table(get_vertex_attribute(sim[[1]]$nw[[1]], "race")))

  expect_equal(dat.attr, nw.attr)

  # Numeric attribute
  dat.attr <- prop.table(table(sim[[1]]$attr$region))
  nw.attr <- prop.table(table(get_vertex_attribute(sim[[1]]$nw[[1]], "region")))

  expect_equal(dat.attr, nw.attr)

  # Second simulation
  dat.attr <- prop.table(table(sim[[2]]$attr$race))
  nw.attr <- prop.table(table(get_vertex_attribute(sim[[2]]$nw[[1]], "race")))

  expect_equal(dat.attr, nw.attr)

  dat.attr <- prop.table(table(sim[[2]]$attr$region))
  nw.attr <- prop.table(table(get_vertex_attribute(sim[[2]]$nw[[1]], "region")))

  expect_equal(dat.attr, nw.attr)

})
