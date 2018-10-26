context("Utility Functions")

test_that("brewer_ramp", {

  expect_true(length(brewer_ramp(100, plt = "Spectral")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Spectral", delete.lights = FALSE)) == 100)
  expect_false(length(brewer_ramp(100, plt = "Spectral")) == 50)
  expect_false(length(brewer_ramp(100, plt = "Spectral", delete.lights = FALSE)) == 50)

  expect_true(length(brewer_ramp(100, plt = "Accent")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Accent", delete.lights = FALSE)) == 100)

  expect_true(length(brewer_ramp(100, plt = "Oranges")) == 100)
  expect_true(length(brewer_ramp(100, plt = "Oranges", delete.lights = FALSE)) == 100)

  expect_true(length(brewer_ramp(100, plt = "Set1")) == 100)

  expect_error(brewer_ramp(100, plt = "Jimmy"))
  expect_error(brewer_ramp(100, plt = "Jimmy", delete.lights = FALSE))
  expect_error(brewer_ramp(-1, plt = "Spectral"))

})

test_that("transco", {

  cols <- transco(c("steelblue", "black"), 0.5)
  cols2 <- transco(1, c(0.5, 1))

  expect_is(class(cols), "character")
  expect_true(length(cols) == 2)
  expect_is(class(cols2), "character")
  expect_true(length(cols2) == 2)
  expect_error(transco(1:2, c(0.2, 0.3)))
  expect_is(class(transco(1, 1)), "character")
  expect_is(transco(1, 1, invisible = TRUE), "character")
  expect_error(transco("bob", 1))

})

test_that("color_tea", {
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  param <- param.net(inf.prob = 0.3, rec.rate = 0.001)
  init <- init.net(i.num = 25, r.num = 25)
  control <- control.net(type = "SIR", nsteps = 10, nsims = 1, verbose.int = 0)
  mod <- netsim(est, param, init, control)
  nd <- get_network(mod)
  nd <- color_tea(nd)
  expect_is(nd, "networkDynamic")
  expect_true(length(unique(get.vertex.attribute.active(nd, "ndtvcol", at = 1))) == 3)
})

test_that("calc_eql for dcm", {
  param <- param.dcm(inf.prob = 0.2, inf.prob.g2 = 0.1, act.rate = 0.5,
                     balance = "g1", rec.rate = 1 / 50, rec.rate.g2 = 1 / 50,
                     b.rate = 1 / 100, b.rate.g2 = NA, ds.rate = 1 / 100,
                     ds.rate.g2 = 1 / 100, di.rate = 1 / 90, di.rate.g2 = 1 / 90)
  init <- init.dcm(s.num = 500, i.num = 1,
                   s.num.g2 = 500, i.num.g2 = 1)
  control <- control.dcm(type = "SIS", nsteps = 500, verbose = FALSE)
  x <- dcm(param, init, control)
  expect_output(calc_eql(x, nsteps = 100), "Rel. Diff.:    0.0001")
  expect_output(calc_eql(x, nsteps = 100, numer = "i.num.g2", "num.g2"),
                "Start Prev.:   0.5028")
  expect_error(calc_eql(x, nsteps = 100, numer = "foo", denom = "num"),
               "numer must be an output compartment on x")
  expect_error(calc_eql(x, nsteps = 100, numer = "i.num", denom = "bar"),
               "denom must be an output compartment on x")

})

test_that("calc_eql for icm", {
  skip_on_cran()
  set.seed(1)
  param <- param.icm(inf.prob = 0.2, act.rate = 0.25, b.rate = 1/100,
                     ds.rate = 1/100, di.rate = 1/90)
  init <- init.icm(s.num = 500, i.num = 1)
  control <- control.icm(type = "SI", nsteps = 500, nsims = 1, verbose = FALSE)
  x <- icm(param, init, control)
  expect_output(calc_eql(x, nsteps = 100), "Start Prev.:   0.7607")
  expect_output(calc_eql(x, nsteps = 100), "<= Threshold:  FALSE")
})

test_that("calc_eql for netsim", {
  skip_on_cran()
  set.seed(35160)
  nw <- network.initialize(n = 50, directed = FALSE)
  nw <- set.vertex.attribute(nw, "race", rbinom(50, 1, 0.5))
  est <- netest(nw, formation = ~edges + nodematch("race"), target.stats = c(25, 10),
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.5, act.rate = 1)
  init <- init.net(i.num = 25)
  control <- control.net(type = "SI", nsims = 1, nsteps = 100, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_output(calc_eql(mod, nsteps = 20), "Start Prev.:   1")
  expect_output(calc_eql(mod, nsteps = 20), "Rel. Diff.:    0")
})

test_that("deleteAttr", {

  l <- list(a = 1:5, b = 6:10)
  expect_is(deleteAttr(l, 5), "list")
  expect_true(length(unique(sapply(deleteAttr(l, 2:3), length))) == 1)

  l2 <- list(a = 1:3, b = 5:20)
  expect_error(deleteAttr(l2, 2:4))
  expect_error(deleteAttr(as.data.frame(l), 1))

  expect_equal(l, deleteAttr(l, NULL))

})

test_that("ssample", {

  expect_true(length(ssample(1:5, 1)) == 1)
  expect_equal(ssample(5, 1), 5)
  expect_equal(ssample(5, 2), 5)
  expect_null(ssample(5, 0))

})

test_that("bipvals", {
  nw <- network.initialize(n = 10, bipartite = 5)
  nw <- set.vertex.attribute(nw, "male", rep(0:1, each = 5))
  expect_true(all(bipvals(nw, mode = 1, "male")) == 0)
  expect_true(all(bipvals(nw, mode = 2, "male")) == 1)
  expect_error(bipvals(nw, val = "male"))

  nw <- network.initialize(n = 10)
  nw <- set.vertex.attribute(nw, "male", rep(0:1, each = 5))
  expect_error(bipvals(nw, 1, "male"), "nw must be a bipartite network")
})

test_that("check_bip_degdist", {
  expect_output(check_bip_degdist(num.m1 = 500, num.m2 = 500,
                    deg.dist.m2 = c(0.40, 0.55, 0.03, 0.02),
                    deg.dist.m1 = c(0.48, 0.41, 0.08, 0.03)),
                "-0.015 Rel Diff")
  expect_output(check_bip_degdist(num.m1 = 500, num.m2 = 500,
                    deg.dist.m1 = c(0.40, 0.55, 0.04, 0.01),
                    deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03)),
                "Edges balanced")
  expect_output(check_bip_degdist(num.m1 = 500, num.m2 = 500,
                                  deg.dist.m1 = c(0.45, 0.55, 0.04, 0.01),
                                  deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03)),
                "deg.dist.m1 TOTAL != 1")
  expect_output(check_bip_degdist(num.m1 = 500, num.m2 = 500,
                                  deg.dist.m1 = c(0.40, 0.55, 0.04, 0.01),
                                  deg.dist.m2 = c(0.55, 0.41, 0.08, 0.03)),
                "deg.dist.m2 TOTAL != 1")
})

test_that("edgelist_censor", {
  skip_on_cran()
  nw <- network.initialize(n = 100, directed = FALSE)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  sim <- netdx(est, nsims = 1, nsteps = 100, keep.tedgelist = TRUE, verbose = FALSE)
  el <- as.data.frame(sim)
  expect_is(edgelist_censor(el), "matrix")
})

test_that("get_degree", {
  nw <- network.initialize(500, directed = FALSE)

  set.seed(1)
  fit <- ergm(nw ~ edges, target.stats = 250, eval.loglik = FALSE)
  sim <- simulate(fit)

  # Slow ERGM-based method
  ergm.method <- unname(summary(sim ~ sociality(base = 0)))
  ergm.method

  # Fast tabulate method with network object
  deg.net <- get_degree(sim)
  deg.net

  # Even faster if network already transformed into an edgelist
  el <- as.edgelist(sim)
  deg.el <- get_degree(el)
  deg.el

  all.equal(ergm.method, deg.net, deg.el)
})

test_that("dissolution_coefs returns appropriate error for incompatible departure rate",{
  dissolution = ~offset(edges)
  d.rate_ <- round(1-sqrt(59/60), 5)
  duration <- 60
  err.msg <- paste("The competing risk of mortality is too high for the given", 
                   " duration of ", duration[1], "; specify a d.rate lower than ", 
                   d.rate_,".",sep="")
  expect_that(dissolution_coefs(dissolution, duration = 60, d.rate = 1/60), throws_error(err.msg))
  expect_that(dissolution_coefs(dissolution, duration = 60, d.rate = 0.01), throws_error(err.msg))
})
