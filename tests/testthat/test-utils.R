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
  nw <- network_initialize(n = 100)
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

test_that("delete_attr", {
  dat <- list()

  l1 <- list(a = 1:5, b = 6:10)
  dat$attr <- l1
  dat <- delete_attr(dat, 5)
  expect_is(dat, "list")

  dat <- delete_attr(dat, 2:3)
  expect_true(length(unique(sapply(dat$attr, length))) == 1)

  l2 <- dat
  l3 <- delete_attr(dat, NULL)

  expect_equal(l2, l3)

  l4 <- list(a = 1:3, b = 5:20)
  dat$attr <- l4
  expect_error(delete_attr(dat, 2:4))
  expect_error(delete_attr(dat, 1))



})

test_that("ssample", {

  expect_true(length(ssample(1:5, 1)) == 1)
  expect_equal(ssample(5, 1), 5)
  expect_equal(ssample(5, 2), 5)
  expect_null(ssample(5, 0))

})

test_that("check_bip_degdist", {
  expect_output(check_bip_degdist(num.g1 = 500, num.g2 = 500,
                                  deg.dist.g2 = c(0.40, 0.55, 0.03, 0.02),
                                  deg.dist.g1 = c(0.48, 0.41, 0.08, 0.03)),
                "-0.015 Rel Diff")
  expect_output(check_bip_degdist(num.g1 = 500, num.g2 = 500,
                                  deg.dist.g1 = c(0.40, 0.55, 0.04, 0.01),
                                  deg.dist.g2 = c(0.48, 0.41, 0.08, 0.03)),
                "Edges balanced")
  expect_output(check_bip_degdist(num.g1 = 500, num.g2 = 500,
                                  deg.dist.g1 = c(0.45, 0.55, 0.04, 0.01),
                                  deg.dist.g2 = c(0.48, 0.41, 0.08, 0.03)),
                "deg.dist.g1 TOTAL != 1")
  expect_output(check_bip_degdist(num.g1 = 500, num.g2 = 500,
                                  deg.dist.g1 = c(0.40, 0.55, 0.04, 0.01),
                                  deg.dist.g2 = c(0.55, 0.41, 0.08, 0.03)),
                "deg.dist.g2 TOTAL != 1")
})

test_that("edgelist_censor", {
  skip_on_cran()
  nw <- network_initialize(n = 100)
  formation <- ~edges
  target.stats <- 50
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
  est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
  sim <- netdx(est, nsims = 1, nsteps = 100, keep.tedgelist = TRUE, verbose = FALSE)
  el <- as.data.frame(sim)
  expect_is(edgelist_censor(el), "matrix")
})

test_that("get_degree", {
  nw <- network_initialize(n = 500)

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
  #Dividing by zero:
  d.rate_ <- round(1-sqrt(59/60), 5)
  err.msg <- paste("The competing risk of departure is too high for the given",
                   " duration of ", 60, "; specify a d.rate lower than ",
                   d.rate_,".",sep="")
  dissolution = ~offset(edges)
  expect_that(dissolution_coefs(dissolution, duration = 60, d.rate = 1/60), throws_error(err.msg))
  expect_that(dissolution_coefs(dissolution, duration = 60, d.rate = 0.01), throws_error(err.msg))
  dissolution = ~offset(edges)+offset(nodematch("group", diff = TRUE))
  duration <- c(60, 30, 80, 100, 125, 160)
  dissolution = ~offset(edges) + offset(nodematch("age.grp", diff = TRUE))
  err.msg <- paste("The competing risk of departure is too high for the given",
                   "edge duration of 60 in place 1.",
                   "Specify a d.rate lower than 0.00837.")
  expect_that(dissolution_coefs(dissolution, duration = duration, d.rate = 1/60), throws_error(err.msg))
})


test_that("get_formula_term_attr checks", {
  nw <- network_initialize(n = 100)

  expect_null(get_formula_term_attr(~edges, nw))

  nw <- network_initialize(n = 100)
  riskg <- sample(rep(1:2, each = 50))
  race <- sample(rep(0:1, each = 50))
  nw <- set_vertex_attribute(nw, "riskg", riskg)
  nw <- set_vertex_attribute(nw, "race", race)
  expect_null(get_formula_term_attr(~edges, nw))

  expect_equal(get_formula_term_attr(~edges + nodefactor("race"), nw), "race")
  expect_equal(get_formula_term_attr(~edges + nodefactor("race") + nodematch("riskg"), nw),
               c("riskg", "race"))
  expect_equal(get_formula_term_attr(~edges + nodefactor(c("race", "riskg")) + nodematch("riskg"), nw),
               c("riskg", "race"))

})

test_that("Users using birth parameters are informed of change in language",{
  expect_error(param.net(b.rate = 2), paste0("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. See documentation for details."))
  expect_that(param.dcm(b.rate = 2), shows_message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. See documentation for details."))
  expect_that(param.icm(b.rate = 2), shows_message("EpiModel 1.7.0 onward renamed the birth rate parameter b.rate to a.rate. See documentation for details."))
})

test_that("Users using birth and death functions are informed of change in language",{
  temp <- function(x){x = x; return(x)}
  expect_error(control.net(type = NULL, nsims = 1, nsteps = 10,
                          departures.FUN = temp, arrivals.FUN = temp,
                          prevalence.FUN = prevalence.net, infection.FUN = infection.net,
                          recovery.FUN = recovery.net, births.FUN = temp, resimulate.network = FALSE), paste0("EpiModel 1.7.0 onward renamed the birth function births.FUN to arrivals.FUN. See documentation for details."))
  expect_error(control.net(type = NULL, nsims = 1, nsteps = 10,
                          departures.FUN = temp, arrivals.FUN = temp,
                          prevalence.FUN = prevalence.net, infection.FUN = infection.net,
                          recovery.FUN = recovery.net, deaths.FUN = temp, resimulate.network = FALSE), paste0("EpiModel 1.7.0 onward renamed the death function deaths.FUN to departures.FUN. See documentation for details."))
})

test_that("Users using old mode syntax are informed of change to group syntax", {

  err.param <- paste0("EpiModel 2.0 onward has updated parameter suffixes reflecting a move from mode to group networks. ",
                      "All .m2 parameters changed to .g2. See documentation for more details.")
  expect_warning(param <- param.net(inf.prob = 0.1, inf.prob.m2 = 0.1, act.rate = 1),
                 err.param)
  err.init <- paste0("EpiModel 2.0 onward has updated initial condition suffixes reflecting a move from mode to group networks. ",
                    "All .m2 initial conditions changed to .g2. See documentation for more details.")
  expect_warning(init <- init.net(i.num = 1, i.num.m2 = 1), err.init)
})

