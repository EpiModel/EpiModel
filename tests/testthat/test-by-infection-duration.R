context("Parameters by Duration of Infection")

# Shared fixtures ------------------------------------------------------------

nw <- network_initialize(n = 100)
est <- netest(nw, formation = ~edges, target.stats = 60,
              coef.diss = dissolution_coefs(~offset(edges), 10),
              verbose = FALSE)
init <- init.net(i.num = 20)


# Constructor ------------------------------------------------------------------

test_that("by_infection_duration builds, prints, and validates", {
  x <- by_infection_duration(c(0.5, 0.5, 0.1))
  expect_s3_class(x, "by_infection_duration")
  expect_equal(as.numeric(x), c(0.5, 0.5, 0.1))
  expect_equal(attr(x, "durations"), c(1, 1, Inf))
  expect_equal(format(x), "by_infection_duration(c(0.5, 0.5, 0.1))")
  expect_output(print(x), "duration of infection")
  expect_s3_class(by_infection_duration(1:3), "by_infection_duration")
  expect_s3_class(by_infection_duration(0.2), "by_infection_duration")

  expect_error(by_infection_duration("a"), "numeric vector")
  expect_error(by_infection_duration(numeric(0)), "numeric vector")
  expect_error(by_infection_duration(c(0.1, NA)), "numeric vector")

  out <- capture.output(print(param.net(inf.prob = x, act.rate = 1)))
  expect_true(any(grepl("inf.prob = by_infection_duration(c(0.5, 0.5, 0.1))",
                        out, fixed = TRUE)))
})

test_that("by_infection_duration takes values by stage", {
  x <- by_infection_duration(c(acute = 0.25, chronic = 0.02),
                             durations = c(10, Inf))
  expect_equal(as.numeric(x), c(0.25, 0.02))
  expect_equal(names(x), c("acute", "chronic"))
  expect_equal(attr(x, "durations"), c(10, Inf))
  expect_equal(format(x), paste0("by_infection_duration(c(acute = 0.25, ",
                                 "chronic = 0.02), durations = c(10, Inf))"))
  out <- capture.output(print(x))
  expect_true(any(grepl("acute +1-10 +0.25", out)))
  expect_true(any(grepl("chronic +11\\+ +0.02", out)))
  out <- capture.output(print(by_infection_duration(c(0.1, 0.2, 0.3),
                                                    durations = c(1, 5, Inf))))
  expect_true(any(grepl("^ +1 +1 +0.1$", out)))
  expect_true(any(grepl("^ +2 +2-6 +0.2$", out)))
  expect_true(any(grepl("^ +3 +7\\+ +0.3$", out)))

  # the stage form and the per-step form give the same values
  y <- by_infection_duration(c(rep(0.25, 10), 0.02))
  expect_equal(infection_duration_value(x, 1:40, "inf.prob"),
               infection_duration_value(y, 1:40, "inf.prob"))
  expect_equal(infection_duration_value(x, c(1, 10, 11, 100), "inf.prob"),
               c(0.25, 0.25, 0.02, 0.02))
  expect_equal(as_duration_vector(x), c(rep(0.25, 10), 0.02))
  one <- by_infection_duration(0.3, durations = Inf)
  expect_equal(infection_duration_value(one, 1:5, "inf.prob"), rep(0.3, 5))

  # arithmetic scales the values and keeps the stages
  x2 <- x * 2
  expect_s3_class(x2, "by_infection_duration")
  expect_equal(attr(x2, "durations"), c(10, Inf))
  expect_equal(infection_duration_value(x2, c(10, 11), "inf.prob"),
               c(0.5, 0.04))

  out <- capture.output(print(param.net(inf.prob = x, act.rate = 1)))
  expect_true(any(grepl("durations = c(10, Inf))", out, fixed = TRUE)))

  expect_error(by_infection_duration(c(1, 2), durations = c(10, 5)),
               "must be Inf")
  expect_error(by_infection_duration(c(1, 2), durations = c(0, Inf)),
               "whole number of time steps")
  expect_error(by_infection_duration(c(1, 2), durations = c(1.5, Inf)),
               "whole number of time steps")
  expect_error(by_infection_duration(c(1, 2), durations = c(Inf, Inf)),
               "whole number of time steps")
  expect_error(by_infection_duration(c(1, 2), durations = 10),
               "one duration per value")
  expect_error(by_infection_duration(c(1, 2), durations = c(NA, Inf)),
               "one duration per value")
})


# Input checks -----------------------------------------------------------------

test_that("the built-in model types refuse a plain vector", {
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, verbose = FALSE)
  expect_error(netsim(est, param.net(inf.prob = c(0.5, 0.1), act.rate = 1),
                      init, control),
               "`inf.prob` is a plain vector of length 2")
  expect_error(netsim(est, param.net(inf.prob = 0.3, act.rate = c(2, 1)),
                      init, control),
               "`act.rate` is a plain vector of length 2")
  expect_error(netsim(list(est, est),
                      param.net(inf.prob = multilayer(c(0.5, 0.1), 0.1),
                                act.rate = 1),
                      init, control),
               "`inf.prob` \\(layer 1\\) is a plain vector")

  control_sis <- control.net(type = "SIS", nsteps = 5, nsims = 1,
                             verbose = FALSE)
  expect_error(netsim(est, param.net(inf.prob = 0.3, act.rate = 1,
                                     rec.rate = c(0.1, 0.2)),
                      init, control_sis),
               "`rec.rate` is a plain vector")
  expect_error(netsim(est, param.net(inf.prob = 0.3, act.rate = 1,
                                     rec.rate = multilayer(0.1)),
                      init, control_sis),
               "cannot be a multilayer")
})

test_that("custom modules keep their own reading of plain vectors", {
  skip_on_cran()
  # a custom infection module that reads inf.prob by layer, as EpiModelCOVID
  # does, is not affected
  by_layer <- function(dat, at) {
    dat <- set_epi(dat, "p2", at, get_param(dat, "inf.prob")[2])
    return(dat)
  }
  control <- control.net(type = NULL, nsteps = 3, nsims = 1,
                         infection.FUN = by_layer, verbose = FALSE)
  sim <- netsim(est, param.net(inf.prob = c(0.3, 0.1), act.rate = 1), init,
                control)
  expect_equal(sim$epi$p2[3, 1], 0.1)

  # the built-in infection module used from an extension model still stops
  control <- control.net(type = NULL, nsteps = 3, nsims = 1,
                         infection.FUN = infection.net, verbose = FALSE)
  expect_error(netsim(est, param.net(inf.prob = c(0.3, 0.1), act.rate = 1),
                      init, control),
               "`inf.prob` is a plain vector")
})

test_that("DCM and ICM refuse by_infection_duration parameters", {
  expect_error(dcm(param.dcm(inf.prob = by_infection_duration(c(0.2, 0.1)),
                             act.rate = 1),
                   init.dcm(s.num = 100, i.num = 1),
                   control.dcm(type = "SI", nsteps = 5)),
               "netsim\\(\\) only")
  expect_error(icm(param.icm(inf.prob = by_infection_duration(c(0.2, 0.1)),
                             act.rate = 1),
                   init.icm(s.num = 100, i.num = 1),
                   control.icm(type = "SI", nsteps = 5, nsims = 1)),
               "netsim\\(\\) only")
})


# Simulation -------------------------------------------------------------------

test_that("transmission follows the infected partner's duration of infection", {
  skip_on_cran()
  # seed infections at t = 1, so that every seed starts inside the window in
  # which inf.prob is 1 (with i.num, infTime is drawn over past steps, and about
  # one fixture network in five gave no transmission)
  status <- rep("s", 100)
  status[seq(1, 100, by = 5)] <- "i"
  init_t1 <- init.net(status.vector = status,
                      infTime.vector = ifelse(status == "i", 1, NA))
  for (tergmLite in c(TRUE, FALSE)) {
    # certain in the first two steps of infection, impossible after
    set.seed(31)
    param <- param.net(inf.prob = by_infection_duration(c(1, 0),
                                                        durations = c(2, Inf)),
                       act.rate = 1)
    control <- control.net(type = "SI", nsteps = 15, nsims = 1,
                           tergmLite = tergmLite, resimulate.network = TRUE,
                           verbose = FALSE)
    sim <- netsim(est, param, init_t1, control)
    tm <- get_transmat(sim)
    expect_gt(nrow(tm), 0)
    expect_true(all(tm$infDur <= 2))
    expect_true(all(tm$transProb == 1))

    # no acts in the first step of infection, two acts after
    set.seed(32)
    param <- param.net(inf.prob = 0.5,
                       act.rate = by_infection_duration(c(0, 2)))
    sim <- netsim(est, param, init, control)
    tm <- get_transmat(sim)
    expect_gt(nrow(tm), 0)
    expect_true(all(tm$infDur >= 2))
    expect_true(all(tm$actRate == 2))
  }
})

test_that("recovery follows the duration of infection", {
  skip_on_cran()
  # recovery is certain in the third step of infection and impossible before,
  # so at the end every node still infected is in its first or second step
  set.seed(33)
  param <- param.net(inf.prob = 0.5, act.rate = 1,
                     rec.rate = by_infection_duration(c(0, 1),
                                                      durations = c(2, Inf)))
  control <- control.net(type = "SIR", nsteps = 15, nsims = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         verbose = FALSE, save.run = TRUE)
  sim <- netsim(est, param, init.net(i.num = 20, r.num = 0), control)
  attr <- sim$run[[1]]$attr
  inf_dur <- pmax(15 - attr$infTime[attr$status == "i"], 1)
  expect_gt(sum(attr$status == "r"), 0)
  expect_true(all(inf_dur <= 2))

  # two groups: group 1 never recovers, group 2 recovers in its second step
  nw2 <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
  est2 <- netest(nw2, formation = ~edges, target.stats = 60,
                 coef.diss = dissolution_coefs(~offset(edges), 10),
                 verbose = FALSE)
  set.seed(34)
  param2 <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.5, act.rate = 1,
                      rec.rate = 0,
                      rec.rate.g2 = by_infection_duration(c(0, 1)))
  sim2 <- netsim(est2, param2, init.net(i.num = 10, i.num.g2 = 10, r.num = 0, r.num.g2 = 0), control)
  attr2 <- sim2$run[[1]]$attr
  grp <- attr2$group
  expect_false(any(attr2$status[grp == 1] == "r"))
  expect_gt(sum(attr2$status[grp == 2] == "r"), 0)
  still_inf_g2 <- grp == 2 & attr2$status == "i"
  expect_true(all(pmax(15 - attr2$infTime[still_inf_g2], 1) <= 1))
})

test_that("an updater can switch a parameter to vary by duration of infection", {
  skip_on_cran()
  set.seed(35)
  updater <- list(list(at = 6, verbose = FALSE,
                       param = list(inf.prob = by_infection_duration(c(0.9, 0.1)))))
  param <- param.net(inf.prob = 0.3, act.rate = 1,
                     .param.updater.list = updater)
  control <- control.net(type = "SI", nsteps = 15, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE, verbose = FALSE)
  sim <- netsim(est, param, init, control)
  tm <- get_transmat(sim)
  before <- tm[tm$at < 6, ]
  after <- tm[tm$at >= 6, ]
  expect_gt(nrow(after), 0)
  expect_true(all(before$transProb == 0.3))
  expect_true(all(after$transProb == ifelse(after$infDur == 1, 0.9, 0.1)))
})
