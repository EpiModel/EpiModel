context("comp_plot.dcm")

run_dcm <- function(type = "SI", vital = FALSE, nsteps = 12) {
  if (vital) {
    if (type == "SIR") {
      p <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 0.05,
                     a.rate = 0.001, ds.rate = 0.001,
                     di.rate = 0.001, dr.rate = 0.001)
      i <- init.dcm(s.num = 100, i.num = 1, r.num = 0)
    } else if (type == "SIS") {
      p <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 0.05,
                     a.rate = 0.001, ds.rate = 0.001, di.rate = 0.001)
      i <- init.dcm(s.num = 100, i.num = 1)
    } else {
      p <- param.dcm(inf.prob = 0.2, act.rate = 1,
                     a.rate = 0.001, ds.rate = 0.001, di.rate = 0.001)
      i <- init.dcm(s.num = 100, i.num = 1)
    }
  } else {
    if (type == "SIR") {
      p <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 0.05)
      i <- init.dcm(s.num = 100, i.num = 1, r.num = 0)
    } else if (type == "SIS") {
      p <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 0.05)
      i <- init.dcm(s.num = 100, i.num = 1)
    } else {
      p <- param.dcm(inf.prob = 0.2, act.rate = 1)
      i <- init.dcm(s.num = 100, i.num = 1)
    }
  }
  control <- control.dcm(type = type, nsteps = nsteps)
  dcm(p, i, control)
}

test_that("comp_plot.dcm draws SI/SIR/SIS with and without vital dynamics", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  for (type in c("SI", "SIR", "SIS")) {
    for (vital in c(FALSE, TRUE)) {
      mod <- run_dcm(type = type, vital = vital)
      expect_silent(comp_plot(mod, at = 10))
    }
  }
})

test_that("comp_plot.dcm validates at and rejects two-group models", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  mod <- run_dcm("SI")
  expect_error(comp_plot(mod, at = 0),
               "Specify a time step between 1 and 12")
  expect_error(comp_plot(mod, at = 999),
               "Specify a time step between 1 and 12")

  # Two-group DCM should hit the 1-group guard.
  param <- param.dcm(inf.prob = 0.3, inf.prob.g2 = 0.2, act.rate = 1,
                     balance = "g1")
  init <- init.dcm(s.num = 100, i.num = 1,
                   s.num.g2 = 100, i.num.g2 = 1)
  control <- control.dcm(type = "SI", nsteps = 10)
  mod2g <- dcm(param, init, control)
  expect_error(comp_plot(mod2g, at = 5),
               "Only 1-group")
})

context("comp_plot.icm and comp_plot.netsim")

run_icm <- function(type = "SI", vital = FALSE, nsteps = 8) {
  if (vital) {
    if (type == "SIR") {
      p <- param.icm(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.1,
                     a.rate = 0.001, ds.rate = 0.001,
                     di.rate = 0.001, dr.rate = 0.001)
      i <- init.icm(s.num = 50, i.num = 2, r.num = 0)
    } else if (type == "SIS") {
      p <- param.icm(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.1,
                     a.rate = 0.001, ds.rate = 0.001, di.rate = 0.001)
      i <- init.icm(s.num = 50, i.num = 2)
    } else {
      p <- param.icm(inf.prob = 0.3, act.rate = 0.5,
                     a.rate = 0.001, ds.rate = 0.001, di.rate = 0.001)
      i <- init.icm(s.num = 50, i.num = 2)
    }
  } else {
    if (type == "SIR") {
      p <- param.icm(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.1)
      i <- init.icm(s.num = 50, i.num = 2, r.num = 0)
    } else if (type == "SIS") {
      p <- param.icm(inf.prob = 0.3, act.rate = 0.5, rec.rate = 0.1)
      i <- init.icm(s.num = 50, i.num = 2)
    } else {
      p <- param.icm(inf.prob = 0.3, act.rate = 0.5)
      i <- init.icm(s.num = 50, i.num = 2)
    }
  }
  control <- control.icm(type = type, nsteps = nsteps, nsims = 2,
                         verbose = FALSE)
  icm(p, i, control)
}

test_that("comp_plot.icm draws SI/SIR/SIS with and without vital dynamics", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  for (type in c("SI", "SIR", "SIS")) {
    for (vital in c(FALSE, TRUE)) {
      mod <- run_icm(type = type, vital = vital)
      expect_silent(comp_plot(mod, at = 5))
    }
  }
})

test_that("comp_plot.icm validates at and rejects two-group models", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  mod <- run_icm("SI")
  expect_error(comp_plot(mod, at = 0),
               "Specify a timestep between 1 and 8")
  expect_error(comp_plot(mod, at = 999),
               "Specify a timestep between 1 and 8")

  param2 <- param.icm(inf.prob = 0.3, inf.prob.g2 = 0.2,
                      act.rate = 0.5, balance = "g1")
  init2 <- init.icm(s.num = 50, i.num = 2, s.num.g2 = 50, i.num.g2 = 1)
  control2 <- control.icm(type = "SI", nsteps = 5, nsims = 1, verbose = FALSE)
  mod2g <- icm(param2, init2, control2)
  expect_error(comp_plot(mod2g, at = 3),
               "Only 1-group")
})

test_that("comp_plot.netsim delegates to comp_plot.icm without error", {
  skip_on_cran()
  nw <- network_initialize(n = 30)
  est <- netest(nw, formation = ~edges, target.stats = 10,
                coef.diss = dissolution_coefs(~offset(edges), 10, 0),
                verbose = FALSE)
  control <- control.net(type = "SIS", nsteps = 5, nsims = 2,
                         verbose = FALSE)
  mod <- netsim(est,
                param.net(inf.prob = 0.3, rec.rate = 0.05),
                init.net(i.num = 5),
                control)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(comp_plot(mod, at = 3))
})
