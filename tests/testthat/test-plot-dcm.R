context("plot.dcm")

# Cheap multi-run model for exercising legend / run / lcomp branches.
build_dcm_multirun <- function(nsteps = 12) {
  param <- param.dcm(inf.prob = 0.2, act.rate = 1:4, rec.rate = 0.05)
  init <- init.dcm(s.num = 100, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = nsteps)
  dcm(param, init, control)
}

build_dcm_singlerun <- function(nsteps = 12) {
  param <- param.dcm(inf.prob = 0.2, act.rate = 1, rec.rate = 0.05)
  init <- init.dcm(s.num = 100, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = nsteps)
  dcm(param, init, control)
}

test_that("plot.dcm validates run and y arguments", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  mod <- build_dcm_multirun()
  expect_error(plot(mod, run = 99), "Specify run between 1 and")
  expect_error(plot(mod, y = "not_a_var"), "Specified y is unavailable")
})

test_that("plot.dcm draws single-run model with default and explicit y", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  mod <- build_dcm_singlerun()
  expect_silent(plot(mod))
  expect_silent(plot(mod, y = "i.num"))
  expect_silent(plot(mod, y = c("s.num", "i.num", "r.num"),
                     col = "Set1", legend = "full"))
  expect_silent(plot(mod, y = "i.num", popfrac = TRUE,
                     grid = TRUE, leg.name = "Infected"))
})

test_that("plot.dcm draws multi-run model with default, single-run, and multi-run subsets", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  mod <- build_dcm_multirun()  # nruns = 4

  # Default y, default run (all runs of i.num): triggers lcomp=1 / nruns>1
  # multi-run loop and the "lim" legend (since nruns >= 3).
  expect_silent(plot(mod))

  # Default y with brewer-ramp palette path (nruns > 5 forces Spectral).
  param_many <- param.dcm(inf.prob = 0.2, act.rate = 1:6, rec.rate = 0.05)
  init <- init.dcm(s.num = 100, i.num = 1, r.num = 0)
  control <- control.dcm(type = "SIR", nsteps = 8)
  mod_many <- dcm(param_many, init, control)
  expect_silent(plot(mod_many))

  # Single run selected: triggers length(run)==1 branch.
  expect_silent(plot(mod, run = 2))

  # Multiple runs selected: triggers length(run)>1 branch + leg.names from
  # the named runs.
  expect_silent(plot(mod, run = c(1, 3)))

  # Multi-y on a multi-run model: must NOT specify multiple runs, but a
  # single-run pick exercises the lcomp>1 / nruns>1 / norun=FALSE branch.
  expect_silent(plot(mod, y = c("s.num", "i.num"), run = 1, legend = "full"))
})

test_that("plot.dcm errors on multiple runs + multiple y, warns on ignored leg.name", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  mod <- build_dcm_multirun()
  expect_error(plot(mod, y = c("s.num", "i.num"), run = c(1, 2)),
               "Plotting multiple runs of multiple y is not supported")

  # Multi-run + multi-y + leg.name triggers the warning at line 401.
  expect_warning(plot(mod, y = c("s.num", "i.num"), leg.name = "anything"),
                 "Legend names ignored")
})

test_that("plot.dcm respects add = TRUE and custom lwd/lty/xlim/ylim", {
  skip_on_cran()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  mod <- build_dcm_singlerun()
  # First plot establishes the device, second overlays.
  plot(mod, y = "i.num")
  expect_silent(plot(mod, y = "s.num", add = TRUE,
                     lwd = 1.5, lty = 2, xlim = c(0, 10),
                     ylim = c(0, 110)))
})
