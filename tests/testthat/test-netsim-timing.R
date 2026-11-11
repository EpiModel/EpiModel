context("Module timing in netsim")

timing_test_sim <- function(timing) {
  set.seed(1)
  nw <- network_initialize(n = 100)
  est <- netest(nw, formation = ~edges, target.stats = 50,
                coef.diss = dissolution_coefs(~offset(edges), duration = 20),
                verbose = FALSE)
  param <- param.net(inf.prob = 0.3, act.rate = 1, rec.rate = 0.05)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SIS", nsteps = 10, nsims = 1, ncores = 1,
                         tergmLite = TRUE, resimulate.network = TRUE,
                         verbose = FALSE, .timing = timing)
  set.seed(2)
  netsim(est, param, init, control)
}

test_that("`.timing = TRUE` records one timing output per module", {
  skip_on_cran()
  sim <- timing_test_sim(timing = TRUE)
  epi <- as.data.frame(sim)
  module_names <- c("resim_nets.FUN", "summary_nets.FUN", "infection.FUN",
                    "recovery.FUN", "departures.FUN", "arrivals.FUN",
                    "nwupdate.FUN", "prevalence.FUN")
  timing_names <- c(paste0("timing.", module_names), "timing.step")
  expect_setequal(grep("^timing\\.", names(epi), value = TRUE), timing_names)

  # No module runs at the first time step, which is the initialization
  simulated <- epi[epi$time > 1, timing_names]
  expect_true(all(is.na(epi[epi$time == 1, timing_names])))
  expect_false(anyNA(simulated))
  expect_true(all(simulated >= 0))
  module_total <- rowSums(simulated[, paste0("timing.", module_names)])
  expect_true(all(simulated$timing.step >= module_total))
})

test_that("timing leaves every other output unchanged", {
  skip_on_cran()
  timed <- as.data.frame(timing_test_sim(timing = TRUE))
  untimed <- as.data.frame(timing_test_sim(timing = FALSE))
  expect_length(grep("^timing\\.", names(untimed)), 0)
  expect_identical(timed[, names(untimed)], untimed)
})

test_that("a module passed through `module.order` is timed under its name", {
  skip_on_cran()
  set.seed(1)
  nw <- network_initialize(n = 50)
  est <- netest(nw, formation = ~edges, target.stats = 25,
                coef.diss = dissolution_coefs(~offset(edges), duration = 20),
                verbose = FALSE)
  count_steps <- function(dat, at) {
    set_epi(dat, "steps.run", at - 1)
  }
  control <- control.net(
    type = NULL, nsteps = 5, nsims = 1, ncores = 1, verbose = FALSE,
    infection.FUN = infection.net, counter.FUN = count_steps,
    module.order = c("resim_nets.FUN", "summary_nets.FUN", "infection.FUN",
                     "counter.FUN", "nwupdate.FUN", "prevalence.FUN"),
    resimulate.network = TRUE, .timing = TRUE
  )
  sim <- netsim(est, param.net(inf.prob = 0.3, act.rate = 1),
                init.net(i.num = 5), control)
  epi <- as.data.frame(sim)
  expect_true("timing.counter.FUN" %in% names(epi))
  expect_false(anyNA(epi$timing.counter.FUN[epi$time > 1]))
})
