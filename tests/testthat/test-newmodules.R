context("New Network Models")

test_that("New network models vignette example", {

  ## New Aging Module
  aging <- function(dat, at) {

    age <- get_attr(dat, "age", override.null.error = TRUE)
    if (is.null(age)) {
      active <- get_attr(dat, "active")
      n <- sum(active == 1)
      age <- sample(18:49, n, replace = TRUE)
    } else {
      age <- get_attr(dat, "age") + 1 / 12
    }
    dat <- set_attr(dat, "age", age)

    return(dat)
  }


  ## Replacement Departure Module
  dfunc <- function(dat, at) {
    active <- get_attr(dat, "active")
    exitTime <- get_attr(dat, "exitTime")
    idsElig <- which(active == 1)
    nElig <- length(idsElig)

    nDepartures <- 0

    if (nElig > 0) {
      ages <- get_attr(dat, "age")[idsElig]
      life.expt <- get_param(dat, "life.expt")
      departure.rates <- pmin(1, 1 / (life.expt * 12 - ages * 12))
      vecDepartures <- which(rbinom(nElig, 1, departure.rates) == 1)
      idsDepartures <- idsElig[vecDepartures]
      nDepartures <- length(idsDepartures)
      if (nDepartures > 0) {
        active[idsDepartures] <- 0
        exitTime[idsDepartures] <- at
        dat <- set_attr(dat, "active", active)
        dat <- set_attr(dat, "exitTime", exitTime)
      }
    }

    # Output
    dat <- set_epi(dat, "d.flow", at, nDepartures)
    return(dat)
  }


  ## Replacement Arrival Module
  afunc <- function(dat, at) {

    # Variables
    growth.rate <- get_param(dat, "growth.rate")
    exptPopSize <- get_epi(dat, "num", 1) * (1 + growth.rate * at)
    active <- get_attr(dat, "active")
    numNeeded <- exptPopSize - sum(active == 1)

    if (numNeeded > 0) {
      nArrivals <- rpois(1, numNeeded)
    } else {
      nArrivals <- 0
    }

    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", 0, nArrivals)

    # Output
    dat <- set_epi(dat, "a.flow", at, nArrivals)

    return(dat)
  }


  ## Network Model
  nw <- network.initialize(50, directed = FALSE)
  est <- netest(nw, formation = ~edges, target.stats = 15,
                coef.diss = dissolution_coefs(~offset(edges), 60, 0.000274),
                verbose = FALSE)


  ## EpiModel
  param <- param.net(inf.prob = 0.35, growth.rate = 0.00083, life.expt = 70)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = 5,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         tergmLite = FALSE, resimulate.network = TRUE, verbose = FALSE)
  mod1 <- netsim(est, param, init, control)
  capture_output(
    mod1
  )

  expect_is(mod1, "netsim")
  expect_output(print(mod1), "resim_nets.FUN")
  expect_output(print(mod1), "infection.FUN")
  expect_output(print(mod1), "departures.FUN")
  expect_output(print(mod1), "arrivals.FUN")
  expect_output(print(mod1), "aging.FUN")

  ## Test module reordering
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         module.order = c("resim_nets.FUN", "infection.FUN",
                                          "aging.FUN", "arrivals.FUN",
                                          "departures.FUN", "prevalence.FUN"),
                         tergmLite = FALSE, resimulate.network = TRUE, verbose = FALSE)
  mod2 <- netsim(est, param, init, control)
  expect_is(mod2, "netsim")

  ### tergmLite replication
  param <- param.net(inf.prob = 0.35, growth.rate = 0.00083, life.expt = 70)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         infection.FUN = infection.net,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         tergmLite = TRUE, verbose = FALSE,
                         resimulate.network = TRUE)
  mod3 <- netsim(est, param, init, control)
  expect_is(mod3, "netsim")

  ## Test module reordering
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         module.order = c("resim_nets.FUN", "infection.FUN",
                                          "aging.FUN", "arrivals.FUN",
                                          "departures.FUN", "prevalence.FUN"),
                         tergmLite = TRUE, resimulate.network = TRUE, verbose = FALSE)
  mod4 <- netsim(est, param, init, control)
  expect_is(mod4, "netsim")

  ## "updated" infection module
  infect <- infection.net
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infect,
                         module.order = c("resim_nets.FUN", "infection.FUN",
                                          "aging.FUN", "arrivals.FUN",
                                          "departures.FUN", "prevalence.FUN"),
                         tergmLite = TRUE, resimulate.network = TRUE, verbose = FALSE)
  mod5 <- netsim(est, param, init, control)
  expect_is(mod5, "netsim")

  expect_output(print(mod5), "resim_nets.FUN")
  expect_output(print(mod5), "infection.FUN")
  expect_output(print(mod5), "departures.FUN")
  expect_output(print(mod5), "arrivals.FUN")
  expect_output(print(mod5), "aging.FUN")

})
