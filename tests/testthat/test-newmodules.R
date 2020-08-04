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
      age <- get_attr(dat, "age") + 1/12
    }
    dat <- set_attr(dat, "age", age)

    return(dat)
  }


  ## Replacement Departure Module
  dfunc <- function(dat, at) {
    active <- get_attr(dat, "active")
    idsElig <- which(active == 1)
    nElig <- length(idsElig)

    nDepartures <- 0

    if (nElig > 0) {
      ages <- get_attr(dat, "age")[idsElig]
      life.expt <- get_param(dat, "life.expt")
      departure.rates <- pmin(1, 1/(life.expt*12 - ages*12))
      vecDepartures <- which(rbinom(nElig, 1, departure.rates) == 1)
      idsDepartures <- idsElig[vecDepartures]
      nDepartures <- length(idsDepartures)
      if (nDepartures > 0) {
        active[idsDepartures] <- 0
        dat <- set_attr(dat, "active", active)
      }
    }

    # Output ----------------------------------
    dat <- set_epi(dat, "d.flow", at, nDepartures)
    return(dat)
  }


  ## Replacement Arrival Module
  afunc <- function(dat, at) {

    # Variables ---------------------------------------------------------------
    growth.rate <- get_param(dat, "growth.rate")
    exptPopSize <- get_epi(dat, "num", 1)*(1 + growth.rate*at)
    n <- sum(get_attr(dat, "active") == 1)
        active <- get_attr(dat, "active")
    numNeeded <- exptPopSize - sum(active == 1)

    if (numNeeded > 0) {
      nArrivals <- rpois(1, numNeeded)
    } else {
      nArrivals <- 0
    }

    # Output ------------------------------------------------------------------
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
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         tergmLite = FALSE, resimulate.network = TRUE)
  mod1 <- netsim(est, param, init, control)
  expect_is(mod1, "netsim")

  ## Test module reordering
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         module.order = c("resim_nets.FUN", "infection.FUN",
                                          "aging.FUN", "arrivals.FUN", "departures.FUN",
                                          "prevalence.FUN"),
                         tergmLite = FALSE, resimulate.network = TRUE)
  undebug(EpiModel:::netsim_loop)
  mod2 <- netsim(est, param, init, control)
  expect_is(mod2, "netsim")

  ### tergmLite replication
  param <- param.net(inf.prob = 0.35, growth.rate = 0.00083, life.expt = 70)
  init <- init.net(i.num = 10)
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         infection.FUN = infection.net,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         tergmLite = TRUE, verbose = FALSE, resimulate.network = TRUE)
  mod3 <- netsim(est, param, init, control)
  expect_is(mod3, "netsim")

  ## Test module reordering
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc,
                         arrivals.FUN = afunc, aging.FUN = aging,
                         infection.FUN = infection.net,
                         module.order = c("resim_nets.FUN", "infection.FUN",
                                          "aging.FUN", "arrivals.FUN", "departures.FUN",
                                          "prevalence.FUN"),
                         tergmLite = TRUE, resimulate.network = TRUE)
  mod4 <- netsim(est, param, init, control)
  expect_is(mod4, "netsim")

})
