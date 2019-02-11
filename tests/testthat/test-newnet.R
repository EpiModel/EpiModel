context("New Network Models")

test_that("New network models vignette example", {


  ## New Aging Module
  aging <- function(dat, at) {

    if (at == 2) {
      n <- sum(dat$attr$active == 1)
      dat$attr$age <- sample(18:49, n, replace = TRUE)
    } else {
      dat$attr$age <- dat$attr$age + 1/12
    }

    return(dat)
  }


  ## Replacement Departure Module
  dfunc <- function(dat, at) {

    idsElig <- which(dat$attr$active == 1)
    nElig <- length(idsElig)

    nDepartures <- 0

    if (nElig > 0) {
      ages <- dat$attr$age[idsElig]
      life.expt <- dat$param$life.expt
      departure.rates <- pmin(1, 1/(life.expt*12 - ages*12))
      vecDepartures <- which(rbinom(nElig, 1, departure.rates) == 1)
      idsDepartures <- idsElig[vecDepartures]
      nDepartures <- length(idsDepartures)
      if (nDepartures > 0) {
        dat$attr$active[idsDepartures] <- 0
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDepartures, deactivate.edges = TRUE)
      }
    }

    if (at == 2) {
      dat$epi$d.flow <- c(0, nDepartures)
    } else {
      dat$epi$d.flow[at] <- nDepartures
    }

    return(dat)
  }


  ## Replacement Arrival Module
  afunc <- function(dat, at) {

    # Variables ---------------------------------------------------------------
    growth.rate <- dat$param$growth.rate
    exptPopSize <- dat$epi$num[1] * (1 + growth.rate*at)
    n <- network.size(dat$nw)
    tea.status <- dat$control$tea.status

    numNeeded <- exptPopSize - sum(dat$attr$active == 1)
    if (numNeeded > 0) {
      nArrivals <- rpois(1, numNeeded)
    } else {
      nArrivals <- 0
    }
    if (nArrivals > 0) {
      dat$nw <- add.vertices(dat$nw, nv = nArrivals)
      newNodes <- (n + 1):(n + nArrivals)
      dat$nw <- activate.vertices(dat$nw, onset = at,
                                  terminus = Inf, v = newNodes)
    }


    # Update Nodal Attributes -------------------------------------------------
    if (nArrivals > 0) {
      dat$attr$active <- c(dat$attr$active, rep(1, nArrivals))
      if (tea.status == TRUE) {
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = 0, onset = at,
                                            terminus = Inf, v = newNodes)
      }
      dat$attr$status <- c(dat$attr$status, rep("s", nArrivals))
      dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nArrivals))
      dat$attr$age <- c(dat$attr$age, rep(18, nArrivals))
    }


    # Output ------------------------------------------------------------------
    if (at == 2) {
      dat$epi$a.flow <- c(0, nArrivals)
    } else {
      dat$epi$a.flow[at] <- nArrivals
    }

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
                         departures.FUN = dfunc, arrivals.FUN = afunc, aging.FUN = aging,
                         get_prev.FUN = get_prev.net, infection.FUN = infection.net,
                         recovery.FUN = recovery.net, depend = TRUE, save.network = FALSE, 
                         verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")

  ## Test module reordering
  control <- control.net(type = NULL, nsims = 1, nsteps = 10,
                         departures.FUN = dfunc, arrivals.FUN = afunc, aging.FUN = aging,
                         get_prev.FUN = get_prev.net, infection.FUN = infection.net,
                         recovery.FUN = recovery.net, 
                         module.order = c("aging.FUN", "arrivals.FUN", "departures.FUN", 
                                          "get_prev.FUN", "infection.FUN", "recovery.FUN"),
                         depend = TRUE, save.network = FALSE, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")

})
