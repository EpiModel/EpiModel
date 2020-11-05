context("New ICM models")

test_that("ICM new modules", {

  ## New Aging Module
  aging <- function(dat, at) {

    if (at == 2) {
      n <- sum(dat$attr$active == 1)
      dat$attr$age <- sample(18:49, n, replace = TRUE)
    } else {
      dat$attr$age <- dat$attr$age + 1 / 12
    }

    return(dat)
  }

  ## Replacement Departure Module
  dfunc <- function(dat, at) {

    idsElig <- which(dat$attr$active == 1)
    nElig <- length(idsElig)
    if (nElig > 0) {
      ages <- dat$attr$age[idsElig]
      life.expt <- dat$param$life.expt
      departure.rates <- pmin(1, 1 / (life.expt * 12 - ages * 12))
      vecDepartures <- which(rbinom(nElig, 1, departure.rates) == 1)
      idsDepartures <- idsElig[vecDepartures]
      nDepartures <- length(idsDepartures)
      if (nDepartures > 0) {
        dat$attr$active[idsDepartures] <- 0
      }
    } else {
      nDepartures <- 0
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

    growth.rate <- dat$param$growth.rate
    exptPopSize <- dat$epi$num[1] * (1 + growth.rate * at)

    numNeeded <- exptPopSize - sum(dat$attr$active == 1)
    if (numNeeded > 0) {
      nArrivals <- rpois(1, numNeeded)
    } else {
      nArrivals <- 0
    }

    dat$attr$active <- c(dat$attr$active, rep(1, nArrivals))
    dat$attr$group <- c(dat$attr$group, rep(1, nArrivals))
    dat$attr$status <- c(dat$attr$status, rep("s", nArrivals))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nArrivals))
    dat$attr$age <- c(dat$attr$age, rep(18, nArrivals))

    if (at == 2) {
      dat$epi$a.flow <- c(0, nArrivals)
    } else {
      dat$epi$a.flow[at] <- nArrivals
    }

    return(dat)
  }


  param <- param.icm(inf.prob = 0.15,
                     act.rate = 0.15,
                     growth.rate = 0.01 / 12,
                     life.expt = 70)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsims = 2, nsteps = 250,
                         departures.FUN = dfunc, arrivals.FUN = afunc,
                         aging.FUN = aging,
                         verbose = FALSE)



  mod <- icm(param, init, control)
  expect_is(mod, "icm")

})