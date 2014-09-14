context("New ICM models")

test_that("ICM new modules", {

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

  ## Replacement Death Module
  dfunc <- function(dat, at) {

    idsElig <- which(dat$attr$active == 1)
    nElig <- length(idsElig)
    if (nElig > 0) {
      ages <- dat$attr$age[idsElig]
      life.expt <- dat$param$life.expt
      death.rates <- pmin(1, 1/(life.expt*12 - ages*12))
      vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
      idsDeaths <- idsElig[vecDeaths]
      nDeaths <- length(idsDeaths)
      if (nDeaths > 0) {
        dat$attr$active[idsDeaths] <- 0
      }
    } else {
      nDeaths <- 0
    }

    if (at == 2) {
      dat$epi$d.flow <- c(0, nDeaths)
    } else {
      dat$epi$d.flow[at] <- nDeaths
    }

    return(dat)
  }


  ## Replacement Birth Module
  bfunc <- function(dat, at) {

    growth.rate <- dat$param$growth.rate
    exptPopSize <- dat$epi$num[1] * (1 + growth.rate*at)

    numNeeded <- exptPopSize - sum(dat$attr$active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    } else {
      nBirths <- 0
    }

    dat$attr$active <- c(dat$attr$active, rep(1, nBirths))
    dat$attr$group <- c(dat$attr$group, rep(1, nBirths))
    dat$attr$status <- c(dat$attr$status, rep("s", nBirths))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
    dat$attr$age <- c(dat$attr$age, rep(18, nBirths))

    if (at == 2) {
      dat$epi$b.flow <- c(0, nBirths)
    } else {
      dat$epi$b.flow[at] <- nBirths
    }

    return(dat)
  }


  param <- param.icm(inf.prob = 0.15,
                     act.rate = 0.15,
                     growth.rate = 0.01/12,
                     life.expt = 70)
  init <- init.icm(s.num = 1000, i.num = 100)
  control <- control.icm(type = "SI", nsims = 2, nsteps = 250,
                         deaths.FUN = dfunc, births.FUN = bfunc,
                         aging.FUN = aging,
                         verbose = FALSE)



  mod <- icm(param, init, control)
  expect_is(mod, "icm")

})




