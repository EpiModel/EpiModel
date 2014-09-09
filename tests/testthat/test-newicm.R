context("ICM new models")

test_that("ICM new modules", {

  ## New Aging Module
  aging <- function(all, at) {

    if (at == 2) {
      n <- sum(all$attr$active == 1)
      all$attr$age <- sample(18:49, n, replace = TRUE)
    } else {
      all$attr$age <- all$attr$age + 1/12
    }

    return(all)
  }

  ## Replacement Death Module
  dfunc <- function(all, at) {

    idsElig <- which(all$attr$active == 1)
    nElig <- length(idsElig)
    if (nElig > 0) {
      ages <- all$attr$age[idsElig]
      life.expt <- all$param$life.expt
      death.rates <- pmin(1, 1/(life.expt*12 - ages*12))
      vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
      idsDeaths <- idsElig[vecDeaths]
      nDeaths <- length(idsDeaths)
      if (nDeaths > 0) {
        all$attr$active[idsDeaths] <- 0
      }
    } else {
      nDeaths <- 0
    }

    if (at == 2) {
      all$out$d.flow <- c(0, nDeaths)
    } else {
      all$out$d.flow[at] <- nDeaths
    }

    return(all)
  }


  ## Replacement Birth Module
  bfunc <- function(all, at) {

    growth.rate <- all$param$growth.rate
    exptPopSize <- all$out$num[1] * (1 + growth.rate*at)

    numNeeded <- exptPopSize - sum(all$attr$active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    } else {
      nBirths <- 0
    }

    all$attr$active <- c(all$attr$active, rep(1, nBirths))
    all$attr$group <- c(all$attr$group, rep(1, nBirths))
    all$attr$status <- c(all$attr$status, rep("s", nBirths))
    all$attr$infTime <- c(all$attr$infTime, rep(NA, nBirths))
    all$attr$age <- c(all$attr$age, rep(18, nBirths))

    if (at == 2) {
      all$out$b.flow <- c(0, nBirths)
    } else {
      all$out$b.flow[at] <- nBirths
    }

    return(all)
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




