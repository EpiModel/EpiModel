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


  ## Replacement Death Module
  dfunc <- function(dat, at) {

    idsElig <- which(dat$attr$active == 1)
    nElig <- length(idsElig)

    nDeaths <- 0

    if (nElig > 0) {
      ages <- dat$attr$age[idsElig]
      life.expt <- dat$param$life.expt
      death.rates <- pmin(1, 1/(life.expt*12 - ages*12))
      vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
      idsDeaths <- idsElig[vecDeaths]
      nDeaths <- length(idsDeaths)
      if (nDeaths > 0) {
        dat$attr$active[idsDeaths] <- 0
        dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                      v = idsDeaths, deactivate.edges = TRUE)
      }
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

    # Variables ---------------------------------------------------------------
    growth.rate <- dat$param$growth.rate
    exptPopSize <- dat$epi$num[1] * (1 + growth.rate*at)
    n <- network.size(dat$nw)
    tea.status <- dat$control$tea.status

    numNeeded <- exptPopSize - sum(dat$attr$active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    } else {
      nBirths <- 0
    }
    if (nBirths > 0) {
      dat$nw <- add.vertices(dat$nw, nv = nBirths)
      newNodes <- (n + 1):(n + nBirths)
      dat$nw <- activate.vertices(dat$nw, onset = at,
                                  terminus = Inf, v = newNodes)
    }


    # Update Nodal Attributes -------------------------------------------------
    if (nBirths > 0) {
      dat$attr$active <- c(dat$attr$active, rep(1, nBirths))
      if (tea.status == TRUE) {
        dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus",
                                            value = 0, onset = at,
                                            terminus = Inf, v = newNodes)
      }
      dat$attr$status <- c(dat$attr$status, rep("s", nBirths))
      dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
      dat$attr$age <- c(dat$attr$age, rep(18, nBirths))
    }


    # Output ------------------------------------------------------------------
    if (at == 2) {
      dat$epi$b.flow <- c(0, nBirths)
    } else {
      dat$epi$b.flow[at] <- nBirths
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
  control <- control.net(type = "SI", nsims = 1, nsteps = 10,
                         deaths.FUN = dfunc, births.FUN = bfunc, aging.FUN = aging,
                         depend = TRUE, save.network = FALSE, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")

  ## Test module reordering
  control <- control.net(type = "SI", nsims = 1, nsteps = 10,
                         deaths.FUN = dfunc, births.FUN = bfunc, aging.FUN = aging,
                         module.order = c("aging.FUN", "births.FUN", "deaths.FUN"),
                         depend = TRUE, save.network = FALSE, verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")

})
