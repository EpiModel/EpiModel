context("New Network Models")

test_that("New network models vignette example", {


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

    nDeaths <- 0

    if (nElig > 0) {
      ages <- all$attr$age[idsElig]
      life.expt <- all$param$life.expt
      death.rates <- pmin(1, 1/(life.expt*12 - ages*12))
      vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
      idsDeaths <- idsElig[vecDeaths]
      nDeaths <- length(idsDeaths)
      if (nDeaths > 0) {
        all$attr$active[idsDeaths] <- 0
        all$nw <- deactivate.vertices(all$nw,
                                      onset = at,
                                      terminus = Inf,
                                      v = idsDeaths,
                                      deactivate.edges = TRUE)
      }
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

    # Variables ---------------------------------------------------------------
    growth.rate <- all$param$growth.rate
    exptPopSize <- all$out$num[1] * (1 + growth.rate*at)
    n <- network.size(all$nw)
    tea.status <- all$control$tea.status

    numNeeded <- exptPopSize - sum(all$attr$active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    } else {
      nBirths <- 0
    }
    if (nBirths > 0) {
      all$nw <- add.vertices(all$nw,
                             nv = nBirths)
      newNodes <- (n + 1):(n + nBirths)
      all$nw <- activate.vertices(all$nw,
                                  onset = at,
                                  terminus = Inf,
                                  v = newNodes)
    }


    # Update Nodal Attributes -------------------------------------------------
    if (nBirths > 0) {
      all$attr$active <- c(all$attr$active, rep(1, nBirths))
      if (tea.status == TRUE) {
        all$nw <- activate.vertex.attribute(all$nw,
                                            prefix = "testatus",
                                            value = 0,
                                            onset = at,
                                            terminus = Inf,
                                            v = newNodes)
      }
      all$attr$status <- c(all$attr$status, rep("s", nBirths))
      all$attr$infTime <- c(all$attr$infTime, rep(NA, nBirths))
      all$attr$age <- c(all$attr$age, rep(18, nBirths))
    }


    # Output ------------------------------------------------------------------
    if (at == 2) {
      all$out$b.flow <- c(0, nBirths)
    } else {
      all$out$b.flow[at] <- nBirths
    }

    return(all)
  }


  ## Network Model
  nw <- network.initialize(50, directed = FALSE)
  est <- netest(nw,
                formation = ~ edges,
                dissolution = ~offset(edges),
                target.stats = 15,
                coef.diss = dissolution_coefs(~offset(edges), 60, 0.0002747253),
                verbose = FALSE)


  ## EpiModel
  param <- param.net(inf.prob = 0.35,
                     growth.rate = 0.00083,
                     life.expt = 70)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI",
                         nsims = 1,
                         nsteps = 10,
                         deaths.FUN = dfunc,
                         births.FUN = bfunc,
                         aging.FUN = aging,
                         depend = TRUE,
                         save.network = FALSE,
                         verbose = FALSE)
  mod <- netsim(est, param, init, control)
  expect_is(mod, "netsim")

})
