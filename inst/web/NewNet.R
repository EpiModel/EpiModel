## ----setup, include=FALSE------------------------------------------------
library(EpiModel)

## ----controlArgs, echo = FALSE-------------------------------------------
args(control.net)

## ----Ex1agingFunc--------------------------------------------------------
aging <- function(dat, at) {
  
  ## Attributes
  if (at == 2) {
    n <- sum(dat$attr$active == 1)
    dat$attr$age <- sample(18:49, n, replace = TRUE)
  } else {
    dat$attr$age <- dat$attr$age + 1/12
  }
  
  ## Summary statistics
  if (at == 2) {
    dat$epi$meanAge <- rep(mean(dat$attr$age, na.rm = TRUE), 2)
  } else {
    dat$epi$meanAge[at] <- mean(dat$attr$age, na.rm = TRUE)
  }
  
  return(dat)
}

## ----Ex1deathEx----------------------------------------------------------
ages <- 18:49
death.rates <- 1/(70*12 - ages*12)
par(mar = c(3.2, 3.2, 1, 1), mgp = c(2, 1, 0))
plot(ages, death.rates, pch = 20, xlab = "age", ylab = "Death Risk")

## ----deathFunc-----------------------------------------------------------
dfunc <- function(dat, at) {
  
  # Parameters
  idsElig <- which(dat$attr$active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0
  
  # Processes
  if (nElig > 0) {
    ages <- dat$attr$age[idsElig]
    life.expt <- dat$param$life.expt
    death.rates <- pmin(1, 1/(life.expt*12 - ages*12))
    vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)
    
    # Update nodal attributes on attr and networkDynamic object
    if (nDeaths > 0) {
      dat$attr$active[idsDeaths] <- 0
      dat$attr$exitTime[idsDeaths] <- at
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf, 
                                    v = idsDeaths, deactivate.edges = TRUE)
    }
  }
  
  # Summary statistics
  if (at == 2) {
    dat$epi$d.flow <- c(0, nDeaths)
  } else {
    dat$epi$d.flow[at] <- nDeaths
  }
  
  return(dat)
}

## ----Ex1birthFunc--------------------------------------------------------
bfunc <- function(dat, at) {
  
  # Variables
  growth.rate <- dat$param$growth.rate
  exptPopSize <- dat$epi$num[1] * (1 + growth.rate * at)
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
    dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf, v = newNodes)
  }
  
  # Update attributes
  if (nBirths > 0) {
    dat$attr$active <- c(dat$attr$active, rep(1, nBirths))
    dat$attr$status <- c(dat$attr$status, rep("s", nBirths))  
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nBirths))
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nBirths))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nBirths))
    dat$attr$age <- c(dat$attr$age, rep(18, nBirths))
    if (tea.status == TRUE) {
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus", 
                                          value = 0, onset = at, 
                                          terminus = Inf, v = newNodes)
    }
  }

  # Summary statistics
  if (at == 2) {
    dat$epi$b.flow <- c(0, nBirths)
  } else {
    dat$epi$b.flow[at] <- nBirths
  }
  
  return(dat)
}

## ----Ex1nwParam----------------------------------------------------------
nw <- network.initialize(500, directed = FALSE)
est <- netest(nw, formation = ~edges, target.stats = 150, 
              coef.diss = dissolution_coefs(~offset(edges), 60, mean(death.rates)))

## ----Ex1epiParam---------------------------------------------------------
param <- param.net(inf.prob = 0.15, growth.rate = 0.00083, life.expt = 70)
init <- init.net(i.num = 50)

## ----Ex1epiControl-------------------------------------------------------
control <- control.net(type = "SI", nsims = 4, ncores = 4, nsteps = 250, 
                       departures.FUN = dfunc, arrivals.FUN = bfunc, aging.FUN = aging, 
                       depend = TRUE, save.network = FALSE)

## ----Ex1runNetsim, cache = TRUE, results = "hide"------------------------
mod <- netsim(est, param, init, control)

## ----Ex1printMod---------------------------------------------------------
mod

## ----Ex1plotMod1---------------------------------------------------------
par(mfrow = c(1,2))
plot(mod, main = "State Prevalences")
plot(mod, popfrac = FALSE, main = "State Sizes", sim.lines = TRUE, 
     qnts = FALSE, mean.smooth = FALSE)

## ----Ex1plotMod2---------------------------------------------------------
par(mfrow = c(1, 2))
plot(mod, y = "num", popfrac = FALSE, main = "Population Size", ylim = c(0, 1000))
plot(mod, y = "meanAge", main = "Mean Age", ylim = c(18, 70))

## ----Ex1plotMod3---------------------------------------------------------
par(mfrow = c(1, 2))
plot(mod, y = "d.flow", popfrac = FALSE, mean.smooth = TRUE, qnts = 1, main = "Deaths")
plot(mod, y = "b.flow", popfrac = FALSE, mean.smooth = TRUE, qnts = 1, main = "Births")

## ----Ex2omigFunc---------------------------------------------------------
outmigrate <- function(dat, at) {
  
  # Variables
  active <- dat$attr$active
  status <- dat$attr$status
  rates <- dat$param$omig.rates
  
  # Process
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  
  nMig <- 0
  if (nElig > 0) {
    ratesElig <- rates[as.numeric(status == "i") + 1]
    vecMig <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecMig) > 0) {
      idsMig <- idsElig[vecMig]
      nMig <- length(idsMig)
      
      dat$attr$active[idsMig] <- 0
      dat$attr$exitTime[idsMig] <- at
      dat$nw <- deactivate.vertices(dat$nw, onset = at, terminus = Inf, 
                                    v = idsMig, deactivate.edges = TRUE)
    }
  }
  
  # Summary statistics
  if (at == 2) {
    dat$epi$omig.flow <- c(0, nMig)
  } else {
    dat$epi$omig.flow[at] <- nMig
  }
  
  return(dat)
}

## ----Ex2imigFunc---------------------------------------------------------
inmigrate <- function(dat, at) {
  
  # Variables
  nw <- dat$nw
  n <- network.size(nw)
  active <- dat$attr$active
  
  exp.inmig <- dat$param$exp.inmig
  risk <- dat$param$imig.risk
  tea.status <- dat$control$tea.status
  
  # Add Nodes
  nMig <- 0
  idsElig <- which(active == 1)
  nElig <- length(idsElig)
  if (nElig > 0) {
    nMig <- rpois(1, exp.inmig)
    if (nMig > 0) {
      dat$nw <- add.vertices(dat$nw, nv = nMig)
      newNodes <- (n + 1):(n + nMig)
      dat$nw <- activate.vertices(dat$nw, onset = at, terminus = Inf,
                                  v = newNodes)
    }
  }
  
  # Update attributes
  if (nMig > 0) {
    dat$attr$active <- c(dat$attr$active, rep(1, nMig))
    newStatus <- rbinom(nMig, 1, risk)
    newStatus <- ifelse(newStatus == 1, "i", "s")
    if (tea.status == TRUE) {
      dat$nw <- activate.vertex.attribute(dat$nw, prefix = "testatus", value = newStatus, 
                                          onset = at, terminus = Inf, v = newNodes)
    }
    dat$attr$status <- c(dat$attr$status, newStatus)
    infTime <- ifelse(newStatus == "i", at, NA)    
    dat$attr$infTime <- c(dat$attr$infTime, infTime)
    dat$attr$entrTime <- c(dat$attr$entrTime, rep(at, nMig))
    dat$attr$exitTime <- c(dat$attr$exitTime, rep(NA, nMig))
  }
  
  # Summary statistics
  if (at == 2) {
    dat$epi$imig.flow <- c(0, nMig)
  } else {
    dat$epi$imig.flow[at] <- nMig
  }
  
  return(dat)
}

## ----Ex2nwParam----------------------------------------------------------
nw <- network.initialize(500, directed = FALSE)
est <- netest(nw, formation = ~edges, target.stats = 150, 
              coef.diss = dissolution_coefs(~offset(edges), 60, 0.0046))

## ----Ex2epiParam---------------------------------------------------------
param <- param.net(inf.prob = 0.2, rec.rate = 0.02, 
                   omig.rates = c(0.005, 0.001), exp.inmig = 3, imig.risk = 0.1)
init <- init.net(i.num = 50)

## ----Ex2epiControl-------------------------------------------------------
control <- control.net(type = "SIS", nsims = 4, ncores = 4, nsteps = 500, 
                       omig.FUN = outmigrate, imig.FUN = inmigrate, 
                       depend = TRUE, save.network = FALSE)

## ----Ex2runNetsim, cache = TRUE, results = "hide"------------------------
mod <- netsim(est, param, init, control)

## ----Ex2print------------------------------------------------------------
mod
plot(mod)

## ----Ex2plot-------------------------------------------------------------
par(mfrow = c(1, 2))
plot(mod, y = "omig.flow", main = "Out Migration")
plot(mod, y = "imig.flow", main = "In Migration")

