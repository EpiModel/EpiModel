##
## R Script File for Journal of Statistical Software Manuscript
## EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks
##

## Section 2. Orientation --------------------------------------------------

## Install the packages
install.packages("EpiModel", dependencies = TRUE)
library("EpiModel")

## Access the help files
help(package = "EpiModel")

## Running the Shiny web applications
if (interactive()) {
  epiweb("dcm")
  epiweb("icm")
  epiweb("net")
}

## Section 4. Base models --------------------------------------------------

## Example 1: Independent SIS Model

set.seed(12345)

## Initialize the network
nw <- network::network.initialize(n = 1000, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "risk", rep(0:1, each = 500))

## ERGM formation formula
formation <- ~ edges + nodefactor("risk") + nodematch("risk") + concurrent

## Target statistics for formula
target.stats <- c(250, 375, 225, 100)

## ERGM dissolution coefficients
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 80)
coef.diss

## Estimate the model
est1 <- netest(nw, formation, target.stats, coef.diss)

## Run diagnostics on the model
dx <- netdx(est1, nsims = 10, nsteps = 1000)
dx

## Plot the diagnostics for the formation formula
par(mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(dx)

## Plot the diagnostics for the dissolution formula
par(mfrow = c(1, 2))
plot(dx, type = "duration")
plot(dx, type = "dissolution")

## Set the initial conditions for the epidemic model
init <- init.net(i.num = 50)

## Set the epidemic parameters
param <- param.net(inf.prob = 0.1, act.rate = 5, rec.rate = 0.02)

## Set the controls
control <- control.net(type = "SIS", nsteps = 500,
                       nsims = 10, epi.by = "risk")

## Simulate the epidemic model
sim1 <- netsim(est1, param, init, control)

## Print the model results
sim1

## Summarize the model results at time step 500
summary(sim1, at = 500)

## Plot the model results (default prevalence plot)
par(mfrow = c(1, 1))
plot(sim1)

## Plot the incidence and recoveries
plot(sim1, y = c("si.flow", "is.flow"), legend = TRUE)

## Plot prevalence stratified by risk group
plot(sim1, y = c("i.num.risk0", "i.num.risk1"), legend = TRUE)

## Plot the static network structure at time steps 1 and 500
par(mfrow = c(1, 2), mar = c(2, 0, 2, 0), mgp = c(1, 1, 0))
plot(sim1, type = "network", at = 1, sims = "mean",
     col.status = TRUE, main = "Prevalence at t1")
plot(sim1, type = "network", at = 500, sims = "mean",
     col.status = TRUE, main = "Prevalence at t500")

## Example 2: Dependent SI Model

set.seed(12345)

## Initial the network
num.m1 <- num.m2 <- 500
nw <- network::network.initialize(num.m1 + num.m2,
                                  bipartite = num.m1, directed = FALSE)

## Enter the sex-specific degree distributions
deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)

## Check the balancing of the degree distributions
check_bip_degdist(num.m1, num.m2, deg.dist.m1, deg.dist.m2)

## Enter the formation model for the ERGM
formation <- ~ edges + b1degree(0:1) + b2degree(0:1)

## Target statistics for the formation model
target.stats <- c(330, 200, 275, 240, 205)

## Calculate coefficients for the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges),
                               duration = 25, d.rate = 0.006)
coef.diss

## Estimate the ERGM
est2 <- netest(nw, formation, target.stats, coef.diss)

## Run diagnostics on the model fit
dx <- netdx(est2, nsims = 10, nsteps = 1000)
dx

## Plot the diagnostics
par(mfrow = c(1, 1))
plot(dx, plots.joined = TRUE)

## Initial conditions for the epidemic model
init <- init.net(i.num = 50, i.num.m2 = 50)

## Epidemic model parameters
param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.1,
                   b.rate = 0.006, b.rate.m2 = NA,
                   ds.rate = 0.005, ds.rate.m2 = 0.006,
                   di.rate = 0.008, di.rate.m2 = 0.009)

## Control settings
control <- control.net(type = "SI", nsims = 10, nsteps = 500,
                       nwstats.formula = ~ edges + meandeg, delete.nodes = TRUE)

## Simulate the epidemic model
sim2 <- netsim(est2, param, init, control)

## Post-simulation diagnostics for ERGM target statistics
plot(sim2, type = "formation", plots.joined = FALSE)
abline(h = 0.66, lty = 2, lwd = 2)

## Plot the epidemic output from the model
par(mfrow = c(1, 2))
plot(sim2, popfrac = TRUE)
plot(sim2, popfrac = FALSE)


## Section 5. Model extension API ------------------------------------------

## Example 1: Age-Dependent Mortality Extension ##

set.seed(12345)

## New aging module
aging <- function(dat, at) {

  ## Attributes
  if (at == 2) {
    n <- sum(dat$attr$active == 1)
    dat$attr$age <- sample(seq(from = 18, to = 69 + 11/12, by = 1/12),
                           n, replace = TRUE)
  } else {
    dat$attr$age <- dat$attr$age + 1/12
  }

  ## Summary statistics
  if (at == 2) {
    dat$epi$meanAge <- c(NA_real_, mean(dat$attr$age, na.rm = TRUE))
  } else {
    dat$epi$meanAge[at] <- mean(dat$attr$age, na.rm = TRUE)
  }

  return(dat)
}

## Express mortality rate as a function of proximity to upper age of 70
ages <- 18:69
death.rates <- 1 / (70 * 12 - ages * 12)
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(ages, death.rates, type = "b", xlab = "age", ylab = "Death Risk")

## Death function
dfunc <- function(dat, at) {

  idsElig <- which(dat$attr$active == 1)
  nElig <- length(idsElig)
  nDeaths <- 0

  if (nElig > 0) {
    ages <- dat$attr$age[idsElig]
    max.age <- dat$param$max.age
    death.rates <- pmin(1, 1 / (max.age * 12 - ages * 12))
    vecDeaths <- which(rbinom(nElig, 1, death.rates) == 1)
    idsDeaths <- idsElig[vecDeaths]
    nDeaths <- length(idsDeaths)

    if (nDeaths > 0) {
      dat$attr$active[idsDeaths] <- 0
      dat$attr$exitTime[idsDeaths] <- at
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

## Arrivals function
bfunc <- function(dat, at) {

  growth.rate <- dat$param$growth.rate
  exptPopSize <- dat$epi$num[1]*(1 + growth.rate*at)
  n <- network.size(dat$nw)

  numNeeded <- exptPopSize - sum(dat$attr$active == 1)
  if (numNeeded > 0) {
    nArrvls <- rpois(1, numNeeded)
  } else {
    nArrvls <- 0
  }
  if (nArrvls > 0) {
    dat$nw <- add.vertices(dat$nw, nv = nArrvls)
    newNodes <- (n + 1):(n + nArrvls)
    dat$nw <- activate.vertices(dat$nw, onset = at,
                                terminus = Inf, v = newNodes)

    dat$attr$active <- c(dat$attr$active, rep(1, nArrvls))
    dat$attr$status <- c(dat$attr$status, rep("s", nArrvls))
    dat$attr$infTime <- c(dat$attr$infTime, rep(NA, nArrvls))
    dat$attr$age <- c(dat$attr$age, rep(18, nArrvls))
  }

  if (at == 2) {
    dat$epi$b.flow <- c(0, nArrvls)
  } else {
    dat$epi$b.flow[at] <- nArrvls
  }

  return(dat)
}

## Initialize the network and estimate the ERGM
nw <- network::network.initialize(500, directed = FALSE)
est3 <- netest(nw, formation = ~ edges, target.stats = 150,
               coef.diss = dissolution_coefs(~ offset(edges), 60,
                                             mean(death.rates)))

## Epidemic model parameterization
param <- param.net(inf.prob = 0.15, growth.rate = 0.01/12, max.age = 70)
init <- init.net(i.num = 50)
control <- control.net(type = NULL, nsims = 5, nsteps = 500,
                       deaths.FUN = dfunc, births.FUN = bfunc,
                       get_prev.FUN = get_prev.net, infection.FUN = infection.net,
                       aging.FUN = aging, depend = TRUE)

## Simulate the epidemic model
sim3 <- netsim(est3, param, init, control)
sim3

## Plot the results of the simulation
plot(sim3, popfrac = TRUE, main = "State Prevalences")
plot(sim3, popfrac = FALSE, main = "State Sizes", sim.lines = TRUE,
     qnts = FALSE, mean.smooth = FALSE)

## Example 2: SEIR Epidemic Model ##

set.seed(12345)

## Modified infection module
infect <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  nw <- dat$nw

  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)

  nElig <- length(idsInf)
  nInf <- 0

  if (nElig > 0 && nElig < nActive) {
    del <- discord_edgelist(dat, at)
    if (!(is.null(del))) {
      del$transProb <- dat$param$inf.prob
      del$actRate <- dat$param$act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)
      if (nInf > 0) {
        dat$attr$status[idsNewInf] <- "e"
        dat$attr$infTime[idsNewInf] <- at
      }
    }
  }

  if (at == 2) {
    dat$epi$se.flow <- c(0, nInf)
  }
  else {
    dat$epi$se.flow[at] <- nInf
  }
  dat$nw <- nw
  return(dat)
}


## New disease progression module
progress <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status

  ei.rate <- dat$param$ei.rate
  ir.rate <- dat$param$ir.rate

  ## E to I progression
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "e")
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }

  ## I to R progression
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)

  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }

  dat$attr$status <- status

  if (at == 2) {
    dat$epi$ei.flow <- c(0, nInf)
    dat$epi$ir.flow <- c(0, nRec)
    dat$epi$e.num <- c(0, sum(active == 1 & status == "e"))
    dat$epi$r.num <- c(0, sum(active == 1 & status == "r"))
  }
  else {
    dat$epi$ei.flow[at] <- nInf
    dat$epi$ir.flow[at] <- nRec
    dat$epi$e.num[at] <- sum(active == 1 & status == "e")
    dat$epi$r.num[at] <- sum(active == 1 & status == "r")
  }

  return(dat)
}

## Initialize the network and estimate the ERGM
nw <- network::network.initialize(500, directed = FALSE)
est4 <- netest(nw, formation = ~ edges, target.stats = 150,
               coef.diss = dissolution_coefs(~ offset(edges), 10))

## Epidemic model parameterization
param <- param.net(inf.prob = 0.5, act.rate = 2, ei.rate = 0.01, ir.rate = 0.005)
init <- init.net(i.num = 10)
control <- control.net(nsteps = 500, nsims = 5, infection.FUN = infect,
                       progress.FUN = progress, recovery.FUN = NULL,
                       skip.check = TRUE, depend = FALSE, verbose.int = 0)

## Simulate the epidemic model
sim4 <- netsim(est4, param, init, control)

## Plot the results
par(mfrow = c(1, 1))
plot(sim4, y = c("s.num", "i.num", "e.num", "r.num"),
     mean.col = 1:4, qnts = 1, qnts.col = 1:4, legend = TRUE)
