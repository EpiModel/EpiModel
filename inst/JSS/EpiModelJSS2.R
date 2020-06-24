##
## R Script File for Journal of Statistical Software Manuscript
## EpiModel: An R Package for Mathematical Modeling of Infectious Disease over Networks
##

## Section 2. Orientation --------------------------------------------------

## Install the packages
#install.packages("EpiModel", dependencies = TRUE)
#library("EpiModel")

## Access the help files
#help(package = "EpiModel")

## Running the Shiny web applications
#if (interactive()) {
#  epiweb("dcm")
#  epiweb("icm")
#  epiweb("net")
#}

## Section 4. Base models --------------------------------------------------

## Example 1: Independent SIS Model

set.seed(12345)

## Initialize the network
nw <- network_initialize(n = 1000)
nw <- set_vertex_attribute(nw, "risk", rep(0:1, each = 500))

## ERGM formation formula
formation <- ~edges + nodefactor("risk") + nodematch("risk") + concurrent

## Target statistics for formula
target.stats <- c(250, 375, 225, 100)

## ERGM dissolution coefficients
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 80)
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
nw <- network_initialize(n = num.m1 + num.m2)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 500))
## Enter the sex-specific degree distributions
deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)

## Check the balancing of the degree distributions
check_bip_degdist(num.m1, num.m2, deg.dist.m1, deg.dist.m2)

## Enter the formation model for the ERGM
formation <- ~ edges + degree(0:1, by = "group") + nodematch("group")

## Target statistics for the formation model
target.stats <- c(330, 200, 275, 240, 205, 0)

## Calculate coefficients for the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
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
init <- init.net(i.num = 50, i.num.g2 = 50)

## Epidemic model parameters
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.1,
                   a.rate = 0.006, a.rate.g2 = NA,
                   ds.rate = 0.005, ds.rate.g2 = 0.006,
                   di.rate = 0.008, di.rate.g2 = 0.009)

## Control settings
control <- control.net(type = "SI", nsims = 10, nsteps = 500,
                       nwstats.formula = ~edges + meandeg)

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

## New Aging Module
aging <- function(dat, at) {

  if (at == 2) {
    active <- get_attr(dat, "active")
    n <- sum(active == 1)
    age <- sample(18:69, n, replace = TRUE)
    dat <- set_attr(dat, "age", age)
  } else {
    age <- get_attr(dat, "age") + 1/12
    dat <- set_attr(dat, "age", age)
  }

  return(dat)
}

## Express departure rate as a function of proximity to upper age of 70
ages <- 18:69
departure.rates <- 1 / (70 * 12 - ages * 12)
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(ages, departure.rates, type = "b", xlab = "age", ylab = "Departure Risk")

## Departure function
dfunc <- function(dat, at) {
  active <- get_attr(dat, "active")
  exitTime <- get_attr(dat, "exitTime")
  idsElig <- which(active == 1)
  nElig <- length(idsElig)

  nDepartures <- 0

  if (nElig > 0) {
    ages <- get_attr(dat, "age")[idsElig]
    max.age <- get_param(dat, "max.age")
    departure.rates <- pmin(1, 1/(max.age*12 - ages*12))
    vecDepartures <- which(rbinom(nElig, 1, departure.rates) == 1)
    idsDepartures <- idsElig[vecDepartures]
    nDepartures <- length(idsDepartures)
    if (nDepartures > 0) {
      active[idsDepartures] <- 0
      dat <- set_attr(dat, "active", active)
      exitTime[idsDepartures] <- at
      dat <- set_attr(dat, "exitTime",exitTime)
      if (get_control(dat, "tergmLite") == FALSE) {
        dat$nw[[1]] <- deactivate.vertices(dat$nw[[1]], onset = at, terminus = Inf,
                                           v = idsDepartures, deactivate.edges = TRUE)
      } else {
        dat <- delete_attr(dat, idsDepartures)
        dat$el[[1]] <- delete_vertices(dat$el[[1]], idsDepartures)
      }
    }
  }
  # Output ----------------------------------
  dat <- set_epi(dat, "d.flow", at, nDepartures)
  return(dat)
}

## Arrivals function
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

  if (nArrivals > 0) {
    newNodes <- (n + 1):(n + nArrivals)
    if (get_control(dat, "tergmLite") == FALSE) {
      dat$nw[[1]] <- add.vertices(dat$nw[[1]], nv = nArrivals)
      dat$nw[[1]] <- activate.vertices(dat$nw[[1]], onset = at,
                                       terminus = Inf, v = newNodes)
    } else {
      dat$el[[1]] <- add_vertices(dat$el[[1]], nv = sum(nArrivals))
    }

    dat <- append_attr(dat, "active", 1, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "entrTime", at, nArrivals)
    dat <- append_attr(dat, "exitTime", NA, nArrivals)
    dat <- append_attr(dat, "age", 18, nArrivals)
  }

  # Output ------------------------------------------------------------------
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}

## Initialize the network and estimate the ERGM
nw <- network_initialize(500)
est3 <- netest(nw, formation = ~ edges, target.stats = 150,
               coef.diss = dissolution_coefs(~ offset(edges), 60,
                                             mean(departure.rates)))

## Epidemic model parameterization
param <- param.net(inf.prob = 0.15, growth.rate = 0.01/12, max.age = 70)
init <- init.net(i.num = 50)
control <- control.net(type = NULL, nsims = 5, nsteps = 500,
                       departures.FUN = dfunc, arrivals.FUN = afunc,
                       prevalence.FUN = prevalence.net, infection.FUN = infection.net,
                       resim_nets.FUN = resim_nets, aging.FUN = aging)

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

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")

  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)

  nElig <- length(idsInf)
  nInf <- 0

  if (nElig > 0 && nElig < nActive) {
    del <- discord_edgelist(dat, at)
    if (!(is.null(del))) {
      del$transProb <- inf.prob
      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
      }
    }
  }

  # Output ---------------------------------
  dat <- set_epi(dat, "se.flow", at, nInf)

  return(dat)
}


## New disease progression module
progress <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")

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

  dat <- set_attr(dat, "status", status)

  dat <- set_epi(dat, "ei.flow", at, nInf)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "e.num", at, sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "r.num", at, sum(active == 1 & status == "r"))

  return(dat)
}

## Initialize the network and estimate the ERGM
nw <- network_initialize(500, directed = FALSE)
est4 <- netest(nw, formation = ~ edges, target.stats = 150,
               coef.diss = dissolution_coefs(~ offset(edges), 10))

## Epidemic model parameterization
param <- param.net(inf.prob = 0.5, act.rate = 2, ei.rate = 0.01, ir.rate = 0.005)
init <- init.net(i.num = 10)
control <- control.net(type = NULL, nsteps = 500, nsims = 5,
                       skip.check = TRUE, verbose.int = 0,
                       infection.FUN = infect, progress.FUN = progress,
                       resim_net.FUN = resim_nets, prevalence.FUN = prevalence.net,
                       resimulate.network = FALSE)

## Simulate the epidemic model
sim4 <- netsim(est4, param, init, control)

## Plot the results
par(mfrow = c(1, 1))
plot(sim4, y = c("s.num", "i.num", "e.num", "r.num"),
     mean.col = 1:4, qnts = 1, qnts.col = 1:4, legend = TRUE)
