## EpiModel 2.0 Worked Example - SEIR Module Extension:
## Description: An closed population of size 1000 in which an epidemic
## of type SIER is introduced. The new exposed state is modeled by defining a new
## infection.FUN module infect2; progression from an exposed state to infectivity
## is handled by the user created function progress.net.
## -----------------------------------------------------------------------------

## Modified infection module using new accessor functions
infect2 <- function(dat, at) {

  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  act.rate <- get_param(dat, "act.rate")
  inf.prob <- get_param(dat, "inf.prob")

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

  dat <- set_epi(dat, "se.flow", at, nInf)

  return(dat)
}

## New disease progression module using new accessor functions.
progress2 <- function(dat, at) {

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

nw <- network_initialize(n = 500)
est <- netest(nw, formation = ~edges, target.stats = 150,
              coef.diss = dissolution_coefs(~offset(edges), 10))

## Epidemic model parameterization
param <- param.net(inf.prob = 0.5, act.rate = 2, ei.rate = 0.01, ir.rate = 0.005)
init <- init.net(i.num = 10)
## If any user defined functions, type must equal NULL
control <- control.net(type = NULL, nsteps = 100, nsims = 1, infection.FUN = infect2,
                       progress.FUN = progress2)

## Simulate the epidemic model
sim <- netsim(est, param, init, control)

## Print the results
print(sim)

## Plot the results
par(mfrow = c(1, 1))
plot(sim, y = c("s.num", "i.num", "e.num", "r.num"),
     mean.col = 1:4, qnts = 1, qnts.col = 1:4, legend = TRUE)
