
## Demo of EpiModel Overhaul Functionality

suppressMessages(library(EpiModel))


#EpiModel Overhaul - Version 2----

# Example - One Group

# Network with no attribute
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est1 <- netest(nw, formation, target.stats, coef.diss)
dx1 <- netdx(est1, nsims = 5, nsteps = 250)
print(dx1)

# Parameters
init <- init.net(i.num = 50)
param <- param.net(inf.prob = 0.1, act.rate = 5, rec.rate = 0.02)
control <- control.net(type = "SIS", nsteps = 100, nsims = 5)

# Epidemic simulation
sim1 <- netsim(est1, param, init, control)
summary(sim1, at = 50)
par(mfrow = c(1, 1), mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim1)


# Example - Two Group

# Network with 'group' (.g2) attribute
num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodematch("group") + nodefactor("group")
target.stats <- c(150, 125, 62.5)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est2 <- netest(nw, formation, target.stats, coef.diss)
dx2 <- netdx(est2, nsims = 5, nsteps = 250)
print(dx2)

#Parameters
init <- init.net(i.num = 5, i.num.g2 = 5)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   rec.rate = 0.02, rec.rate.g2 = 0.02,
                   act.rate = 5)
control <- control.net(type = "SIS", nsteps = 100, nsims = 5)
control

##Epidemic simulation
sim2 <- netsim(est2, param, init, control)

summary(sim2, at = 50)
plot(sim2)


#Worked Example - User Defined Functions

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

## User-defined departure function
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

## User-defined arrivals function
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
#User must specify all base functions when "type" is NULL
control <- control.net(type = NULL, nsims = 1, nsteps = 100,
                       departures.FUN = dfunc, arrivals.FUN = bfunc,
                       get_prev.FUN = get_prev.net, infection.FUN = infection.net,
                       recovery.FUN = recovery.net, aging.FUN = aging, depend = TRUE)

## Simulate the epidemic model
sim3 <- netsim(est3, param, init, control)

#EpiModel Overhaul - Version 3----

#Closed Population; type = "SI"
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50)
param <- param.net(inf.prob = 0.4, act.rate = 5, rec.rate = 0.02)
control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

sim <- netsim(est, param, init, control)

#Open Population; type = "SIR"
nw <- network.initialize(n = 500, directed = FALSE)
formation <- ~edges + concurrent + degrange(from = 4)
target.stats <- c(175, 110, 0)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)

# Parameters
init <- init.net(i.num = 10, r.num = 15)
param <- param.net(inf.prob = 0.1, act.rate = 5,
                   a.rate = 0.01, rec.rate = 0.02,
                   ds.rate = 0.01, di.rate = 0.01, dr.rate = 0.01)
control <- control.net(type = "SIR", nsims = 1, nsteps = 500, depend = TRUE)
sim <- netsim(est, param, init, control)
