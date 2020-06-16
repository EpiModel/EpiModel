
##
## Tutorial: Time-Varying Behavior and Biology
## Day 3 | Network Modeling for Epidemics
##

## Load EpiModel
library(EpiModel)


# 1. Generic Stage-Specific Behavior/Biology ------------------------------

## Infection probabilities per time step of infection
probs <- c(0.5, 0.05)
durs <- c(10, 1)
inf.probs <- rep(probs, durs)

## Variability in act rates over infection
act.rates <- c(rpois(10, lambda = 3), 3)

## Plot the time series
par(mfrow = c(1, 2))
plot(inf.probs, type = "S", lwd = 2, ylim = 0:1)
plot(act.rates, type = "S", lwd = 2)
abline(h = 3, lty = 2)

## Fit a basic ERGM, run a basic SIS
nw <- network.initialize(n = 100, directed = FALSE)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
param <- param.net(inf.prob = inf.probs, act.rate = act.rates, rec.rate = 0.01)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 100, nsims = 1)
sim <- netsim(est, param, init, control)

## Get the transmission matrix
tm <- get_transmat(sim)
head(tm, 10)

## Calculate the percent of infections in "acute" stage infection
mean(tm$infDur <= 10)


# 2. Hollingsworth Model for Stage-Specific HIV Transmission --------------

## Stage-specific transmission *rates* per month
probs <- c(0.2055, 0.0088, 0.0614, 0)
durs <- c(3, 100, 9, 10)
inf.probs <- rep(probs, durs)
inf.probs

## Way of parameterizing this
probs <- c(0.2055, 0.0088, 0.0614, 0.1)
acts <- c(1, 1, 1, 0)
durs <- c(3, 100, 9, 10)
inf.probs <- rep(probs, durs)
act.rates <- rep(acts, durs)

## Plot the time series
par(mfrow = c(1,1))
plot(inf.probs, type = "S", ylim = c(0, 1), lwd = 2)
lines(act.rates, type = "S", col = 2, lwd = 2)
legend(1, 0.8, legend = c("inf.probs", "act.rates"),
       lwd = 3, col = 1:2, bty = "n")
