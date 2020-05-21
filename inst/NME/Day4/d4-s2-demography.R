
##
## Tutorial: Demography and Attributes
## Day 4 | Network Modeling for Epidemics
##

# Load EpiModel
library(EpiModel)

# Initialize the network
nw <- network.initialize(n = 500, directed = FALSE)

# Set up a race attribute on the network;
# This function is the long-form version of the %v% function
nw <- set.vertex.attribute(nw, attrname = "race", value = rep(0:1, each = 250))

# Parameterizing the model
formation <- ~edges + nodefactor("race") + nodematch("race")

# Target stats for model 1 (the null model)
target.stats <- c(125, 125, 62.5)

# Parameterizing the dissolution model (will be the same for model 2)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 40, d.rate = 0.01)
coef.diss

# Fit the model 1
est1 <- netest(nw, formation, target.stats, coef.diss)

# Examine the model coefficients
summary(est1)

# Model diagnostics. We can input other terms to diagnose using nwstats.formula,
# and here we only change the reference category of the nodefactor term (levels = NULL
# means that there is no reference category, so we can see the stats for both race
# 0 and race 1 groups)
dx1 <- netdx(est1, nsims = 10, nsteps = 1000, ncores = 4,
             nwstats.formula = ~edges +
                                nodefactor("race", levels = NULL) +
                                nodematch("race"))
dx1
plot(dx1)

# Target stats for model 2: higher mean degree for race 1 and high assortative mixing
target.stats <- c(125, 187.5, 112.5)

# Fit the second model (everything except target.stats is the same)
est2 <- netest(nw, formation, target.stats, coef.diss)

# Examine the model coefficients here
summary(est2)

# Model diagnostics for model 2
dx2 <- netdx(est2, nsims = 10, nsteps = 1000, ncores = 4,
            nwstats.formula = ~edges +
                               nodefactor("race", levels = NULL) +
                               nodematch("race"))

dx2
plot(dx2)

# Parameterizing the model
param <- param.net(inf.prob = 0.1, act.rate = 5,
                   a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01)

# Initial conditions
init <- init.net(i.num = 50)

# Control settings (reduced number of simulations for computational efficiency)
control <- control.net(type = "SI", nsteps = 500, nsims = 1,
                       epi.by = "race", delete.nodes = TRUE)

# Simulate the two counterfactual models
sim1 <- netsim(est1, param, init, control)
sim2 <- netsim(est2, param, init, control)

# Examine the model results
sim1

# Post-simulation diagnostics
par(mfrow = c(1, 2))
plot(sim1, type = "formation", stats = "edges")
plot(sim2, type = "formation", stats = "edges")

# Model prevalence overall, using the add argument to plot one model on top of
# the other
plot(sim1, y = "i.num", qnts = 1, mean.col = "steelblue",
     qnts.col = "steelblue", main = "Total Prevalence")
plot(sim2, y = "i.num", qnts = 1, mean.col = "firebrick",
     qnts.col = "firebrick", add = TRUE)
legend("topleft", c("Model 1", "Model 2"), lwd = 3,
       col = c("steelblue", "firebrick"), bty = "n")

# Model results stratified by race
par(mfrow = c(1, 2))
plot(sim1, y = c("i.num.race0", "i.num.race1"),  legend = TRUE, qnts = 1,
     ylim = c(0, 200), main = "M1: Disease Prevalence by Race")
plot(sim2, y = c("i.num.race0", "i.num.race1"), legend = TRUE,  qnts = 1,
     ylim = c(0, 200), main = "M2: Disease Prevalence by Race")

