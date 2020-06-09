
##
## Tutorial: SIR Epidemic in a Bipartite Network
## Day 3 | Network Modeling for Epidemics
##

# Load EpiModel
library(EpiModel)


# Network model estimation ------------------------------------------------

# Initialize the network
num.m1 <- num.m2 <- 250
nw <- network.initialize(n = num.m1 + num.m2,
                         bipartite = num.m1,
                         directed = FALSE)
nw

# Fractional degree distributions by mode
deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.54, 0.31, 0.10, 0.05)

# Compare those mode-specific distributions to expected given mean degree
pois.dists <- c(dpois(0:2, lambda = 0.66),
                ppois(2, lambda = 0.66, lower.tail = FALSE))

# Plot the comparison
par(mar = c(3,3,2,1), mgp = c(2,1,0), mfrow = c(1,1))
barplot(cbind(deg.dist.m1, deg.dist.m2, pois.dists),
        beside = TRUE, ylim = c(0, 0.6), col = rainbow(4))
legend("topright", legend = paste0("deg", 3:0),
       pch = 15, col = rev(rainbow(4)),
       cex = 0.9, bg = "white")

# Helper function to check two proportional degree distributions
check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)

# Formation formula and target statistics
formation <- ~edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(165, 100, 137.5, 135, 77.5)

# Dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 5, nsteps = 500,
            nwstats.formula = ~edges + b1degree(0:3) + b2degree(0:3))
dx

# Compare the simulated stats for "out-of-model" fit (degrees 2 and 3)
check_bip_degdist(num.m1, num.m2,
                  deg.dist.m1, deg.dist.m2)

# Plot the formation diagnostics
plot(dx, stats = c("edges", "b1deg1", "b2deg1"))

# Plot the dissolution diagnostics
plot(dx, type = "duration")


# Epidemic model simulation -----------------------------------------------

# Model parameters
param <- param.net(inf.prob = 0.2, inf.prob.m2 = 0.2,
                   rec.rate = 0.02, rec.rate.m2 = 0.02)

# Initial conditions
init <- init.net(i.num = 10, i.num.m2 = 10,
                 r.num = 0, r.num.m2 = 0)

# Control settings
control <- control.net(type = "SIR", nsims = 5, nsteps = 500)

# Simulate with netsim
sim <- netsim(est, param, init, control)

# Plot the results
plot(sim)

# Other ways to examine the epi data
par(mfrow = c(1,2))
plot(sim, y = c("i.num", "i.num.m2"), popfrac = TRUE,
     qnts = 0.5, ylim = c(0, 0.4), legend = TRUE)
plot(sim, y = c("si.flow", "si.flow.m2"),
     qnts = 0.5, ylim = c(0, 3), legend = TRUE)

# Static network plot
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot(sim, type = "network", col.status = TRUE, at = 50,
     sims = "mean", shp.bip = "square")


# Add new variables to the object with mutate_epi
sim <- mutate_epi(sim, ir.m1 = (si.flow / s.num) * 100 * 52)

# Plot the new variable
par(mar = c(3,3,2,1))
plot(sim, y = "ir.m1", qnts = 0.5, mean.lwd = 1.25, ylim = c(0, 200))

# Extract the data.frame for analysis
df <- as.data.frame(sim, out = "mean")

# Cumulative incidence
sum(df$si.flow)
sum(df$si.flow.m2)

# Time step and value of peak standardized incidence
which.max(df$ir.m1)
df$ir.m1[which.max(df$ir.m1)]

# Extract the transmission matrix
tm <- get_transmat(sim)
head(tm, 15)

# The proportion of transmissions that occured within the first 10 weeks of infection
mean(tm$infDur < 10)
