
##
## Tutorial: SIS Epidemic in a One-Mode Network
## Day 3 | Network Modeling for Epidemics
##

library("EpiModel")

# Network model estimation ------------------------------------------------

# Initialize the network
nw <- network.initialize(n = 500, directed = FALSE)

# Define the formation model
formation <- ~edges + concurrent + degrange(from = 4)

# Input the appropriate target statistics for each term
target.stats <- c(175, 110, 0)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

# Review the arguments for the netest function
args(netest)

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 10, nsteps = 1000,
            nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent)
dx

# Plot the formation diagnostics
plot(dx)

nwstats1 <- get_nwstats(dx, sim = 1)
head(nwstats1, 20)

# Plot the dissolution diagnostics
par(mfrow = c(1, 2), mar = c(3,3,2,1))
plot(dx, type = "duration", mean.col = "black")
plot(dx, type = "dissolution", qnts = 0.5,
     mean.lines = FALSE, sim.lines = FALSE)

# Static (cross-sectional ERGM) model diagnostics (if the dynamic dx are not looking good)
dx.static <- netdx(est, nsims = 10000, dynamic = FALSE)
dx.static

par(mfrow = c(1,1))
plot(dx.static, sim.lines = TRUE, sim.lwd = 0.1)

nwstats2 <- get_nwstats(dx.static)
head(nwstats2, 20)


# Epidemic model simulation -----------------------------------------------

# Parameterizing an SIS epidemic
param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.01)

# Initial conditions
init <- init.net(i.num = 10)

# Control settings
control <- control.net(type = "SIS", nsims = 5, nsteps = 500)

# Run the network model simulation with netsim
sim <- netsim(est, param, init, control)

# Print the output to show the model contents
sim

# Plot the output to epi stats
par(mfrow = c(1, 1), mar = c(3,3,2,1), mgp = c(2,1,0))
plot(sim)

# Many different arguments for plot.netsim
par(mfrow = c(1, 2))
plot(sim, sim.lines = TRUE, mean.line = FALSE, qnts = FALSE, popfrac = TRUE)
plot(sim, mean.smooth = FALSE, qnts = 1, qnts.smooth = FALSE, popfrac = TRUE)

# Use the y argument to pull out non-default stats, such as incidence
par(mfrow = c(1,1))
plot(sim, y = c("si.flow", "is.flow"), qnts = FALSE,
     ylim = c(0, 10), legend = TRUE, main = "Flow Sizes")

# Static network plot from one sim at two time points
par(mar = c(0,0,0,0), mfrow = c(1, 2))
plot(sim, type = "network", col.status = TRUE, at = 1, sims = 1)
plot(sim, type = "network", col.status = TRUE, at = 500, sims = 1)

# Summary stats
summary(sim, at = 500)

# Convert model to a data frame for further analysis
df <- as.data.frame(sim)
head(df, 10)
tail(df, 10)

# Extracting mean values also possible
df <- as.data.frame(sim, out = "mean")
head(df, 10)
tail(df, 10)

# Extract the full dynamic network for further analysis
nw1 <- get_network(sim, sim = 1)
nw1

# Temporal edgelist
nwdf <- as.data.frame(nw1)
head(nwdf, 25)

# A transmission matrix contains the time-ordered chain of transmissions
tm1 <- get_transmat(sim, sim = 1)
head(tm1, 10)

tmPhylo <- as.phylo.transmat(tm1)
plot(tmPhylo, show.node.label = TRUE,
     root.edge = TRUE,
     cex = 0.5)

# Plotting with ggplot
df <- as.data.frame(sim)
df.mean <- as.data.frame(sim, out = "mean")

library(ggplot2)
ggplot() +
  geom_line(data = df, mapping = aes(time, i.num, group = sim), alpha = 0.25,
            lwd = 0.25, color = "firebrick") +
  geom_bands(data = df, mapping = aes(time, i.num),
             lower = 0.1, upper = 0.9, fill = "firebrick") +
  geom_line(data = df.mean, mapping = aes(time, i.num)) +
  theme_minimal()
