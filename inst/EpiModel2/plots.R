
## SI Model without Network Feedback
# Initialize network and set network model parameters
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)

# Estimate the network model
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Simulate the epidemic model
param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.15, rec.rate = 0.01, rec.rate.g2 = 0.01)
init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 0, r.num.g2 = 0)
control <- control.net(type = "SIR", nsteps = 20, nsims = 3,
                       verbose = FALSE, save.nwstats = TRUE,
                       nwstats.formula = ~edges + meandeg + concurrent)
mod <- netsim(est, param, init, control)

# Plot epidemic trajectory
plot(mod)
plot(mod, y = c("i.num", "i.num.g2", "r.num", "r.num.g2"), legend = TRUE)
plot(mod, y = c("i.num", "i.num.g2"), mean.col = 1:2, qnts.col = 1:2, legend = TRUE)
plot(mod, type = "epi", grid = TRUE)
plot(mod, type = "epi", popfrac = TRUE)
plot(mod, type = "epi", y = "si.flow", qnts = 1, ylim = c(0, 4))

## SI Model without Network Feedback
# Initialize network and set network model parameters
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)

# Estimate the network model
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Simulate the epidemic model
param <- param.net(inf.prob = 0.3, rec.rate = 0.05)
init <- init.net(i.num = 10, r.num = 0)
control <- control.net(type = "SIR", nsteps = 20, nsims = 3,
                       verbose = FALSE, save.nwstats = TRUE,
                       nwstats.formula = ~edges + meandeg + concurrent)
mod <- netsim(est, param, init, control)

# Plot epidemic trajectory
plot(mod)
plot(mod, type = "epi", grid = TRUE)
plot(mod, type = "epi", popfrac = TRUE)
plot(mod, type = "epi", y = "si.flow", qnts = 1, ylim = c(0, 4), legend = TRUE)
