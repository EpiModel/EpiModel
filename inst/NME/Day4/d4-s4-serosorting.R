
##
## Tutorial: Serosorting
## Day 4 | Network Modeling for Epidemics
##

library(EpiModel)

n <- 500
nw <- network.initialize(n, directed = FALSE)

prev <- 0.2
infIds <- sample(1:n, n*prev)
nw <- set.vertex.attribute(nw, "status", "s")
nw <- set.vertex.attribute(nw, "status", "i", infIds)

mean.deg.inf <- 0.3
inedges.inf <- mean.deg.inf * n * prev

mean.deg.sus <- 0.8
inedges.sus <- mean.deg.sus * n * (1 - prev)

edges <- (inedges.inf + inedges.sus)/2

p <- inedges.sus/(edges*2)
q <- 1 - p
nn <- p^2
np <- 2*p*q
pp <- q^2
round(nn + pp, 3)

fit <- netest(nw,
              formation = ~edges + nodefactor("status"),
              target.stats = c(edges, inedges.sus),
              coef.diss = dissolution_coefs(~offset(edges), duration = 1))
sim <- netdx(fit, dynamic = FALSE, nsims = 1e4,
             nwstats.formula = ~edges + nodematch("status"))
stats <- get_nwstats(sim)
head(stats)
mean(stats$nodematch.status/stats$edges)

nmatch <- edges * 0.91

formation <- ~edges + nodefactor("status") + nodematch("status")
target.stats <- c(edges, inedges.sus, nmatch)

coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 50)

est <- netest(nw, formation, target.stats, coef.diss)

dx <- netdx(est, nsims = 5, nsteps = 500,
            nwstats.formula = ~edges +
                               meandeg +
                               nodefactor("status", levels = NULL) +
                               nodematch("status"), verbose = FALSE)
dx
plot(dx, plots.joined = FALSE)
plot(dx, type = "duration")

nw <- network.initialize(n, directed = FALSE)
est2 <- netest(nw, formation = ~edges, target.stats = edges, coef.diss)

dx2 <- netdx(est2, nsims = 5, nsteps = 1000,
             nwstats.formula = ~edges + meandeg)
dx2
plot(dx2, plots.joined = FALSE)

param <- param.net(inf.prob = 0.03)
init <- init.net()
control <- control.net(type = "SI", nsteps = 500, nsims = 5,
                       nwstats.formula = ~edges +
                                          meandeg +
                                          nodefactor("status", levels = NULL) +
                                          nodematch("status"),
                       delete.nodes = TRUE)

sim <- netsim(est, param, init, control)

param <- param.net(inf.prob = 0.03)
init <- init.net(i.num = n*prev, status.rand = FALSE)
control <- control.net(type = "SI", nsteps = 500, nsims = 5,
                       nwstats.formula = ~edges + meandeg,
                       delete.nodes = TRUE)
sim2 <- netsim(est2, param, init, control)

par(mfrow = c(1,2))
plot(sim, main = "Serosorting")
plot(sim2, main = "No Serosorting")

par(mfrow = c(1,1))
plot(sim, y = "i.num", popfrac = TRUE, sim.lines = FALSE, qnts = 1)
plot(sim2, y = "i.num", popfrac = TRUE, sim.lines = FALSE, qnts = 1,
     mean.col = "firebrick", qnts.col = "firebrick", add = TRUE)
legend("topleft", c("Serosort", "Non-Serosort"), lty = 1, lwd = 3,
       col = c("steelblue", "firebrick"), cex = 0.9, bty = "n")

plot(sim, type = "formation", plots.joined = FALSE)

plot(sim, type = "formation", stats = c("nodefactor.status.s",
                                        "nodefactor.status.i"))

plot(sim, type = "formation",
     stats = c("edges", "nodematch.status"),
     sims = 1, sim.lwd = 2)

