context("Network models with tergmLite")

################################################################################

test_that("tergmLite: 1G, Closed", {

# SI
num <- 1000
nw <- network.initialize(num, directed = FALSE)
formation <- ~edges
target.stats <- 400
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
est <- netest(nw, formation, target.stats, coef.diss)

param <- param.net(inf.prob = 0.1, act.rate = 5)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 20, nsims = 1, ncores = 1,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = 1)
summary(sim, at = 20)


# SIS
param <- param.net(inf.prob = 0.5, act.rate = 1, rec.rate = 0.02)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = 1)
summary(sim, at = 25)


# SIR
param <- param.net(inf.prob = 0.5, act.rate = 1, rec.rate = 0.02)
init <- init.net(i.num = 10, r.num = 5)
control <- control.net(type = "SIR", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = 1)
summary(sim, at = 25)

})

test_that("tergmLite: 2G, Closed", {

# SI
num1 <- num2 <- 500
nw <- network.initialize(num1 + num2, directed = FALSE)
nw <- set.vertex.attribute(nw, "group", rep(1:2, each = num1))
formation <- ~edges + nodematch("group")
target.stats <- c(400, 0)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2)
init <- init.net(i.num = 20, i.num.g2 = 20)
control <- control.net(type = "SI", nsteps = 20, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim)
summary(sim, at = 20)

# SIS
param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
init <- init.net(i.num = 20, i.num.g2 = 20)
control <- control.net(type = "SIS", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim)
summary(sim, at = 25)

# SIR
param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
control <- control.net(type = "SIR", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, mean.col = c(1,2,3,1,2,3))
summary(sim, at = 25)

})

test_that("tergmLite: 1G, Open", {

# SI
num <- 1000
nw <- network.initialize(num, directed = FALSE)
formation <- ~edges
target.stats <- 400
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25,
                               d.rate = 0.005)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
param <- param.net(inf.prob = 0.4, act.rate = 1,
                   a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 25, nsims = 1, ncores = 1,
                       tergmLite = TRUE, verbose = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, y = "si.flow", ylim = c(0, 20))
summary(sim, at = 25)

# SIS
param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                   a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, y = c("si.flow", "is.flow"), legend = TRUE)
summary(sim, at = 25)

# SIR
param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                   a.rate = 0.005, di.rate = 0.005, ds.rate = 0.005,
                   dr.rate = 0.005)
init <- init.net(i.num = 10, r.num = 0)
control <- control.net(type = "SIR", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, qnts = 1)
plot(sim, y = "num", ylim = c(800, 1200))

})

test_that("tergmLite: 2G, Open", {

# SI
num1 <- num2 <- 500
nw <- network.initialize(num1 + num2, directed = FALSE)
nw <- set.vertex.attribute(nw, "group", rep(1:2, each = num1))
formation <- ~ edges + nodematch("group")
target.stats <- c(400, 0)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 25, d.rate = 0.005)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                   act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                   di.rate = 0.005, ds.rate = 0.005,
                   di.rate.g2 = 0.005, ds.rate.g2 = 0.005)
init <- init.net(i.num = 50, i.num.g2 = 50)
control <- control.net(type = "SI", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, qnts = 1, ylim = c(0, 500))

# SIS
param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                   act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                   di.rate = 0.005, ds.rate = 0.005,
                   di.rate.g2 = 0.005, ds.rate.g2 = 0.005,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
init <- init.net(i.num = 50, i.num.g2 = 50)
control <- control.net(type = "SIS", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, qnts = 1, ylim = c(0, 500))

# SIR
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   act.rate = 5, a.rate = 0.005, a.rate.g2 = 0.005,
                   di.rate = 0.005, ds.rate = 0.005,
                   di.rate.g2 = 0.005, ds.rate.g2 = 0.005,
                   dr.rate = 0.005, dr.rate.g2 = 0.005,
                   rec.rate = 0.005, rec.rate.g2 = 0.005)
init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
control <- control.net(type = "SIR", nsteps = 25, nsims = 1, ncores = 1, tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, qnts = 1, ylim = c(0, 500))

})
