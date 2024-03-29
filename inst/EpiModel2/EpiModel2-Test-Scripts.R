
###
### EpiModel Overhaul
### Worked Examples
###

suppressMessages(library("EpiModel"))
par(mar = c(3, 3, 2, 1), mgp = c(2, 1, 0))

# No Vital Dynamics -------------------------------------------------------

## One Group ##

# SI
num <- 1000
nw <- network.initialize(num, directed = FALSE)
formation <- ~edges
target.stats <- 400
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
est <- netest(nw, formation, target.stats, coef.diss)

param <- param.net(inf.prob = 0.1, act.rate = 5)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = FALSE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = 1)
summary(sim, at = 50)

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net(dat, at)
#   cat("*")
# }

# SIS
param <- param.net(inf.prob = 0.5, act.rate = 1, rec.rate = 0.02)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = 1)
summary(sim, at = 250)

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net(dat, at)
#   dat <- recovery.net(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net(dat, at)
#   cat("*")
# }

# SIR
param <- param.net(inf.prob = 0.5, act.rate = 1, rec.rate = 0.02)
init <- init.net(i.num = 10, r.num = 5)
control <- control.net(type = "SIR", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = FALSE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = 1)
summary(sim, at = 50)

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net(dat, at)
#   dat <- recovery.net(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net(dat, at)
#   cat("*")
# }


## Two Group ##

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
control <- control.net(type = "SI", nsteps = 200, nsims = 1, ncores = 1,
                       tergmLite = FALSE)

sim <- netsim(est, param, init, control)
plot(sim)
summary(sim, at = 50)

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net.grp(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net.grp(dat, at)
#   cat("*")
# }

# SIS
param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
init <- init.net(i.num = 20, i.num.g2 = 20)
control <- control.net(type = "SIS", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = FALSE)

sim <- netsim(est, param, init, control)
plot(sim)
summary(sim, at = 225)

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net(dat, at)
#   dat <- recovery.net(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net(dat, at)
#   cat("*")
# }

# SIR
param <- param.net(inf.prob = 0.4, inf.prob.g2 = 0.2,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
control <- control.net(type = "SIR", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, mean.col = c(1, 2, 3, 1, 2, 3)) # TODO: fix color issue with defaults
summary(sim, at = 225)

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net.grp(dat, at)
#   dat <- recovery.net.grp(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net.grp(dat, at)
#   cat("*")
# }


# Vital Dynamics -----------------------------------------------------------

## One Group ##

# SI
num <- 1000
nw <- network.initialize(num, directed = FALSE)
formation <- ~edges
target.stats <- 400
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 25, d.rate = 0.005)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
param <- param.net(inf.prob = 0.4, act.rate = 1,
                   a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 100, nsims = 1, ncores = 1,
                       tergmLite = FALSE, verbose = TRUE,
                       resimulate.network = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, y = "si.flow", ylim = c(0, 20))
summary(sim, at = 250)

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net(dat, at)
#   dat <- departures.net(dat, at)
#   dat <- arrivals.net(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net(dat, at)
#   cat("\n", at)
# }

# SIS
param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                   a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 250, nsims = 1, ncores = 1,
                       tergmLite = FALSE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, y = c("si.flow", "is.flow"), mean.smooth = FALSE, qnts = FALSE,
     legend = TRUE)
summary(sim, at = 250)

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net(dat, at)
#   dat <- recovery.net(dat, at)
#   dat <- departures.net(dat, at)
#   dat <- arrivals.net(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net(dat, at)
#   cat("*")
# }

# SIR
param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                   a.rate = 0.005, di.rate = 0.005, ds.rate = 0.005,
                   dr.rate = 0.005)
init <- init.net(i.num = 10, r.num = 10)
control <- control.net(type = "SIR", nsteps = 300, nsims = 10, ncores = 5,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, qnts = 1)
plot(sim, y = "num", ylim = c(800, 1200))

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net(dat, at)
#   dat <- recovery.net(dat, at)
#   dat <- departures.net(dat, at)
#   dat <- arrivals.net(dat,at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net(dat, at)
#   cat("*")
# }


## Two Group ##

# SI
num1 <- num2 <- 500
nw <- network.initialize(num1 + num2, directed = FALSE)
nw <- set.vertex.attribute(nw, "group", rep(1:2, each = num1))
formation <- ~ edges + nodematch("group")
target.stats <- c(400, 0)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges),
                               duration = 25, d.rate = 0.005)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                   act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                   di.rate = 0.005, ds.rate = 0.005,
                   di.rate.g2 = 0.005, ds.rate.g2 = 0.005)
init <- init.net(i.num = 50, i.num.g2 = 50)
control <- control.net(type = "SI", nsteps = 250, nsims = 1, ncores = 1,
                       tergmLite = FALSE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)

df <- as.data.frame(sim, out = "mean")
df

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net.grp(dat, at)
#   dat <- departures.net.grp(dat, at)
#   dat <- arrivals.net.grp(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net.grp(dat, at)
#   cat("*")
# }

# SIS
param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                   act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                   di.rate = 0.005, ds.rate = 0.005,
                   di.rate.g2 = 0.005, ds.rate.g2 = 0.005,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
init <- init.net(i.num = 50, i.num.g2 = 50)
control <- control.net(type = "SIS", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, qnts = 1, ylim = c(0, 500))

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net.grp(dat, at)
#   dat <- recovery.net.grp(dat, at)
#   dat <- departures.net.grp(dat, at)
#   dat <- arrivals.net.grp(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net.grp(dat, at)
#   cat("*")
# }

# SIR
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   act.rate = 5, a.rate = 0.005, a.rate.g2 = 0.005,
                   di.rate = 0.005, ds.rate = 0.005,
                   di.rate.g2 = 0.005, ds.rate.g2 = 0.005,
                   dr.rate = 0.005, dr.rate.g2 = 0.005,
                   rec.rate = 0.005, rec.rate.g2 = 0.005)
init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
control <- control.net(type = "SIR", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)
plot(sim, qnts = 1, ylim = c(0, 500))

# set.seed(123)
# crosscheck.net(est, param, init, control)
# dat <- initialize.net(est, param, init, control, s = 1)
#
# for (at in 2:control$nsteps) {
#   dat <- resim_nets(dat, at)
#   dat <- infection.net.grp(dat, at)
#   dat <- recovery.net.grp(dat, at)
#   dat <- departures.net.grp(dat, at)
#   dat <- arrivals.net.grp(dat, at)
#   dat <- nwupdate.net(dat, at)
#   dat <- prevalence.net.grp(dat, at)
#   cat("*")
# }


# Attr Copying Tests ------------------------------------------------------

## Closed Pop, 1G
num <- 1000
nw <- network.initialize(num, directed = FALSE)
nw <- set.vertex.attribute(nw, "race", sample(c("B", "W"), num, TRUE))
nw <- set.vertex.attribute(nw, "age", sample(15:65, num, TRUE))
formation <- ~ edges + nodematch("race") + absdiff("age")
target.stats <- c(400, 0, 800)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 25,
                               d.rate = 0.005)
est <- netest(nw, formation, target.stats, coef.diss)
summary(est)

param <- param.net(inf.prob = 0.1, act.rate = 5)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = 1)
summary(sim, at = 50)


## Open Pop, 1G
param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                   a.rate = 0.005, ds.rate = 0.005, di.rate = 0.005,
                   dr.rate = 0.005)
init <- init.net(i.num = 50, r.num = 0)
control <- control.net(type = "SIS", nsteps = 100, nsims = 1, ncores = 1,
                       tergmLite = FALSE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)



## Closed Pop, 2G
num1 <- num2 <- 500
nw <- network.initialize(num1 + num2, directed = FALSE)
nw <- set.vertex.attribute(nw, "group", rep(1:2, each = num1))
nw <- set.vertex.attribute(nw, "race", sample(c("B", "W"), num1 + num2,
                                              replace = TRUE))
nw <- set.vertex.attribute(nw, "age", sample(15:65, num1 + num2, TRUE))
formation <- ~ edges + nodematch("group") + absdiff("age")
target.stats <- c(400, 0, 800)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 25,
                               d.rate = 0.005)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3, act.rate = 1,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
init <- init.net(i.num = 50, i.num.g2 = 50)
control <- control.net(type = "SIS", nsteps = 100, nsims = 1, ncores = 1,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
plot(sim, qnts = FALSE, sim.lines = TRUE)



## Open Pop, 2G
param <- param.net(inf.prob = 0.5, inf.prob.g2 = 0.3,
                   act.rate = 1, a.rate = 0.005, a.rate.g2 = NA,
                   rec.rate = 0.01, rec.rate.g2 = 0.01,
                   ds.rate = 0.005, di.rate = 0.005,
                   ds.rate.g2 = 0.005, di.rate.g2 = 0.005,
                   dr.rate = 0.005, dr.rate.g2 = 0.005)
init <- init.net(i.num = 50, i.num.g2 = 50, r.num = 0, r.num.g2 = 0)
control <- control.net(type = "SI", nsteps = 100, nsims = 1, ncores = 1,
                       tergmLite = FALSE)

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
control <- control.net(type = "SIR", nsteps = 250, nsims = 5, ncores = 5,
                       tergmLite = TRUE)

sim <- netsim(est, param, init, control)
