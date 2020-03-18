
###
### EpiModel Overhaul
### Worked Examples
###

# remotes::install_github("statnet/EpiModel", ref = "EpiModel_Overhaul_s4")
suppressMessages(library("EpiModel"))
par(mar = c(3,3,2,1), mgp = c(2,1,0))

# No Vital Dynamics -------------------------------------------------------

## One Group ##

# SI
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 10)
param <- param.net(inf.prob = 0.1, act.rate = 5)
control <- control.net(type = "SI", nsteps = 100, nsims = 5, ncores = 5, tgl = TRUE, depend = FALSE,
                       save.nwstats = FALSE)

sim <- netsim(est, param, init, control)
plot(sim)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net(dat, at)
  cat("*")
}

# SIS
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50)
param <- param.net(inf.prob = 0.02, act.rate = 5, rec.rate = 0.02)
control <- control.net(type = "SIS", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net(dat, at)
  dat <- recovery.net(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net(dat, at)
  cat("*")
}

# SIR
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 10, r.num = 5)
param <- param.net(inf.prob = 0.02, act.rate = 5, rec.rate = 0.02)
control <- control.net(type = "SIR", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net(dat, at)
  dat <- recovery.net(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net(dat, at)
  cat("*")
}


## Two Group ##

# SI
num1 <- num2 <- 250
nw <- network::network.initialize(num1 + num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(1:2, each = num1))
formation <- ~edges + nodefactor("group") + nodematch("group")
target.stats <- c(150, 150, 75)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 20, i.num.g2 = 20)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.1, act.rate = 5)
control <- control.net(type = "SI", nsteps = 100, nsims =5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net.grp(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net.grp(dat, at)
  cat("*")
}

# SIS
num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(1:2, each = num1))
formation <- ~ edges + nodefactor("group") + nodematch("group")
target.stats <- c(150, 150, 75)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 20, i.num.g2 = 20)
param <- param.net(inf.prob = 0.02, inf.prob.g2 = 0.02,
                   rec.rate = 0.02, rec.rate.g2 = 0.02,
                   act.rate = 5)
control <- control.net(type = "SIS", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net(dat, at)
  dat <- recovery.net(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net(dat, at)
  cat("*")
}

# SIR
num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodefactor("group") + nodematch("group")
target.stats <- c(150, 150, 75)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
param <- param.net(inf.prob = 0.02, inf.prob.g2 = 0.02,
                   rec.rate = 0.02, rec.rate.g2 = 0.02,
                   act.rate = 5)
control <- control.net(type = "SIR", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net.grp(dat, at)
  dat <- recovery.net.grp(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net.grp(dat, at)
  cat("*")
}


# Vital Dynamics -----------------------------------------------------------

## One Group ##

# SI
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 10)
param <- param.net(inf.prob = 0.4, act.rate = 5,
                   a.rate = 0.02, ds.rate = 0.02,
                   di.rate = 0.02)
control <- control.net(type = "SI", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net(dat, at)
  dat <- departures.net(dat, at)
  dat <- arrivals.net(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net(dat, at)
  cat("*")
}

# SIS
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 10)
param <- param.net(inf.prob = 0.02, act.rate = 5, rec.rate = 0.02,
                   a.rate = 0.02, ds.rate = 0.02, di.rate = 0.02)
control <- control.net(type = "SIS", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net(dat, at)
  dat <- recovery.net(dat, at)
  dat <- departures.net(dat, at)
  dat <- arrivals.net(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net(dat, at)
  cat("*")
}

# SIR
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 10, r.num = 10)
param <- param.net(inf.prob = 0.4, act.rate = 5, rec.rate = 0.02,
                   a.rate = 0.02, di.rate = 0.02, ds.rate = 0.02,
                   dr.rate = 0.02)
control <- control.net(type = "SIR", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net(dat, at)
  dat <- recovery.net(dat, at)
  dat <- departures.net(dat, at)
  dat <- arrivals.net(dat,at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net(dat, at)
  cat("*")
}


## Two Group ##

# SI
num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodefactor("group") + nodematch("group")
target.stats <- c(150, 150, 75)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50, i.num.g2 = 50)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   act.rate = 5, a.rate = 0.02, a.rate.g2 = 0.02,
                   di.rate = 0.03, ds.rate = 0.03,
                   di.rate.g2 = 0.03, ds.rate.g2 = 0.03)
control <- control.net(type = "SI", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net.grp(dat, at)
  dat <- departures.net.grp(dat, at)
  dat <- arrivals.net.grp(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net.grp(dat, at)
  cat("*")
}

# SIS
num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodefactor("group") + nodematch("group")
target.stats <- c(150, 150, 75)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50, i.num.g2 = 50)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   act.rate = 5, a.rate = 0.02, a.rate.g2 = 0.02,
                   di.rate = 0.03, ds.rate = 0.03,
                   di.rate.g2 = 0.03, ds.rate.g2 = 0.03,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
control <- control.net(type = "SIS", nsteps = 100, nsims =5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net.grp(dat, at)
  dat <- recovery.net.grp(dat, at)
  dat <- departures.net.grp(dat, at)
  dat <- arrivals.net.grp(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net.grp(dat, at)
  cat("*")
}

# SIR
num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodefactor("group") + nodematch("group")
target.stats <- c(150, 150, 75)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 10, i.num.g2 = 10, r.num = 5, r.num.g2 = 5)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   act.rate = 5, a.rate = 0.02, a.rate.g2 = 0.02,
                   di.rate = 0.03, ds.rate = 0.03,
                   di.rate.g2 = 0.03, ds.rate.g2 = 0.03,
                   dr.rate = 0.02, dr.rate.g2 = 0.02,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)
control <- control.net(type = "SIR", nsteps = 100, nsims = 5, tgl = TRUE,
                       save.nwstats = FALSE)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:control$nsteps) {
  dat <- resim_nets(dat, at)
  dat <- infection.net.grp(dat, at)
  dat <- recovery.net.grp(dat, at)
  dat <- departures.net.grp(dat, at)
  dat <- arrivals.net.grp(dat, at)
  dat <- nwupdate.net(dat, at)
  dat <- prevalence.net.grp(dat, at)
  cat("*")
}
