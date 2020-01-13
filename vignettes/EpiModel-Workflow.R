#EpiModel Overhaul v3: Worked Examples

remotes::install_github("statnet/EpiModel", ref = "overhaul_s3_b")

##One Group

#SI

#Closed population; type = "SI"
num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 10)
param <- param.net(inf.prob = 0.4, act.rate = 5)
control <- control.net(type = "SI", nsteps = 100, nsims = 5)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (i in 2:100) {
  dat <- resim_nets(dat, at = i)
  dat <- infection.net(dat, at = i)
  dat <- nw.update.net(dat, at = i)
  dat <- get_prev.net(dat, at = i)
}

#EpiMOdel (Master) and EpiModel Overhaul v3 Plots

#SIR/S

num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50)
param <- param.net(inf.prob = 0.4, act.rate = 5, rec.rate = 0.2)
control <- control.net(type = "SIS", nsteps = 100, nsims = 5)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- resim_nets(dat, at = i)
  dat <- recovery.net(dat, at = i)
  dat <- infection.net(dat, at = i)
  dat <- nw.update.net(dat, at = i)
  dat <- get_prev.net(dat, at = i)
}

num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50)
param <- param.net(inf.prob = 0.4, act.rate = 5, rec.rate = 0.2)
control <- control.net(type = "SIR", nsteps = 100, nsims = 5)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- resim_nets(dat, at = i)
  dat <- recovery.net(dat, at = i)
  dat <- infection.net(dat, at = i)
  dat <- nw.update.net(dat, at = i)
  dat <- get_prev.net(dat, at = i)
}

#Arrivals/Departures

num <- 200
nw <- network::network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 60
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50)
param <- param.net(inf.prob = 0.4, act.rate = 5, rec.rate = 0.2,
                   a.rate = 0.025, ds.rate = 0.02, di.rate = 0.01,
                   dr.rate = 0.02)
control <- control.net(type = "SIR", nsteps = 100, nsims = 5)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- resim_nets(dat, at = i)
  dat <- departures.net(dat, at = i)
  dat <- arrivals.net(dat, at = i)
  dat <- recovery.net(dat, at = i)
  dat <- infection.net(dat, at = i)
  dat <- nw.update.net(dat, at = i)
  dat <- get_prev.net(dat, at = i)
}

##Two Group

#SI

#Closed population; type = "SI"
num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodematch("group") + nodefactor("group")
target.stats <- c(150, 125, 62.5)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50, i.num.g2 = 50)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   act.rate = 5)
control <- control.net(type = "SI", nsteps = 100)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- resim_nets(dat, at = i)
  dat <- infection.net(dat, at = i)
  dat <- nw.update.net(dat, at = i)
  dat <- get_prev.net(dat, at = i)
}

#EpiMOdel (Master) and EpiModel Overhaul v3 Plots

#SIR/S

num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodematch("group") + nodefactor("group")
target.stats <- c(150, 125, 62.5)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50, i.num.g2 = 50)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   rec.rate = 0.02, rec.rate.g2 = 0.02,
                   act.rate = 5)
control <- control.net(type = "SIS", nsteps = 100)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- resim_nets(dat, at = i)
  dat <- recovery.net(dat, at = i)
  dat <- infection.net(dat, at = i)
  dat <- nw.update.net(dat, at = i)
  dat <- get_prev.net(dat, at = i)
}

num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodematch("group") + nodefactor("group")
target.stats <- c(150, 125, 62.5)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50, i.num.g2 = 50)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   rec.rate = 0.02, rec.rate.g2 = 0.02,
                   act.rate = 5)
control <- control.net(type = "SIR", nsteps = 100)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- resim_nets(dat, at = i)
  dat <- recovery.net(dat, at = i)
  dat <- infection.net(dat, at = i)
  dat <- nw.update.net(dat, at = i)
  dat <- get_prev.net(dat, at = i)
}

#Arrivals/Departures


num1 <- num2 <- 250
nw <- network::network.initialize(num1+num2, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "group", rep(c(1,2), each = num1))
formation <- ~ edges + nodematch("group") + nodefactor("group")
target.stats <- c(150, 125, 62.5)
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 50)
est <- netest(nw, formation, target.stats, coef.diss)

# Parameters
init <- init.net(i.num = 50, i.num.g2 = 50)
param <- param.net(inf.prob = 0.1, inf.prob.g2 = 0.2,
                   rec.rate = 0.02, rec.rate.g2 = 0.02,
                   a.rate = 0.025, a.rate.g2 = 0.03,
                   ds.rate = 0.01, di.rate = 0.01, dr.rate = 0.01,
                   act.rate = 5)
control <- control.net(type = "SIR", nsteps = 100)

set.seed(123)
crosscheck.net(est, param, init, control)
dat <- initialize.net(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- resim_nets(dat, at = i)
  dat <- departures.net(dat, at = i)
  dat <- arrivals.net(dat, at = i)
  dat <- recovery.net(dat, at = i)
  dat <- infection.net(dat, at = i)
  dat <- nw.update.net(dat, at = i)
  dat <- get_prev.net(dat, at = i)
}
