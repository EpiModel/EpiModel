## ----setup, include = FALSE----------------------------------------------
require(knitr)
require(EpiModel)
opts_chunk$set(comment = NA, message = FALSE, tidy = FALSE)

## ----netInit-------------------------------------------------------------
nw <- network.initialize(n = 1000, directed = FALSE)
nw <- set.vertex.attribute(nw, "race", rep(0:1, each = 500))

## ----netForm1------------------------------------------------------------
formation <- ~edges + nodefactor("race") + nodematch("race") + concurrent

## ----netTs1--------------------------------------------------------------
target.stats <- c(250, 375, 225, 100)

## ----netDiss1------------------------------------------------------------
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
coef.diss

## ----netestArgs, echo = FALSE--------------------------------------------
args(netest)

## ----netEst1, results = "hide"-------------------------------------------
est1 <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE)

## ----netDx1, results = "hide"--------------------------------------------
dx <- netdx(est1, nsims = 5, nsteps = 500,
            nwstats.formula = ~edges + nodefactor("race", base = 0) +
                               nodematch("race") + concurrent)

## ----printNetDx1---------------------------------------------------------
dx

## ----plotNetDx1----------------------------------------------------------
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(dx)

## ----plotNetDx2----------------------------------------------------------
par(mfrow = c(1, 2))
plot(dx, type = "duration")
plot(dx, type = "dissolution")

## ----netParam1-----------------------------------------------------------
param <- param.net(inf.prob = 0.1, act.rate = 5, rec.rate = 0.02)

## ----netInit1------------------------------------------------------------
status.vector <- c(rbinom(500, 1, 0.1), rep(0, 500))
status.vector <- ifelse(status.vector == 1, "i", "s")
init <- init.net(status.vector = status.vector)

## ----netControl1, results = "hide"---------------------------------------
control <- control.net(type = "SIS", nsteps = 500, nsims = 10, epi.by = "race")

## ----netsim1, results = "hide", cache = TRUE-----------------------------
sim1 <- netsim(est1, param, init, control)

## ----printSim1-----------------------------------------------------------
sim1

## ----summSim1------------------------------------------------------------
summary(sim1, at = 500)

## ----netAdf1-------------------------------------------------------------
head(as.data.frame(sim1), 10)

## ----netAdf2-------------------------------------------------------------
head(as.data.frame(sim1, out = "mean"), 10)

## ----netGnet1------------------------------------------------------------
nw <- get_network(sim1, sim = 1)
nw

## ----netGetTm------------------------------------------------------------
head(get_transmat(sim1, sim = 1), 10)

## ----netPlot1------------------------------------------------------------
par(mfrow = c(1,1), mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim1)

## ----netPlot2------------------------------------------------------------
plot(sim1, mean.line = FALSE, qnts = FALSE, sim.lines = TRUE)

## ----netPlot3------------------------------------------------------------
plot(sim1, y = c("si.flow", "is.flow"), qnts = 1, legend = TRUE)

## ----netPlot4------------------------------------------------------------
plot(sim1, y = c("i.num.race0", "i.num.race1"), legend = TRUE)

## ----net19---------------------------------------------------------------
par(mfrow = c(1,2), mar = c(0,0,1,0))
plot(sim1, type = "network", at = 1, col.status = TRUE,
     main = "Prevalence at t1")
plot(sim1, type = "network", at = 500, col.status = TRUE,
     main = "Prevalence at t500")

## ----net2Init------------------------------------------------------------
num.m1 <- 500
num.m2 <- 500
nw <- network.initialize(num.m1 + num.m2, bipartite = num.m1, directed = FALSE)

## ----net2DegDist---------------------------------------------------------
deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)

## ----net2PoisDist--------------------------------------------------------
pois.dists <- c(dpois(0:2, lambda = 0.66), ppois(2, lambda = 0.66, lower = FALSE))

## ----net2PoisPlot--------------------------------------------------------
par(mar = c(3, 3, 2, 1), mfrow = c(1, 1))
cols <- transco(RColorBrewer::brewer.pal(4, "Set1"), 0.8)
barplot(cbind(deg.dist.m1, deg.dist.m2, pois.dists),
        beside = FALSE, ylim = c(0, 1), col = cols)
legend("topright", legend = paste0("deg", 3:0),
       pch = 15, col = rev(cols), bg = "white")

## ----net2CheckDegDist----------------------------------------------------
check_bip_degdist(num.m1, num.m2, deg.dist.m1, deg.dist.m2)

## ----net2FormMod---------------------------------------------------------
formation <- ~edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(330, 200, 275, 240, 205)

## ----net2DissMod---------------------------------------------------------
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25, 
                               d.rate = 0.005)
coef.diss

## ----net2Est, results = "hide"-------------------------------------------
est2 <- netest(nw, formation, target.stats, coef.diss)

## ----net2Dx, results = "hide"--------------------------------------------
dx <- netdx(est2, nsims = 5, nsteps = 500)
dx

## ----net2DxPlot1---------------------------------------------------------
plot(dx, stats = "edges")

## ----net2DxPlot2---------------------------------------------------------
plot(dx, stats = c("b1deg1", "b2deg1"), sim.lines = TRUE, sim.lwd = 0.4,
     qnts = FALSE, mean.lwd = 4, mean.smooth = FALSE)

## ----net2Param-----------------------------------------------------------
param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.1,
                   a.rate = 0.005, a.rate.m2 = NA,
                   ds.rate = 0.005, ds.rate.m2 = 0.005,
                   di.rate = 0.005, di.rate.m2 = 0.005)

## ----net2RunSim, results = "hide", cache = TRUE--------------------------
init <- init.net(i.num = 50, i.num.m2 = 50)
control <- control.net(type = "SI", nsims = 5, nsteps = 500,
                       nwstats.formula = ~edges + meandeg, delete.nodes = TRUE)
sim2 <- netsim(est2, param, init, control)

## ----net2Print-----------------------------------------------------------
sim2

## ----net2Plot1-----------------------------------------------------------
plot(sim2, type = "formation", plots.joined = FALSE)

## ----net2Plot2-----------------------------------------------------------
plot(sim2, popfrac = FALSE)

