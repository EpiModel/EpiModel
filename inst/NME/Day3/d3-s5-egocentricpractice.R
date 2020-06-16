
##
## PRACTICE WITH EGOCENTRIC DATA AND TARGET STATISTICS
## Day 3 | Network Modeling for Epidemics
##

library(EpiModel)
mynet <- network.initialize(2000, directed = FALSE, bipartite = 1000)
cmty <- c(rep(1, 400), rep(2, 600), rep(1, 400), rep(2, 600))
mynet <- set.vertex.attribute(mynet, "cmty", cmty)

par(mar = c(3,3,2,1))
sides <- cmty + 2
plot(mynet, vertex.col = rep(1:2, each = 1000), vertex.sides = sides)

plot(get.vertex.attribute(mynet, "cmty"))







# female degree
fd <- c(1, 0, 2, 1, 1, 1, 1, 0, 0, 1)

# female frac deg dist
prop.table(table(fd))

# male degree
md <- c(2, 0, 0, 1, 1, 0, 2, 1, 2, 0)

# male frac deg dist
prop.table(table(md))

# empirical data
check_bip_degdist(1000, 1000,
                  c(0.3, 0.6, 0.1),
                  c(0.4, 0.3, 0.3))

# one solution (change first mode only)
check_bip_degdist(1000, 1000,
                  c(0.3, 0.5, 0.2),
                  c(0.4, 0.3, 0.3))

# solution by averaging
check_bip_degdist(1000, 1000,
                  c(0.3, 0.55, 0.15),
                  c(0.4, 0.35, 0.25))

total.edges <- 850
prop.same.cmty.edges <- 12/17
nodematch.stat <- total.edges * prop.same.cmty.edges





























formation <- ~edges + nodematch('cmty') + degrange(from=3) + b1degree(1) + b2degree(1)

target.stats <- c(850, 600, 0, 550, 350)

myfit <- netest(mynet,
                formation = formation,
                target.stats = target.stats,
                coef.diss = dissolution_coefs(~offset(edges), 60))

mydx <- netdx(myfit, nsims = 10, nsteps = 100)
mydx
get_nwstats(mydx)
plot(mydx)

mycontrol <- control.net(type = "SIS", nsteps = 1000, nsims = 1,
                         nwstats.formula = ~edges + nodematch('cmty') + b1degree(0:5) + b2degree(0:5),
                         verbose = TRUE)
myinit <- init.net(i.num = 10, i.num.m2 = 10)
myparam <- param.net(inf.prob = 0.5, act.rate = 0.6,
                     rec.rate = 0.05, rec.rate.m2 = 0.05)

mySIS <- netsim(myfit,
                param = myparam,
                control = mycontrol,
                init = myinit)
plot(mySIS)

get_nwstats(mySIS)
plot(mySIS, type = "formation", sim.lines = TRUE)
