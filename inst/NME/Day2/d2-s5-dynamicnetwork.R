
##
## Dynamic Network Modeling in EpiModel
## Day 2 | Network Modeling for Epidemics
##

## Set up

set.seed(0)
library(EpiModel)

# Scenario 1: MSM Networks

net1 <- network.initialize(100, directed = FALSE)

coef.diss.1 <- dissolution_coefs(~offset(edges), 90)
coef.diss.1

fit1 <- netest(net1,
               formation = ~edges,
               target.stats = 20,
               coef.diss = coef.diss.1)

summary(fit1)

sim1 <- netdx(fit1, nsteps = 1000, nsims = 10,
              keep.tedgelist = TRUE)
sim1

plot(sim1, type = "formation")
plot(sim1, type = "duration")
plot(sim1, type = "dissolution")

tel <- as.data.frame(sim1, sim = 1)

hist(tel$duration)
mean(tel$duration[tel$onset < 100])
sum(tel$terminus.censored == TRUE)
plot(tel$onset, tel$terminus)
table(c(tel$head,tel$tail))
hist(table(c(tel$head,tel$tail)))


# Scenario 2: Modifying Network Size

net2 <- network.initialize(1000, directed = FALSE)

fit2 <- netest(net2,
   formation = ~edges,
   target.stats = 200,
   coef.diss = dissolution_coefs(~offset(edges), 90))

sim2 <- netdx(fit2, nsteps = 1000, nsims = 10, keep.tedgelist = TRUE)

plot(sim2, type = "formation")
plot(sim2, type = "duration")
plot(sim2, type = "dissolution")


# Scenario 3: Adding ERGM Terms

n <- 500
net3 <- network.initialize(n, directed = FALSE)
set.vertex.attribute(net3, "race", c(rep("B", n/2), rep("W", n/2)))
net3

form.formula.3 <- ~edges + nodematch("race") + degree(0) + concurrent
target.stats.3 <- c(0.9*n/2, (0.9*n/2)*(5/6), 0.36*n, 0.18*n)

diss.formula.3 <- ~offset(edges) + offset(nodematch("race"))

?dissolution_coefs

fit3 <- netest(net3,
             formation = form.formula.3,
             target.stats = target.stats.3,
             coef.diss = dissolution_coefs(~offset(edges) + offset(nodematch("race")),
                                           c(200, 100)))

sim3 <- netdx(fit3, nsteps = 1000, nsims = 10, keep.tedgelist = TRUE)
sim3

plot(sim3, type = "formation")
plot(sim3, type = "duration")

race <- get.vertex.attribute(net3, "race")

tel3 <- as.data.frame(sim3, sim = 1)
mean(tel3$duration[(race[tel3$tail] != race[tel3$head]) & tel3$onset < 100])
mean(tel3$duration[(race[tel3$tail] == race[tel3$head]) & tel3$onset < 100])


# Scenario 4: Full STERGM

fit4 <- netest(net3,
    formation = form.formula.3,
    target.stats = target.stats.3,
    coef.diss = dissolution_coefs(~offset(edges)+offset(nodematch("race")),
                                  c(20, 10)))

sim4 <- netdx(fit4, nsteps = 1000, nsims = 10, keep.tedgelist = TRUE)
sim4

plot(sim4, type = "formation")

# Without approximation
fit5 <- netest(net3,
         formation = form.formula.3,
         target.stats = target.stats.3,
         coef.diss = dissolution_coefs(~offset(edges) + offset(nodematch("race")),
                                       c(20, 10)),
         edapprox = FALSE)

sim5 <- netdx(fit5, nsteps = 1000, nsims = 10, keep.tedgelist = TRUE)
sim5

plot(sim5, type = "formation")

race <- get.vertex.attribute(net3, "race")

tel5 <- as.data.frame(sim5, sim = 1)
mean(tel5$duration[(race[tel5$tail] != race[tel5$head]) & tel5$onset < 100])
mean(tel5$duration[(race[tel5$tail] == race[tel5$head]) & tel5$onset < 100])
