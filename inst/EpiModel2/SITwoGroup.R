## EpiModel 2.0 Worked Example - SI Model with no Vital Dynamics on a Two-Group
## Network:
## Description: A two-group closed population of size 1000 (500 individuals in
## each group) in which an epidemic of type SI is introduced.
## -----------------------------------------------------------------------------

## update to new initialization function network_initialize.
num1 <- num2 <- 500
nw <- network_initialize(n = num1 + num2)

##Update to new vertex assignment function set_vertex_attribute
nw <- set_vertex_attribute(nw, "group", rep(1:2, c(num1, num2)))

## If looking to replicate 'bipartite' network:
formation <- ~edges + nodematch("group")
target.stats <- c(400, 0)

## If interested instead in within group mixing; here 100 of the total edges
## match on group:
#target.stats <- c(400, 100)

coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss)

## Epidemic model parameterization
param <- param.net(inf.prob = 0.2, inf.prob.g2 = 0.2)
init <- init.net(i.num = 20, i.num.g2 = 20)
control <- control.net(type = "SI", nsteps = 100, resimulate.network = FALSE,
                       tergmLite = FALSE)

## Simulate the epidemic model
sim <- netsim(est, param, init, control)

## Print the results
print(sim)

## Plot the results
plot(sim)
