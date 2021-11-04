## EpiModel 2.0 Worked Example - SIR Model with Vital Dynamics:
## Description: An open population of size 1000 in which an epidemic
## of type SIR is introduced.
## -----------------------------------------------------------------------------

## update to new initialization function network_initialize.
num <- 1000
nw <- network_initialize(num)
formation <- ~edges
target.stats <- 400
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 20, d.rate = 0.005)
est <- netest(nw, formation, target.stats, coef.diss)

## Epidemic model parameterization
param <- param.net(inf.prob = 0.4, act.rate = 1, rec.rate = 0.02,
                   a.rate = 0.005, di.rate = 0.005, ds.rate = 0.005,
                   dr.rate = 0.005)
init <- init.net(i.num = 30, r.num = 5)

## control.net setting depend has been replaced by resimulate.network.
## control.net settings save.network/save.transmat have been folded in tergmLite
## functionality.
control <- control.net(type = "SIR", nsteps = 100, resimulate.network = TRUE,
                       tergmLite = TRUE)

## Simulate the epidemic model
sim <- netsim(est, param, init, control)

## Print the results
print(sim)

## Plot the results
plot(sim)

