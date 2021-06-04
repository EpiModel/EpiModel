
num=20
toy <- network_initialize(n=num, directed=F)

ft <- netest(toy,  formation = ~edges ,
             target.stats = c(100),
             coef.diss = dissolution_coefs(dissolution = ~offset(edges), duration = 1), edapprox = T)

nws <- list()
for (i in 1:10) {
  nws[[i]] <- simulate(ft$fit)
}

nwd <- networkDynamic(network.list = nws[1], onsets = 1, termini = 1, end = 10)
as.data.frame(nwd)



param <- param.net(inf.prob = 0.2, act.rate = 1, rec.rate = 1/3)
init <- init.net(i.num = 10, r.num = 0)
control <- control.net(type = "SIR", nsteps = 20, nsims = 1, ncores = 1,  tergmLite = F,  resimulate.network = T)

simz <- netsim(ft, param, init, control)

dat <- as.data.frame(simz$network$sim1)
table(dat$duration)  # now all = 1
