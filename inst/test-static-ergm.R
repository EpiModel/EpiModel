
## Extended tests for reimplementation of network simulations in EpiModel v2.1.0

num = 20
toy <- network_initialize(num)
param <- param.net(inf.prob = 0.2, act.rate = 1, rec.rate = 1/3)
init <- init.net(i.num = 10, r.num = 0)

# TERGM resimulate.network = TRUE tergmLite = FALSE
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 10),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = FALSE tergmLite = FALSE
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 10),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = FALSE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# ERGM resimulate.network = TRUE tergmLite = FALSE
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 1),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# ERGM resimulate.network = FALSE tergmLite = FALSE
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 1),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = FALSE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = TRUE tergmLite = TRUE
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 10),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = TRUE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
# nwd <- as.data.frame(net)
# head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = FALSE tergmLite = TRUE
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 10),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = TRUE,  resimulate.network = FALSE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# ERGM resimulate.network = TRUE tergmLite = TRUE
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 1),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = TRUE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# ERGM resimulate.network = FALSE tergmLite = TRUE (resets to resim = TRUE)
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 1),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = TRUE,  resimulate.network = FALSE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)


# TERGM resimulate.network = TRUE tergmLite = FALSE full STERGM = TRUE
ft <- netest(toy,  formation = ~edges ,
             target.stats = 100,
             coef.diss = dissolution_coefs(dissolution = ~offset(edges),
                                           duration = 10),
             edapprox = FALSE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = FALSE tergmLite = FALSE full STERGM = TRUE
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = FALSE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = TRUE tergmLite = TRUE full STERGM = TRUE
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = TRUE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = FALSE tergmLite = TRUE full STERGM = TRUE  (resets to resim = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = TRUE,  resimulate.network = FALSE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = FALSE tergmLite = FALSE full STERGM = FALSE multi dur = TRUE
num = 20
toy <- network_initialize(num)
toy <- set_vertex_attribute(toy, "age.grp", rbinom(num, 1, 0.5))
get_vertex_attribute(toy, "age.grp")

ft <- netest(toy,  formation = ~edges + nodematch("age.grp"),
             target.stats = c(100, 50),
             coef.diss = dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("age.grp")),
                                           duration = c(10, 20)),
             edapprox = TRUE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = FALSE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = TRUE tergmLite = FALSE full STERGM = FALSE multi dur = TRUE
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = TRUE tergmLite = TRUE full STERGM = FALSE multi dur = TRUE
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = TRUE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)

# TERGM resimulate.network = TRUE tergmLite = FALSE full STERGM = TRUE multi dur = TRUE
ft <- netest(toy,  formation = ~edges + nodematch("age.grp"),
             target.stats = c(100, 50),
             coef.diss = dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("age.grp")),
                                           duration = c(10, 20)),
             edapprox = FALSE)
control <- control.net(type = "SIR", nsteps = 10, nsims = 2, ncores = 1,
                       tergmLite = FALSE,  resimulate.network = TRUE, verbose = FALSE)
simz <- netsim(ft, param, init, control)
net <- get_network(simz, sim = 2)
nwd <- as.data.frame(net)
head(nwd, 10)
get_nwstats(simz)
