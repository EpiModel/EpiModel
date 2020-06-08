
##
## Tutorial: Dynamic Network Visualization
## Day 4 | Network Modeling for Epidemics
##

library("EpiModel")
library("ndtv")

nw <- network.initialize(n = 100, directed = FALSE)
nw <- set.vertex.attribute(nw, "race", rbinom(100, 1, 0.5))
nw <- set.vertex.attribute(nw, "age", sample(18:50, 100, TRUE))

formation <- ~edges + nodematch("race") + nodefactor("race") +
  absdiff("age") + concurrent
target.stats <- c(40, 30, 20, 80, 20)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges) +
                                 offset(nodematch("race")),
                               duration = c(20, 10))
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

param <- param.net(inf.prob = 1)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 25, nsims = 1, verbose = FALSE)
sim <- netsim(est, param, init, control)

nw <- get_network(sim)

nw <- color_tea(nw, verbose = FALSE)

slice.par <- list(start = 1, end = 25, interval = 1,
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))

compute.animation(nw, slice.par = slice.par, verbose = TRUE)

race <- get.vertex.attribute(nw, "race")
race.shape <- ifelse(race == 1, 4, 50)

age <- get.vertex.attribute(nw, "age")
age.size <- age/25

render.d3movie(
    nw,
    render.par = render.par,
    plot.par = plot.par,
    vertex.cex = age.size,
    vertex.sides = race.shape,
    vertex.col = "ndtvcol",
    edge.col = "darkgrey",
    vertex.border = "lightgrey",
    displaylabels = FALSE,
    filename = paste0(getwd(), "/movie.html"))
