n <- 1e2
nw <- network_initialize(n = n)
formation <- ~edges
target.stats <- 0.4 * n
coef.diss <- dissolution_coefs(dissolution = ~ offset(edges), duration = 10)
est <- netest(nw, formation, target.stats, coef.diss)

# 2. Epidemic simulation
param <- param.net(inf.prob = 1)
init <- init.net(i.num = 0.1 * n)
control <- control.net(
  type = "SI",
  nsteps = 100,
  nsims = 1,
  tergmLite = TRUE,
  tracked.attributes = c('active', 'status'),
  cumulative.edgelist = TRUE,
  truncate.el.cuml = Inf,
  save.cumulative.edgelist = TRUE,
  save.run = TRUE,
  verbose = FALSE
)
sim <- netsim(est, param, init, control)

nw <- EpiModel:::make_networkDynamic(sim)
nw <- color_tea(nw, old.var = "status", verbose = FALSE)

# 4. Compute layout + render
slice.par <- list(
  start = 1,
  end = 25,
  interval = 1,
  aggregate.dur = 1,
  rule = "any"
)
compute.animation(nw, slice.par = slice.par)
render.d3movie(
  nw,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html")
)
