
# Setup
n <- 100
nw <- network_initialize(n)

prev <- 0.2
infIds <- sample(1:n, n*prev)
nw <- set_vertex_attribute(nw, "status", "s")
nw <- set_vertex_attribute(nw, "status", "i", infIds)

mean.deg.inf <- 0.3
inedges.inf <- mean.deg.inf * n * prev

mean.deg.sus <- 0.8
inedges.sus <- mean.deg.sus * n * (1 - prev)

edges <- (inedges.inf + inedges.sus)/2
nmatch <- edges * 0.91

formation <- ~edges + nodefactor("status") + nodematch("status")
target.stats <- c(edges, inedges.sus, nmatch)

coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 50)

suppressMessages(
  est <- netest(nw, formation, target.stats, coef.diss)
)

param <- param.net(inf.prob = 0.03)
init <- init.net()

test_that("nw stats in full tergm", {
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = FALSE,
                         nwstats.formula = ~edges +
                           meandeg +
                           nodefactor("status", levels = NULL) +
                           nodematch("status"), verbose = FALSE)
  sim <- netsim(est, param, init, control)

  expect_s3_class(sim, "netsim")
  nws <- get_nwstats(sim)
  expect_equal(names(nws), c("time", "sim", "edges", "meandeg",
                             "nodefactor.status.i", "nodefactor.status.s",
                             "nodematch.status"))
  expect_equal(dim(nws), c(5, 7))

  # default = 'formation'
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = FALSE,
                         verbose = FALSE)
  sim <- netsim(est, param, init, control)

  expect_s3_class(sim, "netsim")
  nws <- get_nwstats(sim)

  expect_equal(names(nws), c("time", "sim", "edges",
                             "nodefactor.status.s", "nodematch.status"))
  expect_equal(dim(nws), c(5, 5))

})

test_that("nw stats in tergmLite", {
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = TRUE,
                         nwstats.formula = ~edges +
                           meandeg +
                           nodefactor("status", levels = NULL) +
                           nodematch("status"), verbose = FALSE)
  sim <- netsim(est, param, init, control)

  expect_s3_class(sim, "netsim")
  nws <- get_nwstats(sim)
  expect_equal(names(nws), c("time", "sim", "edges", "meandeg",
                             "nodefactor.status.i", "nodefactor.status.s",
                             "nodematch.status"))
  expect_equal(dim(nws), c(5, 7))

  # default = 'formation'
  control <- control.net(type = "SI", nsteps = 5, nsims = 1, ncores = 1,
                         resimulate.network = TRUE, tergmLite = FALSE,
                         verbose = FALSE)
  sim <- netsim(est, param, init, control)

  expect_s3_class(sim, "netsim")
  nws <- get_nwstats(sim)

  expect_equal(names(nws), c("time", "sim", "edges",
                             "nodefactor.status.s", "nodematch.status"))
  expect_equal(dim(nws), c(5, 5))
})
