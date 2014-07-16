context("Diagnostics")


# Setup -------------------------------------------------------------------

## Model 1: Edges Only
num <- 50
nw <- network.initialize(num, directed = FALSE)
formation <- ~ edges
target.stats <- 15
dissolution <- ~ offset(edges)
duration <- 20
coef.diss <- dissolution_coefs(dissolution, duration)
est1 <- netest(nw,
               formation,
               dissolution,
               target.stats,
               coef.diss,
               verbose = FALSE)


# Tests -------------------------------------------------------------------

test_that("netdx works for single simulations of edges only model", {
  dx <- netdx(est1, nsims = 1, nsteps = 10, verbose = FALSE)
  expect_is(dx, "netdx")
})

test_that("netdx works for multiple simulations of edges only model", {
  dx <- netdx(est1, nsims = 2, nsteps = 10, verbose = FALSE)
  expect_is(dx, "netdx")
})

