devtools::load_all()

nw <- network_initialize(n = 50)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(nw, formation = ~edges + nodematch("race"),
              target.stats = c(25, 10),
              coef.diss = dissolution_coefs(~offset(edges), 10, 0),
              verbose = FALSE)

param <- param.net(inf.prob = 0.3, act.rate = 0.5)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
mod <- netsim(est, param, init, control)
expect_is(mod, "netsim")
plot(mod)
plot(mod, type = "formation")
plot(mod, type = "network")
test_net(mod)

param_random <- function(n, values, weights = NULL) {
  f <- function() {
    return(sample(x = values, size = n, prob = weights, replace = TRUE))
  }
  return(f)
}

param <- param.net(
  inf.prob = 0.3,
  act.rate = param_random(1, c(0.25, 0.5, 0.75)),
  tx.halt.part.prob = function() rbeta(1, 1, 2)
)

generate_rng_param <- function(param) {
  rng_names <- names(Filter(is.function, param))
  if (!is.null(rng_names)) {
    param[rng_names] <- lapply(param[rng_names], , do.call, args = list())
  }

  return(param)
}

generate_rng_param(param)

param <- param.net(inf.prob = 0.3, act.rate = 0.5)
param <- param.net(
  inf.prob = 0.3,
  .rng_params = list(
    act.rate = param_random(1, c(0.25, 0.5, 0.75)),
    tx.halt.part.prob = function() rbeta(1, 1, 2)
  )
)

generate_rng_param <- function(param) {
  rng_names <- names(param$.rng_params)
  param[rng_names] <- lapply(param$.rng_params, do.call, args = list())

  return(param)
}

generate_rng_param(param)

