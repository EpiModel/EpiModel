context("Cumulative edgelist (All SOC)")

nw <- network_initialize(n = 100)
est <- netest(
  nw,
  formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~ offset(edges), 10, 0),
  verbose = FALSE
)

param <- param.net(
  inf.prob = 0.3,
  act.rate = 0.1
)

init <- init.net(i.num = 10)

test_that("netsim, SI, Cumulative Edgelist", {
  skip_on_cran()
  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 100,
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    truncate.el.cuml = NULL,
    cumulative.edgelist = TRUE,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)
  d <- get_cumulative_edgelists_df(mod[[1]])

  expect_is(d, "data.frame")
  expect_equal(colnames(d), c("head", "tail", "start", "stop", "network"))

  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 100,
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    truncate.el.cuml = 40,
    cumulative.edgelist = TRUE,
    save.cumulative.edgelist = TRUE,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)
  d <- get_cumulative_edgelists_df(mod[[1]])
  expect_gte(min(d$stop, na.rm = TRUE), 40)

  d <- get_partners(dat = mod[[1]], index_posit_ids = 1:10, networks = 1)
  expect_equal(colnames(d), c("index", "partner", "start", "stop", "network"))

  # expect errors when trying to access networks that don't exist
  expect_error(
    get_partners(dat = mod[[1]], index_posit_ids = 1:10, networks = 2)
  )

  expect_error(
    update_cumulative_edgelist(dat = mod[[1]], network = 2)
  )

  # Cumulative edgelist extraction and reachable set
  sim <- process_out.net(mod)
  el_cuml <- sim$cumulative.edgelist[[1]]
  el_cuml <- dedup_cumulative_edgelist(el_cuml)

  nnodes <- max(el_cuml$head, el_cuml$tail)

  from_step <- min(el_cuml$start)
  to_step <- 100
  nodes <- sample(nnodes, 10)
  # Only test that the functions run without error. See the functions examples
  # to check correctness
  el_tp <- get_forward_reachable(el_cuml, from_step, to_step, nodes)
  el_tp <- get_backward_reachable(el_cuml, from_step, to_step, nodes, "yes")
})

test_that("netsim, SI, Cumulative Edgelist - missing args", {
  skip_on_cran()
  control <- control.net(
    type = "SI",
    nsims = 1,
    nsteps = 10,
    resimulate.network = TRUE,
    tergmLite = TRUE,
    verbose = FALSE,
    raw.output = TRUE
  )

  mod <- netsim(est, param, init, control)

  expect_error(get_cumulative_edgelists_df(mod[[1]]))

  expect_equal(mod[[1]]$control$cumulative.edgelist, FALSE)
  expect_equal(mod[[1]]$control$truncate.el.cuml, 0)

})

test_that("netsim, SI, Cumulative Edgelist with arrivals and departures", {
  skip_on_cran()
  nw <- network_initialize(n = 100)

  formation <- ~edges
  target.stats <- 70
  coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 38, d.rate = 0.01)
  est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  param <- param.net(inf.prob = 0.1, act.rate = 1,
                     a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01)

  init <- init.net(i.num = 10)
  control <- control.net(
    type = "SI", nsteps = 100, nsims = 1,
    tergmLite = TRUE, # required as we want the removal of inactive nodes
    resimulate.network = TRUE, verbose = FALSE,
    truncate.el.cuml = Inf,
    cumulative.edgelist = TRUE,
    raw.output = TRUE
  )

  sim1 <- netsim(est1, param, init, control)

  dat <- sim1[[1]]
  uid <- get_unique_ids(dat)
  pid <- get_posit_ids(dat)
  is_active_posit_ids(dat, pid)

  el_cuml <- get_cumulative_edgelists_df(dat)
  spids <- sample(pid, 20)
  partners <- get_partners(dat, spids)
  degree_cuml <- get_cumulative_degree(dat, spids)

  # checks the unicity of UIDs
  expect_true(all(table(uid) == 1))
  # checks that all attributed UID are not in the currently attributed UIDs
  expect_false(setequal(seq_len(dat$run$last_unique_id), get_unique_ids(dat)))
  # checks that all attributed UID are not in the currently attributed PIDs
  expect_false(setequal(seq_len(dat$run$last_unique_id), get_posit_ids(dat)))

  # checks the translation between UID and PID in get_partners
  expect_true(all(get_posit_ids(dat, unique(partners$index)) %in% spids))

})
