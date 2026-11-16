# Toy model: 2 scenarios ("a" and "b") of 2 batch files of 2 SI simulations, as
# written by `netsim_scenarios`
make_toy_scenario_batches <- function(scen_dir) {
  set.seed(1)
  nw <- network_initialize(n = 100)
  est <- netest(nw,
    formation = ~edges, target.stats = 50,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )
  control <- control.net(
    type = "SI", nsims = 2, nsteps = 5,
    tergmLite = TRUE, resimulate.network = TRUE, save.run = TRUE,
    verbose = FALSE
  )
  param <- param.net(inf.prob = 0.3, act.rate = 1)
  init <- init.net(i.num = 10)
  dir.create(scen_dir, recursive = TRUE)
  for (scen in c("a", "b")) {
    for (b in 1:2) {
      sim <- netsim(est, param, init, control)
      saveRDS(sim, file.path(scen_dir, paste0("sim__", scen, "__", b, ".rds")))
    }
  }
}

test_that("netsim_scenarios runs the batches of each scenario in the current session", {
  skip_on_cran()
  output_dir <- file.path(tempdir(), "scen_run_out")
  on.exit(unlink(output_dir, recursive = TRUE))
  dir.create(output_dir, recursive = TRUE)

  set.seed(10)
  nw <- network_initialize(n = 200)
  est <- netest(nw,
    formation = ~edges, target.stats = 60,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )
  param <- param.net(inf.prob = 0.9, rec.rate = 0.01, act.rate = 2)
  control <- control.net(type = "SIS", nsims = 1, nsteps = 2, verbose = FALSE)
  init <- init.net(i.num = 10)
  # run the simulations one after the other rather than in `future` workers
  control$future.use.plan <- future::tweak(future::sequential)

  est_path <- file.path(output_dir, "est.rds")
  saveRDS(est, est_path)

  scenarios.df <- dplyr::tribble(
    ~.scenario.id, ~.at, ~inf.prob, ~rec.rate,
    "base", 0, 0.9, 0.01,
    "multiple_changes", 0, 0.1, 0.04,
    "multiple_changes", 20, 0.9, 0.01,
    "multiple_changes", 40, 0.1, 0.1
  )
  scenarios.list <- create_scenario_list(scenarios.df)

  n_rep <- 3
  n_cores <- 2
  n_scen <- length(scenarios.list)
  res <- withVisible(
    netsim_scenarios(
      path_to_x = est_path,
      param, init, control,
      scenarios_list = scenarios.list,
      n_rep = n_rep, n_cores = n_cores,
      output_dir = output_dir,
      libraries = NULL
    )
  )
  expect_false(res$visible)
  expect_null(res$value)

  expect_s3_class(readRDS(file.path(output_dir, "sim__base__1.rds")), "netsim")
  expect_length(
    list.files(output_dir),
    n_scen * ceiling(n_rep / n_cores) + 1 # +1 for the est file
  )
})

test_that("get_scenarios_batches_infos sorts by scenario and batch number", {
  scen_dir <- file.path(tempdir(), "batches_order")
  on.exit(unlink(scen_dir, recursive = TRUE))
  dir.create(scen_dir, recursive = TRUE)
  file.create(file.path(scen_dir, paste0("sim__b__", 1:12, ".rds")))
  file.create(file.path(scen_dir, paste0("sim__a__", c(10, 2, 1), ".rds")))

  infos <- get_scenarios_batches_infos(scen_dir)
  expect_equal(infos$scenario_name, c(rep("a", 3), rep("b", 12)))
  expect_equal(infos$batch_number, c(1L, 2L, 10L, 1:12))
  expect_type(infos$file_path, "character")
  expect_equal(
    basename(infos$file_path)[1:3],
    c("sim__a__1.rds", "sim__a__2.rds", "sim__a__10.rds")
  )
})

test_that("the scenario file listings only match the batch file names", {
  scen_dir <- file.path(tempdir(), "batches_strict")
  on.exit(unlink(scen_dir, recursive = TRUE))
  dir.create(scen_dir, recursive = TRUE)
  dir.create(file.path(scen_dir, "sim__dir__1.rds"))
  dir.create(file.path(scen_dir, "df__dir.rds"))
  file.create(file.path(
    scen_dir,
    c("est.rds", "sim__x.rds", "merged__a.rds", "df__a", "sim__ok__1.rds", "df__ok.rds")
  ))

  b_infos <- get_scenarios_batches_infos(scen_dir)
  expect_equal(b_infos$scenario_name, "ok")
  expect_equal(b_infos$batch_number, 1L)

  d_infos <- get_scenarios_tibble_infos(scen_dir)
  expect_equal(d_infos$scenario_name, "ok")
  expect_equal(basename(d_infos$file_path), "df__ok.rds")
})

test_that("get_batch_plan keeps a tweaked plan, else tweaks the default", {
  user_plan <- future::tweak(future::multisession, workers = 3)
  expect_identical(
    get_batch_plan(list(future.use.plan = user_plan), future::multicore, 2),
    user_plan
  )

  default_plan <- get_batch_plan(
    list(future.use.plan = FALSE), future::multicore, 2
  )
  expect_s3_class(default_plan, "tweaked")
  expect_s3_class(default_plan, "multicore")
  expect_equal(attr(default_plan, "tweaks")$workers, 2)

  # the tweak holds the value of `n_cores`, not a promise on the caller
  n <- 4
  default_plan <- get_batch_plan(list(), future::multisession, n)
  n <- 1
  expect_equal(attr(default_plan, "tweaks")$workers, 4)
})

test_that("batch files keep the cumulative edgelists of all their simulations", {
  skip_on_cran()
  out_dir <- file.path(tempdir(), "scen_cel_out")
  est_path <- file.path(tempdir(), "scen_cel_est.rds")
  on.exit(unlink(c(out_dir, est_path), recursive = TRUE))

  set.seed(1)
  nw <- network_initialize(n = 100)
  est <- netest(nw,
    formation = ~edges, target.stats = 50,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )
  saveRDS(est, est_path)
  control <- control.net(
    type = "SI", nsteps = 3, tergmLite = TRUE, resimulate.network = TRUE,
    cumulative.edgelist = TRUE, save.cumulative.edgelist = TRUE,
    verbose = FALSE
  )
  scenario <- create_scenario_list(
    data.frame(.at = 0, .scenario.id = "base")
  )[[1]]

  # 3 replications over 2 batches: a batch of 2 merged sims, and one of 1
  for (b in 1:2) {
    netsim_run_one_scenario(
      scenario, b, est_path,
      param.net(inf.prob = 0.3, act.rate = 1), init.net(i.num = 10), control,
      libraries = NULL, output_dir = out_dir,
      n_batch = 2, n_rep = 3, n_cores = 2
    )
  }
  batches <- lapply(1:2, function(b) {
    readRDS(file.path(out_dir, paste0("sim__base__", b, ".rds")))
  })
  expect_length(batches[[1]]$cumulative.edgelist, 2)
  expect_length(batches[[2]]$cumulative.edgelist, 1)

  merged <- merge(batches[[1]], batches[[2]], keep.cumulative.edgelist = TRUE)
  expect_length(merged$cumulative.edgelist, 3)
})

test_that("merge_netsim_scenarios merges the batches of each scenario once", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "scen_merge_sim")
  merged_dir <- file.path(tempdir(), "scen_merge_out")
  on.exit(unlink(c(scen_dir, merged_dir), recursive = TRUE))
  make_toy_scenario_batches(scen_dir)

  with(future::plan("sequential"), local = TRUE)
  merge_netsim_scenarios(scen_dir, merged_dir)

  expect_setequal(list.files(merged_dir), c("merged__a.rds", "merged__b.rds"))
  merged <- readRDS(file.path(merged_dir, "merged__a.rds"))
  expect_s3_class(merged, "netsim")
  expect_equal(merged$control$nsims, 4)
  expect_equal(ncol(merged$epi$i.num), 4)
})

test_that("merge_netsim_scenarios_tibble writes one tibble per scenario", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "scen_tibble_sim")
  df_dir <- file.path(tempdir(), "scen_tibble_out")
  df_sel_dir <- file.path(tempdir(), "scen_tibble_sel_out")
  on.exit(unlink(c(scen_dir, df_dir, df_sel_dir), recursive = TRUE))
  make_toy_scenario_batches(scen_dir)

  with(future::plan("sequential"), local = TRUE)
  merge_netsim_scenarios_tibble(scen_dir, df_dir, steps_to_keep = 1)

  expect_setequal(list.files(df_dir), c("df__a.rds", "df__b.rds"))
  d <- readRDS(file.path(df_dir, "df__a.rds"))
  expect_s3_class(d, "tbl_df")
  expect_equal(names(d)[1:3], c("batch_number", "sim", "time"))
  # 2 batches x 2 sims x 2 time steps kept
  expect_equal(nrow(d), 8)
  expect_equal(sort(unique(d$batch_number)), 1:2)
  # `sim` is renumbered over the batches, `sim_number` keeps the in-batch number
  expect_equal(sort(unique(d$sim)), 1:4)
  expect_equal(sort(unique(d$sim_number)), 1:2)
  expect_true(all(c("i.num", "num") %in% names(d)))

  # `cols` selects among the columns beyond `batch_number`, `sim` and `time`,
  # and can rename them
  merge_netsim_scenarios_tibble(
    scen_dir, df_sel_dir,
    steps_to_keep = 1, cols = c(cases = "i.num")
  )
  d_sel <- readRDS(file.path(df_sel_dir, "df__a.rds"))
  expect_equal(
    names(d_sel),
    c("batch_number", "sim", "time", "cases", "sim_number")
  )

  infos <- get_scenarios_tibble_infos(df_dir)
  expect_setequal(infos$scenario_name, c("a", "b"))
  expect_type(infos$file_path, "character")
  expect_true(all(file.exists(infos$file_path)))
})
