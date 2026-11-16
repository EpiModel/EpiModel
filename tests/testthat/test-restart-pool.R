# Copy the content of directory `from` into `to`, removing `to` first
copy_dir <- function(from, to) {
  if (dir.exists(to)) unlink(to, recursive = TRUE)
  dir.create(to, recursive = TRUE)
  file.copy(
    list.files(from, full.names = TRUE, all.files = TRUE, no.. = TRUE),
    to,
    recursive = TRUE,
    copy.date = TRUE
  )
  to
}

# Toy model: 2 batch files of 3 SI simulations, as written by
# `netsim_scenarios` for a scenario named "base"
make_toy_batches <- function(scen_dir) {
  set.seed(1)
  nw <- network_initialize(n = 100)
  est <- netest(nw,
    formation = ~edges, target.stats = 50,
    coef.diss = dissolution_coefs(~offset(edges), 10, 0),
    verbose = FALSE
  )
  control <- control.net(
    type = "SI", nsims = 3, nsteps = 5,
    tergmLite = TRUE, resimulate.network = TRUE, save.run = TRUE,
    verbose = FALSE
  )
  dir.create(scen_dir, recursive = TRUE)
  for (b in 1:2) {
    sim <- netsim(est, toy_param, toy_init, control)
    saveRDS(sim, file.path(scen_dir, paste0("sim__base__", b, ".rds")))
  }
}

toy_param <- param.net(inf.prob = 0.3, act.rate = 1)
toy_init <- init.net(i.num = 10)
toy_control_restart <- control.net(
  type = "SI", start = 2, nsteps = 4,
  tergmLite = TRUE, resimulate.network = TRUE, save.run = TRUE,
  verbose = FALSE
)

# Copy a pool and tag each element with its pool number in an extra epi
# column, copied as is into time step 1 of the simulations restarted from it
tag_pool <- function(pool_dir, tagged_dir) {
  copy_dir(pool_dir, tagged_dir)
  for (f in list.files(tagged_dir, pattern = "^[0-9]+\\.rds$", full.names = TRUE)) {
    x <- readRDS(f)
    x$epi$pool_tag <- data.frame(sim1 = as.numeric(sub("\\.rds$", "", basename(f))))
    saveRDS(x, f)
  }
  tagged_dir
}

test_that("make_restart_pool makes pool element i from row i of keep_sims", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_make")
  pool_dir <- file.path(tempdir(), "rp_pool_make")
  on.exit(unlink(c(scen_dir, pool_dir), recursive = TRUE))
  make_toy_batches(scen_dir)

  keep_sims <- data.frame(batch_number = c(2, 1, 1, 2), sim_number = c(3, 1, 3, 1))
  pool_index <- make_restart_pool(pool_dir, scen_dir, "base", keep_sims, "infTime")

  expect_setequal(list.files(pool_dir), c(paste0(1:4, ".rds"), "pool_index.csv"))
  pool_index_csv <- utils::read.csv(file.path(pool_dir, "pool_index.csv"))
  expect_equal(pool_index_csv$pool_num, 1:4)
  expect_equal(pool_index_csv$batch_number, keep_sims$batch_number)
  expect_equal(pool_index_csv$sim_number, keep_sims$sim_number)
  expect_equal(
    pool_index_csv$epimodel_version,
    rep(as.character(utils::packageVersion("EpiModel")), 4)
  )
  expect_equal(pool_index$batch_number, keep_sims$batch_number)

  for (i in seq_len(nrow(keep_sims))) {
    src <- readRDS(file.path(scen_dir, paste0("sim__base__", keep_sims$batch_number[i], ".rds")))
    expected <- make_restart_point(src, "infTime", sim_num = keep_sims$sim_number[i])
    elt <- readRDS(file.path(pool_dir, paste0(i, ".rds")))
    expect_identical(elt$run[[1]]$attr, expected$run[[1]]$attr)
    expect_equal(elt$epi, expected$epi)
    expect_equal(elt$control$nsims, 1)
  }

  expect_error(
    make_restart_pool(pool_dir, scen_dir, "base", keep_sims, "infTime"),
    "already exists"
  )
})

test_that("make_restart_pool rejects bad keep_sims without leaving a pool", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_bad")
  pool_dir <- file.path(tempdir(), "rp_pool_bad")
  on.exit(unlink(c(scen_dir, pool_dir), recursive = TRUE))
  make_toy_batches(scen_dir)

  expect_error(
    make_restart_pool(pool_dir, scen_dir, "base",
                      data.frame(batch_number = 3, sim_number = 1), "infTime"),
    "No simulation file"
  )
  expect_error(
    make_restart_pool(pool_dir, scen_dir, "other",
                      data.frame(batch_number = 1, sim_number = 1), "infTime"),
    "No simulation file"
  )
  expect_error(
    make_restart_pool(pool_dir, scen_dir, "base",
                      data.frame(batch_number = 1, sim_number = 4), "infTime"),
    "out of range"
  )
  expect_error(
    make_restart_pool(pool_dir, scen_dir, "base",
                      data.frame(batch_number = c(1, 1), sim_number = c(2, 2)), "infTime"),
    "more than once"
  )
  expect_error(
    make_restart_pool(pool_dir, scen_dir, "base",
                      data.frame(batch_number = 1, sim_number = 1), "not_an_attr"),
    "not present"
  )
  expect_false(dir.exists(pool_dir))
})

test_that("validate_restart_pool checks pools and single files", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_valid")
  pool_dir <- file.path(tempdir(), "rp_pool_valid")
  copy_dir_path <- file.path(tempdir(), "rp_pool_valid_copy")
  on.exit(unlink(c(scen_dir, pool_dir, copy_dir_path), recursive = TRUE))
  make_toy_batches(scen_dir)
  keep_sims <- data.frame(batch_number = c(1, 2, 2), sim_number = c(2, 1, 3))
  make_restart_pool(pool_dir, scen_dir, "base", keep_sims, "infTime")
  epimodel_ver <- utils::packageVersion("EpiModel")

  # each case below starts from a fresh copy of the pool
  reset_copy <- function() copy_dir(pool_dir, copy_dir_path)
  edit_copy_index <- function(f) {
    index_path <- file.path(copy_dir_path, "pool_index.csv")
    utils::write.csv(f(utils::read.csv(index_path)), index_path, row.names = FALSE)
  }

  # a valid pool
  expect_warning(infos <- validate_restart_pool(pool_dir), NA)
  expect_equal(infos$pool_num, 1:3)
  expect_equal(infos$n_runs, rep(1, 3))
  expect_equal(infos$nsteps, rep(1, 3))
  expect_equal(infos$epimodel_version, rep(as.character(epimodel_ver), 3))

  # a single restart point file
  expect_warning(validate_restart_pool(file.path(pool_dir, "1.rds")), NA)
  # a single file with several simulations
  expect_warning(
    validate_restart_pool(file.path(scen_dir, "sim__base__1.rds")),
    "holds 3 simulations"
  )
  expect_error(validate_restart_pool(file.path(tempdir(), "no_such_pool")), "neither")

  # a directory without `pool_index.csv` is not a restart pool
  reset_copy()
  file.remove(file.path(copy_dir_path, "pool_index.csv"))
  expect_error(validate_restart_pool(copy_dir_path), "not a restart pool")

  # an index without the version column, or not numbering the elements 1 to N
  reset_copy()
  edit_copy_index(function(idx) {
    idx$epimodel_version <- NULL
    idx
  })
  expect_error(validate_restart_pool(copy_dir_path), "no `epimodel_version` column")
  reset_copy()
  edit_copy_index(function(idx) idx[-2, ])
  expect_error(validate_restart_pool(copy_dir_path), "`1` to `N`")

  # an element listed in the index without its file
  reset_copy()
  file.remove(file.path(copy_dir_path, "2.rds"))
  expect_error(validate_restart_pool(copy_dir_path), "Missing restart pool element.*2\\.rds")

  # an element file not listed in the index
  reset_copy()
  file.copy(file.path(pool_dir, "1.rds"), file.path(copy_dir_path, "4.rds"))
  expect_warning(validate_restart_pool(copy_dir_path), "not listed.*4\\.rds")

  # an element holding several simulations
  reset_copy()
  file.copy(
    file.path(scen_dir, "sim__base__1.rds"), file.path(copy_dir_path, "2.rds"),
    overwrite = TRUE
  )
  expect_error(validate_restart_pool(copy_dir_path), "holds 3 simulations")

  # an element without `control`
  x <- readRDS(file.path(pool_dir, "1.rds"))
  reset_copy()
  x_bad <- x
  x_bad$control <- NULL
  saveRDS(x_bad, file.path(copy_dir_path, "2.rds"))
  expect_error(validate_restart_pool(copy_dir_path), "restart element\\(s\\): control")

  # another major.minor EpiModel version, or none recorded (before 2.7), in
  # the index (checked by `netsim_path_wrapper`) or in the element itself
  # (checked by `netsim`). Only the patch version may differ
  patch_ver <- paste0(epimodel_ver[, 1:2], ".99")
  reset_copy()
  edit_copy_index(function(idx) {
    idx$epimodel_version[2] <- "2.6.0"
    idx
  })
  expect_error(
    validate_restart_pool(copy_dir_path),
    "Element 2 of the restart pool .* made with EpiModel v2.6.0"
  )
  expect_error(validate_restart_pool(copy_dir_path), "Please rebuild the restart object")
  edit_copy_index(function(idx) {
    idx$epimodel_version[2] <- NA
    idx
  })
  expect_error(
    validate_restart_pool(copy_dir_path),
    "Element 2 .* made with EpiModel < 2.7"
  )
  edit_copy_index(function(idx) {
    idx$epimodel_version[2] <- patch_ver
    idx
  })
  expect_warning(validate_restart_pool(copy_dir_path), NA)

  x_bad <- x
  attr(x_bad, "version") <- numeric_version("2.6.0")
  saveRDS(x_bad, file.path(copy_dir_path, "2.rds"))
  expect_error(validate_restart_pool(copy_dir_path), "2\\.rds\" was made with EpiModel v2.6.0")
  attr(x_bad, "version") <- NULL
  saveRDS(x_bad, file.path(copy_dir_path, "2.rds"))
  expect_error(validate_restart_pool(copy_dir_path), "2\\.rds\" was made with EpiModel < 2.7")
  attr(x_bad, "version") <- numeric_version(patch_ver)
  saveRDS(x_bad, file.path(copy_dir_path, "2.rds"))
  expect_warning(validate_restart_pool(copy_dir_path), NA)

  # elements holding different numbers of time steps
  reset_copy()
  x_bad <- x
  x_bad$control$nsteps <- 5
  saveRDS(x_bad, file.path(copy_dir_path, "2.rds"))
  expect_error(
    validate_restart_pool(copy_dir_path),
    "1 step\\(s\\) for element\\(s\\) 1, 3; 5 step\\(s\\) for element\\(s\\) 2"
  )
})

test_that("internal_run_batch restarts simulation k from pool element (k - 1) %% N + 1", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_run")
  pool_dir <- file.path(tempdir(), "rp_pool_run")
  tagged_dir <- file.path(tempdir(), "rp_pool_run_tagged")
  out_dir <- file.path(tempdir(), "rp_out_run")
  on.exit(unlink(c(scen_dir, pool_dir, tagged_dir, out_dir), recursive = TRUE))
  make_toy_batches(scen_dir)
  keep_sims <- data.frame(batch_number = c(2, 1, 2, 1), sim_number = c(1, 3, 2, 2))
  make_restart_pool(pool_dir, scen_dir, "base", keep_sims, "infTime")
  tag_pool(pool_dir, tagged_dir)

  scenario <- create_scenario_list(
    data.frame(.at = 0, .scenario.id = "base")
  )[[1]]
  # 6 replications over 3 batches of 2: sims 1 to 6 use elements 1 2 3 4 1 2
  for (b in 1:3) {
    internal_run_batch(
      scenario = scenario,
      batch_num = b,
      path_to_x = tagged_dir,
      param = toy_param,
      init = toy_init,
      control = toy_control_restart,
      libraries = NULL,
      output_dir = out_dir,
      n_batch = 3,
      n_rep = 6,
      n_cores = 2
    )
  }
  expect_setequal(
    list.files(out_dir),
    paste0("sim__base__", 1:3, ".rds")
  )

  pool_epi <- lapply(1:4, function(i) readRDS(file.path(pool_dir, paste0(i, ".rds")))$epi)
  tags <- c()
  for (b in 1:3) {
    sim <- readRDS(file.path(out_dir, paste0("sim__base__", b, ".rds")))
    expect_equal(sim$control$nsims, 2)
    expect_equal(nrow(sim$epi$i.num), 4)
    for (j in 1:2) {
      k <- (b - 1) * 2 + j
      expected_elt <- (k - 1) %% 4 + 1
      tags <- c(tags, sim$epi$pool_tag[1, j])
      expect_equal(sim$epi$i.num[1, j], pool_epi[[expected_elt]]$i.num[1, 1])
      expect_equal(sim$epi$num[1, j], pool_epi[[expected_elt]]$num[1, 1])
    }
  }
  expect_equal(tags, c(1, 2, 3, 4, 1, 2))
})

test_that("internal_run_batch gives each simulation its own checkpoint directory", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_ckpt")
  pool_dir <- file.path(tempdir(), "rp_pool_ckpt")
  out_dir <- file.path(tempdir(), "rp_out_ckpt")
  ckpt_dir <- file.path(tempdir(), "rp_ckpt")
  on.exit(unlink(c(scen_dir, pool_dir, out_dir, ckpt_dir), recursive = TRUE))
  make_toy_batches(scen_dir)
  make_restart_pool(
    pool_dir, scen_dir, "base",
    data.frame(batch_number = 1, sim_number = 1:2), "infTime"
  )

  control <- toy_control_restart
  control$.checkpoint.dir <- ckpt_dir
  control$.checkpoint.steps <- 1
  control$.checkpoint.keep <- TRUE
  scenario <- create_scenario_list(
    data.frame(.at = 0, .scenario.id = "base")
  )[[1]]
  internal_run_batch(
    scenario = scenario,
    batch_num = 1,
    path_to_x = pool_dir,
    param = toy_param,
    init = toy_init,
    control = control,
    libraries = NULL,
    output_dir = out_dir,
    n_batch = 1,
    n_rep = 2,
    n_cores = 2
  )

  sim <- readRDS(file.path(out_dir, "sim__base__1.rds"))
  batch_ckpt_dir <- paste0(ckpt_dir, "/sim__base__1")
  expect_equal(sim$control$nsims, 2)
  expect_equal(sim$control$.checkpoint.dir, batch_ckpt_dir)
  expect_setequal(
    list.files(batch_ckpt_dir),
    c("sim_1", "sim_2")
  )

  # without `.checkpoint.keep`, the batch directory goes once the batch is saved
  unlink(c(out_dir, ckpt_dir), recursive = TRUE)
  control$.checkpoint.keep <- FALSE
  internal_run_batch(
    scenario = scenario,
    batch_num = 1,
    path_to_x = pool_dir,
    param = toy_param,
    init = toy_init,
    control = control,
    libraries = NULL,
    output_dir = out_dir,
    n_batch = 1,
    n_rep = 2,
    n_cores = 2
  )
  expect_true(file.exists(file.path(out_dir, "sim__base__1.rds")))
  expect_false(dir.exists(batch_ckpt_dir))
})

test_that("a failure while writing the pool removes the partial directory", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_partial")
  pool_dir <- file.path(tempdir(), "rp_pool_partial")
  on.exit(unlink(c(scen_dir, pool_dir), recursive = TRUE))
  make_toy_batches(scen_dir)

  # batch 1 is written before batch 2 fails, so this exercises the cleanup of
  # an already partially written pool
  keep_sims <- data.frame(batch_number = c(1, 2), sim_number = c(1, 99))
  expect_error(
    make_restart_pool(pool_dir, scen_dir, "base", keep_sims, "infTime"),
    "out of range"
  )
  expect_false(dir.exists(pool_dir))
  expect_false(dir.exists(get_partial_pool_dir(pool_dir)))
})

test_that("netsim_path_wrapper checks versions from the pool index, netsim checks the elements", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_wrapper")
  pool_dir <- file.path(tempdir(), "rp_pool_wrapper")
  on.exit(unlink(c(scen_dir, pool_dir), recursive = TRUE))
  make_toy_batches(scen_dir)
  make_restart_pool(
    pool_dir, scen_dir, "base",
    data.frame(batch_number = 1, sim_number = 1:3), "infTime"
  )
  index_path <- file.path(pool_dir, "pool_index.csv")
  pool_index <- utils::read.csv(index_path)
  run <- function(sim_nums) {
    netsim_path_wrapper(
      path_to_x = pool_dir,
      param = toy_param,
      init = toy_init,
      control = toy_control_restart,
      sim_nums = sim_nums
    )
  }
  with(future::plan("sequential"), local = TRUE)

  # first check: the index records element 2 as made with another EpiModel
  # version. Only the elements used are checked
  bad_index <- pool_index
  bad_index$epimodel_version[2] <- "2.6.0"
  utils::write.csv(bad_index, index_path, row.names = FALSE)
  expect_error(run(1:2), "Element 2 of the restart pool .* made with EpiModel v2.6.0")
  expect_length(run(c(1, 3)), 2)
  utils::write.csv(pool_index, index_path, row.names = FALSE)

  # second check, by `netsim`: the element itself was made with another version
  x <- readRDS(file.path(pool_dir, "2.rds"))
  x_bad <- x
  attr(x_bad, "version") <- numeric_version("2.6.0")
  saveRDS(x_bad, file.path(pool_dir, "2.rds"))
  # `future.apply` warns that it cancels the other iterations
  expect_error(suppressWarnings(run(2)), "Please rebuild the restart object")
  saveRDS(x, file.path(pool_dir, "2.rds"))

  # the pool size comes from the index (3 elements), not from the files:
  # simulation 4 uses element 1 even when element 3 is missing, and the
  # missing element only fails when read
  file.rename(file.path(pool_dir, "3.rds"), file.path(pool_dir, "3.rds.bak"))
  expect_equal(
    run(4)[[1]]$restart_point$sim1,
    as.character(file.path(pool_dir, "1.rds"))
  )
  expect_error(suppressWarnings(run(3)), "cannot open")

  # a directory without `pool_index.csv` is not a restart pool
  file.remove(index_path)
  expect_error(suppressWarnings(run(1)), "not a restart pool")
})

test_that("the scenario runners reject raw.output = TRUE at setup", {
  control <- control.net(
    type = "SI", start = 2, nsteps = 4, tergmLite = TRUE,
    resimulate.network = TRUE, raw.output = TRUE, verbose = FALSE
  )

  # `internal_scenarios_setup` is called by `netsim_scenarios` and by the
  # EpiModelHPC `slurmworkflow` step templates
  expect_error(
    internal_scenarios_setup(
      path_to_x = file.path(tempdir(), "rp_pool_raw"),
      param = toy_param,
      init = toy_init,
      control = control,
      scenarios_list = NULL,
      n_rep = 1,
      n_cores = 1,
      output_dir = file.path(tempdir(), "rp_out_raw"),
      libraries = NULL
    ),
    "raw.output = TRUE"
  )
})

test_that("netsim_path_wrapper turns verbose off for several simulations", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_verbose")
  on.exit(unlink(scen_dir, recursive = TRUE))
  make_toy_batches(scen_dir)
  x_path <- file.path(scen_dir, "restart.rds")
  saveRDS(
    make_restart_point(
      readRDS(file.path(scen_dir, "sim__base__1.rds")),
      time_attrs = "infTime", sim_num = 1
    ),
    x_path
  )
  control <- toy_control_restart
  control$verbose <- TRUE

  with(future::plan("sequential"), local = TRUE)
  expect_output(
    netsim_path_wrapper(
      path_to_x = x_path,
      param = toy_param,
      init = toy_init,
      control = control,
      sim_nums = 1
    ),
    "Timestep"
  )
  expect_silent(
    netsim_path_wrapper(
      path_to_x = x_path,
      param = toy_param,
      init = toy_init,
      control = control,
      sim_nums = 1:2
    )
  )
})

test_that("the partial pool left by a killed make_restart_pool is removed", {
  skip_on_cran()
  scen_dir <- file.path(tempdir(), "rp_scen_killed")
  pool_dir <- file.path(tempdir(), "rp_pool_killed")
  partial_dir <- get_partial_pool_dir(pool_dir)
  on.exit(unlink(c(scen_dir, pool_dir), recursive = TRUE))
  make_toy_batches(scen_dir)

  # what a `make_restart_pool` killed while writing leaves behind
  expect_equal(basename(partial_dir), ".rp_pool_killed.partial")
  dir.create(partial_dir, recursive = TRUE)
  saveRDS("stale", file.path(partial_dir, "1.rds"))
  saveRDS("stale", file.path(partial_dir, "3.rds"))

  expect_message(
    make_restart_pool(
      pool_dir, scen_dir, "base",
      data.frame(batch_number = 1, sim_number = 1:2), "infTime"
    ),
    "partial restart pool"
  )
  expect_false(dir.exists(partial_dir))
  expect_setequal(
    list.files(pool_dir, all.files = TRUE, no.. = TRUE),
    c("1.rds", "2.rds", "pool_index.csv")
  )
  expect_s3_class(readRDS(file.path(pool_dir, "1.rds")), "netsim")
})
