
library("EpiModel")

n <- 1000
nw <- network.initialize(n = n, directed = FALSE)
formation <- ~edges + concurrent
target.stats <- c(1*(n/2), 250)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 10)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

dx <- netdx(est, nsims = 10, nsteps = 500, ncores = 4, verbose = FALSE)
dx

plot(dx)

library("EasyABC")
netest_refit_abc <- function(est, nsims, ncores, nsteps,
                             prior.min = 0, prior.max = 0,
                             p_acc_min = 0.1) {

  est_orig <- est
  save(est_orig, nsteps, file = "temp-refit-abc.rda")

  myfunc <- function(x) {
    set.seed(x[1])
    require(EpiModel)
    load("temp-refit-abc.rda")
    est_temp <- est_orig
    est_temp$coef.form <- est_temp$coef.form + x[2:length(x)]
    dx <- netdx(est_temp, nsims = 1, nsteps = nsteps, verbose = FALSE)
    out <- get_nwstats(dx)
    out <- out[, which(!names(out) %in% c("time", "sim")), drop = FALSE]
    out <- colMeans(out)
    return(out)
  }

  targets <- est_orig$target.stats
  n_targets <- length(targets)
  priors <- list()
  for (ii in seq_len(n_targets)) {
    priors[[ii]] <- c("unif", prior.min, prior.max)
  }

  refit <- ABC_sequential(
    method = "Lenormand",
    model = myfunc,
    prior = priors,
    nb_simul = nsims,
    summary_stat_target = targets,
    p_acc_min = p_acc_min,
    progress_bar = TRUE,
    verbose = FALSE,
    n_cluster = ncores,
    use_seed = TRUE
  )

  if (n_targets == 1) {
    coef.adj <- sum(refit$param * refit$weights)
  } else {
    coef.adj <- rep(NA, length(est$coef.form))
    for (jj in seq_len(ncol(refit$param))) {
      coef.adj[jj] <- sum(refit$param[, jj] * refit$weights)
    }
  }

  est_new <- est
  est_new$coef.form <- est_new$coef.form + coef.adj
  est_new$refit <- refit

  unlink("temp-refit-abc.rda")

  return(est_new)
}

est_new <- netest_refit_abc(est, nsims = 25, ncores = 5, nsteps = 300,
                            prior.min = -0.5, prior.max = 0)

est_new$refit
est$coef.form
est_new$coef.form

dx_new <- netdx(est_new, nsims = 10, nsteps = 300, ncores = 4, verbose = FALSE)
dx_new

par(mfrow = c(1,1))
plot(dx_new)

est_full <- netest(nw, formation, target.stats, coef.diss, edapprox = FALSE)

dx_full <- netdx(est, nsims = 10, nsteps = 500, ncores = 4, verbose = FALSE)
dx_full
