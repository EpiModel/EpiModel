sim <- readRDS("../EpiModelHIV-Template/data/run/calibration/sim__bad_calib__1.rds")

sim1 <- get_sims(sim, 1)
rp <- make_restart_point(sim, c("hiv.dx.time"), sim_num = 1, keep_steps = 1)

merge(rp, sim1, param.error = FALSE)
debug(merge.netsim)

out <- readRDS("nwstats_running_done.rds")
running <- out$running
done <- out$done
rm(out)

running$stats
done$stats


running$stats$nwstats[[1]]

#TODO: use this to finish initialize of $stats
#  - make the padding in the summary net + transmat creator
processed <- done$stats$nwstats$sim1[[1]] |>
  as.matrix() |>
  apply(1, c, simplify = FALSE) |>
  unname()

identical(processed, running$stats$nwstats[[1]])

dat <- list()

x <- running

dat$stats <- list()

for (nw in seq_len(1)) {
  x$stats$nwstats$sim1[[nw]] |>
    as.matrix() |>
    apply(1, c, simplify = FALSE) |>
    unname()
}
