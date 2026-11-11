#!/usr/bin/env Rscript
# Reference case for the run time of netsim(): a built-in SIS epidemic on
# 10,000 nodes with arrivals and departures, tergmLite, 100 time steps.
#
#   Rscript inst/bench/sis-10k.R
#
# Prints the milliseconds per time step of each module and of the whole step.
# To compare two versions of EpiModel, run it with each one back to back, at
# least twice, on the same machine and with the same R_VSIZE and R_NSIZE. For
# a source tree: Rscript -e 'devtools::load_all(); source("inst/bench/sis-10k.R")'

library(EpiModel)

network_size <- 10000
nsteps <- 100

# Three groups holding 40%, 20% and 40% of the nodes, with mean degrees of
# 0.5, 0.6 and 0.7, and 70% of the 3,000 partnerships within a group
nw <- network_initialize(n = network_size)
nw <- set_vertex_attribute(nw, "race", rep(1:3, times = network_size * c(0.4, 0.2, 0.4)))
formation <- ~edges + nodematch("race", diff = TRUE) + nodefactor("race") +
  concurrent + degrange(from = 4)
target_stats <- c(edges = 3000, same_race_1 = 700, same_race_2 = 300, same_race_3 = 1100,
                  degree_race_2 = 1200, degree_race_3 = 2800, concurrent = 900, degree_4_up = 0)

set.seed(1)
est <- suppressMessages(
  netest(nw, formation = formation, target.stats = unname(target_stats),
         coef.diss = dissolution_coefs(~offset(edges), duration = 50, d.rate = 0.001),
         verbose = FALSE)
)

param <- param.net(inf.prob = 0.1, act.rate = 2, rec.rate = 0.02,
                   a.rate = 0.001, ds.rate = 0.001, di.rate = 0.001)
init <- init.net(i.num = network_size / 10)
control <- control.net(type = "SIS", nsteps = nsteps, nsims = 1, ncores = 1,
                       tergmLite = TRUE, resimulate.network = TRUE,
                       verbose = FALSE, .timing = TRUE)

set.seed(2)
gc_seconds_before <- gc.time()[[3]]
sim <- netsim(est, param, init, control)
gc_ms <- (gc.time()[[3]] - gc_seconds_before) * 1000 / (nsteps - 1)

# The first simulated step (time 2) is left out: it carries one-time costs
epi <- as.data.frame(sim)
timing_columns <- grep("^timing[.]", names(epi), value = TRUE)
ms_per_step <- colMeans(epi[epi$time > 2, timing_columns]) * 1000
names(ms_per_step) <- sub("^timing[.]", "", timing_columns)
step_ms <- ms_per_step[["step"]]
module_ms <- ms_per_step[names(ms_per_step) != "step"]

# A fixed workload timed in this process (the fastest of five runs), so that
# step times from machines of different speed compare as a ratio to it
cpu_reference_ms <- min(replicate(5, {
  reference_start <- Sys.time()
  values <- (seq_len(1e6) * 7919) %% 1000003
  sorted_values <- sort(values)
  positions <- match(values, sorted_values[seq_len(1e5)])
  running <- 0
  for (index in seq_len(5e5)) running <- running + values[index] %% 7
  as.numeric(Sys.time() - reference_start, units = "secs") * 1000
}))

floor_setting <- function(variable) {
  value <- Sys.getenv(variable)
  if (nzchar(value)) value else "not set"
}
cat(
  "EpiModel ", format(packageVersion("EpiModel")), ", ", R.version.string, ", ",
  R.version$platform, "\n",
  "R_VSIZE ", floor_setting("R_VSIZE"), ", R_NSIZE ", floor_setting("R_NSIZE"), "\n",
  format(network_size, big.mark = ","), " nodes, ", nsteps, " steps\n\n",
  sep = ""
)
print(data.frame(
  ms.per.step = round(c(module_ms, other = step_ms - sum(module_ms), step = step_ms), 2)
))
cat(
  "\nGarbage collection, inside the times above: ", round(gc_ms, 2), " ms per step\n",
  "CPU reference: ", round(cpu_reference_ms), " ms; step over reference: ",
  round(step_ms / cpu_reference_ms, 4), "\n",
  sep = ""
)
