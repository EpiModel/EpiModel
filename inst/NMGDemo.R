
##
## EpiModel NMG Demo
## February 2021
##


# 1. Accessor Functions  --------------------------------------------------

# General demonstration
dat <- list(
  attr = list(
    active = rbinom(100, 1, 0.9)
  ),
  epi = list(),
  param = list(),
  init = list(),
  control = list(
    nsteps = 150
  )
)

dat <- add_attr(dat, "age")
dat <- set_attr(dat, "age", runif(100))
dat <- set_attr(dat, "status", rbinom(100, 1, 0.9))
dat <- set_attr(dat, "status", rep(1, 150), override.length.check = TRUE)
dat <- append_attr(dat, "status", 1, 10)
dat <- append_attr(dat, "age", NA, 10)
get_attr_list(dat)
get_attr_list(dat, c("age", "active"))
get_attr(dat, "status")
get_attr(dat, "status", c(1, 4))

dat <- add_epi(dat, "i.num")
dat <- set_epi(dat, "i.num", 150, 10)
dat <- set_epi(dat, "s.num", 150, 90)
get_epi_list(dat)
get_epi_list(dat, c("i.num", "s.num"))
get_epi(dat, "i.num")
get_epi(dat, "i.num", c(1, 4))
get_epi(dat, "i.num", rbinom(150, 1, 0.2) == 1)

dat <- add_param(dat, "x")
dat <- set_param(dat, "x", 0.4)
dat <- set_param(dat, "y", 0.8)
get_param_list(dat)
get_param_list(dat, c("x", "y"))
get_param(dat, "x")

dat <- add_init(dat, "x")
dat <- set_init(dat, "x", 0.4)
dat <- set_init(dat, "y", 0.8)
get_init_list(dat)
get_init_list(dat, c("x", "y"))
get_init(dat, "x")

dat <- add_control(dat, "x")
dat <- set_control(dat, "x", 0.4)
dat <- set_control(dat, "y", 0.8)
get_control_list(dat)
get_control_list(dat, c("x", "y"))
get_control(dat, "x")


# Specific example of use
recovery.net <- function(dat, at) {

  ## Only run with SIR/SIS
  type <- get_control(dat, "type")
  if (!(type %in% c("SIR", "SIS")) && !is.null(type)) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  recovState <- ifelse(type == "SIR", "r", "s")

  rec.rate <- get_param(dat, "rec.rate")

  nRecov <- 0
  idsElig <- which(active == 1 & status == "i")
  nElig <- length(idsElig)


  # Time-Varying Recovery Rate ----------------------------------------------
  infDur <- at - infTime[active == 1 & status == "i"]
  infDur[infDur == 0] <- 1
  lrec.rate <- length(rec.rate)
  if (lrec.rate == 1) {
    ratesElig <- rec.rate
  } else {
    ratesElig <- ifelse(infDur <= lrec.rate, rec.rate[infDur],
                        rec.rate[lrec.rate])
  }


  # Process -----------------------------------------------------------------
  if (nElig > 0) {
    vecRecov <- which(rbinom(nElig, 1, ratesElig) == 1)
    if (length(vecRecov) > 0) {
      idsRecov <- idsElig[vecRecov]
      nRecov <- length(idsRecov)
      status[idsRecov] <- recovState
    }
  }
  dat <- set_attr(dat, "status", status)

  # Output ------------------------------------------------------------------
  outName <- ifelse(type == "SIR", "ir.flow", "is.flow")
  dat <- set_epi(dat, outName, at, nRecov)

  return(dat)
}



# 2. UIDs and Core Attributes ---------------------------------------------

initialize.net <- function(x, param, init, control, s) {

  ### ...

  # Nodal Attributes

  # Standard attributes
  num <- network.size(nw)
  dat <- append_core_attr(dat, 1, num)

  ## Pull attr on nw to dat$attr
  dat <- copy_nwattr_to_datattr(dat)

  ## Store current proportions of attr
  nwterms <- get_network_term_attr(nw)
  if (!is.null(nwterms)) {
    dat$temp$nwterms <- nwterms
    dat$temp$t1.tab <- get_attr_prop(dat, nwterms)
  }

  ## Infection Status and Time
  dat <- init_status.net(dat)

  # Conversions for tergmLite
  if (control$tergmLite == TRUE) {
    dat <- tergmLite::init_tergmLite(dat)
  }

  ### ...

  return(dat)
}

# Inside the function
append_core_attr <- function(dat, at, n.new) {
  dat <- append_attr(dat, "active", 1, n.new)
  dat <- append_attr(dat, "entrTime", at, n.new)
  dat <- append_attr(dat, "exitTime", NA, n.new)

  dat <- update_uids(dat, n.new)

  return(dat)
}

update_uids <- function(dat, n.new) {
  last_uid <- if (is.null(dat[["_last_uid"]])) 0L else dat[["_last_uid"]]
  next_uids <- seq_len(n.new) + last_uid
  dat[["_last_uid"]] <- last_uid + as.integer(n.new)
  dat <- append_attr(dat, "uid", next_uids, n.new)

  return(dat)
}

# Use in arrival module
arrivals.net <- function(dat, at) {

  # Conditions --------------------------------------------------------------
  vital <- get_param(dat, "vital")
  if (vital == FALSE) {
    return(dat)
  }

  # Variables ---------------------------------------------------------------
  a.rate <- get_param(dat, "a.rate")
  index <- at - 1
  nOld <- get_epi(dat, "num", index)
  nArrivals <- 0

  # Add Nodes ---------------------------------------------------------------
  if (nOld > 0) {
    nArrivals <- rbinom(1, nOld, a.rate)
    if (nArrivals > 0) {
      dat <- append_core_attr(dat, at, nArrivals)
      dat <- append_attr(dat, "status", "s", nArrivals)
      dat <- append_attr(dat, "infTime", NA, nArrivals)
    }
  }

  # Output ------------------------------------------------------------------
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  return(dat)
}


# 3. Saving Transmission Matrix -------------------------------------------

# http://statnet.org/nme/d5-s4-COVIDDemog.html


# 4. Revised EDA ----------------------------------------------------------

# The NME TERGM example that is supposed to break the EDA
n <- 500
net3 <- network_initialize(n)
net3 <- set_vertex_attribute(net3, "race", rep(c("B", "W"), each = n/2))
net3

form.formula.3 <- ~edges + nodematch("race") + degree(0) + concurrent
target.stats.3 <- c(0.9*n/2, (0.9*n/2)*(5/6), 0.36*n, 0.18*n)

diss.formula.3 <- ~offset(edges) + offset(nodematch("race"))

fit4 <- netest(net3,
               formation = form.formula.3,
               target.stats = target.stats.3,
               coef.diss = dissolution_coefs(diss.formula.3, c(20, 10)))
sim4 <- netdx(fit4, nsteps = 1000, nsims = 10, ncores = 5, keep.tedgelist = TRUE)
sim4

plot(sim4)


# 5. Random Parameters ----------------------------------------------------

## Example SIR model parameterization with fixed and random parameters
# Network model estimation
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

# Random parameter list
my_randoms <- list(
  act.rate = param_random(1:3),
  inf.prob = function() rbeta(1, 1, 2),
)

# Parameters, initial conditions, and control settings
param <- param.net(rec.rate = 0.02, random.params = my_randoms)
init <- init.net(i.num = 10, r.num = 0)
control <- control.net(type = "SIR", nsteps = 100, nsims = 1,
                       resimulate.network = TRUE, tergmLite = TRUE)

# Simulate the model
sim <- netsim(est, param, init, control)

# Print and plot
sim


# Under the hood of netsim

# Define random parameter list
my_randoms <- list(
  act.rate = param_random(c(0.25, 0.5, 0.75)),
  tx.prob = function() rbeta(1, 1, 2),
  stratified.test.rate = function() c(
    rnorm(1, 0.05, 0.01),
    rnorm(1, 0.15, 0.03),
    rnorm(1, 0.25, 0.05)
  )
)

# Parameter model with deterministic and random parameters
param <- param.net(inf.prob = 0.3, random.params = my_randoms)

# Parameters are drawn automatically in netsim by calling the function
# within netsim_loop. Demonstrating draws here but this is not used by
# end user.
paramDraw <- generate_random_params(param, verbose = TRUE)


# 6. Update Parameters ----------------------------------------------------

## Set parameters during calibration
param <- param_msm(netstats = netstats,
                   epistats = epistats,
                   hiv.test.rate = c(0.00385, 0.00380, 0.00690),
                   tx.init.prob = c(0.1775, 0.190, 0.2521),
                   tx.halt.partial.prob = c(0.0062, 0.0055, 0.0031),
                   tx.reinit.partial.prob = c(0.00255, 0.00255, 0.00255),
                   trans.scale = c(2.44, 0.424, 0.270),
                   riskh.start = 52 * 59,
                   prep.start = (52 * 60) + 1,
                   prep.start.prob = 0.66)

saveRDS(param, "data/input/param.burnin1.rds")


## Ready to run intervention scenarios

param <- readRDS("data/input/param.burnin1.rds")

init <- init_msm()
control <- control_msm(simno = fsimno,
                       start = (52 * 60) + 1,
                       nsteps = 52 * 65,
                       nsims = ncores,
                       ncores = ncores,
                       initialize.FUN = reinit_msm)

# Intervention parameters
interv_params <- list(
  part.ident.start = Inf,
  part.index.window = 0,
  part.index.degree = 1,
  part.ident.main.window = 12,
  part.ident.casl.window = 12,
  part.ident.ooff.window = 12,
  part.ident.main.prob = 1,
  part.ident.casl.prob = 1,
  part.ident.ooff.prob = 1,
  part.hiv.test.rate = c(1, 1, 1),
  part.prep.start.prob = 0.5,
  part.tx.init.prob = c(0.6, 0.6, 0.8),
  part.tx.halt.prob = c(0.00102, 0.00102, 0.00071),
  part.tx.reinit.prob = c(0.5, 0.5, 0.5)
)
param <- update_params(param, interv_params)
param



# 7. Smart Error Messages -------------------------------------------------

nw <- network_initialize(n = 100)
est <- netest(nw, formation = ~edges, target.stats = 25,
              coef.diss = dissolution_coefs(~offset(edges), 10, 0),
              verbose = FALSE)

## "SI, 1M, CL: 1 sim"
param <- param.net(inf.prob = 0.5)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsims = 1, nsteps = 25, verbose = FALSE)
x <- netsim(est, param, init, control)

# Addition of a module generating errors and warning
error_module <- function(dat, at) {
  if (at > 19) {
    warning("this is a warning")
  }

  if (at > 20) {
    stop("This is an error")
  }

  dat
}

control <- control.net(type = NULL, nsims = 1, nsteps = 25, verbose = FALSE,
                       error_module.FUN = error_module)
x <- netsim(est, param, init, control)


# 8. TODO: Time-varying attribute storage ---------------------------------

# https://github.com/statnet/EpiModel/issues/386
# https://github.com/statnet/EpiModel/issues/478


# 9. TODO: Flexible Prevalence Module with Factory Functions --------------

epi_sup_dx_this_year <- function(r_ind) { # r_ind is the race(s) of interest
  function(dat, at) {
    with(dat$attr, {

      cond <- at - vl.last.supp <= 52
      pop <- (at - diag.time <= 52 + 52 & at - diag.time >= 52 &
                race %in% r_ind)

      sum(cond & pop, na.rm = TRUE) / sum(pop, na.rm = TRUE)
    })
  }
}

# for the same tracker but without race distinction, I could then use

param$epi_funs <- c(
  "sup_dx_this_year" = epi_sup_dx_this_year(1:3),
  "sup_dx_this_year.B" = epi_sup_dx_this_year(1),
  "sup_dx_this_year.H" = epi_sup_dx_this_year(2),
  "sup_dx_this_year.W" = epi_sup_dx_this_year(3)
)

prevalence.net <- function(dat, at) {
  epi_funs <- get_params(dat, epi_funs)

  for (fn in names(epi_funs)) {
    dat <- set_epi(dat, fn, at, epi_funs[[fn]](dat, at))
  }
}


# 10. TODO: Move Initialization Function Outside netsim -------------------

# https://github.com/statnet/EpiModel/issues/385

# Example for EpiModelHIV

initialize_msm <- function(x, param, init, control, s) {

  ## Master Data List Setup ##
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control


  ## Network Setup ##
  # Initial network simulations
  dat$nw <- list()
  for (i in 1:3) {
    dat$nw[[i]] <- simulate(x[[i]]$fit, basis = x[[i]]$fit$newnetwork)
  }
  nw <- dat$nw

  # Pull Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Nodal Attributes Setup ##
  dat$attr <- param$netstats$attr

  num <- network.size(nw[[1]])
  # dat$attr$active <- rep(1, num)
  # dat$attr$uid <- 1:num
  dat <- append_core_attr(dat, 1, num)
  # dat$attr$arrival.time <- rep(1, num)

  # Circumcision
  rates <- param$circ.prob[dat$attr$race]
  dat$attr$circ <- rbinom(length(rates), 1, rates)

  # Insertivity Quotient
  ins.quot <- rep(NA, num)
  role.class <- dat$attr$role.class
  ins.quot[role.class == 0]  <- 1
  ins.quot[role.class == 1]  <- 0
  ins.quot[role.class == 2]  <- runif(sum(role.class == 2))
  dat$attr$ins.quot <- ins.quot

  # HIV-related attributes
  dat <- init_status_msm(dat)

  # STI Status
  dat <- init_sti_msm(dat)

  # PrEP-related attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)
  dat$attr$prepStartTime <- rep(NA, num)
  dat$attr$prepLastRisk <- rep(NA, num)
  dat$attr$prepLastStiScreen <- rep(NA, num)

  # Partner Identification attributes
  dat$attr$part.scrnd <- rep(NA, num)
  dat$attr$part.ident <- rep(NA, num)
  dat$attr$part.art <- rep(NA, num)

  ## Other Setup ##
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  # Prevalence Tracking
  # dat$temp$max.uid <- num
  dat <- prevalence_msm(dat, at = 1)

  # Setup Partner List
  dat <- init_plist(dat)

  # Network statistics
  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_nwstats(dat, at = 1)
  }

  # dat$param$netstats <- NULL
  class(dat) <- "dat"
  return(dat)
}
