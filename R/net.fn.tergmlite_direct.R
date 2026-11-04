# tergmLite direct engine (prototype) ------------------------------------------
#
# Under `control.net(tergmLite = TRUE, tergmLite.direct = TRUE)`, simulate_dat()
# hands each model-based layer to simulate_dat_direct() instead of
# stats::simulate(). Per layer, the ergm model and the Metropolis-Hastings
# proposal are built once and kept in `dat$tergmLite.cache` (derived state: it
# is rebuilt whenever it is missing or stale, so it never needs to be saved).
# Each time step then
#
#   1. recomputes the inputs of the cached model's nodal-attribute terms from
#      the current attributes and network size (nodematch, nodefactor,
#      degrange), against the levels frozen at build time; if the set of
#      observed attribute values has changed, the model is rebuilt instead,
#   2. computes eta from the current coefficients (edges_correct() changes the
#      edges coefficient every step),
#   3. runs tergm::tergm_MCMC_sample() once on a networkLite basis, and
#   4. reads the new edgelist straight from the returned ergm_state, with no
#      as.network() / add.edges() round trip.
#
# The model, proposal, eta and sampler are those simulate() would use, so the
# Markov chain is the same one, draw for draw, given the same RNG state.
#
# Terms are recomputed only when the patch can be verified: at build time each
# patch is applied to the build-time network and must reproduce the inputs
# ergm computed. A term with no verified patch makes the layer fall back to a
# full model rebuild every step (still without the conversion round trip).

simulate_dat_direct <- function(dat, at, network) {
  nwparam <- get_nwparam(dat, network = network)
  control <- get_network_control(dat, network, "set.control.tergm")
  dur1 <- !(nwparam$coef.diss$duration[1] > 1)
  if (dur1) {
    coef <- nwparam$coef.form
    control$MCMC.prop.args <- list(discordance_fraction = 0)
  } else {
    coef <- c(nwparam$coef.form, nwparam$coef.diss$coef.adj)
  }
  key <- list(formation = nwparam$formation,
              dissolution = if (!dur1) nwparam$coef.diss$dissolution,
              constraints = nwparam$constraints,
              control = control[c("MCMC.prop", "MCMC.prop.weights",
                                  "MCMC.prop.args", "term.options")])

  ## the basis: networkLite on the current edgelist and attributes
  nw <- get_network(dat, network)
  n <- network.size(nw)

  cache <- if (length(dat$tergmLite.cache) >= network) dat$tergmLite.cache[[network]]
  fresh <- is.null(cache) || !identical(cache$key, key)
  if (!fresh) {
    cache <- tld_update(cache, dat, nw, n)
    fresh <- is.null(cache)
  }
  if (fresh) {
    cache <- tld_build(key, nwparam, dur1, control, nw, n)
    b <- dat$tergmLite.stats$builds
    b <- c(b, rep(0, max(0, network - length(b))))
    b[network] <- b[network] + 1
    dat$tergmLite.stats$builds <- b
  }

  ## time bookkeeping, as tergm's simulate_formula.network() does it:
  ## the chain starts at `at - 1` and its TICK advances the clock to `at`
  nw %n% "time" <- at - 1L
  ctl <- control
  ctl$time.burnin <- 0
  ctl$time.interval <- 1
  ctl$time.samplesize <- 1L
  ctl$collect <- FALSE
  ctl$changes <- TRUE
  if (length(coef) != cache$npar) {
    stop("coef has ", length(coef), " elements, while the model requires ",
         cache$npar, " parameters.")
  }
  eta <- ergm.eta(coef, cache$model$etamap)
  z <- tergm::tergm_MCMC_sample(nw, cache$model, NULL, cache$proposal,
                                control = ctl, eta = eta)

  ## edgelist straight from the ergm_state
  st <- z$newnetwork
  el <- matrix(c(st$el$.tail, st$el$.head), ncol = 2L)
  attr(el, "n") <- as.integer(n)
  attr(el, "vnames") <- network.vertex.names(nw)
  attr(el, "bipartite") <- FALSE
  attr(el, "directed") <- FALSE
  attr(el, "loops") <- FALSE
  class(el) <- c("matrix_edgelist", "edgelist", "matrix", "array")
  dat$run$el[[network]] <- el
  if (get_network_control(dat, network, "tergmLite.track.duration") == TRUE) {
    dat$run$net_attr[[network]][["time"]] <- st$nw0 %n% "time"
    dat$run$net_attr[[network]][["lasttoggle"]] <- st$nw0 %n% "lasttoggle"
  }
  if (is.null(dat$tergmLite.cache)) dat$tergmLite.cache <- list()
  dat$tergmLite.cache[[network]] <- cache
  dat
}

# Build the model and proposal for a layer, as tergm's
# simulate_formula.network() does, and work out the per-step patch plan.
tld_build <- function(key, nwparam, dur1, control, nw, n) {
  if (dur1) {
    formula <- nwparam$formation
  } else {
    formula <- ~Form(nwparam$formation) + Persist(nwparam$coef.diss$dissolution)
  }
  proposal <- ergm_proposal(nwparam$constraints,
                            arguments = control$MCMC.prop.args, nw = nw,
                            weights = control$MCMC.prop.weights,
                            hints = control$MCMC.prop, class = "t")
  model <- ergm_model(formula, nw, dynamic = TRUE,
                      term.options = control$term.options,
                      extra.aux = list(proposal = proposal$auxiliaries,
                                       system = ~.lasttoggle))
  proposal$aux.slots <- model$slots.extra.aux$proposal

  ## A proposal for an unconstrained model depends on neither the attributes
  ## nor the network size; its R-side `arguments` carry the build-time network
  ## (tens of MB at 100k nodes) and are not read by the sampler.
  static_prop <- tld_trivial_constraints(nwparam$constraints)
  if (static_prop) proposal$arguments <- NULL

  plan <- if (static_prop) tld_plan(model, nw, n, environment(formula)) else NULL
  list(key = key, model = model, proposal = proposal, plan = plan,
       reusable = !is.null(plan), npar = nparam(model))
}

tld_trivial_constraints <- function(constraints) {
  rhs <- statnet.common::list_rhs.formula(constraints)
  length(rhs) == 0L || all(vapply(rhs, function(t) identical(t, as.name(".")),
                                  logical(1)))
}

# Recompute the attribute-dependent inputs; NULL means "rebuild".
tld_update <- function(cache, dat, nw, n) {
  if (!cache$reusable) return(NULL)
  attrs <- dat$run$attr
  memo <- list()
  enc <- list()
  model <- cache$model
  for (spec in cache$plan) {
    if (spec$kind == "degrange") {
      trm <- tld_get(model, spec$path)
      to <- ifelse(spec$to == Inf, pmax(spec$from, n) + 1, spec$to)
      trm$inputs <- as.double(c(rbind(spec$from, to)))
      if (any(spec$from == 0)) {
        trm$emptynwstats <- as.double(replace(rep(0, length(spec$from)),
                                              spec$from == 0, n))
      }
      model <- tld_set(model, spec$path, trm)
      next
    }
    ## nodal attribute terms: match the current values to the frozen levels
    akey <- spec$akey
    m <- memo[[akey]]
    if (is.null(m)) {
      raw <- if (spec$fast) attrs[[spec$attr]] else ergm_get_vattr(spec$attr, nw)
      if (length(raw) != n) return(NULL)
      m <- match(raw, spec$U)
      if (anyNA(m) || any(tabulate(m, length(spec$U)) == 0L)) return(NULL)
      memo[[akey]] <- m
    }
    trm <- tld_get(model, spec$path)
    v <- enc[[spec$ekey]]
    if (is.null(v)) {
      if (spec$kind == "nodematch") {
        code <- spec$map[m]
        nd <- which(is.na(code))
        code[nd] <- spec$nu + seq_along(nd)
        v <- as.double(if (spec$diff) c(seq_len(spec$nu), code) else code)
      } else {
        v <- as.integer(spec$map[m])
      }
      enc[[spec$ekey]] <- v
    }
    if (spec$kind == "nodematch") trm$inputs <- v else trm$iinputs <- v
    model <- tld_set(model, spec$path, trm)
  }
  cache$model <- model
  cache
}

tld_get <- function(x, path) {
  for (p in path) x <- x[[p]]
  x
}

tld_set <- function(x, path, value) {
  if (length(path) == 1L) {
    x[[path[[1L]]]] <- value
    return(x)
  }
  x[[path[[1L]]]] <- tld_set(x[[path[[1L]]]], path[-1L], value)
  x
}

# Walk the model (and the submodels of operator terms such as Form() and
# Persist()) and return a list of verified patch specs, or NULL if any term
# cannot be patched.
tld_plan <- function(model, nw, n, env, path = list()) {
  specs <- list()
  for (i in seq_along(model$terms)) {
    trm <- model$terms[[i]]
    p <- c(path, list("terms", i))
    if (inherits(trm$submodel, "ergm_model")) {
      if (length(trm$inputs) || length(trm$iinputs) || !is.null(trm$ext.encode)) {
        return(NULL)
      }
      ## the operator's formula argument, evaluated where the term was
      subf <- tryCatch(eval(trm$call[[2L]], env), error = function(e) NULL)
      if (!inherits(subf, "formula")) return(NULL)
      sub <- tld_plan(trm$submodel, nw, n, environment(subf),
                      c(p, list("submodel")))
      if (is.null(sub)) return(NULL)
      specs <- c(specs, sub)
      next
    }
    spec <- tld_term_spec(trm, nw, n, env)
    if (is.null(spec)) return(NULL)
    if (spec$kind != "static") {
      spec$path <- p
      specs <- c(specs, list(spec))
    }
  }
  specs
}

# Patch spec for one leaf term, verified against the inputs ergm computed.
tld_term_spec <- function(trm, nw, n, env) {
  if (length(trm$coef.names) == 0L) {
    ## auxiliary: static if it carries no inputs
    if (length(trm$inputs) || length(trm$iinputs)) return(NULL)
    return(list(kind = "static"))
  }
  if (!identical(trm$pkgname, "ergm")) return(NULL)
  args <- if (is.call(trm$call)) {
    tryCatch(lapply(as.list(trm$call)[-1L], eval, envir = env),
             error = function(e) NULL)
  } else {
    list()
  }
  if (is.null(args)) return(NULL)
  switch(trm$name,
    edges = list(kind = "static"),
    concurrent = if (!length(trm$inputs) && !length(trm$iinputs)) list(kind = "static"),
    degree = list(kind = "static"),
    degrange = tld_spec_degrange(trm, args, n),
    nodematch = tld_spec_nodal(trm, args, nw, "nodematch"),
    nodefactor = tld_spec_nodal(trm, args, nw, "nodefactor"),
    NULL
  )
}

tld_match_args <- function(args, stub) {
  tryCatch(do.call(stub, args), error = function(e) NULL)
}

tld_spec_degrange <- function(trm, args, n) {
  a <- tld_match_args(args, function(from, to = Inf, by = NULL,
                                     homophily = FALSE, levels = NULL) {
    list(from = from, to = to, by = by)
  })
  if (is.null(a) || !is.null(a$by)) return(NULL)
  from <- a$from
  to <- a$to
  if (length(to) == 1 && length(from) > 1) to <- rep(to, length(from))
  if (length(from) == 1 && length(to) > 1) from <- rep(from, length(to))
  spec <- list(kind = "degrange", from = from, to = to)
  to2 <- ifelse(to == Inf, pmax(from, n) + 1, to)
  if (!identical(as.double(c(rbind(from, to2))), trm$inputs)) return(NULL)
  spec
}

tld_spec_nodal <- function(trm, args, nw, kind) {
  if (kind == "nodematch") {
    stub <- function(attr, diff = FALSE, keep = NULL, levels = NULL) {
      list(attr = attr, diff = diff, keep = keep, levels = levels,
           levels_missing = missing(levels))
    }
  } else {
    stub <- function(attr, base = 1, levels = NULL) {
      list(attr = attr, base = base, levels = levels,
           levels_missing = missing(levels))
    }
  }
  a <- tld_match_args(args, stub)
  if (is.null(a)) return(NULL)
  ## level specs that depend only on the observed values
  if (is.function(a$levels) || inherits(a$levels, "formula")) return(NULL)
  nodecov <- tryCatch(ergm_get_vattr(a$attr, nw), error = function(e) NULL)
  if (is.null(nodecov)) return(NULL)
  U <- sort(unique(nodecov))
  u <- tryCatch(ergm_attr_levels(a$levels, nodecov, nw, levels = U),
                error = function(e) NULL)
  if (is.null(u)) return(NULL)
  fast <- is.character(a$attr) && length(a$attr) == 1L && !inherits(a$attr, "AsIs")
  spec <- list(kind = kind, attr = a$attr, fast = fast, U = U,
               akey = paste(deparse(a$attr), collapse = ""))
  m <- match(nodecov, U)
  if (kind == "nodematch") {
    if (a$levels_missing && !is.null(a$keep)) u <- u[a$keep]
    spec$diff <- isTRUE(a$diff)
    spec$nu <- length(u)
    spec$map <- match(U, u)
    code <- spec$map[m]
    nd <- which(is.na(code))
    code[nd] <- spec$nu + seq_along(nd)
    got <- as.double(if (spec$diff) c(seq_len(spec$nu), code) else code)
    if (!identical(got, trm$inputs)) return(NULL)
    spec$ekey <- paste(kind, spec$akey, spec$diff, paste(spec$map, collapse = ","))
  } else {
    if (a$levels_missing && any(NVL(a$base, 0) != 0)) u <- u[-a$base]
    if (length(u) == 0) return(NULL)
    spec$map <- match(U, u, nomatch = 0) - 1L
    got <- as.integer(spec$map[m])
    if (!identical(got, trm$iinputs)) return(NULL)
    spec$ekey <- paste(kind, spec$akey, paste(spec$map, collapse = ","))
  }
  spec
}
