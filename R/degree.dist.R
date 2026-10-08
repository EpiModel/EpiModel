#' @title Extract Degree Distributions from Network Diagnostics
#'
#' @description Calculates the cumulative or momentary degree distribution from
#'              the timed edgelist of a dynamic [netdx()] simulation.
#'
#' @param x An `EpiModel` object of class `netdx`, run with
#'        `dynamic = TRUE` and `keep.tedgelist = TRUE`.
#' @param type Distribution to calculate, either `"cumulative"` for the number
#'        of partners each node accumulates over the observation window, or
#'        `"momentary"` for the number of partners each node has at a single
#'        point in time, averaged over the time steps in the window.
#' @param sims A vector of simulation numbers to include. The default is all
#'        simulations in `x`.
#' @param window A vector of length 2 giving the first and last time step of
#'        the observation window, on the time scale of the timed edgelist
#'        (which starts at 0). The default is the full simulation,
#'        `c(0, x$nsteps)`.
#' @param unique.partners If `TRUE`, a dyad that forms, dissolves, and then
#'        re-forms within the window contributes one partner to the cumulative
#'        degree of each node; if `FALSE`, it contributes one per partnership.
#'        Only used when `type = "cumulative"`.
#'
#' @details
#' The momentary degree distribution is the distribution plotted by
#' [plot.netdx()] when `degree(k)` terms are included in the `nwstats.formula`
#' of [netdx()]: it counts the partners a node has at one moment in time. The
#' cumulative degree distribution counts the partners a node accumulates across
#' the whole observation window, so a node with one partner at a time and a
#' node with several simultaneous partners may have the same cumulative degree.
#'
#' Partnerships are counted if they are active at any point in the window: a
#' partnership contributes to the cumulative degree of both of its nodes when
#' its spell overlaps `window`, and to the momentary degree of both nodes at
#' each time step of the window in which it is active. Partnerships that are
#' ongoing when the simulation ends are counted for the steps they are
#' observed.
#'
#' Because the cumulative degree grows with the length of the observation
#' window, `window` should be set to a substantively meaningful period (for
#' example, one year of time steps) when comparing across models or against
#' data.
#'
#' @return
#' A `data.frame` with one row per simulation and degree value, and the
#' columns:
#'
#'  * `sim`: the simulation number.
#'  * `degree`: the degree value, from 0 to the maximum observed across the
#'    selected simulations.
#'  * `count`: the number of nodes with that degree. For
#'    `type = "momentary"`, this is the mean number of nodes per time step,
#'    and so need not be a whole number.
#'  * `prop`: `count` divided by the network size.
#'
#' @seealso [plot.netdx()] with `type = "cumldeg"` to plot these distributions,
#'   and [netdx()] to generate the input object.
#'
#' @keywords extract
#' @export
#'
#' @examples
#' \donttest{
#' nw <- network_initialize(n = 100)
#' est <- netest(nw, formation = ~edges, target.stats = 50,
#'               coef.diss = dissolution_coefs(~offset(edges), duration = 10),
#'               verbose = FALSE)
#' dx <- netdx(est, nsims = 2, nsteps = 50, keep.tedgelist = TRUE,
#'             verbose = FALSE)
#'
#' # Cumulative degree over the full simulation
#' cd <- get_degree_dist(dx, type = "cumulative")
#' head(cd)
#'
#' # Mean cumulative degree per simulation
#' tapply(cd$degree * cd$prop, cd$sim, sum)
#'
#' # Momentary degree over the last 25 steps
#' get_degree_dist(dx, type = "momentary", window = c(25, 50))
#' }
#'
get_degree_dist <- function(x, type = c("cumulative", "momentary"), sims = NULL,
                            window = NULL, unique.partners = TRUE) {

  if (!inherits(x, "netdx")) {
    stop("x must be an object of class netdx")
  }
  type <- match.arg(type)

  if (x$dynamic == FALSE || is.null(x$tedgelist)) {
    stop("Degree distributions require the timed edgelist: rerun netdx with ",
         "dynamic = TRUE and keep.tedgelist = TRUE")
  }

  sims <- if (is.null(sims)) seq_len(x$nsims) else sims
  if (max(sims) > x$nsims) {
    stop("Maximum sim number is ", x$nsims)
  }

  window <- validate_degree_window(window, x$nsteps)
  n <- network.size(x$nw)

  dists <- lapply(x$tedgelist[sims], function(tel) {
    if (type == "cumulative") {
      degs <- tedgelist_cumulative_degree(tel, n, window, unique.partners)
      tabulate(degs + 1, nbins = max(degs) + 1)
    } else {
      tedgelist_momentary_degree(tel, n, window)
    }
  })

  ## put all sims on a common degree axis, running to the highest degree
  ## observed in any of them
  maxdeg <- max(0, vapply(dists, function(z) max(which(z > 0), 0) - 1, numeric(1)))

  out <- do.call(
    rbind,
    Map(function(sim, dist) {
      count <- numeric(maxdeg + 1)
      obs <- seq_len(min(length(dist), maxdeg + 1))
      count[obs] <- dist[obs]
      data.frame(sim = sim, degree = 0:maxdeg, count = count)
    }, sims, dists)
  )
  out$prop <- out$count / n
  row.names(out) <- NULL

  return(out)
}

## Validates the observation window for a degree distribution, returning a
## length-2 numeric of the first and last time step to observe. Timed edgelists
## from netdx start at time 0, so the full window is c(0, nsteps).
validate_degree_window <- function(window, nsteps) {
  window <- NVL(window, c(0, nsteps))
  if (length(window) != 2 || !is.numeric(window) || anyNA(window)) {
    stop("window must be a numeric vector of length 2, c(start, end)")
  }
  if (window[1] >= window[2]) {
    stop("window start must be less than window end")
  }
  if (window[1] < 0 || window[2] > nsteps) {
    stop("window must fall within the simulated time range, c(0, ", nsteps, ")")
  }
  window
}

#' @title Cumulative Degree of Each Node in a Timed Edgelist
#'
#' @param tedgelist A timed edgelist, as produced by
#'   [`networkDynamic::as.data.frame.networkDynamic`].
#' @param n Size of the network the edgelist was simulated on.
#' @param window Length-2 numeric of the first and last time step to observe.
#' @param unique.partners If `TRUE`, count distinct partners rather than
#'   distinct partnerships.
#'
#' @return An integer vector of length `n`, giving the number of partners each
#'   node has over `window`.
#' @keywords internal
tedgelist_cumulative_degree <- function(tedgelist, n, window,
                                        unique.partners = TRUE) {
  ## spells are half open, [onset, terminus), so a spell is observed if it is
  ## active at any point in [window[1], window[2]]
  active <- tedgelist$onset <= window[2] & tedgelist$terminus > window[1]
  tail <- tedgelist$tail[active]
  head <- tedgelist$head[active]

  if (unique.partners == TRUE) {
    dyads <- cbind(pmin(tail, head), pmax(tail, head))
    keep <- !duplicated(dyads)
    tail <- tail[keep]
    head <- head[keep]
  }

  tabulate(c(tail, head), nbins = n)
}

#' @title Momentary Degree Distribution of a Timed Edgelist
#'
#' @inheritParams tedgelist_cumulative_degree
#'
#' @return A numeric vector of length `n`, giving the mean number of nodes per
#'   time step with degree 0 through `n - 1` over `window`.
#' @keywords internal
tedgelist_momentary_degree <- function(tedgelist, n, window) {
  nsteps <- window[2] - window[1] + 1

  ## time indices of the steps at which each partnership starts and stops
  ## contributing to degree, clamped to the window and shifted to start at 1
  starts <- pmax(tedgelist$onset, window[1]) - window[1] + 1
  stops <- pmin(tedgelist$terminus, window[2] + 1) - window[1] + 1
  active <- starts < stops
  starts <- starts[active]
  stops <- stops[active]
  nodes <- c(tedgelist$tail[active], tedgelist$head[active])

  ## one entry per time step, holding the nodes gaining or losing a partner
  ## at that step; stops may fall one step past the window, for partnerships
  ## still active at its end
  gains <- split(nodes, factor(rep(starts, 2), levels = seq_len(nsteps)))
  losses <- split(nodes, factor(rep(stops, 2), levels = seq_len(nsteps + 1)))

  degs <- integer(n)
  dist <- numeric(n)
  for (at in seq_len(nsteps)) {
    degs <- degs + tabulate(gains[[at]], nbins = n) -
      tabulate(losses[[at]], nbins = n)
    dist <- dist + tabulate(degs + 1, nbins = n)
  }

  dist / nsteps
}
