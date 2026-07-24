
#' @title Definition for absdiffnodemix ERGM Term
#'
#' @description This function defines and initializes the absdiffnodemix ERGM
#'              term that allows for targeting homophily based on a non-binary
#'              attribute (e.g., age) by combinations of a binary attribute
#'              (e.g., race).
#'
#' @param nw An object of class `network`.
#' @param arglist A list of arguments as specified in the `ergm.userterms`
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        `ergm.userterms` package framework.
#'
#' @details
#' This ERGM user term was written to allow for age-based homophily in
#' partnership formation that is heterogeneous by race. The `absdiff`
#' component targets the distribution of age mixing on that continuous
#' variable, and the `nodemix` component differentiates this for
#' black-black, black-white, and white-white couples.
#'
#' @aliases absdiffnodemix
#'
InitErgmTerm.absdiffnodemix <- function(nw, arglist, ...) {

  a <- check.ErgmTerm(nw,
                      arglist,
                      directed = FALSE,
                      bipartite = FALSE,
                      varnames = c("attr", "by"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  nodecov <- ergm_get_vattr(a$attr, nw, accept = "numeric")
  nodecovby <- ergm_get_vattr(a$by, nw)
  nodecovbyname <- attr(nodecovby, "name")
  u <- sort(unique(nodecovby))
  if (any(is.na(nodecovby))) {
    u <- c(u, NA)
  }

  nodecovby <- match(nodecovby, u, nomatch = length(u) + 1)
  ui <- seq(along = u)

  uui <- matrix(seq_along(ui) ^ 2, length(ui), length(ui))
  urm <- t(sapply(ui, rep, length(ui)))
  ucm <- sapply(ui, rep, length(ui))
  uun <- outer(u, u, paste, sep = ".")
  uui <- uui[upper.tri(uui, diag = TRUE)]
  urm <- urm[upper.tri(urm, diag = TRUE)]
  ucm <- ucm[upper.tri(ucm, diag = TRUE)]
  uun <- uun[upper.tri(uun, diag = TRUE)]

  inputs <-  c(length(nodecov), length(urm), nodecov, nodecovby, urm, ucm)

  list(name = "absdiffnodemix",
       coef.names = paste("absdiffnodemix", attr(nodecov, "name"),
                          nodecovbyname, uun, sep = "."),
       pkgname = "EpiModel",
       inputs = inputs,
       dependence = FALSE)
}


#' @title Definition for absdiffby ERGM Term
#'
#' @description This function defines and initializes the absdiffby ERGM term
#'              that allows for representing homophily with respect to a
#'              non-binary attribute (e.g., age) differentially by a binary
#'              attribute (e.g., sex).
#'
#' @param nw An object of class `network`.
#' @param arglist A list of arguments as specified in the `ergm.userterms`
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        `ergm.userterms` package framework.
#'
#' @details
#' This ERGM user term was written to allow for age-based homophily in
#' partnership formation that is asymmetric by sex. The `absdiff` component
#' targets age-based homophily while the `by` component allows that to be
#' structured by a binary attribute such as "male", in order to enforce an
#' offset in the average difference. This allows, for example, a average age
#' difference in partnerships, but with males (on average) older than females.
#'
#' @aliases absdiffby
#'
InitErgmTerm.absdiffby <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw,
                      arglist,
                      directed = FALSE,
                      bipartite = FALSE,
                      varnames = c("attr", "by", "assym"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_VATTR_SPEC, "numeric"),
                      required = c(TRUE, TRUE, TRUE),
                      defaultvalues = list(NULL, NULL, NULL))

  nodecov <- ergm_get_vattr(a$attr, nw, accept = "numeric")
  nodeby <- ergm_get_vattr(a$by, nw)
  coef.names <- paste("absdiffby", attr(nodecov, "name"),
                      attr(nodeby, "name"), sep = ".")

  list(name = "absdiffby",
    coef.names = coef.names,
    pkgname = "EpiModel",
    inputs = c(a$assym, nodecov, nodeby),
    dependence = FALSE,
    emptynwstats = 0
  )
}


#' @title Definition for fuzzynodematch ERGM Term
#'
#' @description This function defines and initializes the fuzzynodematch ERGM
#'              term that allows for generalized homophily.
#'
#' @param nw An object of class `network`.
#' @param arglist A list of arguments as specified in the `ergm.userterms`
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        `ergm.userterms` package framework.
#'
#' @details
#' This ERGM user term was written to allow for generalized homophily.The
#' `attr` term argument should specify a character vertex attribute
#' encoding the "venues" associated to each node.  The `split` argument
#' should specify a string that separates different "venues" in the attribute
#' value for each node, as handled by `strsplit` with `fixed = TRUE`.
#' For example, if `split` is `"|"` (the default), and the attribute
#' value for a given node is `"a12|b476"`, then the associated venues for
#' this node are `"a12"` and `"b476"`.  The empty string `""` is
#' interpreted as "no venues".
#'
#' If the `binary` term argument is `FALSE` (the default), the change
#' statistic for an on-toggle is the number of unique venues associated to both
#' nodes (informally speaking, this could be described as the number of venues
#' on which the two nodes "match"); if `binary` is `TRUE`, the change
#' statistic for an on-toggle is `1` if any venue is associated to both
#' nodes, and `0` otherwise.
#'
#' @aliases fuzzynodematch
#'
InitErgmTerm.fuzzynodematch <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "split", "binary"),
                      vartypes = c(ERGM_VATTR_SPEC, "character", "logical"),
                      defaultvalues = list(NULL, "|", FALSE),
                      required = c(TRUE, FALSE, FALSE))

  nodecov <- ergm_get_vattr(a$attr, nw, accept = "character")
  venues <- strsplit(nodecov, split = a$split, fixed = TRUE)

  ## drop "" from venues and enforce uniqueness of venues for each node
  venues <- lapply(venues, function(x) unique(x[nchar(x) > 0L]))

  ## record number of venues and offset in position for each node
  lengths <- unlist(lapply(venues, length))
  positions <- cumsum(lengths) - lengths

  ## convert venues to vector
  venues <- unlist(venues)

  ## convert venues from strings to integers
  levels <- sort(unique(venues))
  venues <- match(venues, levels)

  ## sort venues for each node
  venues <- unlist(lapply(seq_len(network.size(nw)), function(i) sort(venues[positions[i] + seq_len(lengths[i])])))

  binary <- a$binary

  list(name = "fuzzynodematch",
       coef.names = paste("fuzzynodematch", attr(nodecov, "name"),
                          binary, sep = "."),
       binary = as.integer(binary),
       venues = as.integer(venues),
       lengths = as.integer(lengths),
       positions = as.integer(positions),
       dependence = FALSE,
       minval = 0)
}


#' @title Definition for edist ERGM Term
#'
#' @description This function defines and initializes the edist ERGM term that
#'              adds one statistic per edge equal to the Euclidean distance
#'              between the two incident nodes in a space defined by two or more
#'              numeric nodal coordinates (e.g., latitude and longitude), raised
#'              to one or more powers.
#'
#' @param nw An object of class `network`.
#' @param arglist A list of arguments as specified in the `ergm.userterms`
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        `ergm.userterms` package framework.
#'
#' @details
#' For an edge between nodes \eqn{i} and \eqn{j} with numeric coordinate
#' attributes \eqn{x_{i}} (each a vector over the supplied `attr` names), the
#' per-edge statistic is the Euclidean distance
#' \eqn{d_{ij} = \sqrt{\sum_k (x_{ik} - x_{jk})^2}} raised to the power `pow`,
#' i.e., \eqn{d_{ij}^{\mathrm{pow}}}. The term statistic is the sum of
#' \eqn{d_{ij}^{\mathrm{pow}}} over all edges.
#'
#' Unlike an equivalent `edgecov` specification, `edist` computes the distance
#' on the fly from the length-\eqn{n} coordinate vectors and never materializes
#' the \eqn{n \times n} distance matrix, so its memory cost is linear rather
#' than quadratic in network size. Because the statistic is dyad-independent,
#' the model contribution to the log-odds of an edge is
#' \eqn{\theta \, d_{ij}^{\mathrm{pow}}}.
#'
#' The `pow` argument controls the functional form of that dependence. With
#' `pow = 1` (the default) the log-odds are linear in Euclidean distance. Other
#' values make the log-odds non-linear in distance: for example `pow = 2`
#' yields squared Euclidean distance (a Gaussian-type decay in tie probability),
#' and larger powers produce sharper declines in tie probability as distance
#' increases. Supplying a vector of powers adds one statistic per power, so a
#' single term can capture a polynomial in distance (e.g., `pow = c(1, 2)` fits
#' \eqn{\theta_1 d_{ij} + \theta_2 d_{ij}^2}).
#'
#' @aliases edist
#'
InitErgmTerm.edist <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      directed = FALSE,
                      bipartite = FALSE,
                      varnames = c("attr", "pow"),
                      vartypes = c(ERGM_VATTR_SPEC, "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))

  coords <- ergm_get_vattr(a$attr, nw, accept = "numeric", multiple = "matrix")
  coords <- as.matrix(coords)
  ndim <- ncol(coords)
  if (ndim < 2) {
    ergm_Init_abort("The 'edist' term requires two or more numeric nodal ",
                    "attributes (e.g., 'lat' and 'long'); use 'absdiff' for a ",
                    "single dimension.")
  }
  cname <- attr(coords, "name")
  pow <- a$pow

  suffix <- ifelse(pow == 1, "", paste0(".pow", pow))
  coef.names <- paste0("edist.", cname, suffix)

  list(name = "edist",
       coef.names = coef.names,
       pkgname = "EpiModel",
       inputs = c(ndim, length(pow), pow, as.vector(coords)),
       dependence = FALSE,
       emptynwstats = numeric(length(pow)))
}
