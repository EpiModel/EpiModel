
#' @title Definition for absdiffnodemix ERGM Term
#'
#' @description This function defines and initialize the absdiffnodemix ERGM term
#'              that allows for targeting age homophily by race.
#'
#' @param nw An object of class `network`.
#' @param arglist A list of arguments as specified in the `ergm.userterms`
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        `ergm.userterms` package framework.
#'
#' @details
#' This ERGM user term was written to allow for age-based homophily in partnership
#' formation that is heterogenous by race. The absdiff component allows targets
#' the distribution of age mixing on that continuous variable, and the nodemix
#' component differentiates this for black-black, black-white, and white-white
#' couples.
#'
#' @author Steven M. Goodreau
#'
#' @export
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

  uui <- matrix(1:length(ui) ^ 2, length(ui), length(ui))
  urm <- t(sapply(ui, rep, length(ui)))
  ucm <- sapply(ui, rep, length(ui))
  uun <- outer(u, u, paste, sep = ".")
  uui <- uui[upper.tri(uui, diag = TRUE)]
  urm <- urm[upper.tri(urm, diag = TRUE)]
  ucm <- ucm[upper.tri(ucm, diag = TRUE)]
  uun <- uun[upper.tri(uun, diag = TRUE)]

  inputs = c(length(nodecov), length(urm), nodecov, nodecovby, urm, ucm)

  list(name = "absdiffnodemix",
       coef.names = paste("absdiffnodemix", attr(nodecov, "name"), nodecovbyname, uun, sep = "."),
       pkgname = "EpiModel",
       inputs = inputs,
       dependence = FALSE)
}


#' @title Definition for absdiffby ERGM Term
#'
#' @description This function defines and initialize the absdiffby ERGM term
#'              that allows for targeting age homophily by sex.
#'
#' @param nw An object of class `network`.
#' @param arglist A list of arguments as specified in the `ergm.userterms`
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        `ergm.userterms` package framework.
#'
#' @details
#' This ERGM user term was written to allow for age-based homophily in partnership
#' formation that is asymetric by sex. The absdiff component targets age homophily
#' while the by component allows that to be structed by a binary attribute such
#' as "male", in order to enforce an offset in the average difference.
#'
#' @author Samuel M. Jenness
#'
#' @export
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
  coef.names <- paste("absdiffby", attr(nodecov, "name"), attr(nodeby, "name"), sep = ".")

  list(name = "absdiffby",
       coef.names = coef.names,
       pkgname = "EpiModel",
       inputs = c(a$assym, nodecov, nodeby),
       dependence = FALSE,
       emptynwstats = 0
  )
}
