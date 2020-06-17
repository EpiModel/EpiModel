
#' @title Deterministic Compartmental Model Functions
#'
#' @description These functions parameterize the base deterministic
#'              compartmental models solved using the \code{dcm} function.
#'
#' @param t Time vector, passed into model function internally through
#'        \code{\link{dcm}} via the control settings in \code{\link{control.dcm}}.
#' @param t0 Initial conditions for model, passed into model function internally
#'        through \code{\link{dcm}} via the initial conditions in
#'        \code{\link{init.dcm}}.
#' @param parms Model parameters, passed into model function internally through
#'        \code{\link{dcm}} via the parameter settings in \code{\link{param.dcm}}.
#'
#' @details
#' This help page shows the names of all the base deterministic compartmental
#' model functions supported in EpiModel. Base models are those already
#' programmed interally within the software. The model functions may be printed
#' to see their internal structure, either directly on the console or by using
#' the \code{print.mod} argument in \code{\link{control.dcm}}.
#'
#' The naming convention for the models listed here follows the format:
#' \code{mod_<disease type>_<number of groups>_<vital dynamics>}. The supported
#' disease types are SI, SIS, and SIR; the number of groups are 1 or 2; and the
#' vital dynamic options are closed (fixed population composition) or open (with
#' arrivals and departures).
#' @name dcm.mods
#'
NULL


# SI, 1 group, closed pop -------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SI_1g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num <- s.num + i.num
    
    # Parameters
    lambda <- inf.prob * act.rate * i.num / num
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda <- lambda * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda * s.num
    
    # ODEs
    dS <- -si.flow
    dI <- si.flow
    
    # Output
    list(c(dS, dI,
           si.flow),
         num = num)
  })
}


# SI, 1 group, open pop ---------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SI_1g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num <- s.num + i.num
    
    # Parameters
    lambda <- inf.prob * act.rate * i.num / num
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda <- lambda * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda * s.num
    a.flow <- a.rate * num
    ds.flow <- ds.rate * s.num
    di.flow <- di.rate * i.num
    
    # ODEs
    dS <- -si.flow + a.flow - ds.flow
    dI <- si.flow - di.flow
    
    # Output
    list(c(dS, dI,
           si.flow, a.flow, ds.flow, di.flow),
         num = num)
  })
}


# SI, 2 group, closed pop -------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SI_2g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num.g1 <- s.num + i.num
    num.g2 <- s.num.g2 + i.num.g2
    
    # Act Balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }
    
    # Group Lambdas
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda.g1 <- lambda.g1 * (1 - inter.eff)
      lambda.g2 <- lambda.g2 * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda.g1 * s.num
    si.flow.g2 <- lambda.g2 * s.num.g2
    
    # ODEs
    dSm1 <- -si.flow
    dIm1 <-  si.flow
    dSm2 <- -si.flow.g2
    dIm2 <-  si.flow.g2
    
    # Output
    list(c(dSm1, dIm1, dSm2, dIm2,
           si.flow, si.flow.g2),
         num = num.g1, num.g2 = num.g2)
  })
}


# SI, 2 group, open pop ---------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SI_2g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num.g1 <- s.num + i.num
    num.g2 <- s.num.g2 + i.num.g2
    
    # Act Balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }
    
    # Group Lambdas
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda.g1 <- lambda.g1 * (1 - inter.eff)
      lambda.g2 <- lambda.g2 * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda.g1 * s.num
    si.flow.g2 <- lambda.g2 * s.num.g2
    if (is.na(a.rate.g2)) {
      a.flow <- a.rate * num.g1
      a.flow.g2 <- a.rate * num.g1
    } else {
      a.flow <- a.rate * num.g1
      a.flow.g2 <- a.rate.g2 * num.g2
    }
    ds.flow <- ds.rate * s.num
    ds.flow.g2 <- ds.rate.g2 * s.num.g2
    di.flow <- di.rate * i.num
    di.flow.g2 <- di.rate.g2 * i.num.g2
    
    # ODEs
    dSm1 <- -si.flow + a.flow - ds.flow
    dIm1 <-  si.flow - di.flow
    dSm2 <- -si.flow.g2 + a.flow.g2 - ds.flow.g2
    dIm2 <-  si.flow.g2 - di.flow.g2
    
    # Output
    list(c(dSm1, dIm1, dSm2, dIm2,
           si.flow, a.flow, ds.flow, di.flow,
           si.flow.g2, a.flow.g2, ds.flow.g2, di.flow.g2),
         num = num.g1,
         num.g2 = num.g2)
  })
}


# SIR, 1 group, closed pop ------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIR_1g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num <- s.num + i.num + r.num
    
    # Parameters
    lambda <- inf.prob * act.rate * i.num / num
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda <- lambda * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda * s.num
    ir.flow <- rec.rate * i.num
    
    # ODEs
    dS <- -si.flow
    dI <- si.flow - ir.flow
    dR <- ir.flow
    
    # Output
    list(c(dS, dI, dR,
           si.flow, ir.flow),
         num = num)
  })
}


# SIR, 1 group, open pop --------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIR_1g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num <- s.num + i.num + r.num
    
    # Parameters
    lambda <- inf.prob * act.rate * i.num / num
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda <- lambda * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda * s.num
    ir.flow <- rec.rate * i.num
    a.flow <- a.rate * num
    ds.flow <- ds.rate * s.num
    di.flow <- di.rate * i.num
    dr.flow <- dr.rate * r.num
    
    # ODEs
    dS <- -si.flow + a.flow - ds.flow
    dI <- si.flow - ir.flow - di.flow
    dR <- ir.flow - dr.flow
    
    # Output
    list(c(dS, dI, dR,
           si.flow, ir.flow, a.flow,
           ds.flow, di.flow, dr.flow),
         num = num)
  })
}


# SIR, 2 group, closed pop ------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIR_2g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num.g1 <- s.num + i.num + r.num
    num.g2 <- s.num.g2 + i.num.g2 + r.num.g2
    
    # Act Balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }
    
    # Group Lambdas
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda.g1 <- lambda.g1 * (1 - inter.eff)
      lambda.g2 <- lambda.g2 * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda.g1 * s.num
    si.flow.g2 <- lambda.g2 * s.num.g2
    ir.flow <- rec.rate * i.num
    ir.flow.g2 <- rec.rate.g2 * i.num.g2
    
    # ODEs
    dSm1 <- -si.flow
    dIm1 <- si.flow - ir.flow
    dRm1 <- ir.flow
    dSm2 <- -si.flow.g2
    dIm2 <- si.flow.g2 - ir.flow.g2
    dRm2 <- ir.flow.g2
    
    # Output
    list(c(dSm1, dIm1, dRm1, dSm2, dIm2, dRm2,
           si.flow, ir.flow, si.flow.g2, ir.flow.g2),
         num = num.g1,
         num.g2 = num.g2)
  })
}


# SIR, 2 group, open pop --------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIR_2g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num.g1 <- s.num + i.num + r.num
    num.g2 <- s.num.g2 + i.num.g2 + r.num.g2
    
    # Act Balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }
    
    # Group Lambdas
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda.g1 <- lambda.g1 * (1 - inter.eff)
      lambda.g2 <- lambda.g2 * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda.g1 * s.num
    si.flow.g2 <- lambda.g2 * s.num.g2
    ir.flow <- rec.rate * i.num
    ir.flow.g2 <- rec.rate.g2 * i.num.g2
    if (is.na(a.rate.g2)) {
      a.flow <- a.rate * num.g1
      a.flow.g2 <- a.rate * num.g1
    } else {
      a.flow <- a.rate * num.g1
      a.flow.g2 <- a.rate.g2 * num.g2
    }
    ds.flow <- ds.rate * s.num
    ds.flow.g2 <- ds.rate.g2 * s.num.g2
    di.flow <- di.rate * i.num
    di.flow.g2 <- di.rate.g2 * i.num.g2
    dr.flow <- dr.rate * r.num
    dr.flow.g2 <- dr.rate.g2 * r.num.g2
    
    # ODEs
    dSm1 <- -si.flow + a.flow - ds.flow
    dIm1 <- si.flow - ir.flow - di.flow
    dRm1 <- ir.flow - dr.flow
    dSm2 <- -si.flow.g2 + a.flow.g2 - ds.flow.g2
    dIm2 <- si.flow.g2 - ir.flow.g2 - di.flow.g2
    dRm2 <- ir.flow.g2 - dr.flow.g2
    
    # Output
    list(c(dSm1, dIm1, dRm1, dSm2, dIm2, dRm2,
           si.flow, ir.flow, a.flow, ds.flow, di.flow, dr.flow,
           si.flow.g2, ir.flow.g2, a.flow.g2, ds.flow.g2,
           di.flow.g2, dr.flow.g2),
         num = num.g1, num.g2 = num.g2)
  })
}


# SIS, 1 group, closed pop ------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIS_1g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num <- s.num + i.num
    
    # Parameters
    lambda <- inf.prob * act.rate * i.num / num
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda <- lambda * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda * s.num
    is.flow <- rec.rate * i.num
    
    # ODEs
    dS <- -si.flow + is.flow
    dI <- si.flow - is.flow
    
    # Output
    list(c(dS, dI, si.flow, is.flow),
         num = num)
  })
}


# SIS, 1 group, open pop --------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIS_1g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num <- s.num + i.num
    
    # Parameters
    lambda <- inf.prob * act.rate * i.num / num
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda <- lambda * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda * s.num
    is.flow <- rec.rate * i.num
    a.flow <- a.rate * num
    ds.flow <- ds.rate * s.num
    di.flow <- di.rate * i.num
    
    # ODEs
    dS <- -si.flow + is.flow + a.flow - ds.flow
    dI <- si.flow - is.flow - di.flow
    
    # Output
    list(c(dS, dI, si.flow, is.flow, a.flow, ds.flow, di.flow),
         num = num)
  })
}


# SIS, 2 group, closed pop ------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIS_2g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num.g1 <- s.num + i.num
    num.g2 <- s.num.g2 + i.num.g2
    
    # Act Balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }
    
    # Group Lambdas
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda.g1 <- lambda.g1 * (1 - inter.eff)
      lambda.g2 <- lambda.g2 * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda.g1 * s.num
    is.flow <- rec.rate * i.num
    si.flow.g2 <- lambda.g2 * s.num.g2
    is.flow.g2 <- rec.rate.g2 * i.num.g2
    
    # ODEs
    dSm1 <- -si.flow + is.flow
    dIm1 <-  si.flow - is.flow
    dSm2 <- -si.flow.g2 + is.flow.g2
    dIm2 <-  si.flow.g2 - is.flow.g2
    
    # Output
    list(c(dSm1, dIm1, dSm2, dIm2,
           si.flow, is.flow, si.flow.g2, is.flow.g2),
         num = num.g1, num.g2 = num.g2)
  })
}


# SIS, 2 group, open pop --------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIS_2g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Derivations
    num.g1 <- s.num + i.num
    num.g2 <- s.num.g2 + i.num.g2
    
    # Act Balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }
    
    # Group Lambdas
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda.g1 <- lambda.g1 * (1 - inter.eff)
      lambda.g2 <- lambda.g2 * (1 - inter.eff)
    }
    
    # Flows
    si.flow <- lambda.g1 * s.num
    si.flow.g2 <- lambda.g2 * s.num.g2
    is.flow <- rec.rate * i.num
    is.flow.g2 <- rec.rate.g2 * i.num.g2
    if (is.na(a.rate.g2)) {
      a.flow <- a.rate * num.g1
      a.flow.g2 <- a.rate * num.g1
    } else {
      a.flow <- a.rate * num.g1
      a.flow.g2 <- a.rate.g2 * num.g2
    }
    ds.flow <- ds.rate * s.num
    ds.flow.g2 <- ds.rate.g2 * s.num.g2
    di.flow <- di.rate * i.num
    di.flow.g2 <- di.rate.g2 * i.num.g2
    
    # ODEs
    dSm1 <- -si.flow + is.flow + a.flow - ds.flow
    dIm1 <-  si.flow - is.flow - di.flow
    dSm2 <- -si.flow.g2 + is.flow.g2 + a.flow.g2 - ds.flow.g2
    dIm2 <-  si.flow.g2 - is.flow.g2 - di.flow.g2
    
    # Output
    list(c(dSm1, dIm1, dSm2, dIm2,
           si.flow, is.flow, a.flow, ds.flow, di.flow,
           si.flow.g2, is.flow.g2, a.flow.g2, ds.flow.g2, di.flow.g2),
         num = num.g1, num.g2 = num.g2)
  })
}
