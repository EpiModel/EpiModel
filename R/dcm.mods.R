
#' @title Deterministic Compartmental Model Functions
#'
#' @description These functions parameterize the built-in deterministic
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
#' This help page shows the names of all the integrated deterministic compartmental
#' model functions supported in EpiModel. Integrated models are those already
#' programmed interally within the software. The model functions may be printed
#' to see their internal structure, either directly on the console or by using
#' the \code{print.mod} argument in \code{\link{control.dcm}}.
#'
#' The naming convention for the models listed here follows the format:
#' \code{mod_<disease type>_<number of groups>_<vital dynamics>}. The supported
#' disease types are SI, SIS, and SIR; the number of groups are 1 or 2; and the
#' vital dynamic options are closed (fixed population composition) or open (with
#' births and deaths).

#' @name dcm.mods
#'
NULL


# SI, 1 group, closed pop -------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SI_1g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num <- s.num + i.num

    # varying parameters
    lambda <- inf.prob * act.rate * i.num / num

    # main ODEs
    dS <- -lambda*s.num
    dI <- lambda*s.num

    list(c(dS, dI),
         num = num,
         si.flow = lambda*s.num)
  })
}


# SI, 1 group, open pop ---------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SI_1g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num <- s.num + i.num

    # varying parameters
    lambda <- inf.prob * act.rate * i.num / num

    # main ODEs
    dS <- -lambda*s.num + b.rate*num - ds.rate*s.num
    dI <- lambda*s.num - di.rate*i.num

    list(c(dS, dI),
         num = num,
         si.flow = lambda*s.num,
         b.flow = b.rate*num,
         ds.flow = ds.rate*s.num,
         di.flow = di.rate*i.num)
  })
}


# SI, 2 group, closed pop -------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SI_2g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num.g1 <- s.num + i.num
    num.g2 <- s.num.g2 + i.num.g2
    num <- num.g1 + num.g2

    # act rate balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }

    # group-specific foi
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1

    # main ODEs
    dSm1 <- -lambda.g1*s.num
    dIm1 <- lambda.g1*s.num

    dSm2 <- -lambda.g2*s.num.g2
    dIm2 <- lambda.g2*s.num.g2

    list(c(dSm1, dIm1,
           dSm2, dIm2),
         num = num.g1,
         num.g2 = num.g2,
         si.flow = lambda.g1*s.num,
         si.flow.g2 = lambda.g2*s.num.g2)
  })
}


# SI, 2 group, open pop ---------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SI_2g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num.g1 <- s.num + i.num
    num.g2 <- s.num.g2 + i.num.g2
    num <- num.g1 + num.g2

    # act rate balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }

    # group-specific foi
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1

    # birth rates
    if (is.na(b.rate.g2)) {
      br.g1 <- b.rate*num.g1
      br.g2 <- b.rate*num.g1
    } else {
      br.g1 <- b.rate*num.g1
      br.g2 <- b.rate*num.g2
    }

    # main ODEs
    dSm1 <- -lambda.g1*s.num + br.g1 - ds.rate*s.num
    dIm1 <- lambda.g1*s.num - di.rate*i.num

    dSm2 <- -lambda.g2*s.num.g2 + br.g2 - ds.rate.g2*s.num.g2
    dIm2 <- lambda.g2*s.num.g2 - di.rate.g2*i.num.g2

    list(c(dSm1, dIm1,
           dSm2, dIm2),
         num = num.g1,
         num.g2 = num.g2,
         si.flow = lambda.g1 * s.num,
         si.flow.g2 = lambda.g2 * s.num.g2,
         b.flow = br.g1,
         ds.flow = ds.rate * s.num,
         di.flow = di.rate * i.num,
         b.flow.g2 = br.g2,
         ds.flow.g2 = ds.rate.g2 * s.num.g2,
         di.flow.g2 = di.rate.g2 * i.num.g2)
  })
}


# SIR, 1 group, closed pop ------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIR_1g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num <- s.num + i.num + r.num

    # varying parameters
    lambda <- inf.prob * act.rate * i.num / num

    # main ODEs
    dS <- -lambda*s.num
    dI <- lambda*s.num - rec.rate*i.num
    dR <- rec.rate*i.num

    list(c(dS, dI, dR),
         num = num,
         si.flow = lambda * s.num,
         ir.flow = rec.rate * i.num)
  })
}


# SIR, 1 group, open pop --------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIR_1g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num <- s.num + i.num + r.num

    # varying parameters
    lambda <- inf.prob * act.rate * i.num / num

    # main ODEs
    dS <- -lambda*s.num + b.rate*num - ds.rate*s.num
    dI <- lambda*s.num - rec.rate*i.num - di.rate*i.num
    dR <- rec.rate*i.num - dr.rate*r.num

    list(c(dS, dI, dR),
         num = num,
         si.flow = lambda * s.num,
         ir.flow = rec.rate * i.num,
         b.flow = b.rate * num,
         ds.flow = ds.rate * s.num,
         di.flow = di.rate * i.num,
         dr.flow = dr.rate * r.num)
  })
}


# SIR, 2 group, closed pop ------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIR_2g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num.g1 <- s.num + i.num + r.num
    num.g2 <- s.num.g2 + i.num.g2 + r.num.g2
    num <- num.g1 + num.g2

    # act rate balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }

    # group-specific foi
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1

    # main ODEs
    dSm1 <- -lambda.g1*s.num
    dIm1 <- lambda.g1*s.num - rec.rate*i.num
    dRm1 <- rec.rate*i.num

    dSm2 <- -lambda.g2*s.num.g2
    dIm2 <- lambda.g2*s.num.g2 - rec.rate.g2*i.num.g2
    dRm2 <- rec.rate.g2*i.num.g2

    list(c(dSm1, dIm1, dRm1,
           dSm2, dIm2, dRm2),
         num = num.g1,
         num.g2 = num.g2,
         si.flow = lambda.g1 * s.num,
         ir.flow = rec.rate * i.num,
         si.flow.g2 = lambda.g2 * s.num.g2,
         ir.flow.g2 = rec.rate.g2 * i.num.g2)
  })
}


# SIR, 2 group, open pop --------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIR_2g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num.g1 <- s.num + i.num + r.num
    num.g2 <- s.num.g2 + i.num.g2 + r.num.g2
    num <- num.g1 + num.g2

    # act rate balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }

    # group-specific foi
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2 / num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num / num.g1

    # birth rates
    if (is.na(b.rate.g2)) {
      br.g1 <- b.rate*num.g1
      br.g2 <- b.rate*num.g1
    } else {
      br.g1 <- b.rate*num.g1
      br.g2 <- b.rate*num.g2
    }

    # main ODEs
    dSm1 <- -lambda.g1*s.num + br.g1 - ds.rate*s.num
    dIm1 <- lambda.g1*s.num - rec.rate*i.num - di.rate*i.num
    dRm1 <- rec.rate*i.num - dr.rate*r.num

    dSm2 <- -lambda.g2*s.num.g2 + br.g2 - ds.rate.g2*s.num.g2
    dIm2 <- lambda.g2*s.num.g2 - rec.rate.g2*i.num.g2 - di.rate.g2*i.num.g2
    dRm2 <- rec.rate.g2*i.num.g2 - dr.rate.g2*r.num.g2

    list(c(dSm1, dIm1, dRm1,
           dSm2, dIm2, dRm2),
         num = num.g1,
         num.g2 = num.g2,
         si.flow = lambda.g1 * s.num,
         ir.flow = rec.rate * i.num,
         si.flow.g2 = lambda.g2 * s.num.g2,
         ir.flow.g2 = rec.rate.g2 * i.num.g2,
         b.flow = br.g1,
         ds.flow = ds.rate * s.num,
         di.flow = di.rate * i.num,
         dr.flow = dr.rate * r.num,
         b.flow.g2 = br.g2,
         ds.flow.g2 = ds.rate.g2 * s.num.g2,
         di.flow.g2 = di.rate.g2 * i.num.g2,
         dr.flow.g2 = dr.rate.g2 * r.num.g2)
  })
}


# SIS, 1 group, closed pop ------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIS_1g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num <- s.num + i.num

    # varying parameters
    lambda <- inf.prob * act.rate * i.num / num

    # main ODEs
    dS <- -lambda*s.num + rec.rate*i.num
    dI <- lambda*s.num - rec.rate*i.num

    list(c(dS, dI),
         num = num,
         si.flow = lambda * s.num,
         is.flow = rec.rate * i.num)
  })
}


# SIS, 1 group, open pop --------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIS_1g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num <- s.num + i.num

    # varying parameters
    lambda <- inf.prob * act.rate * i.num / num

    # main ODEs
    dS <- -lambda*s.num + rec.rate*i.num + b.rate*num - ds.rate*s.num
    dI <- lambda*s.num - rec.rate*i.num - di.rate*i.num

    list(c(dS, dI),
         num = num,
         si.flow = lambda * s.num,
         is.flow = rec.rate * i.num,
         b.flow = b.rate * num,
         ds.flow = ds.rate * s.num,
         di.flow = di.rate * i.num)
  })
}


# SIS, 2 group, closed pop ------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIS_2g_cl <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num.g1 <- s.num + i.num
    num.g2 <- s.num.g2 + i.num.g2
    num <- num.g1 + num.g2

    # act rate balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }

    # group-specific foi
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2/num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num/num.g1

    # main ODEs
    dSm1 <- -lambda.g1*s.num + rec.rate*i.num
    dIm1 <- lambda.g1*s.num - rec.rate*i.num

    dSm2 <- -lambda.g2*s.num.g2 + rec.rate.g2*i.num.g2
    dIm2 <- lambda.g2*s.num.g2 - rec.rate.g2*i.num.g2

    list(c(dSm1, dIm1,
           dSm2, dIm2),
         num = num.g1,
         num.g2 = num.g2,
         si.flow = lambda.g1 * s.num,
         is.flow = rec.rate * i.num,
         si.flow.g2 = lambda.g2 * s.num.g2,
         is.flow.g2 = rec.rate.g2 * i.num.g2)
  })
}


# SIS, 2 group, open pop --------------------------------------------------
#' @rdname dcm.mods
#' @keywords internal
#' @export
mod_SIS_2g_op <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {

    # derived totals
    num.g1 <- s.num + i.num
    num.g2 <- s.num.g2 + i.num.g2
    num <- num.g1 + num.g2

    # act rate balancing
    if (balance == "g1") {
      ar.g1 <- act.rate
      ar.g2 <- ar.g1 * num.g1 / num.g2
    }
    if (balance == "g2") {
      ar.g2 <- act.rate.g2
      ar.g1 <- ar.g2 * num.g2 / num.g1
    }

    # group-specific foi
    lambda.g1 <- inf.prob * ar.g1 * i.num.g2/num.g2
    lambda.g2 <- inf.prob.g2 * ar.g2 * i.num/num.g1

    # birth rates
    if (is.na(b.rate.g2)) {
      br.g1 <- b.rate*num.g1
      br.g2 <- b.rate*num.g1
    } else {
      br.g1 <- b.rate*num.g1
      br.g2 <- b.rate*num.g2
    }

    # main ODEs
    dSm1 <- -lambda.g1*s.num + rec.rate*i.num + br.g1 - ds.rate*s.num
    dIm1 <- lambda.g1*s.num - rec.rate*i.num - di.rate*i.num

    dSm2 <- -lambda.g2*s.num.g2 + rec.rate.g2*i.num.g2 + br.g2 - ds.rate.g2*s.num.g2
    dIm2 <- lambda.g2*s.num.g2 - rec.rate.g2*i.num.g2 - di.rate.g2*i.num.g2

    list(c(dSm1, dIm1,
           dSm2, dIm2),
         num = num.g1,
         num.g2 = num.g2,
         si.flow = lambda.g1 * s.num,
         is.flow = rec.rate * i.num,
         si.flow.g2 = lambda.g2 * s.num.g2,
         is.flow.g2 = rec.rate.g2 * i.num.g2,
         b.flow = br.g1,
         ds.flow = ds.rate * s.num,
         di.flow = di.rate * i.num,
         ds.flow.g2 = ds.rate.g2 * s.num.g2,
         di.flow.g2 = di.rate.g2 * i.num.g2)
  })
}
