
##
## Testing EpiModel output
##
##  The main testing functions below check for consistency of output:
##  that the compartment sizes and flows are consistent according to
##  discrete-time difference equations
##


#' @title Write Out Test Progress to Console
#'
#' @description Writes the name of a test and \code{...} to console for
#'              showing testing progress.
#'
#' @param test Character string with the name of a test.
#'
#' @keywords internal
#' @export
#'
mcat <- function(test) { cat("\n", test, " ... ", sep = "") }


#' @title Test the Model Output from a Network Model
#'
#' @description Tests whether the model output from a network model is consistent
#'              with key balancing equations for compartment and flow sizes for
#'              each simulation for each time step.
#'
#' @param x An object of class \code{netsim}.
#'
#' @keywords internal
#' @export
#'
test_net <- function(x) {

  type <- x$control$type
  nsims <- x$control$nsims
  vital <- x$param$vital
  groups <- x$param$groups

  for (s in 1:nsims) {
    df <- as.data.frame(x, out = "vals", sim = s)
    for (i in 2:nrow(df)) {

      if (type == "SI") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i - 1] + a.flow[i] - ds.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i] + a.flow[i] - ds.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- df$num.g2[i] == df$num.g2[i - 1] + df$a.flow.g2[i] -
                                    df$ds.flow.g2[i] - df$di.flow.g2[i]
            if (test == FALSE) stop("\nFailed *num.g2* at SIM", s, "TIME", i)
            test <- df$s.num.g2[i] == df$s.num.g2[i - 1] - df$si.flow.g2[i] +
                                      df$a.flow.g2[i] - df$ds.flow.g2[i]
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
      }

      if (type == "SIR") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - ir.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          test <- with(df, r.num[i] == r.num[i - 1] + ir.flow[i])
          if (test == FALSE) stop("\nFailed *r.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] - ir.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
            test <- with(df, r.num.g2[i] == r.num.g2[i - 1] + ir.flow.g2[i])
            if (test == FALSE) stop("\nFailed *r.num.g2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i - 1] + a.flow[i] - ds.flow[i] -
                         di.flow[i] - dr.flow[i])
          if (test == FALSE) stop("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i] + a.flow[i] - ds.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - ir.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          test <- with(df, r.num[i] == r.num[i - 1] + ir.flow[i] - dr.flow[i])
          if (test == FALSE) stop("\nFailed *r.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, num.g2[i] == num.g2[i - 1] + a.flow.g2[i] - ds.flow.g2[i] -
                                          di.flow.g2[i] - dr.flow.g2[i])
            if (test == FALSE) stop("\nFailed *num* at SIM", s, "TIME", i)
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i] +
                           a.flow.g2[i] - ds.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] -
                                            ir.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
            test <- with(df, r.num.g2[i] == r.num.g2[i - 1] + ir.flow.g2[i] - dr.flow.g2[i])
            if (test == FALSE) stop("\nFailed *r.num* at SIM", s, "TIME", i)
          }
        }
      }

      if (type == "SIS") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i] + is.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - is.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i] + is.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] - is.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i - 1] + a.flow[i] - ds.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i] + is.flow[i] +
                         a.flow[i] - ds.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - is.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, num.g2[i] == num.g2[i - 1] + a.flow.g2[i] -
                           ds.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) stop("\nFailed *num* at SIM", s, "TIME", i)
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i] +
                           is.flow.g2[i] + a.flow.g2[i] - ds.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] -
                           is.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          }
        }
      }

    }
  }

}


#' @title Test the Model Output from a Stochastic Individual Contact Model
#'
#' @description Tests whether the model output from an individual contact model
#'              is consistent with key balancing equations for compartment and
#'              flow sizes for each simulation for each time step.
#'
#' @param x An object of class \code{icm}.
#'
#' @keywords internal
#' @export
#'
test_icm <- function(x) {

  type <- x$control$type
  nsims <- x$control$nsims
  vital <- x$param$vital
  groups <- x$param$groups

  for (s in 1:nsims) {
    df <- as.data.frame(x, out = "vals", sim = s)
    for (i in 2:nrow(df)) {

      if (type == "SI") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i - 1] + a.flow[i] - ds.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i] + a.flow[i] - ds.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- df$num.g2[i] == df$num.g2[i - 1] + df$a.flow.g2[i] -
                                    df$ds.flow.g2[i] - df$di.flow.g2[i]
            if (test == FALSE) stop("\nFailed *num.g2* at SIM", s, "TIME", i)
            test <- df$s.num.g2[i] == df$s.num.g2[i - 1] - df$si.flow.g2[i] +
                                      df$a.flow.g2[i] - df$ds.flow.g2[i]
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
      }

      if (type == "SIR") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - ir.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          test <- with(df, r.num[i] == r.num[i - 1] + ir.flow[i])
          if (test == FALSE) stop("\nFailed *r.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] - ir.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
            test <- with(df, r.num.g2[i] == r.num.g2[i - 1] + ir.flow.g2[i])
            if (test == FALSE) stop("\nFailed *r.num.g2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i - 1] + a.flow[i] - ds.flow[i] -
                         di.flow[i] - dr.flow[i])
          if (test == FALSE) stop("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i] + a.flow[i] - ds.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - ir.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          test <- with(df, r.num[i] == r.num[i - 1] + ir.flow[i] - dr.flow[i])
          if (test == FALSE) stop("\nFailed *r.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, num.g2[i] == num.g2[i - 1] + a.flow.g2[i] - ds.flow.g2[i] -
                           di.flow.g2[i] - dr.flow.g2[i])
            if (test == FALSE) stop("\nFailed *num.g2* at SIM", s, "TIME", i, sep = " ")
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i] +
                           a.flow.g2[i] - ds.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] -
                           ir.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
            test <- with(df, r.num.g2[i] == r.num.g2[i - 1] + ir.flow.g2[i] - dr.flow.g2[i])
            if (test == FALSE) stop("\nFailed *r.num.g2* at SIM", s, "TIME", i)
          }
        }
      }

      if (type == "SIS") {
        if (vital == FALSE) {
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i] + is.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - is.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i] + is.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] - is.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
        if (vital == TRUE) {
          test <- with(df, num[i] == num[i - 1] + a.flow[i] - ds.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *num* at SIM", s, "TIME", i)
          test <- with(df, s.num[i] == s.num[i - 1] - si.flow[i] + is.flow[i] +
                         a.flow[i] - ds.flow[i])
          if (test == FALSE) stop("\nFailed *s.num* at SIM", s, "TIME", i)
          test <- with(df, i.num[i] == i.num[i - 1] + si.flow[i] - is.flow[i] - di.flow[i])
          if (test == FALSE) stop("\nFailed *i.num* at SIM", s, "TIME", i)
          if (groups == 2) {
            test <- with(df, num.g2[i] == num.g2[i - 1] + a.flow.g2[i] -
                           ds.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) stop("\nFailed *num.g2* at SIM", s, "TIME", i)
            test <- with(df, s.num.g2[i] == s.num.g2[i - 1] - si.flow.g2[i] +
                           is.flow.g2[i] + a.flow.g2[i] - ds.flow.g2[i])
            if (test == FALSE) stop("\nFailed *s.num.g2* at SIM", s, "TIME", i)
            test <- with(df, i.num.g2[i] == i.num.g2[i - 1] + si.flow.g2[i] -
                           is.flow.g2[i] - di.flow.g2[i])
            if (test == FALSE) stop("\nFailed *i.num.g2* at SIM", s, "TIME", i)
          }
        }
      }

    }
  }

}
