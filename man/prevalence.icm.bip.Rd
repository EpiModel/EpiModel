% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/icm.utils.R
\name{prevalence.icm.bip}
\alias{prevalence.icm.bip}
\title{Get Epidemic Output from icm Model}
\usage{
prevalence.icm.bip(dat, at)
}
\arguments{
\item{dat}{Main \code{icm_dat} class data object passed through \code{icm}
simulations.}

\item{at}{Current time step.}
}
\value{
The updated \code{icm_dat} class main data object.
}
\description{
This function provides all active model state sizes from
the network at the specified time step, output to a list of
vectors.
}
\keyword{internal}
