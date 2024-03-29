% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/icm.mod.status.R
\name{infection.icm}
\alias{infection.icm}
\title{Primary Infection Module for icm}
\usage{
infection.icm(dat, at)
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
This function simulates the main infection process given the
current state of the actors in the system.
}
\keyword{internal}
