% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net.mod.vital.R
\name{arrivals.net}
\alias{arrivals.net}
\title{Arrivals: netsim Module}
\usage{
arrivals.net(dat, at)
}
\arguments{
\item{dat}{Main \code{netsim_dat} object containing a \code{networkDynamic}
object and other initialization information passed from
\code{\link{netsim}}.}

\item{at}{Current time step.}
}
\value{
The updated \code{netsim_dat} main list object.
}
\description{
This function simulates new arrivals into the network
for use in \code{\link{netsim}} simulations.
}
\seealso{
\code{\link{netsim}}
}
\keyword{internal}
\keyword{netMod}
