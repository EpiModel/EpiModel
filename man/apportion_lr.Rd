% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{apportion_lr}
\alias{apportion_lr}
\title{Apportion Using the Largest Remainder Method}
\usage{
apportion_lr(vector.length, values, proportions, shuffled = FALSE)
}
\arguments{
\item{vector.length}{Length for the output vector.}

\item{values}{Values for the output vector.}

\item{proportions}{Proportion distribution with one number for each value.
This must sum to 1.}

\item{shuffled}{If \code{TRUE}, randomly shuffle the order of the vector.}
}
\value{
A vector of length \code{vector.length} containing the apportioned
values from \code{values}.
}
\description{
Apportions a vector of values given a specified frequency
distribution of those values such that the length of the output
is robust to rounding and other instabilities.
}
\examples{
\dontrun{
## Example 1: Without rounding
apportioned_vec_1 <- apportion_lr(4, c(1, 2, 3, 4, 5),
                                     c(0.25, 0, 0.25, 0.25, 0.25))

## Example 2: With rounding
apportioned_vec_2 <- apportion_lr(5, c(1, 2, 3, 4, 5),
                                     c(0.21, 0, 0.29, 0.25, 0.25))
}

}
