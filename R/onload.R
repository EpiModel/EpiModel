.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste(c(strwrap(paste("EpiModel 2.0 NOTE: With the release of EpiModel 2.0, significant changes have been made to the EpiModel workflow, requiring updates to user code written using EpiModel version 1.x. See https://epimodel.org/ for migration assistance.",sep="")),""),collapse="\n"))
}
