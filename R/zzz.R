.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste(c(strwrap(
      paste("EpiModel NOTE: EpiModel 2.0 implements significant changes to the
            EpiModel application programming interface (API), which may require
            updates to model code written under EpiModel 1.x.
            See https://epimodel.org/ for migration guidance.", sep = "")), ""),
      collapse = "\n"))
}
