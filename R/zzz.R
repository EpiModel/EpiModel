.onAttach <- function(...) {
  
  if (!interactive()) return()

  vcurr <- utils::packageVersion("EpiModel")
  
  msg <- c(
    
    "\n",
    
    paste0("============ Loading EpiModel ", vcurr, " ============"),
    
    "\n"
    
  )

  packageStartupMessage(msg)
}
