#' Logging values
#'
#'
log_val <- function(dat, at, attr, indexes, values) {
  if (is.null(dat[["logs"]])) {
    dat[["logs"]] <- list()
  }

  if ( length(values) != 1 && length(values) != length(indexes) ) {
    stop(
      "When trying to log a value for `", attr, "` at time step ", at,
      "The size of the `values` vector is not equal to the number of node ",
      "selected by the `indexes` vector nor of length 1. \n",
      "Expected: ", length(indexes), " or 1 \n",
      "Given: ", length(values)
    )
  }

  element <- list(at, attr, get_unique_ids(dat, indexes), values)
  names(element) <- c("step", "attribute", "uids", "values")

  dat[["logs"]] <- append(
    dat[["logs"]],
    list(element)
  )

  return(dat)
}

logs2dfs <- function(dat) {
  attrs <- lapply(dat$logs, function(x) x$attribute)
  attrs_names <- unique(attrs)

 dfs <- vector(mode = "list", length(attrs_names))
  names(dfs) <- attrs_names

  for (a in attrs_names) {
    parts <- dat$logs[attrs == a]
    # parts <- Filter(function(x) x$attribute == a, dat$logs)
    parts <- lapply(parts, as.data.frame)
    dfs[[a]] <- do.call("rbind", parts)
  }

  return(dfs)
}


## I will need a way to concat the DFs when there is multiple sims
##
## I will need a `get_uid` function

#dat <- list()
#dat <- log_val(dat, 10, "attra", 1:5, "a")
#dat <- log_val(dat, 11, "attra", 1:5, "b")
#dat <- log_val(dat, 11, "attrb", 6:10, 1:5)
#dat <- log_val(dat, 12, "attrb", 4:8, 1:5)
#dat <- log_val(dat, 12, "attrc", 4:8, TRUE)

#dat$logs

#logs2dfs(dat)

#rng_logs <- function() {
#  n <- sample(1:300, 1)
#  step <- sample(1:1000, 1)
#  uids <- sample(1:300, n)
#  attribute <- sample(LETTERS, 1)
#  values <- rbinom(1, 1, 0.5)
#  values <- if (values == 1) sample(LETTERS, 1) else sample(LETTERS, n, TRUE)

#  list(step, attribute, uids, values)
#}

#log_list <- replicate(1e4, list(rng_logs()))

#dat <- list()
#for (ll in log_list) {
#  dat <- log_val(dat, ll[[1]], ll[[2]], ll[[3]], ll[[4]])
#}
#df_logs <- logs2dfs(dat)

#print(object.size(dat$logs), units = "MiB")
#print(object.size(df_logs), units = "MiB")
