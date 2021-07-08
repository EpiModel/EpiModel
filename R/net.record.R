# placeholder function
get_unique_ids <- function(dat, indexes) {
  return(indexes + 1000)
}

#' Logging values
#'
#'
record_node_value <- function(dat, at, measure, indexes, values) {
  if (is.null(dat[["node_records"]])) {
    dat[["node_records"]] <- list()
  }

  if ( length(values) != 1 && length(values) != length(indexes) ) {
    stop(
      "When trying to record a value for `", measure, "` at time step ", at,
      "The size of the `values` vector is not equal to the number of node ",
      "selected by the `indexes` vector nor of length 1. \n",
      "Expected: ", length(indexes), " or 1 \n",
      "Given: ", length(values)
    )
  }

  element <- list(at, measure, get_unique_ids(dat, indexes), values)
  names(element) <- c("step", "measure", "uids", "values")

  dat[["node_records"]] <- append(
    dat[["node_records"]],
    list(element)
  )

  return(dat)
}

record_model_value <- function(dat, at, measure, value) {
  if (is.null(dat[["model_records"]])) {
    dat[["model_records"]] <- list()
  }

  if ( length(values) != 1) {
    stop(
      "When trying to record a value for `", measure, "` at time step ", at,
      "The size of the `value` vector is not equal to 1 \n",
      "Expected: 1 \n",
      "Given: ", length(value)
    )
  }

  element <- list(at, measure, value)
  names(element) <- c("step", "measure", "value")

  dat[["model_records"]] <- append(
    dat[["model_records"]],
    list(element)
  )

  return(dat)
}

# make one for model and one for node
records2data.frame <- function(dat) {
  measures <- lapply(dat[["records"]], function(x) x[["measure"]])
  measure.names <- unique(measures)

  dfs <- vector(mode = "list", length(measure.names))
  names(dfs) <- measure.names

  for (m in measure.names) {
    # parts <- dat[["records"]][measure == m]
    parts <- Filter(function(x) x[["measure"]] == m, dat[["records"]])
    parts <- lapply(parts, as.data.frame)
    dfs[[m]] <- do.call("rbind", parts)
  }

  return(dfs)
}


## I will need a way to concat the DFs when there is multiple sims
##  - simply add simno ? like saveout?
##
## I will need a `get_uid` function

# dat <- list()
# dat <- record_node_value(dat, 10, "attra", 1:5, "a")
# dat <- record_node_value(dat, 11, "attra", 1:5, "b")
# dat <- record_node_value(dat, 11, "attrb", 6:10, 1:5)
# dat <- record_node_value(dat, 12, "attrb", 4:8, 1:5)
# dat <- record_node_value(dat, 12, "attrc", 4:8, TRUE)

# dat$records

# records2data.frame(dat)

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
