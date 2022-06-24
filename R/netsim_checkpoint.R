# Check if the `.checkpoint` controls are present and well formed
netsim_is_checkpointed <- function(control) {
  ckpt_file_check <-
    is.character(control[[".checkpoint.dir"]]) &&
    length(control[[".checkpoint.dir"]]) == 1

  ckpt_steps_check <-
    is.numeric(control[[".checkpoint.steps"]]) &&
    length(control[[".checkpoint.steps"]]) == 1 &&
    control[[".checkpoint.steps"]] > 0

  return(ckpt_file_check && ckpt_steps_check)
}

# Check the existence of a checkpoint file
netsim_is_resume_checkpoint <- function(control, s) {
  resume_checkpoint <- FALSE
  if (control[[".checkpointed"]])
    resume_checkpoint <- file.exists(netsim_get_checkpoint_filename(control, s))

  return(resume_checkpoint)
}

netsim_load_checkpoint <- function(control, s) {
  readRDS(netsim_get_checkpoint_filename(control, s))
}

netsim_save_checkpoint <- function(dat, s) {
  if (!dir.exists(get_control(dat, ".checkpoint.dir")))
    dir.create(get_control(dat, ".checkpoint.dir"), recursive = TRUE)

  checkpoint_file_name <- netsim_get_checkpoint_filename(dat[["control"]], s)
  tmp_file_name <- paste0(checkpoint_file_name, "_tmp")

  saveRDS(dat, tmp_file_name,
          compress = get_control(dat, ".checkpoint.compress"))

  if (file.exists(checkpoint_file_name)) file.remove(checkpoint_file_name)
  file.rename(tmp_file_name, checkpoint_file_name)

  invisible(dat)
}


netsim_get_checkpoint_filename <- function(control, s) {
  paste0(control[[".checkpoint.dir"]], "/", s, ".rds")
}

# Remove all checkpoint files and directories
#
netsim_clear_checkpoint <- function(control) {
  if (control[[".checkpointed"]])
      unlink(control[[".checkpoint.dir"]], recursive = TRUE)
  invisible(TRUE)
}
