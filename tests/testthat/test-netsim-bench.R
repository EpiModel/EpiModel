context("netsim timing guard")

test_that("the reference case stays within 1.5 times its baseline step time", {
  skip_on_cran()
  # The baseline is from Apple Silicon, and CI runners are not measured yet
  skip_on_ci()
  skip_on_os(c("windows", "linux"))

  # Time per step of `inst/bench/sis-10k.R` over its CPU reference, on an M5
  # Max laptop on 11/11/26. Lower it when a change makes `netsim()` faster.
  baseline <- 0.40

  bench <- new.env()
  bench_script <- system.file("bench", "sis-10k.R", package = "EpiModel")
  bench_output <- capture.output(source(bench_script, local = bench))

  expect_match(bench_output, "resim_nets.FUN", fixed = TRUE, all = FALSE)
  expect_lt(bench$step_ms / bench$cpu_reference_ms, 1.5 * baseline)
})
