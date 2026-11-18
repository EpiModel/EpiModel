# Create a Single Sim File per Scenarios Using the Files From `netsim_scenarios`

Merges the batch files of each scenario into a single `netsim` object,
saved as `merged__<scenario>.rds` in `output_dir`.

## Usage

``` r
merge_netsim_scenarios(
  sim_dir,
  output_dir,
  keep.transmat = TRUE,
  keep.network = TRUE,
  keep.nwstats = TRUE,
  keep.other = TRUE,
  param.error = FALSE,
  keep.diss.stats = TRUE,
  truncate.at = NULL,
  keep.run = TRUE,
  keep.cumulative.edgelist = FALSE,
  keep.attr.history = TRUE
)
```

## Arguments

- sim_dir:

  The folder where the simulation files are to be stored.

- output_dir:

  The folder where the merged files will be stored.

- keep.transmat:

  If `TRUE`, keep the transmission matrices from the original `x` and
  `y` elements. Note: transmission matrices only saved when
  (`save.transmat == TRUE`).

- keep.network:

  If `TRUE`, keep the `networkDynamic` objects from the original `x` and
  `y` elements. Note: network only saved when (`tergmLite == FALSE`).

- keep.nwstats:

  If `TRUE`, keep the network statistics (as set by the
  `nwstats.formula` parameter in `control.net`) from the original `x`
  and `y` elements.

- keep.other:

  If `TRUE`, keep the other simulation elements (as set by the
  `save.other` parameter in `control.net`) from the original `x` and `y`
  elements. Elements with their own `keep.*` argument (e.g. `run`)
  follow that argument instead.

- param.error:

  If `TRUE`, if `x` and `y` have different params (in
  [`param.net()`](https://epimodel.github.io/EpiModel/reference/param.net.md))
  or controls (passed in
  [`control.net()`](https://epimodel.github.io/EpiModel/reference/control.net.md))
  an error will prevent the merge. Use `FALSE` to override that check.
  The controls for how the runs were executed rather than what they
  simulated (`nsims`, `ncores`, `verbose`, `verbose.int`, and the
  checkpoint settings) may differ.

- keep.diss.stats:

  If `TRUE`, keep `diss.stats` from the original `x` and `y` objects.

- truncate.at:

  Time step at which to left-truncate the time series.

- keep.run:

  If `TRUE`, keep the `run` sublists (as set by the `save.run` parameter
  in `control.net`) from the original `x` and `y` elements. These are
  required to restart from a simulation of the merged object in
  [`netsim()`](https://epimodel.github.io/EpiModel/reference/netsim.md),
  selected with
  [`get_sims()`](https://epimodel.github.io/EpiModel/reference/get_sims.md).

- keep.cumulative.edgelist:

  If `TRUE`, keep the cumulative edgelists (as set by the
  `save.cumulative.edgelist` parameter in `control.net`) from the
  original `x` and `y` elements. `FALSE` by default, as these grow with
  the length of the simulation.

- keep.attr.history:

  If `TRUE`, keep the attribute histories (as set by
  [`record_attr_history()`](https://epimodel.github.io/EpiModel/reference/record_attr_history.md))
  from the original `x` and `y` elements.

## Value

`NULL`, invisibly. The merged files are saved in `output_dir`.
