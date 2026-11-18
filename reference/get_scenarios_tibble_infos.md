# Helper function to access the infos on merged scenarios `data.frame`

This function returns the list of scenario tibble files and the
corresponding scenario name present in a given directory. It is meant to
be used after
[`merge_netsim_scenarios_tibble()`](https://epimodel.github.io/EpiModel/reference/merge_netsim_scenarios_tibble.md).
Only files named `df__<scenario>.rds` are listed: other files, and
directories, are ignored.

## Usage

``` r
get_scenarios_tibble_infos(scenario_dir)
```

## Arguments

- scenario_dir:

  the directory where `merge_netsim_scenarios_tibble` saved the merged
  tibbles.

## Value

a `tibble` with two columns: `file_path` - the full path of the scenario
tibble file and `scenario_name` the associated scenario name.
