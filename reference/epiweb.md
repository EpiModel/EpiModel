# EpiModel Web

Launches interactive Shiny applications for exploring deterministic
compartmental models (DCMs) and stochastic individual contact models
(ICMs) of infectious disease transmission.

## Usage

``` r
epiweb(class, ...)
```

## Arguments

- class:

  Model class, with options of `"dcm"` and `"icm"`.

- ...:

  Additional arguments passed to
  [shiny::runApp](https://rdrr.io/pkg/shiny/man/runApp.html), such as
  `port`, `host`, or `launch.browser`.

## Details

`epiweb` launches a web-based graphical interface for configuring,
running, and analyzing epidemic models built on EpiModel. Each
application provides interactive controls for model parameters,
real-time (or on-demand) visualization with plotly, narrative model
summaries, and downloadable data tables. Both applications include a
**Guide** tab with a comprehensive user manual covering the underlying
theory, parameter definitions, and instructions for interpreting output.

## DCM Application (`class = "dcm"`)

The DCM app solves systems of ordinary differential equations for
one-group SI, SIR, and SIS models. Models run instantly and the
interface is fully reactive (no run button). Features include:

- **Scenario presets:** Flu-like (SIR), STI-like (SIS), and Measles-like
  (SIR) configurations that set all parameters to epidemiologically
  plausible values.

- **Interventions:** An optional mid-epidemic intervention that reduces
  the transmission probability by a user-specified efficacy, shown as a
  vertical dashed line on the plot.

- **Sensitivity analysis:** Vary any epidemic parameter across a range
  of values to compare trajectories side by side with a blue-to-red
  color ramp.

- **Vital dynamics:** Optional births and deaths to model endemic
  equilibria over longer time horizons.

- **Summary tab:** Displays the basic reproduction number (R₀), peak
  timing, cumulative infections, attack rate, and contextual
  interpretation.

- **Data tab:** Searchable table of model output with CSV download.

## ICM Application (`class = "icm"`)

The ICM app runs agent-based stochastic microsimulations for one-group
SI, SIR, and SIS models. Because ICMs are computationally more
expensive, a **Run Model** button controls when simulations execute.
Features include:

- **Stochastic visualization:** Multiple simulations are summarized with
  mean trajectory lines and interquartile range (IQR) ribbons showing
  the middle 50\\ individual trajectories.

- **Scenario presets:** The same Flu-like, STI-like, and Measles-like
  configurations as the DCM app, adapted for smaller populations
  suitable for agent-based simulation.

- **Interventions:** Same intervention system as DCM, natively supported
  by
  [param.icm](http://epimodel.github.io/EpiModel/reference/param.icm.md).

- **Vital dynamics:** Optional births and deaths, as in DCM.

- **Summary tab:** Displays R₀, epidemic timeline from the mean
  trajectory, and a stochastic variation section showing the range and
  standard deviation of peak infections and final prevalence across
  simulations.

- **Data tab:** View simulation output as means, standard deviations, or
  individual simulation values, with CSV download.

## Required Packages

Both applications require the `shiny`, `bslib`, `DT`, and `plotly`
packages. Install any missing packages with
`install.packages(c("shiny", "bslib", "DT", "plotly"))`.

## See also

[dcm](http://epimodel.github.io/EpiModel/reference/dcm.md) and
[param.dcm](http://epimodel.github.io/EpiModel/reference/param.dcm.md)
for the DCM modeling API;
[icm](http://epimodel.github.io/EpiModel/reference/icm.md) and
[param.icm](http://epimodel.github.io/EpiModel/reference/param.icm.md)
for the ICM modeling API.

## Examples

``` r
if (FALSE) { # \dontrun{
## Deterministic compartmental models
epiweb(class = "dcm")

## Stochastic individual contact models
epiweb(class = "icm")
} # }
```
