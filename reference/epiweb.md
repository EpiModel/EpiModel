# EpiModel Web

Runs a web browser-based GUI of deterministic compartmental models and
stochastic individual contact models.

## Usage

``` r
epiweb(class, ...)
```

## Arguments

- class:

  Model class, with options of `"dcm"` and `"icm"`.

- ...:

  Additional arguments passed to
  [`shiny::runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Details

`epiweb` runs a web-based GUI of one-group deterministic compartmental
models and stochastic individual contact models with user input on model
type, state sizes, and parameters. Model output may be plotted,
summarized, and saved as raw data using the core `EpiModel`
functionality for these model classes. These applications are built
using the `shiny` package framework.

## References

RStudio. shiny: Web Application Framework for R. R package version
1.0.5. 2015. <https://shiny.posit.co/>.

## See also

[`dcm`](http://epimodel.github.io/EpiModel/reference/dcm.md),
[`icm`](http://epimodel.github.io/EpiModel/reference/icm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Deterministic compartmental models
epiweb(class = "dcm")

## Stochastic individual contact models
epiweb(class = "icm")
} # }
```
